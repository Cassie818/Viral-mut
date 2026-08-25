#!/usr/bin/env python3
"""Aggregate CaLM codon probabilities to amino-acid probabilities.

This addresses the probability-space concern that codon-level probabilities
can split an amino-acid probability mass across synonymous codons.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import matplotlib.pyplot as plt
import statsmodels.api as sm
from calm import CaLM
from calm.sequence import CodonSequence
from scipy.stats import fisher_exact
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests


DEFAULT_INPUT = Path(
    "Results/Revision/len1022_aa_aggregation/"
    "len1022_model_control_score_table_complete_cases_for_aa_aggregation.csv"
)
DEFAULT_WEIGHTS = Path("/Users/cassie/Desktop/genome-protein-mut/CaLM/calm/calm_weights/calm_weights.ckpt")
DEFAULT_GENE_DIR = Path("/Users/cassie/Desktop/Gene")
DEFAULT_OUT_DIR = Path("Results/Revision/len1022_aa_aggregation")
FIG_DIR = Path("Figure")

AA_TO_CODONS = {
    "K": ["AAA", "AAG"],
    "N": ["AAU", "AAC"],
    "I": ["AUA", "AUU", "AUC"],
    "M": ["AUG"],
    "T": ["ACA", "ACU", "ACC", "ACG"],
    "R": ["AGA", "AGG", "CGA", "CGU", "CGC", "CGG"],
    "S": ["AGU", "AGC", "UCA", "UCU", "UCC", "UCG"],
    "*": ["UAA", "UAG", "UGA"],
    "Y": ["UAU", "UAC"],
    "L": ["UUA", "UUG", "CUA", "CUU", "CUC", "CUG"],
    "F": ["UUU", "UUC"],
    "C": ["UGU", "UGC"],
    "W": ["UGG"],
    "Q": ["CAA", "CAG"],
    "H": ["CAU", "CAC"],
    "P": ["CCA", "CCU", "CCC", "CCG"],
    "E": ["GAA", "GAG"],
    "D": ["GAU", "GAC"],
    "V": ["GUA", "GUU", "GUC", "GUG"],
    "A": ["GCA", "GCU", "GCC", "GCG"],
    "G": ["GGA", "GGU", "GGC", "GGG"],
}
CODON_TO_AA = {codon: aa for aa, codons in AA_TO_CODONS.items() for codon in codons}
AA_DEGENERACY = {aa: len(codons) for aa, codons in AA_TO_CODONS.items()}
PATHOGENIC_LABELS = {"pathogenic", "likely_pathogenic"}
BENIGN_LABELS = {"benign", "likely_benign"}


class CaLMProb(CaLM):
    def codon_probabilities(self, sequence: str) -> np.ndarray:
        seq = CodonSequence(sequence)
        tokens = self.tokenize(seq)
        with torch.no_grad():
            logits = self.model(tokens)["logits"]
            probs = F.softmax(logits, dim=-1)
        return probs.detach().cpu().numpy()[0, 1:-1, :]


def read_fasta(path: Path) -> str:
    seq_parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line.upper().replace("T", "U"))
    return "".join(seq_parts)


def standardize_codon(codon: str) -> str:
    return str(codon).upper().replace("T", "U")


def score_gene(calm: CaLMProb, gene: str, group: pd.DataFrame, gene_dir: Path) -> list[dict[str, object]]:
    sequence = read_fasta(gene_dir / f"{gene}.fasta")
    probs = calm.codon_probabilities(sequence)
    tok_to_idx = calm.alphabet.tok_to_idx
    aa_indices = {
        aa: [tok_to_idx[codon] for codon in codons]
        for aa, codons in AA_TO_CODONS.items()
    }
    aa_probs = {
        aa: probs[:, idxs].sum(axis=1)
        for aa, idxs in aa_indices.items()
    }

    rows: list[dict[str, object]] = []
    eps = np.finfo(np.float64).tiny
    for _, row in group.iterrows():
        site = int(row["Site_gene"]) - 1
        ref_codon = standardize_codon(row["Ref_gene"])
        mut_codon = standardize_codon(row["Mut_gene"])
        ref_aa = str(row["Ref_prot"])
        mut_aa = str(row["Mut_prot"])

        if site < 0 or site >= probs.shape[0]:
            continue
        if ref_codon not in tok_to_idx or mut_codon not in tok_to_idx:
            continue
        if ref_aa not in AA_TO_CODONS or mut_aa not in AA_TO_CODONS:
            continue

        codon_llr = np.log(max(float(probs[site, tok_to_idx[mut_codon]]), eps)) - np.log(
            max(float(probs[site, tok_to_idx[ref_codon]]), eps)
        )
        aa_agg_llr = np.log(max(float(aa_probs[mut_aa][site]), eps)) - np.log(
            max(float(aa_probs[ref_aa][site]), eps)
        )

        out = row.to_dict()
        out.update(
            {
                "calm_codon_llr_recomputed": codon_llr,
                "calm_aa_agg_llr": aa_agg_llr,
                "calm_codon_minus_aa_agg_llr": codon_llr - aa_agg_llr,
                "ref_aa_degeneracy": AA_DEGENERACY[ref_aa],
                "mut_aa_degeneracy": AA_DEGENERACY[mut_aa],
                "log_ref_over_mut_degeneracy": np.log(AA_DEGENERACY[ref_aa] / AA_DEGENERACY[mut_aa]),
            }
        )
        rows.append(out)
    return rows


def append_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def compute_scores(args: argparse.Namespace, output: Path) -> pd.DataFrame:
    df = pd.read_csv(args.input)
    if output.exists() and not args.force:
        existing = pd.read_csv(output)
        done = set(existing["Gene_gene"].dropna().unique())
    else:
        if output.exists():
            output.unlink()
        existing = None
        done = set()

    grouped = {gene: group for gene, group in df.groupby("Gene_gene", sort=True)}
    remaining = [gene for gene in grouped if gene not in done]
    if args.sort_by_length:
        remaining = sorted(remaining, key=lambda gene: int(grouped[gene]["Site_gene"].max()))
    if args.max_remaining_genes is not None:
        remaining = remaining[: args.max_remaining_genes]
    if remaining:
        calm = CaLMProb(weights_file=str(args.weights))
        calm.model.eval()
        for idx, gene in enumerate(remaining, start=1):
            group = grouped[gene]
            try:
                rows = score_gene(calm, gene, group, args.gene_dir)
                append_rows(output, rows)
            except Exception as exc:
                print(f"FAILED {gene}: {exc}")
                continue
            if idx % args.report_every == 0 or idx == len(remaining):
                print(f"Scored {idx}/{len(remaining)} remaining genes; latest={gene}")

    return pd.read_csv(output)


def add_labels_and_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "label" not in out:
        out["label"] = out["Label_prot"].map(lambda x: 1 if x in PATHOGENIC_LABELS else 0)
    prot_col = "esm2_650m_llr" if "esm2_650m_llr" in out.columns else "LLR_prot"
    calm_codon_col = "calm_codon_llr_recomputed" if "calm_codon_llr_recomputed" in out.columns else "LLR_gene"
    out["protein_llr_for_fig5"] = out[prot_col]
    out["protein_model_for_fig5"] = "ESM-2 650M" if prot_col == "esm2_650m_llr" else "ESM-2 150M"
    out["calm_codon_llr_for_fig5"] = out[calm_codon_col]
    out["diff_codon"] = out["protein_llr_for_fig5"] - out["calm_codon_llr_for_fig5"]
    out["diff_aa_agg"] = out["protein_llr_for_fig5"] - out["calm_aa_agg_llr"]
    out["pathogenic_binary"] = out["Label_prot"].isin(PATHOGENIC_LABELS).astype(int)
    out["benign_binary"] = out["Label_prot"].isin(BENIGN_LABELS).astype(int)
    out["Pair"] = list(zip(out["Ref_prot"], out["Mut_prot"]))
    return out


def extreme_sets(df: pd.DataFrame, diff_col: str) -> tuple[pd.DataFrame, pd.DataFrame, float, float]:
    lower = float(df[diff_col].quantile(0.01))
    upper = float(df[diff_col].quantile(0.99))
    clm = pd.concat(
        [
            df[(df[diff_col] > upper) & (df["Label_prot"].isin(PATHOGENIC_LABELS))],
            df[(df[diff_col] < lower) & (df["Label_prot"].isin(BENIGN_LABELS))],
        ]
    )
    plm = pd.concat(
        [
            df[(df[diff_col] < lower) & (df["Label_prot"].isin(PATHOGENIC_LABELS))],
            df[(df[diff_col] > upper) & (df["Label_prot"].isin(BENIGN_LABELS))],
        ]
    )
    return clm, plm, lower, upper


def pair_enrichment(df: pd.DataFrame, subset: pd.DataFrame, prefix: str) -> pd.DataFrame:
    overall = Counter(df["Pair"])
    subset_counts = Counter(subset["Pair"])
    total = sum(overall.values())
    subset_total = sum(subset_counts.values())
    rest_total = total - subset_total
    rows = []
    pvals = []
    for pair in sorted(overall):
        n_total = overall[pair]
        n_subset = subset_counts[pair]
        freq_overall = n_total / total
        freq_subset = n_subset / subset_total if subset_total else np.nan
        br = freq_subset / freq_overall if freq_overall and subset_total else np.nan
        a = n_subset
        b = subset_total - n_subset
        c = n_total - n_subset
        d = rest_total - c
        _, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        pvals.append(pval)
        rows.append(
            {
                "Pair": pair,
                "Ref_prot": pair[0],
                "Mut_prot": pair[1],
                f"Count_{prefix}": n_subset,
                "Count_overall": n_total,
                f"Freq_{prefix}": freq_subset,
                "Freq_overall": freq_overall,
                f"BR_{prefix}": br,
                f"pvalue_{prefix}_raw": pval,
            }
        )
    reject, pvals_corr, _, _ = multipletests(pvals, alpha=0.01, method="fdr_bh")
    out = pd.DataFrame(rows)
    out[f"pvalue_{prefix}_corrected"] = pvals_corr
    out[f"significant_{prefix}"] = reject
    out["Ref_degeneracy"] = out["Ref_prot"].map(AA_DEGENERACY)
    out["Mut_degeneracy"] = out["Mut_prot"].map(AA_DEGENERACY)
    out["Degeneracy_pair"] = out["Ref_degeneracy"].astype(str) + "->" + out["Mut_degeneracy"].astype(str)
    return out


def degeneracy_logit(df: pd.DataFrame, flag_col: str) -> pd.Series:
    model_df = df.dropna(subset=[flag_col, "log_ref_over_mut_degeneracy", "Gene_gene"]).copy()
    x = sm.add_constant(model_df[["log_ref_over_mut_degeneracy"]])
    y = model_df[flag_col].astype(float)
    model = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster",
        cov_kwds={"groups": model_df["Gene_gene"]},
        maxiter=200,
    )
    coef = float(model.params["log_ref_over_mut_degeneracy"])
    ci_low, ci_high = model.conf_int().loc["log_ref_over_mut_degeneracy"].tolist()
    return pd.Series(
        {
            "flag": flag_col,
            "n_variants": len(model_df),
            "n_flagged": int(y.sum()),
            "coef_log_ref_over_mut_degeneracy": coef,
            "or_log_ref_over_mut_degeneracy": np.exp(coef),
            "or_95ci_low": np.exp(ci_low),
            "or_95ci_high": np.exp(ci_high),
            "p_value": float(model.pvalues["log_ref_over_mut_degeneracy"]),
        }
    )


def analyze(scored: pd.DataFrame, out_dir: Path, fig_prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = add_labels_and_scores(scored)
    protein_model = str(df["protein_model_for_fig5"].iloc[0])
    aurocs = {
        protein_model: roc_auc_score(df["label"], -df["protein_llr_for_fig5"]),
        "CaLM codon": roc_auc_score(df["label"], -df["calm_codon_llr_for_fig5"]),
        "CaLM AA-aggregated": roc_auc_score(df["label"], -df["calm_aa_agg_llr"]),
    }
    validation = {
        "n_variants": len(df),
        "n_genes": df["Gene_gene"].nunique(),
        "codon_llr_recompute_max_abs_diff": float(
            np.max(np.abs(df["calm_codon_llr_recomputed"] - df["LLR_gene"]))
        ),
        "codon_vs_aa_agg_diff_pearson": float(df["diff_codon"].corr(df["diff_aa_agg"])),
    }

    summaries = []
    pair_tables = []
    for space, diff_col in [("codon", "diff_codon"), ("aa_aggregated", "diff_aa_agg")]:
        clm, plm, lower, upper = extreme_sets(df, diff_col)
        df[f"clm_better_{space}"] = df.index.isin(clm.index).astype(int)
        df[f"plm_better_{space}"] = df.index.isin(plm.index).astype(int)
        clm_pair = pair_enrichment(df, clm, f"clm_{space}")
        plm_pair = pair_enrichment(df, plm, f"plm_{space}")
        pair_table = clm_pair.merge(
            plm_pair.drop(columns=["Ref_prot", "Mut_prot", "Count_overall", "Freq_overall", "Ref_degeneracy", "Mut_degeneracy", "Degeneracy_pair"]),
            on="Pair",
            how="outer",
        )
        pair_table["space"] = space
        pair_tables.append(pair_table)
        summaries.append(
            {
                "space": space,
                "diff_col": diff_col,
                "lower_1pct": lower,
                "upper_99pct": upper,
                "clm_better_n": len(clm),
                "plm_better_n": len(plm),
                "clm_enriched_pairs_fdr01": int(
                    ((clm_pair[f"significant_clm_{space}"]) & (clm_pair[f"BR_clm_{space}"] > 1)).sum()
                ),
                "plm_enriched_pairs_fdr01": int(
                    ((plm_pair[f"significant_plm_{space}"]) & (plm_pair[f"BR_plm_{space}"] > 1)).sum()
                ),
                "protein_model": protein_model,
                "auroc_protein": aurocs[protein_model],
                "auroc_calm_codon": aurocs["CaLM codon"],
                "auroc_calm_aa_aggregated": aurocs["CaLM AA-aggregated"],
                **validation,
            }
        )

    summary = pd.DataFrame(summaries)
    pair_all = pd.concat(pair_tables, ignore_index=True)
    degeneracy = pd.DataFrame(
        [
            degeneracy_logit(df, "clm_better_codon"),
            degeneracy_logit(df, "clm_better_aa_aggregated"),
        ]
    )
    comparison = summary.merge(
        degeneracy.assign(
            space=lambda x: x["flag"].map(
                {
                    "clm_better_codon": "codon",
                    "clm_better_aa_aggregated": "aa_aggregated",
                }
            )
        )[
            [
                "space",
                "n_flagged",
                "coef_log_ref_over_mut_degeneracy",
                "or_log_ref_over_mut_degeneracy",
                "or_95ci_low",
                "or_95ci_high",
                "p_value",
            ]
        ],
        on="space",
        how="left",
    ).rename(
        columns={
            "n_flagged": "degeneracy_n_flagged",
            "coef_log_ref_over_mut_degeneracy": "degeneracy_coef",
            "or_log_ref_over_mut_degeneracy": "degeneracy_or",
            "or_95ci_low": "degeneracy_or_95ci_low",
            "or_95ci_high": "degeneracy_or_95ci_high",
            "p_value": "degeneracy_p",
        }
    )

    scored_with_flags = df.drop(columns=["Pair"])
    scored_with_flags.to_csv(out_dir / "calm_aa_aggregation_variant_scores_with_flags.csv", index=False)
    scored_with_flags.to_csv(out_dir / "len1022_calm_aa_aggregation_variant_scores_with_flags.csv", index=False)
    summary.to_csv(out_dir / "calm_aa_aggregation_control_summary.csv", index=False)
    summary.to_csv(out_dir / "len1022_calm_aa_aggregation_control_summary.csv", index=False)
    pair_all.to_csv(out_dir / "calm_aa_aggregation_pair_enrichment.csv", index=False)
    pair_all.to_csv(out_dir / "len1022_calm_aa_aggregation_pair_enrichment.csv", index=False)
    degeneracy.to_csv(out_dir / "calm_aa_aggregation_degeneracy_logit.csv", index=False)
    degeneracy.to_csv(out_dir / "len1022_calm_aa_aggregation_degeneracy_logit.csv", index=False)
    comparison.to_csv(out_dir / "len1022_calm_aa_aggregation_discordance_model_comparison.csv", index=False)
    plot_control(summary, degeneracy, fig_prefix)
    return summary, degeneracy


def plot_control(summary: pd.DataFrame, degeneracy: pd.DataFrame, fig_prefix: str) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"width_ratios": [1, 1.15]})

    auroc_values = [
        summary["auroc_protein"].iloc[0],
        summary["auroc_calm_codon"].iloc[0],
        summary["auroc_calm_aa_aggregated"].iloc[0],
    ]
    protein_label = str(summary["protein_model"].iloc[0]).replace(" ", "\n", 1)
    labels = [protein_label, "CaLM\ncodon", "CaLM\nAA-agg"]
    colors = ["#9E9E9E", "#4D6F8F", "#B33F49"]
    axes[0].bar(labels, auroc_values, color=colors, width=0.65)
    axes[0].set_ylim(0.74, max(auroc_values) + 0.04)
    axes[0].set_ylabel("AUROC")
    axes[0].set_title("Pathogenicity prediction")
    for idx, val in enumerate(auroc_values):
        axes[0].text(idx, val + 0.006, f"{val:.3f}", ha="center", va="bottom", fontsize=8)

    deg = degeneracy.copy()
    deg["display"] = deg["flag"].map(
        {
            "clm_better_codon": "Codon-level\nCaLM better",
            "clm_better_aa_aggregated": "AA-aggregated\nCaLM better",
        }
    )
    y = np.arange(len(deg))
    axes[1].errorbar(
        deg["or_log_ref_over_mut_degeneracy"],
        y,
        xerr=[
            deg["or_log_ref_over_mut_degeneracy"] - deg["or_95ci_low"],
            deg["or_95ci_high"] - deg["or_log_ref_over_mut_degeneracy"],
        ],
        fmt="o",
        color="#B33F49",
        ecolor="#B33F49",
        capsize=3,
        markersize=5,
    )
    axes[1].axvline(1, color="#555555", linestyle=(0, (3, 2)), linewidth=0.8)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(deg["display"])
    axes[1].set_xlabel("OR for log(ref degeneracy / mut degeneracy)")
    axes[1].set_title("Degeneracy dependence")
    for idx, row in deg.iterrows():
        axes[1].text(
            row["or_95ci_high"] * 1.03,
            idx,
            f"OR={row['or_log_ref_over_mut_degeneracy']:.2f}\np={row['p_value']:.1e}",
            va="center",
            fontsize=7.5,
            color="#1F2933",
        )

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#E8E8E8", linewidth=0.7)

    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{fig_prefix}_probability_space_control.png", dpi=600)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    parser.add_argument("--gene-dir", type=Path, default=DEFAULT_GENE_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--fig-prefix", default="calm_aa_aggregation")
    parser.add_argument("--sort-by-length", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--report-every", type=int, default=25)
    parser.add_argument(
        "--max-remaining-genes",
        type=int,
        default=None,
        help="Score at most this many not-yet-scored genes, then exit normally.",
    )
    parser.add_argument(
        "--score-only",
        action="store_true",
        help="Only append variant scores; skip summaries and figures.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.output is None:
        args.output = args.out_dir / "len1022_calm_aa_aggregation_variant_scores.csv"
    scored = compute_scores(args, args.output)
    if args.score_only:
        print(
            f"Score-only mode complete: {len(scored)} variants across "
            f"{scored['Gene_gene'].nunique()} genes"
        )
        print(f"Wrote {args.output}")
        return
    summary, degeneracy = analyze(scored, args.out_dir, args.fig_prefix)
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(degeneracy.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(f"Wrote {args.output}")
    print(f"Wrote {args.out_dir / 'calm_aa_aggregation_control_summary.csv'}")
    print(f"Wrote {args.out_dir / 'calm_aa_aggregation_pair_enrichment.csv'}")
    print(f"Wrote {args.out_dir / 'calm_aa_aggregation_degeneracy_logit.csv'}")
    print(f"Wrote {FIG_DIR / (args.fig_prefix + '_probability_space_control.png')}")


if __name__ == "__main__":
    main()
