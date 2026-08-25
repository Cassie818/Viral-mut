#!/usr/bin/env python3
"""ClinVar gene-level source tables and sequence-confounder controls.

This script uses the ClinVar cohort to create gene-level summaries for
downstream codon-contribution analyses and to test whether CaLM signal
persists after controlling explicit sequence features.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score


INPUT = Path("Results/Revision/len1022_model_control/len1022_model_control_score_table_complete_cases.csv")
GENE_FASTA_DIR = Path("/Users/cassie/Desktop/Gene")
OUT_DIR = Path("Results/Revision/len1022_core_clinvar")

BASES = ("A", "C", "G", "T")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def best_weight(y: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    best_auc = -np.inf
    best_w = np.nan
    for w in np.linspace(0, 1, 101):
        score = (1 - w) * a + w * b
        auc = roc_auc_score(y, score)
        if auc > best_auc:
            best_auc = auc
            best_w = w
    return float(best_w), float(best_auc)


def per_gene_summary(df: pd.DataFrame, min_pos: int = 5, min_neg: int = 5) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for gene, gene_df in df.groupby("Gene_gene", sort=True):
        y = gene_df["label"].to_numpy()
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        if n_pos < min_pos or n_neg < min_neg:
            continue

        s150 = gene_df["esm2_150m_score"].to_numpy(dtype=float)
        s650 = gene_df["esm2_650m_score"].to_numpy(dtype=float)
        s1b = gene_df["esm1b_650m_score"].to_numpy(dtype=float)
        calm = gene_df["calm_score"].to_numpy(dtype=float)

        w_150_650, auc_150_650 = best_weight(y, s150, s650)
        w_650_calm, auc_650_calm = best_weight(y, s650, calm)
        w_150_calm, auc_150_calm = best_weight(y, s150, calm)
        w_650_1b, auc_650_1b = best_weight(y, s650, s1b)
        w_1b_calm, auc_1b_calm = best_weight(y, s1b, calm)

        auc_150 = roc_auc_score(y, s150)
        auc_650 = roc_auc_score(y, s650)
        auc_1b = roc_auc_score(y, s1b)
        auc_calm = roc_auc_score(y, calm)

        rows.append(
            {
                "gene": gene,
                "n": int(len(gene_df)),
                "n_pathogenic": n_pos,
                "n_benign": n_neg,
                "protein_length": int(gene_df["protein_length"].iloc[0]),
                "ESM-2 150M AUROC": auc_150,
                "ESM-2 650M AUROC": auc_650,
                "ESM-1b 650M AUROC": auc_1b,
                "CaLM AUROC": auc_calm,
                "ESM-2 150M + 650M AUROC": auc_150_650,
                "ESM-2 150M + 650M weight_second": w_150_650,
                "ESM-2 650M + ESM-1b 650M AUROC": auc_650_1b,
                "ESM-2 650M + ESM-1b 650M weight_second": w_650_1b,
                "ESM-2 650M + CaLM AUROC": auc_650_calm,
                "ESM-2 650M + CaLM weight_second": w_650_calm,
                "ESM-2 150M + CaLM AUROC": auc_150_calm,
                "ESM-2 150M + CaLM weight_second": w_150_calm,
                "ESM-1b 650M + CaLM AUROC": auc_1b_calm,
                "ESM-1b 650M + CaLM weight_second": w_1b_calm,
                "650M_minus_150M": auc_650 - auc_150,
                "150M650M_minus_650M": auc_150_650 - auc_650,
                "650M1b_minus_650M": auc_650_1b - auc_650,
                "650MCaLM_minus_650M": auc_650_calm - auc_650,
                "1bCaLM_minus_1b": auc_1b_calm - auc_1b,
                "650MCaLM_minus_150M650M": auc_650_calm - auc_150_650,
                "650MCaLM_minus_650M1b": auc_650_calm - auc_650_1b,
            }
        )
    return pd.DataFrame(rows)


def fit_regression(df: pd.DataFrame, x_col: str, y_col: str, label: str, weighted: bool) -> dict[str, float | str | int]:
    reg_df = df.dropna(subset=[x_col, y_col, "n_pathogenic", "n_benign"]).copy()
    x = reg_df[x_col].to_numpy(dtype=float)
    y = reg_df[y_col].to_numpy(dtype=float)
    X = sm.add_constant(x)
    if weighted:
        weights = (reg_df["n_pathogenic"] * reg_df["n_benign"]) / (reg_df["n_pathogenic"] + reg_df["n_benign"])
        model = sm.WLS(y, X, weights=weights.to_numpy(dtype=float)).fit(cov_type="HC3")
        weight_label = "n_pathogenic*n_benign/(n_pathogenic+n_benign)"
        method = "WLS_HC3"
    else:
        model = sm.OLS(y, X).fit(cov_type="HC3")
        weight_label = ""
        method = "OLS_HC3"
    ci_low, ci_high = [float(v) for v in model.conf_int(alpha=0.05)[1]]
    return {
        "analysis": label,
        "method": method,
        "weighted": int(weighted),
        "n_genes": int(len(reg_df)),
        "x_column": x_col,
        "y_column": y_col,
        "weight_definition": weight_label,
        "intercept": float(model.params[0]),
        "slope": float(model.params[1]),
        "slope_95ci_low": ci_low,
        "slope_95ci_high": ci_high,
        "p_value": float(model.pvalues[1]),
        "r_squared": float(model.rsquared),
    }


def read_fasta(path: Path) -> str | None:
    if not path.exists():
        return None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line.upper().replace("U", "T"))
    return "".join(parts) or None


def gene_composition(gene: str) -> dict[str, float]:
    seq = read_fasta(GENE_FASTA_DIR / f"{gene}.fasta")
    if seq is None:
        return {"gene_gc": np.nan, "gene_cpg_density": np.nan}
    clean = "".join(base for base in seq if base in BASES)
    if not clean:
        return {"gene_gc": np.nan, "gene_cpg_density": np.nan}
    return {
        "gene_gc": (clean.count("G") + clean.count("C")) / len(clean),
        "gene_cpg_density": clean.count("CG") / max(len(clean) - 1, 1),
    }


def mutation_features(ref: str, mut: str) -> dict[str, float]:
    ref = str(ref).upper().replace("U", "T")
    mut = str(mut).upper().replace("U", "T")
    out = {
        "valid_snv_codon": 0.0,
        "ref_codon_gc": np.nan,
        "mut_codon_gc": np.nan,
        "delta_codon_gc": np.nan,
        "ref_cpg": np.nan,
        "mut_cpg": np.nan,
        "delta_cpg": np.nan,
        "is_transition": np.nan,
        "is_transversion": np.nan,
        "mutated_codon_position": np.nan,
        "local_degeneracy": np.nan,
    }
    if len(ref) != 3 or len(mut) != 3 or any(base not in BASES for base in ref + mut):
        return out
    diffs = [idx for idx, (a, b) in enumerate(zip(ref, mut)) if a != b]
    if len(diffs) != 1:
        return out
    pos = diffs[0]
    ref_aa = CODON_TABLE.get(ref)
    if ref_aa is None or ref_aa == "*":
        return out

    syn_count = 0
    for base in BASES:
        candidate = ref[:pos] + base + ref[pos + 1 :]
        if CODON_TABLE.get(candidate) == ref_aa:
            syn_count += 1

    ref_gc = sum(base in {"G", "C"} for base in ref)
    mut_gc = sum(base in {"G", "C"} for base in mut)
    ref_cpg = int("CG" in ref)
    mut_cpg = int("CG" in mut)
    transition = int((ref[pos], mut[pos]) in TRANSITIONS)
    out.update(
        {
            "valid_snv_codon": 1.0,
            "ref_codon_gc": ref_gc / 3,
            "mut_codon_gc": mut_gc / 3,
            "delta_codon_gc": (mut_gc - ref_gc) / 3,
            "ref_cpg": float(ref_cpg),
            "mut_cpg": float(mut_cpg),
            "delta_cpg": float(mut_cpg - ref_cpg),
            "is_transition": float(transition),
            "is_transversion": float(1 - transition),
            "mutated_codon_position": float(pos + 1),
            "local_degeneracy": float(syn_count),
        }
    )
    return out


def zscore(series: pd.Series) -> pd.Series:
    std = series.std(ddof=0)
    if std == 0 or np.isnan(std):
        return series * np.nan
    return (series - series.mean()) / std


def design_matrix(df: pd.DataFrame, columns: list[str], add_constant: bool = True) -> pd.DataFrame:
    x = df[columns].copy()
    for col in ["local_degeneracy", "mutated_codon_position"]:
        if col in x:
            dummies = pd.get_dummies(x[col].astype("Int64").astype(str), prefix=col, drop_first=True)
            x = pd.concat([x.drop(columns=[col]), dummies], axis=1)
    x = x.astype(float)
    if add_constant:
        x = sm.add_constant(x, has_constant="add")
    return x


def fit_glm_clustered(
    df: pd.DataFrame,
    predictors: list[str],
    analysis: str,
    focal_predictor: str | None,
) -> pd.Series:
    model_df = df.dropna(subset=["label", "Gene_gene", *predictors]).copy()
    y = model_df["label"].astype(float)
    x = design_matrix(model_df, predictors)
    model = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster",
        cov_kwds={"groups": model_df["Gene_gene"]},
        maxiter=200,
    )
    pred = pd.Series(model.predict(x), index=model_df.index)
    out = {
        "analysis": analysis,
        "n_variants": int(len(model_df)),
        "n_genes": int(model_df["Gene_gene"].nunique()),
        "focal_predictor": focal_predictor or "",
        "coef": np.nan,
        "coef_95ci_low": np.nan,
        "coef_95ci_high": np.nan,
        "odds_ratio": np.nan,
        "or_95ci_low": np.nan,
        "or_95ci_high": np.nan,
        "p_value": np.nan,
        "auroc": float(roc_auc_score(y, pred)),
        "pseudo_r2_mcfadden": float(1 - model.llf / model.llnull),
    }
    if focal_predictor is not None:
        ci_low, ci_high = model.conf_int().loc[focal_predictor].tolist()
        coef = float(model.params[focal_predictor])
        out.update(
            {
                "coef": coef,
                "coef_95ci_low": float(ci_low),
                "coef_95ci_high": float(ci_high),
                "odds_ratio": float(np.exp(coef)),
                "or_95ci_low": float(np.exp(ci_low)),
                "or_95ci_high": float(np.exp(ci_high)),
                "p_value": float(model.pvalues[focal_predictor]),
            }
        )
    return pd.Series(out)


def sequence_confounder_control(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_df = pd.DataFrame([mutation_features(ref, mut) for ref, mut in zip(df["Ref_gene"], df["Mut_gene"])])
    valid = pd.concat([df.reset_index(drop=True), feature_df], axis=1)
    gene_comp = pd.DataFrame.from_dict(
        {gene: gene_composition(gene) for gene in sorted(valid["Gene_gene"].dropna().unique())},
        orient="index",
    )
    gene_comp.index.name = "Gene_gene"
    valid = valid.merge(gene_comp.reset_index(), on="Gene_gene", how="left")
    valid = valid[valid["valid_snv_codon"] == 1].copy()
    valid = valid.dropna(subset=["calm_score", "esm2_650m_score", "gene_gc", "gene_cpg_density"])

    valid["z_calm_score"] = zscore(valid["calm_score"])
    valid["z_esm2_650m_score"] = zscore(valid["esm2_650m_score"])
    valid["z_gene_gc"] = zscore(valid["gene_gc"])
    valid["z_gene_cpg_density"] = zscore(valid["gene_cpg_density"])

    confounders = [
        "z_esm2_650m_score",
        "ref_codon_gc",
        "delta_codon_gc",
        "ref_cpg",
        "delta_cpg",
        "is_transition",
        "local_degeneracy",
        "mutated_codon_position",
        "z_gene_gc",
        "z_gene_cpg_density",
    ]
    residual_x = design_matrix(valid, confounders)
    residual_model = sm.OLS(valid["z_calm_score"], residual_x).fit()
    valid["calm_incremental_residual"] = residual_model.resid
    valid["z_calm_incremental_residual"] = zscore(valid["calm_incremental_residual"])

    models = [
        fit_glm_clustered(valid, confounders, "base_esm2_plus_sequence_controls", None),
        fit_glm_clustered(valid, confounders + ["z_calm_score"], "base_plus_calm", "z_calm_score"),
        fit_glm_clustered(
            valid,
            ["z_calm_incremental_residual"],
            "calm_incremental_residual_only",
            "z_calm_incremental_residual",
        ),
    ]
    summary = pd.DataFrame(models)
    base_auroc = float(summary.loc[summary["analysis"] == "base_esm2_plus_sequence_controls", "auroc"].iloc[0])
    summary["delta_auroc_vs_base"] = summary["auroc"] - base_auroc

    residual_meta = pd.DataFrame(
        [
            {
                "analysis": "calm_residualization",
                "outcome": "z_calm_score",
                "n_variants": int(len(valid)),
                "n_genes": int(valid["Gene_gene"].nunique()),
                "r_squared": float(residual_model.rsquared),
                "covariates": ";".join(confounders),
            }
        ]
    )

    keep_cols = [
        "label",
        "Gene_gene",
        "variant_id",
        "calm_score",
        "esm2_650m_score",
        "valid_snv_codon",
        "ref_codon_gc",
        "mut_codon_gc",
        "delta_codon_gc",
        "ref_cpg",
        "mut_cpg",
        "delta_cpg",
        "gene_gc",
        "gene_cpg_density",
        "is_transition",
        "is_transversion",
        "mutated_codon_position",
        "local_degeneracy",
        "z_calm_score",
        "z_esm2_650m_score",
        "z_calm_incremental_residual",
    ]
    return summary, residual_meta, valid[keep_cols]


def summarize_core_analyses(new_per_gene: pd.DataFrame, new_regs: pd.DataFrame, seq_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    rows.append(
        {
            "analysis": "per_gene_panel",
            "cohort": "len1022_complete",
            "n_genes": int(len(new_per_gene)),
            "mean_delta_delta": float(new_per_gene["650MCaLM_minus_150M650M"].mean()),
            "median_delta_delta": float(new_per_gene["650MCaLM_minus_150M650M"].median()),
            "fraction_delta_delta_positive": float((new_per_gene["650MCaLM_minus_150M650M"] > 0).mean()),
        }
    )
    for _, row in new_regs.iterrows():
        rows.append(
            {
                "analysis": row["analysis"],
                "cohort": "len1022_complete",
                "n_genes": int(row["n_genes"]),
                "slope": float(row["slope"]),
                "slope_95ci_low": float(row["slope_95ci_low"]),
                "slope_95ci_high": float(row["slope_95ci_high"]),
                "p_value": float(row["p_value"]),
                "r_squared": float(row["r_squared"]),
            }
        )
    for _, row in seq_summary.iterrows():
        if row["focal_predictor"]:
            rows.append(
                {
                    "analysis": row["analysis"],
                    "cohort": "len1022_complete",
                    "n_genes": int(row["n_genes"]),
                    "n_variants": int(row["n_variants"]),
                    "odds_ratio": float(row["odds_ratio"]),
                    "or_95ci_low": float(row["or_95ci_low"]),
                    "or_95ci_high": float(row["or_95ci_high"]),
                    "p_value": float(row["p_value"]),
                    "delta_auroc_vs_base": float(row["delta_auroc_vs_base"]),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT)
    df = df.drop_duplicates(subset=["variant_id", "Gene_gene"]).copy()

    cohort_summary = pd.DataFrame(
        [
            {
                "cohort": "len1022_complete_cases",
                "n_variants": int(len(df)),
                "n_genes": int(df["Gene_gene"].nunique()),
                "n_pathogenic": int((df["label"] == 1).sum()),
                "n_benign": int((df["label"] == 0).sum()),
                "min_protein_length": int(df["protein_length"].min()),
                "max_protein_length": int(df["protein_length"].max()),
            }
        ]
    )
    cohort_summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)

    per_gene = per_gene_summary(df)
    per_gene["codon_dependency_index"] = per_gene["ESM-2 650M + CaLM weight_second"]
    per_gene["delta_delta_auroc"] = per_gene["650MCaLM_minus_150M650M"]
    per_gene["independent_calm_minus_650m"] = per_gene["CaLM AUROC"] - per_gene["ESM-2 650M AUROC"]
    per_gene["inverse_variance_weight"] = (
        per_gene["n_pathogenic"] * per_gene["n_benign"] / (per_gene["n_pathogenic"] + per_gene["n_benign"])
    )
    per_gene.to_csv(OUT_DIR / "per_gene_ensemble_summary.csv", index=False)

    regs = pd.DataFrame(
        [
            fit_regression(
                per_gene,
                "codon_dependency_index",
                "delta_delta_auroc",
                "codon_dependency_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "codon_dependency_index",
                "delta_delta_auroc",
                "weighted_codon_dependency_predicts_cross_modal_gain",
                True,
            ),
            fit_regression(
                per_gene,
                "independent_calm_minus_650m",
                "delta_delta_auroc",
                "independent_calm_signal_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "independent_calm_minus_650m",
                "delta_delta_auroc",
                "weighted_independent_calm_signal_predicts_cross_modal_gain",
                True,
            ),
            fit_regression(
                per_gene,
                "CaLM AUROC",
                "delta_delta_auroc",
                "calm_standalone_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "CaLM AUROC",
                "delta_delta_auroc",
                "weighted_calm_standalone_predicts_cross_modal_gain",
                True,
            ),
        ]
    )
    regs.to_csv(OUT_DIR / "per_gene_delta_delta_regressions.csv", index=False)

    seq_summary, residual_meta, seq_features = sequence_confounder_control(df)
    seq_summary.to_csv(OUT_DIR / "calm_sequence_confounder_regression_summary.csv", index=False)
    residual_meta.to_csv(OUT_DIR / "calm_sequence_confounder_residualization_summary.csv", index=False)
    seq_features.to_csv(OUT_DIR / "calm_sequence_confounder_variant_features.csv", index=False)

    comparison = summarize_core_analyses(per_gene, regs, seq_summary)
    comparison.to_csv(OUT_DIR / "core_analysis_conclusion_check.csv", index=False)

    print("Cohort")
    print(cohort_summary.to_string(index=False))
    print("\nPer-gene panel")
    print(
        per_gene[
            [
                "gene",
                "n",
                "n_pathogenic",
                "n_benign",
                "ESM-2 650M AUROC",
                "CaLM AUROC",
                "ESM-2 150M + 650M AUROC",
                "ESM-2 650M + CaLM AUROC",
                "delta_delta_auroc",
                "codon_dependency_index",
                "independent_calm_minus_650m",
            ]
        ]
        .describe(include="all")
        .to_string()
    )
    print("\nRegressions")
    print(regs.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nSequence confounder control")
    print(seq_summary.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nResidualization")
    print(residual_meta.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
