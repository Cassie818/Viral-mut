#!/usr/bin/env python3
"""Fig. 5: probability-space control for CaLM/PLM discordance."""

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


BASE = Path("Results/Revision/len1022_aa_aggregation")
FIG_DIR = Path("Figure")
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "mpl_cache"))

import matplotlib.pyplot as plt


AA_ORDER = list("WFYCMHNDQEIKPTVAGLRS")
EDGE = "#A5A5A0"
DARK = "#3A3A3A"
BLUE = "#BFE3D8"
ROSE = "#F0B9AD"
POINT_BLUE = "#A9DCCD"
POINT_SLATE = "#8298A1"
SMOKE = "#EEF0EF"
GRID = "#F8F8F7"
MORANDI_DIVERGING = LinearSegmentedColormap.from_list(
    "clean_teal_coral",
    ["#6FA8B8", "#CFE8E7", "#FFFFFF", "#F6D8CE", "#DD8F83"],
)
PATHOGENIC_LABELS = {"pathogenic", "likely_pathogenic"}
BENIGN_LABELS = {"benign", "likely_benign"}


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        clip_on=False,
    )


def format_axes(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(0.9)
    ax.tick_params(axis="both", width=0.85, length=3.0, color=EDGE, labelcolor=DARK)


def plot_discordance_scatter(ax: plt.Axes, variants: pd.DataFrame) -> None:
    calm_codon_col = "calm_codon_llr_for_fig5" if "calm_codon_llr_for_fig5" in variants.columns else "LLR_gene"
    full = variants.dropna(subset=["esm2_650m_llr", calm_codon_col, "calm_aa_agg_llr"]).copy()
    full["diff_codon_650m"] = full["esm2_650m_llr"] - full[calm_codon_col]
    full["diff_aa_agg_650m"] = full["esm2_650m_llr"] - full["calm_aa_agg_llr"]
    r = full["diff_codon_650m"].corr(full["diff_aa_agg_650m"])
    lim = np.nanpercentile(
        np.abs(pd.concat([full["diff_codon_650m"], full["diff_aa_agg_650m"]])),
        99.5,
    )
    lim = float(np.ceil(lim))
    plot_df = full.sample(frac=1.0, random_state=7)
    ax.scatter(
        plot_df["diff_codon_650m"],
        plot_df["diff_aa_agg_650m"],
        s=5.2,
        color="#9FCFB7",
        alpha=0.13,
        edgecolor="none",
        linewidth=0,
        rasterized=True,
    )
    ax.plot([-lim, lim], [-lim, lim], color="#5E5E5A", linewidth=0.65, linestyle=(0, (3, 2)))
    ax.axhline(0, color="#DADAD6", linewidth=0.55)
    ax.axvline(0, color="#DADAD6", linewidth=0.55)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Codon-space discordance")
    ax.set_ylabel("AA-space discordance")
    ax.text(
        0.05,
        0.92,
        f"Pearson r = {r:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.1,
        color=DARK,
    )
    add_panel_label(ax, "A")
    format_axes(ax)


def p_to_stars(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def plot_degeneracy_forest(ax: plt.Axes, comparison: pd.DataFrame) -> None:
    df = comparison[
        (comparison["protein_model"] == "ESM-2 650M")
        & (comparison["space"].isin(["codon", "aa_aggregated"]))
    ].copy()
    df["label"] = df["space"].map({"codon": "CaLM codon", "aa_aggregated": "AA-aggregated"})
    df = df.set_index("space").loc[["codon", "aa_aggregated"]].reset_index()
    x = np.arange(len(df))
    colors = [POINT_BLUE, "#E9AFA7"]
    yerr = [
        df["degeneracy_or"] - df["degeneracy_or_95ci_low"],
        df["degeneracy_or_95ci_high"] - df["degeneracy_or"],
    ]
    ax.bar(
        x,
        df["degeneracy_or"],
        width=0.58,
        color=colors,
        edgecolor=EDGE,
        linewidth=0.85,
        zorder=2,
    )
    ax.errorbar(
        x,
        df["degeneracy_or"],
        yerr=yerr,
        fmt="none",
        ecolor="#5A5A56",
        elinewidth=0.95,
        capsize=4,
        capthick=0.95,
        zorder=3,
    )
    for xi, row in zip(x, df.to_dict("records")):
        stars = p_to_stars(float(row["degeneracy_p"]))
        if stars:
            ax.text(
                xi,
                float(row["degeneracy_or_95ci_high"]) + 0.035,
                stars,
                ha="center",
                va="bottom",
                fontsize=8.4,
                color=DARK,
            )
    ax.axhline(1, color=DARK, linestyle=(0, (3, 2)), linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(["CaLM\ncodon", "AA-\naggregated"])
    ax.set_ylim(0, 1.08)
    ax.set_xlim(-0.55, len(df) - 0.45)
    ax.set_ylabel("Odds ratio")
    add_panel_label(ax, "B")
    format_axes(ax)


def extreme_sets_esm2_650m_aa(variants: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, float, float]:
    df = variants.dropna(subset=["esm2_650m_llr", "calm_aa_agg_llr", "Label_prot"]).copy()
    df["diff_aa_agg_650m"] = df["esm2_650m_llr"] - df["calm_aa_agg_llr"]
    lower = float(df["diff_aa_agg_650m"].quantile(0.01))
    upper = float(df["diff_aa_agg_650m"].quantile(0.99))
    calm = pd.concat(
        [
            df[(df["diff_aa_agg_650m"] > upper) & (df["Label_prot"].isin(PATHOGENIC_LABELS))],
            df[(df["diff_aa_agg_650m"] < lower) & (df["Label_prot"].isin(BENIGN_LABELS))],
        ],
        ignore_index=False,
    )
    plm = pd.concat(
        [
            df[(df["diff_aa_agg_650m"] < lower) & (df["Label_prot"].isin(PATHOGENIC_LABELS))],
            df[(df["diff_aa_agg_650m"] > upper) & (df["Label_prot"].isin(BENIGN_LABELS))],
        ],
        ignore_index=False,
    )
    return calm, plm, lower, upper


def pair_enrichment(variants: pd.DataFrame, subset: pd.DataFrame, prefix: str) -> pd.DataFrame:
    df = variants.dropna(subset=["Ref_prot", "Mut_prot"]).copy()
    df["Pair"] = list(zip(df["Ref_prot"], df["Mut_prot"]))
    sub = subset.dropna(subset=["Ref_prot", "Mut_prot"]).copy()
    sub["Pair"] = list(zip(sub["Ref_prot"], sub["Mut_prot"]))
    overall = Counter(df["Pair"])
    subset_counts = Counter(sub["Pair"])
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
                f"Count_{prefix}_aa_aggregated": n_subset,
                "Count_overall": n_total,
                f"Freq_{prefix}_aa_aggregated": freq_subset,
                "Freq_overall": freq_overall,
                f"BR_{prefix}_aa_aggregated": br,
                f"pvalue_{prefix}_aa_aggregated_raw": pval,
            }
        )
    reject, pvals_corr, _, _ = multipletests(pvals, alpha=0.01, method="fdr_bh")
    out = pd.DataFrame(rows)
    out[f"pvalue_{prefix}_aa_aggregated_corrected"] = pvals_corr
    out[f"significant_{prefix}_aa_aggregated"] = reject
    return out


def build_esm2_650m_pair_enrichment(variants: pd.DataFrame) -> pd.DataFrame:
    complete = variants.dropna(subset=["esm2_650m_llr", "calm_aa_agg_llr", "Ref_prot", "Mut_prot"]).copy()
    calm, plm, lower, upper = extreme_sets_esm2_650m_aa(complete)
    calm_pair = pair_enrichment(complete, calm, "clm")
    plm_pair = pair_enrichment(complete, plm, "plm")
    out = calm_pair.merge(
        plm_pair.drop(columns=["Ref_prot", "Mut_prot", "Count_overall", "Freq_overall"]),
        on="Pair",
        how="outer",
    )
    out["protein_model"] = "ESM-2 650M"
    out["space"] = "aa_aggregated"
    out["lower_1pct"] = lower
    out["upper_99pct"] = upper
    out["clm_leaning_n"] = len(calm)
    out["plm_leaning_n"] = len(plm)
    out.to_csv(BASE / "fig5_pair_enrichment_esm2_650m_aa_aggregated.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "protein_model": "ESM-2 650M",
                "space": "aa_aggregated",
                "n_variants": len(complete),
                "lower_1pct": lower,
                "upper_99pct": upper,
                "clm_leaning_n": len(calm),
                "plm_leaning_n": len(plm),
                "clm_enriched_pairs_fdr01": int(
                    (
                        out["significant_clm_aa_aggregated"].fillna(False)
                        & (out["BR_clm_aa_aggregated"] > 1)
                    ).sum()
                ),
                "plm_enriched_pairs_fdr01": int(
                    (
                        out["significant_plm_aa_aggregated"].fillna(False)
                        & (out["BR_plm_aa_aggregated"] > 1)
                    ).sum()
                ),
            }
        ]
    )
    summary.to_csv(BASE / "fig5_pair_enrichment_esm2_650m_aa_aggregated_summary.csv", index=False)
    return out


def heatmap_matrix(pair_df: pd.DataFrame, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    matrix = pd.DataFrame(np.nan, index=AA_ORDER, columns=AA_ORDER)
    sig = pd.DataFrame(False, index=AA_ORDER, columns=AA_ORDER)
    br_col = f"BR_{prefix}_aa_aggregated"
    q_col = f"pvalue_{prefix}_aa_aggregated_corrected"
    sig_col = f"significant_{prefix}_aa_aggregated"
    for _, row in pair_df.iterrows():
        ref = row["Ref_prot"]
        mut = row["Mut_prot"]
        if ref in matrix.index and mut in matrix.columns and pd.notna(row[br_col]):
            matrix.loc[ref, mut] = np.log2(max(float(row[br_col]), 1e-6))
            sig.loc[ref, mut] = bool(row[sig_col]) and float(row[br_col]) > 1 and float(row[q_col]) < 0.01
    return matrix.to_numpy(float), sig.to_numpy(bool)


def plot_pair_heatmap(
    ax: plt.Axes,
    pair_df: pd.DataFrame,
    prefix: str,
    label: str,
    add_cbar: bool = False,
    cbar_ax: plt.Axes | None = None,
) -> None:
    matrix, sig = heatmap_matrix(pair_df, prefix)
    cmap = MORANDI_DIVERGING.copy()
    cmap.set_bad("#FFFFFF")
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5)
    display = matrix.copy()
    display[np.abs(display) < 0.30] = 0.0
    im = ax.imshow(display, cmap=cmap, norm=norm, aspect="equal", interpolation="nearest")
    ax.set_xticks(np.arange(len(AA_ORDER)))
    ax.set_xticklabels(AA_ORDER, fontsize=7.0)
    ax.set_yticks(np.arange(len(AA_ORDER)))
    ax.set_yticklabels(AA_ORDER, fontsize=7.0)
    ax.set_xlabel("Mutant amino acid")
    ax.set_ylabel("Reference amino acid")
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_facecolor("#FFFFFF")
    for i, j in np.argwhere(sig):
        rect = plt.Rectangle(
            (j - 0.5, i - 0.5),
            1,
            1,
            fill=False,
            edgecolor="#2B2B2B",
            linewidth=0.95,
            joinstyle="miter",
        )
        ax.add_patch(rect)
    ax.set_title(label, fontsize=8.0, color=DARK, pad=7.0)
    if add_cbar:
        cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, fraction=0.038, pad=0.075)
        cbar.set_label(r"log$_2$ enrichment", fontsize=7.4)
        cbar.ax.tick_params(labelsize=6.8, width=0.8, length=2.5)
    format_axes(ax)
    return im


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(BASE / "len1022_calm_aa_aggregation_variant_scores_with_flags.csv")
    comparison = pd.read_csv(BASE / "len1022_calm_aa_aggregation_discordance_model_comparison.csv")
    aa_pair = build_esm2_650m_pair_enrichment(variants)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(7.15, 6.55))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.0],
        height_ratios=[1.0, 1.0],
        hspace=0.34,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_discordance_scatter(ax_a, variants)
    ax_a.set_aspect("equal", adjustable="box")
    plot_degeneracy_forest(ax_b, comparison)
    add_panel_label(ax_c, "C")
    plot_pair_heatmap(ax_c, aa_pair, "clm", "CaLM-leaning discordance", add_cbar=False)
    add_panel_label(ax_d, "D")
    im_d = plot_pair_heatmap(ax_d, aa_pair, "plm", "PLM-leaning discordance", add_cbar=False)

    fig.subplots_adjust(left=0.08, right=0.93, top=0.965, bottom=0.08)
    d_pos = ax_d.get_position()
    cax_d = fig.add_axes([d_pos.x1 + 0.008, d_pos.y0, 0.012, d_pos.height])
    cbar = fig.colorbar(im_d, cax=cax_d)
    cbar.set_label(r"log$_2$ enrichment", fontsize=7.4)
    cbar.ax.tick_params(labelsize=6.8, width=0.8, length=2.5)
    for filename in [
        "fig5.png",
        "fig5_codon_resolution_discordance.png",
    ]:
        path = FIG_DIR / filename
        fig.savefig(path, dpi=600, bbox_inches="tight")
        print(path)


if __name__ == "__main__":
    main()
