#!/usr/bin/env python3
"""Fig. 5: probability-space control for CaLM/PLM discordance."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


BASE = Path("Results/ClinVar/substitution_discordance")
FIG_DIR = Path("Figure")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/viral_mut_mpl_cache")

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
    "muted_blue_coral",
    ["#5E8797", "#B9D0D4", "#FAF9F6", "#EBC2B9", "#CD8178"],
)


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
        s=0.8,
        color="#E7D46A",
        alpha=0.026,
        edgecolor="none",
        linewidth=0,
        rasterized=True,
    )
    ax.plot([-lim, lim], [-lim, lim], color="#5E5E5A", linewidth=0.65, linestyle=(0, (3, 2)))
    ax.axhline(0, color="#DADAD6", linewidth=0.55)
    ax.axvline(0, color="#DADAD6", linewidth=0.55)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Codon-space discordance", fontsize=7.2)
    ax.set_ylabel("AA-space discordance", fontsize=7.2)
    ax.tick_params(axis="both", labelsize=7.0)
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
    colors = ["#A9D5E3", "#E8B8B0"]
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
        linewidth=0.75,
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
    ax.set_xticklabels(["CaLM\ncodon", "AA-\naggregated"], fontsize=6.4)
    ax.set_ylim(0, 1.08)
    ax.set_xlim(-0.55, len(df) - 0.45)
    ax.set_ylabel("OR", fontsize=6.6)
    ax.tick_params(axis="y", labelsize=6.4)
    add_panel_label(ax, "B")
    format_axes(ax)


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
    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5)
    rows, cols = np.where(np.isfinite(matrix))
    values = matrix[rows, cols]
    im = ax.scatter(
        cols,
        rows,
        c=values,
        cmap=cmap,
        norm=norm,
        s=25,
        marker="o",
        edgecolors="none",
        linewidths=0,
    )
    ax.set_xticks(np.arange(len(AA_ORDER)))
    ax.set_xticklabels(AA_ORDER, fontsize=7.0)
    ax.set_yticks(np.arange(len(AA_ORDER)))
    ax.set_yticklabels(AA_ORDER, fontsize=7.0)
    ax.set_xlabel("Mutant amino acid", fontsize=7.2)
    ax.set_ylabel("Reference amino acid", fontsize=7.2)
    ax.set_xlim(-0.6, len(AA_ORDER) - 0.4)
    ax.set_ylim(len(AA_ORDER) - 0.4, -0.6)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_facecolor("#FFFFFF")
    for i, j in np.argwhere(sig):
        ax.scatter(
            j,
            i,
            s=38,
            marker="o",
            facecolors="none",
            edgecolors="#303030",
            linewidths=0.75,
            zorder=3,
        )
    ax.set_title(label, fontsize=8.0, color="#000000", pad=7.0)
    if add_cbar:
        cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, fraction=0.038, pad=0.075)
        cbar.set_label(r"log$_2$ enrichment", fontsize=7.4)
        cbar.ax.tick_params(labelsize=6.8, width=0.8, length=2.5)
    format_axes(ax)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#9C9C98")
        spine.set_linewidth(0.75)
    return im


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(BASE / "variant_discordance_scores.csv.gz")
    comparison = pd.read_csv(BASE / "discordance_model_comparison.csv")
    aa_pair = pd.read_csv(BASE / "pair_enrichment_esm2_650m_aa_aggregated.csv")

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(5.90, 5.30))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.04, 0.96],
        height_ratios=[1.0, 1.0],
        hspace=0.29,
        wspace=0.24,
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

    fig.subplots_adjust(left=0.085, right=0.925, top=0.96, bottom=0.085)
    d_pos = ax_d.get_position()
    cax_d = fig.add_axes([d_pos.x1 + 0.006, d_pos.y0, 0.011, d_pos.height])
    cbar = fig.colorbar(im_d, cax=cax_d)
    cbar.set_label(r"log$_2$ enrichment", fontsize=7.4)
    cbar.ax.tick_params(labelsize=6.8, width=0.8, length=2.5)
    path = FIG_DIR / "fig5.png"
    fig.savefig(path, dpi=600, bbox_inches="tight")
    print(path)


if __name__ == "__main__":
    main()
