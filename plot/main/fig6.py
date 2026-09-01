#!/usr/bin/env python3
"""Gene-level codon contribution figure for the revised ClinVar analyses."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm


BASE = Path("Results/ClinVar/gene_level")
GPROF = BASE / "gprofiler"
FIG_DIR = Path("Figure")

EDGE = "#8E8E8A"
TEXT = "#2F2F2F"
GRID = "#E9E7E2"
SMOKE = "#EEF0EF"
ROSE = "#E9B8BB"
BLUE = "#B9DDE9"
BAR_BLUE = "#C7E6F3"
GENE_POINT = "#A9BEC6"
GENE_LINE = "#5F5F5B"
LOWESS_LINE = "#9FCFB7"
SAGE = "#C7DCC1"
YELLOW = "#F3EBC4"
LILAC = "#D8CAE8"
TEAL = "#A9D8D1"


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "ps.fonttype": 42,
        }
    )


def format_ax(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", color=EDGE, labelcolor=TEXT, width=0.9, length=3)
    ax.grid(False)


def add_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="#111111",
        clip_on=False,
    )


def fit_line(df: pd.DataFrame, x_col: str, y_col: str) -> tuple[float, float, float, float, float]:
    reg = df.dropna(subset=[x_col, y_col]).copy()
    x = sm.add_constant(reg[x_col].to_numpy(float), has_constant="add")
    y = reg[y_col].to_numpy(float)
    model = sm.OLS(y, x).fit(cov_type="HC3")
    ci_low, ci_high = model.conf_int(alpha=0.05)[1]
    return (
        float(model.params[1]),
        float(ci_low),
        float(ci_high),
        float(model.pvalues[1]),
        float(model.rsquared),
    )


def p_math(p: float) -> str:
    if p == 0:
        return r"$p<10^{-300}$"
    exponent = int(np.floor(np.log10(abs(p))))
    mantissa = p / (10**exponent)
    return rf"$p={mantissa:.1f}\times10^{{{exponent}}}$"


def scatter_with_fit(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    color: str,
    stat_mode: str = "ols",
    show_lowess: bool = False,
    show_stats: bool = True,
) -> None:
    reg = df.dropna(subset=[x_col, y_col, "n"]).copy()
    size = np.clip(np.sqrt(reg["n"].to_numpy(float)) * 4.2, 10, 48)
    ax.scatter(
        reg[x_col],
        reg[y_col],
        s=size,
        color=color,
        alpha=0.44,
        edgecolor="none",
        linewidth=0,
        rasterized=True,
    )
    slope, _, _, pval, r2 = fit_line(reg, x_col, y_col)
    xs = np.linspace(float(reg[x_col].min()), float(reg[x_col].max()), 100)
    model = sm.OLS(reg[y_col], sm.add_constant(reg[x_col], has_constant="add")).fit()
    ys = model.params.iloc[0] + model.params.iloc[1] * xs
    ax.plot(xs, ys, color=GENE_LINE, linewidth=0.85, linestyle="-", label="OLS")
    if show_lowess:
        smooth = lowess(reg[y_col], reg[x_col], frac=0.38, return_sorted=True)
        smooth = smooth[smooth[:, 0] <= 0.15]
        ax.plot(
            smooth[:, 0],
            smooth[:, 1],
            color=LOWESS_LINE,
            linewidth=1.05,
            linestyle=(0, (4, 2)),
            label="LOWESS",
        )
    ax.axhline(0, color="#BDBDB8", linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if stat_mode == "spearman":
        rho, rho_p = spearmanr(reg[x_col], reg[y_col])
        note = rf"$\rho={rho:.2f}$" + "\n" + p_math(float(rho_p)) + "\n" + rf"$\beta_{{OLS}}={slope:.3f}$"
    else:
        note = rf"$\beta={slope:.3f}$" + "\n" + p_math(float(pval)) + "\n" + rf"$R^2={r2:.2f}$"
    if show_stats:
        ax.text(
            0.03,
            0.96,
            note,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            color=TEXT,
            bbox=dict(boxstyle="round,pad=0.18", facecolor=(1, 1, 1, 0.78), edgecolor="none"),
        )
    format_ax(ax, "both")


def top_bar(ax: plt.Axes, df: pd.DataFrame, metric: str, label: str, color: str, n: int = 12) -> None:
    top = df.sort_values(metric, ascending=False).head(n).iloc[::-1]
    y = np.arange(len(top))
    ax.barh(y, top[metric], color=color, edgecolor=EDGE, linewidth=0.8, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"])
    ax.set_xlabel(label)
    xmax = max(float(top[metric].max()), 0.01)
    ax.set_xlim(0, xmax * 1.14)
    for yi, val in zip(y, top[metric]):
        ax.text(val + xmax * 0.016, yi, f"{val:.2f}" if metric.endswith("weight") else f"{val:.3f}", va="center", fontsize=6.8, color="#6B6B67")
    format_ax(ax, "x")


def top_gene_table(
    ax: plt.Axes,
    df: pd.DataFrame,
    metric: str,
    label: str,
    n: int = 8,
    metric_label: str = "Weight",
    fmt: str = "{:.2f}",
) -> None:
    top = df.sort_values(metric, ascending=False).head(n).copy()
    ax.axis("off")
    ax.text(0.00, 1.00, label, ha="left", va="top", fontsize=7.5, color=TEXT, transform=ax.transAxes)
    headers = ["Gene", metric_label, "n"]
    xs = [0.00, 0.58, 0.88]
    y = 0.86
    for x, h in zip(xs, headers):
        ax.text(x, y, h, ha="left" if h == "Gene" else "right", va="center", fontsize=6.8, fontweight="bold", color=TEXT, transform=ax.transAxes)
    ax.plot([0.00, 0.95], [0.80, 0.80], color=EDGE, linewidth=0.8, transform=ax.transAxes, clip_on=False)
    for i, (_, row) in enumerate(top.iterrows()):
        yy = 0.72 - i * 0.085
        if i % 2 == 0:
            ax.add_patch(
                plt.Rectangle(
                    (-0.015, yy - 0.035),
                    0.95,
                    0.070,
                    transform=ax.transAxes,
                    facecolor="#F6F6F3",
                    edgecolor="none",
                    zorder=-1,
                )
            )
        ax.text(xs[0], yy, str(row["gene"]), ha="left", va="center", fontsize=6.8, color=TEXT, transform=ax.transAxes)
        ax.text(xs[1], yy, fmt.format(row[metric]), ha="right", va="center", fontsize=6.8, color=TEXT, transform=ax.transAxes)
        ax.text(xs[2], yy, f"{int(row['n'])}", ha="right", va="center", fontsize=6.8, color="#6B6B67", transform=ax.transAxes)


def combined_top_gene_table(ax: plt.Axes, df: pd.DataFrame, n: int = 8) -> None:
    by_weight = df.sort_values("esm2_calm_weight", ascending=False).head(n).copy()
    by_gain = df.sort_values("esm2_cross_modal_gain", ascending=False).head(n).copy()
    ax.axis("off")
    ax.text(0.03, 0.99, "Top genes by CaLM weight", ha="left", va="top", fontsize=8.0, color=TEXT, transform=ax.transAxes)
    ax.text(0.55, 0.99, "Top genes by cross-modal gain", ha="left", va="top", fontsize=8.0, color=TEXT, transform=ax.transAxes)

    left_x = [0.03, 0.29, 0.42]
    right_x = [0.55, 0.80, 0.93]
    y_head = 0.84
    for xs, headers in [
        (left_x, ["Gene", "Weight", "n"]),
        (right_x, ["Gene", "Gain", "n"]),
    ]:
        for x, h in zip(xs, headers):
            ax.text(
                x,
                y_head,
                h,
                ha="left" if h == "Gene" else "right",
                va="center",
                fontsize=7.2,
                fontweight="bold",
                color=TEXT,
                transform=ax.transAxes,
            )
        ax.plot([xs[0], xs[-1]], [0.78, 0.78], color=EDGE, linewidth=0.8, transform=ax.transAxes, clip_on=False)

    for i, (_, row) in enumerate(by_weight.iterrows()):
        yy = 0.70 - i * 0.075
        if i % 2 == 0:
            ax.add_patch(
                plt.Rectangle((0.025, yy - 0.032), 0.42, 0.064, transform=ax.transAxes, facecolor="#F6F6F3", edgecolor="none", zorder=-1)
            )
        ax.text(left_x[0], yy, str(row["gene"]), ha="left", va="center", fontsize=7.0, color=TEXT, transform=ax.transAxes)
        ax.text(left_x[1], yy, f"{row['esm2_calm_weight']:.2f}", ha="right", va="center", fontsize=7.0, color=TEXT, transform=ax.transAxes)
        ax.text(left_x[2], yy, f"{int(row['n'])}", ha="right", va="center", fontsize=7.0, color="#6B6B67", transform=ax.transAxes)

    for i, (_, row) in enumerate(by_gain.iterrows()):
        yy = 0.70 - i * 0.075
        if i % 2 == 0:
            ax.add_patch(
                plt.Rectangle((0.545, yy - 0.032), 0.42, 0.064, transform=ax.transAxes, facecolor="#F6F6F3", edgecolor="none", zorder=-1)
            )
        ax.text(right_x[0], yy, str(row["gene"]), ha="left", va="center", fontsize=7.0, color=TEXT, transform=ax.transAxes)
        ax.text(right_x[1], yy, f"{row['esm2_cross_modal_gain']:.3f}", ha="right", va="center", fontsize=7.0, color=TEXT, transform=ax.transAxes)
        ax.text(right_x[2], yy, f"{int(row['n'])}", ha="right", va="center", fontsize=7.0, color="#6B6B67", transform=ax.transAxes)

    ax.plot([0.50, 0.50], [0.18, 0.91], color="#DAD8D2", linewidth=0.8, transform=ax.transAxes, clip_on=False)


def top_gene_bar_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    metric: str,
    title: str,
    xlabel: str,
    fmt: str,
    n: int = 8,
    xlim: tuple[float, float] | None = None,
    bar_left: float = 0.0,
    label_inside: bool = False,
    show_n: bool = True,
    show_values: bool = True,
) -> None:
    top = df.sort_values(metric, ascending=False).head(n).iloc[::-1].copy()
    y = np.arange(len(top))
    widths = top[metric].to_numpy(float) - bar_left
    ax.barh(
        y,
        widths,
        left=bar_left,
        color=YELLOW,
        edgecolor=EDGE,
        linewidth=0.85,
        height=0.66,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"], fontsize=7.2)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=8.1, color=TEXT, pad=6)
    if xlim is None:
        xmax = float(top[metric].max())
        ax.set_xlim(0, xmax * 1.20)
        label_offset = xmax * 0.018
    else:
        ax.set_xlim(*xlim)
        label_offset = (xlim[1] - xlim[0]) * 0.018
    if show_values:
        for yi, (_, row) in zip(y, top.iterrows()):
            if label_inside:
                x_text = float(row[metric]) - label_offset
                ha = "right"
            else:
                x_text = float(row[metric]) + label_offset
                ha = "left"
            ax.text(
                x_text,
                yi,
                f"{fmt.format(row[metric])}" + (f"  n={int(row['n'])}" if show_n else ""),
                va="center",
                ha=ha,
                fontsize=6.8,
                color="#5F5F5A",
            )
    format_ax(ax, "x")


def top_gene_lollipop_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    metric: str,
    title: str,
    xlabel: str,
    fmt: str,
    n: int = 8,
    xlim: tuple[float, float] = (0.95, 1.018),
) -> None:
    top = df.sort_values(metric, ascending=False).head(n).iloc[::-1].copy()
    y = np.arange(len(top))
    values = top[metric].to_numpy(float)
    ax.scatter(values, y, s=36, color=BLUE, edgecolor=EDGE, linewidth=0.75, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(top["gene"], fontsize=7.2)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=8.1, color=TEXT, pad=6)
    ax.set_xlim(*xlim)
    label_offset = (xlim[1] - xlim[0]) * 0.015
    for yi, (_, row) in zip(y, top.iterrows()):
        ax.text(
            float(row[metric]) + label_offset,
            yi,
            f"{fmt.format(row[metric])}  n={int(row['n'])}",
            va="center",
            ha="left",
            fontsize=6.8,
            color="#5F5F5A",
        )
    format_ax(ax, "x")


def wrap_terms(terms: pd.Series, width: int = 32) -> list[str]:
    return ["\n".join(textwrap.wrap(str(term), width=width)) for term in terms]


def enrichment_panel(ax: plt.Axes) -> None:
    terms = pd.read_csv(GPROF / "all_gene_sets_gprofiler_terms.csv")
    sub = terms[terms["gene_set"] == "esm2_top_codon_weight"].sort_values("p_value").head(6).copy()
    if sub.empty:
        ax.axis("off")
        ax.text(
            0.02,
            0.78,
            "Functional enrichment",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color=TEXT,
        )
        ax.text(
            0.02,
            0.58,
            "No g:Profiler terms remained\nsignificant after using the\n466 evaluable genes as the\ncustom background.",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.2,
            linespacing=1.35,
            color=TEXT,
        )
        ax.text(
            0.02,
            0.18,
            "Exploratory enrichment is\ntherefore shown only as a\nbackground-sensitive trend.",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.7,
            linespacing=1.3,
            color="#6B6B67",
        )
        return
    sub = sub.iloc[::-1]
    y = np.arange(len(sub))
    colors = [SAGE, SAGE, SAGE, YELLOW, YELLOW, YELLOW][: len(sub)]
    ax.barh(y, sub["intersection_size"], color=colors, edgecolor=EDGE, linewidth=0.8, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(wrap_terms(sub["term_name"], 24), fontsize=6.4)
    ax.set_xlabel("Gene count")
    for yi, (_, row) in zip(y, sub.iterrows()):
        ax.text(
            row["intersection_size"] + 0.35,
            yi,
            f"p={row['p_value']:.1e}",
            va="center",
            fontsize=6.5,
            color="#6B6B67",
        )
    ax.set_xlim(0, max(sub["intersection_size"]) + 7)
    format_ax(ax, "x")


def overlap_panel(ax: plt.Axes) -> None:
    esm2_weight = set(pd.read_csv(BASE / "esm2_top_codon_weight_decile_genes.csv")["gene"])
    esm2_gain = set(pd.read_csv(BASE / "esm2_top_cross_modal_gain_decile_genes.csv")["gene"])
    esm1b_weight = set(pd.read_csv(BASE / "esm1b_top_codon_weight_decile_genes.csv")["gene"])
    esm1b_gain = set(pd.read_csv(BASE / "esm1b_top_gain_decile_genes.csv")["gene"])
    labels = [
        "ESM-2 weight vs ESM-1b weight",
        "ESM-2 gain vs ESM-1b gain",
        "ESM-2 weight vs ESM-2 gain",
    ]
    overlaps = [
        len(esm2_weight & esm1b_weight),
        len(esm2_gain & esm1b_gain),
        len(esm2_weight & esm2_gain),
    ]
    jaccards = [
        len(esm2_weight & esm1b_weight) / len(esm2_weight | esm1b_weight),
        len(esm2_gain & esm1b_gain) / len(esm2_gain | esm1b_gain),
        len(esm2_weight & esm2_gain) / len(esm2_weight | esm2_gain),
    ]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, overlaps, color=[TEAL, BLUE, LILAC], edgecolor=EDGE, linewidth=0.8, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.6)
    ax.set_xlabel("Shared genes among top deciles")
    ax.set_xlim(0, 30)
    for yi, val, jac in zip(y, overlaps, jaccards):
        ax.text(val + 0.6, yi, f"{val}/47, J={jac:.2f}", ha="left", va="center", fontsize=6.8, color=TEXT)
    format_ax(ax, "x")


def main() -> None:
    setup()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(BASE / "gene_level_codon_contribution_summary.csv")

    fig = plt.figure(figsize=(5.85, 5.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 0.95], wspace=0.22, hspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
    ax_c1 = fig.add_subplot(gs[1, 0])
    ax_c2 = fig.add_subplot(gs[1, 1])

    scatter_with_fit(
        ax_a,
        df,
        "esm2_calm_weight",
        "esm2_cross_modal_gain",
        r"$\bar{w}_g$",
        r"$\Delta\mathrm{AUROC}_{+\mathrm{CaLM}}$",
        ROSE,
        show_stats=False,
    )
    scatter_with_fit(
        ax_b,
        df,
        "esm2_calm_auroc_residual",
        "esm2_gain_residual",
        r"$r_{\mathrm{CaLM},g}$",
        r"$r_{\mathrm{gain},g}$",
        ROSE,
        show_stats=False,
    )
    ax_b.axvline(0, color="#BDBDB8", linewidth=0.8, linestyle=(0, (3, 2)))
    ax_a.set_ylim(-0.15, 0.30)
    ax_b.set_ylim(-0.15, 0.30)
    top_gene_bar_panel(
        ax_c1,
        df,
        "esm2_calm_weight",
        "Top genes by CaLM weight",
        r"$\bar{w}_g$",
        "{:.2f}",
        n=10,
        xlim=(0, 1.02),
        bar_left=0.0,
        label_inside=True,
        show_n=False,
        show_values=True,
    )
    top_gene_bar_panel(
        ax_c2,
        df,
        "esm2_cross_modal_gain",
        "Top genes by cross-modal gain",
        r"$\Delta\mathrm{AUROC}_{+\mathrm{CaLM}}$",
        "{:.3f}",
        n=10,
        bar_left=0.0,
        show_n=False,
        show_values=True,
    )

    for ax, label in zip([ax_a, ax_b, ax_c1, ax_c2], list("ABCD")):
        add_label(ax, label)

    fig.subplots_adjust(left=0.095, right=0.985, top=0.965, bottom=0.075)

    output_path = FIG_DIR / "fig6.png"
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    print(output_path)


if __name__ == "__main__":
    main()
