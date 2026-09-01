#!/usr/bin/env python3
"""Gene-level PLM-background and generic-ensembling controls."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


DATA = Path("Results/ClinVar/gene_level/gene_level_codon_contribution_summary.csv")
STATS = Path("Results/ClinVar/gene_level/gene_level_baseline_adjusted_regression.csv")
OUT = Path("Figure/figS9.png")

TEXT = "#30302E"
EDGE = "#8D8D88"
POINT = "#B8DCE8"
OLS = "#5D5D59"


def setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def p_text(p: float) -> str:
    exponent = int(np.floor(np.log10(p)))
    mantissa = p / 10**exponent
    return rf"$p={mantissa:.1f}\times10^{{{exponent}}}$"


def panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    background: str,
    residual_axes: bool = False,
    stat_slope: float | None = None,
    stat_p: float | None = None,
) -> None:
    data = df.dropna(subset=[x_col, y_col]).copy()
    ax.scatter(
        data[x_col],
        data[y_col],
        s=15,
        color=POINT,
        alpha=0.48,
        edgecolor="none",
        rasterized=True,
    )

    model = sm.OLS(data[y_col], sm.add_constant(data[x_col], has_constant="add")).fit(cov_type="HC3")
    xs = np.linspace(float(data[x_col].min()), float(data[x_col].max()), 200)
    ys = model.params.iloc[0] + model.params.iloc[1] * xs
    ax.plot(xs, ys, color=OLS, linewidth=0.9, label="OLS fit")

    reported_slope = float(model.params.iloc[1]) if stat_slope is None else stat_slope
    reported_p = float(model.pvalues.iloc[1]) if stat_p is None else stat_p
    stats = rf"$\beta={reported_slope:.3f}$" + "\n" + p_text(reported_p)

    ax.text(
        0.04,
        0.96,
        stats,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        color=TEXT,
        bbox={"facecolor": (1, 1, 1, 0.82), "edgecolor": "none", "pad": 2.0},
    )
    ax.text(0.98, 0.04, background, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.2, color=TEXT)
    ax.axhline(0, color="#C9C9C4", linewidth=0.7, linestyle=(0, (3, 2)))
    if residual_axes:
        ax.axvline(0, color="#C9C9C4", linewidth=0.7, linestyle=(0, (3, 2)))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(0.9)
    ax.tick_params(color=EDGE, labelcolor=TEXT, width=0.8, length=3)


def main() -> None:
    setup()
    df = pd.read_csv(DATA)
    stats = pd.read_csv(STATS).set_index("analysis")
    esm1b_adjusted = stats.loc["esm1b_calm_auroc_predicts_gain_adjusting_for_plm_auroc"]
    plm_control_adjusted = stats.loc[
        "esm2_plm_control_calm_auroc_predicts_gain_adjusting_for_plm_auroc"
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharey=True)

    panel(
        axes[0, 0],
        df,
        "esm1b_calm_weight",
        "esm1b_gain_over_protein",
        r"$\bar{w}_g$",
        r"$\Delta\mathrm{AUROC}_{+\mathrm{CaLM}}$",
        "ESM-1b (650M)",
    )
    panel(
        axes[0, 1],
        df,
        "esm1b_calm_auroc_residual",
        "esm1b_gain_residual",
        r"$r_{\mathrm{CaLM},g}$",
        r"$r_{\mathrm{gain},g}$",
        "ESM-1b (650M)",
        residual_axes=True,
        stat_slope=float(esm1b_adjusted["calm_coefficient"]),
        stat_p=float(esm1b_adjusted["calm_p_value"]),
    )
    panel(
        axes[1, 0],
        df,
        "esm2_calm_weight",
        "esm2_cross_modal_advantage",
        r"$\bar{w}_g$",
        r"$\Delta\mathrm{AUROC}_{\mathrm{adv}}$",
        "PLM+PLM control",
    )
    panel(
        axes[1, 1],
        df,
        "esm2_plm_control_calm_auroc_residual",
        "esm2_plm_control_gain_residual",
        r"$r_{\mathrm{CaLM},g}$",
        r"$r_{\mathrm{adv},g}$",
        "PLM+PLM control",
        residual_axes=True,
        stat_slope=float(plm_control_adjusted["calm_coefficient"]),
        stat_p=float(plm_control_adjusted["calm_p_value"]),
    )

    for ax, label in zip(axes.flat, "ABCD"):
        ax.text(-0.15, 1.04, label, transform=ax.transAxes, fontsize=10.5, fontweight="bold", ha="left", va="bottom")

    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=1, frameon=False, bbox_to_anchor=(0.5, 0.012), fontsize=7.4)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.975, bottom=0.12, wspace=0.25, hspace=0.34)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=600, bbox_inches="tight")
    print(OUT)


if __name__ == "__main__":
    main()
