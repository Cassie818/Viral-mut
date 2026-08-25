#!/usr/bin/env python3
"""Supplementary results for ClinMAVE functional-class analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import t


MODELS = {
    "ESM-2 (650M)": {
        "fold": "clinmave_esm2_650m_calm_cv_fold_metrics.csv",
        "summary": "clinmave_esm2_650m_calm_cv_summary.csv",
        "delta": "delta_combo_vs_esm2_650m",
        "plm_auroc": "auroc_esm2_650m",
        "ensemble_auroc": "auroc_esm2_650m_calm",
        "p": "p_delta_combo_vs_esm2_650m_paired_t",
        "color": "#9BC9DD",
    },
    "ESM-1b (650M)": {
        "fold": "clinmave_esm1b_650m_calm_cv_fold_metrics.csv",
        "summary": "clinmave_esm1b_650m_calm_cv_summary.csv",
        "delta": "delta_combo_vs_esm1b_650m",
        "plm_auroc": "auroc_esm1b_650m",
        "ensemble_auroc": "auroc_esm1b_650m_calm",
        "p": "p_delta_combo_vs_esm1b_650m_paired_t",
        "color": "#BFD8A8",
    },
}
GROUPS = [("DMS", "lof"), ("DMS", "gof"), ("CBGE", "lof"), ("CBGE", "gof")]
EDGE = "#8F8F8A"
TEXT = "#363636"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="Results/ClinMAVE/functional_effects")
    parser.add_argument("--table-dir", default="Results/ClinMAVE/supplementary_tables")
    parser.add_argument("--figure-dir", default="Figure")
    return parser.parse_args()


def load_results(input_dir: Path) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    folds, summaries = {}, {}
    for model, config in MODELS.items():
        folds[model] = pd.read_csv(input_dir / config["fold"])
        summaries[model] = pd.read_csv(input_dir / config["summary"])
    return folds, summaries


def confidence_interval(values: pd.Series) -> tuple[float, float]:
    values = values.dropna().to_numpy(float)
    if len(values) < 2:
        return np.nan, np.nan
    half_width = t.ppf(0.975, len(values) - 1) * values.std(ddof=1) / np.sqrt(len(values))
    return float(values.mean() - half_width), float(values.mean() + half_width)


def p_display(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def build_summary_table(folds: dict[str, pd.DataFrame], summaries: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for assay, case_class in GROUPS:
        for model, config in MODELS.items():
            fold = folds[model]
            subfold = fold[(fold["assay"] == assay) & (fold["case_class"] == case_class)].copy()
            summary = summaries[model]
            subsummary = summary[(summary["assay"] == assay) & (summary["case_class"] == case_class)].iloc[0]
            delta = subfold[config["delta"]]
            ci_low, ci_high = confidence_interval(delta)
            rows.append(
                {
                    "assay": assay,
                    "comparison": f"Functionally normal vs. {'GoF' if case_class == 'gof' else 'LoF'}",
                    "protein_model": model,
                    "n_variants": int(subsummary["n_variants"]),
                    "n_genes": int(subsummary["n_genes"]),
                    "n_evaluable_folds": int(delta.notna().sum()),
                    "n_total_folds": int(len(subfold)),
                    "mean_calm_weight": float(subsummary["mean_calm_weight"]),
                    "sd_calm_weight": float(subsummary["sd_calm_weight"]),
                    "plm_auroc": float(subfold[config["plm_auroc"]].mean()),
                    "ensemble_auroc": float(subfold[config["ensemble_auroc"]].mean()),
                    "delta_auroc": float(delta.mean()),
                    "ci95_delta_low": ci_low,
                    "ci95_delta_high": ci_high,
                    "paired_t_p": float(subsummary[config["p"]]),
                    "paired_t_p_display": p_display(float(subsummary[config["p"]])),
                    "ci95_crosses_zero": bool(ci_low <= 0 <= ci_high),
                }
            )
    return pd.DataFrame(rows)


def style_axis(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(0.9)
    ax.tick_params(color=EDGE, labelcolor=TEXT, width=0.8, length=3)


def plot_gof_deltas(folds: dict[str, pd.DataFrame], figure_path: Path) -> None:
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.65), sharey=True)
    rng = np.random.default_rng(16)
    for ax, assay in zip(axes, ["DMS", "CBGE"]):
        for index, (model, config) in enumerate(MODELS.items()):
            sub = folds[model]
            values = sub[(sub["assay"] == assay) & (sub["case_class"] == "gof")][config["delta"]].dropna().to_numpy(float)
            y = np.full(len(values), index, dtype=float) + rng.normal(0, 0.045, len(values))
            ax.scatter(values, y, s=28, color=config["color"], edgecolor=EDGE, linewidth=0.45, zorder=3)
            if len(values):
                low, high = confidence_interval(pd.Series(values))
                ax.errorbar(values.mean(), index, xerr=[[values.mean() - low], [high - values.mean()]], fmt="o", color=TEXT, markersize=3.8, capsize=3, linewidth=0.9, zorder=4)
        ax.axvline(0, color="#A8A8A3", linestyle=(0, (3, 2)), linewidth=0.8)
        ax.set_title(assay, fontsize=9, fontweight="normal", color=TEXT)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["ESM-2 (650M)", "ESM-1b (650M)"])
        ax.invert_yaxis()
        ax.set_xlabel("Fold-level $\\Delta$AUROC")
        style_axis(ax)
    fig.subplots_adjust(left=0.15, right=0.99, top=0.85, bottom=0.18, wspace=0.15)
    fig.savefig(figure_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_weights(folds: dict[str, pd.DataFrame], figure_path: Path) -> None:
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 9})
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    labels = ["DMS\nLoF", "DMS\nGoF", "CBGE\nLoF", "CBGE\nGoF"]
    rng = np.random.default_rng(17)
    for offset, (model, config) in zip([-0.14, 0.14], MODELS.items()):
        for index, (assay, case_class) in enumerate(GROUPS):
            values = folds[model][(folds[model]["assay"] == assay) & (folds[model]["case_class"] == case_class)]["calm_weight"].to_numpy(float)
            x = np.full(len(values), index + offset) + rng.normal(0, 0.025, len(values))
            ax.scatter(x, values, s=24, color=config["color"], edgecolor=EDGE, linewidth=0.4, alpha=0.95, zorder=3)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Optimised CaLM weight")
    ax.set_ylim(-0.04, 1.06)
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=config["color"],
            markeredgecolor=EDGE,
            markeredgewidth=0.45,
            markersize=5.5,
            label=model,
        )
        for model, config in MODELS.items()
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper right", fontsize=7.4, handletextpad=0.35)
    style_axis(ax)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.95, bottom=0.18)
    fig.savefig(figure_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    table_dir = Path(args.table_dir)
    figure_dir = Path(args.figure_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    folds, summaries = load_results(input_dir)
    build_summary_table(folds, summaries).to_csv(
        table_dir / "clinmave_functional_class_results.csv", index=False
    )
    plot_gof_deltas(folds, figure_dir / "figS5.png")
    plot_weights(folds, figure_dir / "figS6.png")


if __name__ == "__main__":
    main()
