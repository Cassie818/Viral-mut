#!/usr/bin/env python3
"""Fig. 3: ClinMAVE gene-wise modality contribution summary."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


BASE = Path("Results/ClinMAVE/functional_effects")
FIG_DIR = Path("Figure")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/viral_mut_mpl_cache")

import matplotlib.pyplot as plt


SUMMARY_FILES = {
    "ESM-2 (650M)": BASE / "clinmave_esm2_650m_calm_cv_summary.csv",
    "ESM-1b (650M)": BASE / "clinmave_esm1b_650m_calm_cv_summary.csv",
}

FOLD_FILES = {
    "ESM-2 (650M)": BASE / "clinmave_esm2_650m_calm_cv_fold_metrics.csv",
    "ESM-1b (650M)": BASE / "clinmave_esm1b_650m_calm_cv_fold_metrics.csv",
}

MODELS = ["ESM-2 (650M)", "ESM-1b (650M)"]
GROUPS = [
    ("DMS", "lof", "DMS\nLoF"),
    ("DMS", "gof", "DMS\nGoF"),
    ("CBGE", "lof", "CBGE\nLoF"),
    ("CBGE", "gof", "CBGE\nGoF"),
]
MODEL_COLORS = {
    "ESM-2 (650M)": "#A9D8EA",
    "ESM-1b (650M)": "#C8DDBA",
}
EDGE = "#868686"
DARK = "#383838"
LIGHT_EDGE = "#A7A7A7"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    folds = []
    for model, path in SUMMARY_FILES.items():
        df = pd.read_csv(path)
        df["model"] = model
        summaries.append(df)
    for model, path in FOLD_FILES.items():
        df = pd.read_csv(path)
        df["model"] = model
        folds.append(df)
    return pd.concat(summaries, ignore_index=True), pd.concat(folds, ignore_index=True)


def row_for(summary: pd.DataFrame, model: str, assay: str, case_class: str) -> pd.Series:
    hit = summary[
        (summary["model"] == model)
        & (summary["assay"] == assay)
        & (summary["case_class"] == case_class)
    ]
    if len(hit) != 1:
        raise ValueError(f"Expected one row for {model} {assay} {case_class}, got {len(hit)}")
    return hit.iloc[0]


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        clip_on=False,
    )


def format_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", width=1.1, length=3.5, color=EDGE, labelcolor=DARK)


def metric_columns(model: str) -> tuple[str, str, str]:
    if model == "ESM-2 (650M)":
        return (
            "delta_combo_vs_esm2_650m_mean",
            "delta_combo_vs_esm2_650m_sd",
            "p_delta_combo_vs_esm2_650m_paired_t",
        )
    return (
        "delta_combo_vs_esm1b_650m_mean",
        "delta_combo_vs_esm1b_650m_sd",
        "p_delta_combo_vs_esm1b_650m_paired_t",
    )


def fold_delta_ci(folds: pd.DataFrame, model: str, assay: str, case_class: str) -> tuple[float, bool]:
    if model == "ESM-2 (650M)":
        col = "delta_combo_vs_esm2_650m"
    else:
        col = "delta_combo_vs_esm1b_650m"
    sub = folds[
        (folds["model"] == model)
        & (folds["assay"] == assay)
        & (folds["case_class"] == case_class)
    ][col].dropna()
    if len(sub) < 2:
        return np.nan, True
    mean = sub.mean()
    sem = stats.sem(sub)
    ci_half = stats.t.ppf(0.975, len(sub) - 1) * sem
    return ci_half, (mean - ci_half <= 0 <= mean + ci_half)


def collect_plot_rows(summary: pd.DataFrame, folds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_index, (assay, case_class, group_label) in enumerate(GROUPS):
        for model in MODELS:
            row = row_for(summary, model, assay, case_class)
            delta_mean, delta_sd, delta_p = metric_columns(model)
            delta_ci_half, delta_ci_crosses_zero = fold_delta_ci(folds, model, assay, case_class)
            rows.append(
                {
                    "group_index": group_index,
                    "assay": assay,
                    "case_class": case_class,
                    "group_label": group_label.replace("\n", " "),
                    "model": model,
                    "mean_calm_weight": row["mean_calm_weight"],
                    "sd_calm_weight": row["sd_calm_weight"],
                    "delta_auroc_vs_plm": row[delta_mean],
                    "sd_delta_auroc_vs_plm": row[delta_sd],
                    "ci95_delta_auroc_vs_plm": delta_ci_half,
                    "ci95_delta_crosses_zero": delta_ci_crosses_zero,
                    "p_delta_auroc_vs_plm": row[delta_p],
                    "n_folds": int(row["n_folds"]),
                    "n_evaluable_folds": int(row["n_evaluable_folds"]),
                    "n_variants": int(row["n_variants"]),
                    "n_genes": int(row["n_genes"]),
                }
            )
    return pd.DataFrame(rows)


def draw_bars(
    ax: plt.Axes,
    plot_df: pd.DataFrame,
    metric_col: str,
    sd_col: str,
    ylabel: str,
    ylim: tuple[float, float],
    yticks: np.ndarray,
    label: str,
    value_fmt: str,
) -> None:
    x = np.arange(len(GROUPS))
    width = 0.28
    offsets = {"ESM-2 (650M)": -0.20, "ESM-1b (650M)": 0.20}

    for model in MODELS:
        sub = plot_df[plot_df["model"] == model].sort_values("group_index")
        xs = x + offsets[model]
        vals = sub[metric_col].to_numpy(float)
        errs = sub[sd_col].to_numpy(float)
        bars = ax.bar(
            xs,
            vals,
            yerr=errs,
            width=width,
            color=MODEL_COLORS[model],
            edgecolor=LIGHT_EDGE,
            linewidth=1.15,
            capsize=2.7,
            error_kw={"elinewidth": 1.0, "capthick": 1.0, "ecolor": DARK},
            label=model,
            zorder=2,
        )

        for bar, val, err, group_index, n_eval, n_folds in zip(
            bars,
            vals,
            errs,
            sub["group_index"].astype(int),
            sub["n_evaluable_folds"].astype(int),
            sub["n_folds"].astype(int),
        ):
            if metric_col == "delta_auroc_vs_plm":
                va = "bottom" if val >= 0 else "top"
                y = val + err + 0.006 if val >= 0 else val - 0.007
                fontsize = 7.0
                if group_index == 3 and val < 0:
                    y = -0.026 if model == "ESM-2 (650M)" else -0.011
            else:
                va = "bottom"
                y = val + 0.035
                fontsize = 7.4
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                y,
                value_fmt.format(val),
                ha="center",
                va=va,
                fontsize=fontsize,
                color=DARK,
            )
            if n_eval < n_folds:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    ylim[0] + 0.055 * (ylim[1] - ylim[0]),
                    f"{n_eval}/{n_folds}",
                    ha="center",
                    va="bottom",
                    fontsize=7.1,
                    color=DARK,
                )

    if metric_col == "delta_auroc_vs_plm":
        ax.axhline(0, color=DARK, linewidth=0.9, zorder=1)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, _, label in GROUPS])
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    add_panel_label(ax, label)
    format_axes(ax)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    summary, folds = load_data()
    plot_df = collect_plot_rows(summary, folds)
    plot_df.to_csv(BASE / "fig3_plot_data.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.6,
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.65, 3.15), constrained_layout=False)
    draw_bars(
        axes[0],
        plot_df,
        "mean_calm_weight",
        "sd_calm_weight",
        "Optimised CaLM weight",
        (0.0, 1.12),
        np.linspace(0.0, 1.0, 6),
        "A",
        "{:.2f}",
    )
    draw_bars(
        axes[1],
        plot_df,
        "delta_auroc_vs_plm",
        "ci95_delta_auroc_vs_plm",
        r"$\Delta$AUROC vs PLM",
        (-0.09, 0.18),
        np.arange(-0.05, 0.181, 0.05),
        "B",
        "{:+.3f}",
    )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.54, -0.02),
        fontsize=8.2,
        handlelength=1.3,
        columnspacing=1.6,
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.94, bottom=0.24, wspace=0.32)

    for filename in [
        "fig3_clinmave_genewise_modality.png",
        "fig3.png",
    ]:
        path = FIG_DIR / filename
        fig.savefig(path, dpi=600, bbox_inches="tight")
        print(path)


if __name__ == "__main__":
    main()
