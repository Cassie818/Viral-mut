#!/usr/bin/env python3
"""Regenerate bootstrap-aware supplementary figures and summary table."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Figure"
TABLE_OUT = ROOT / "Results/ClinVar/model_control"
MODEL = ROOT / "Results/ClinVar/model_control"
CONTEXT = ROOT / "Results/ClinVar/context_control"
MAVE = ROOT / "Results/ClinMAVE/functional_effects"

EDGE = "#8F918D"
DARK = "#333333"
ZERO = "#AAAAA5"
COLORS = ["#D7DEE2", "#C4B2E6", "#A7DAD2", "#F4D986", "#F2B9A8", "#B9D8E8", "#9DB6C6"]


def style(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(EDGE)
    ax.spines["bottom"].set_color(EDGE)
    ax.tick_params(color=EDGE, labelcolor=DARK)


def fold_deltas(folds, model_a, model_b):
    a = folds[folds.model == model_a].set_index("fold").test_auc
    b = folds[folds.model == model_b].set_index("fold").test_auc
    return (a - b).dropna().to_numpy()


def interval_rows(ax, labels, folds, points, lows, highs, colors, xlim, xlabel):
    y = np.arange(len(labels))[::-1]
    rng = np.random.default_rng(20260908)
    ax.axvline(0, color=ZERO, lw=0.8, ls="--", zorder=0)
    for yi, values, point, low, high, color in zip(y, folds, points, lows, highs, colors):
        jitter = rng.normal(0, 0.055, len(values))
        ax.scatter(values, np.full(len(values), yi) + jitter, s=23, color=color,
                   edgecolor=EDGE, linewidth=0.6, alpha=0.72, zorder=2)
        ax.errorbar(point, yi, xerr=[[point-low], [high-point]], fmt="o", ms=5.5,
                    color=DARK, ecolor=DARK, elinewidth=1.2, capsize=4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
    style(ax)


def fig_s1():
    folds = pd.read_csv(MODEL / "model_control_gene_heldout_fold_results.csv")
    boot = pd.read_csv(MODEL / "model_control_gene_bootstrap.csv").set_index("comparison")
    specs = [
        ("ESM-2 150M + ESM-2 650M vs ESM-2 650M", "ESM-2 150M + ESM-2 650M", "ESM-2 650M",
         "ESM-2 (150M) + ESM-2 (650M)\nvs ESM-2 (650M)"),
        ("ESM-2 650M + ESM-1b 650M vs ESM-2 650M", "ESM-2 650M + ESM-1b 650M", "ESM-2 650M",
         "ESM-2 (650M) + ESM-1b (650M)\nvs ESM-2 (650M)"),
        ("ESM-2 650M + ESM-1b 650M vs ESM-1b 650M", "ESM-2 650M + ESM-1b 650M", "ESM-1b 650M",
         "ESM-2 (650M) + ESM-1b (650M)\nvs ESM-1b (650M)"),
        ("ESM-2 650M + CaLM vs ESM-2 650M", "ESM-2 650M + CaLM", "ESM-2 650M",
         "ESM-2 (650M) + CaLM\nvs ESM-2 (650M)"),
        ("ESM-1b 650M + CaLM vs ESM-1b 650M", "ESM-1b 650M + CaLM", "ESM-1b 650M",
         "ESM-1b (650M) + CaLM\nvs ESM-1b (650M)"),
        ("ESM-2 650M + ESM-1b 650M + CaLM vs ESM-2 650M + ESM-1b 650M",
         "ESM-2 650M + ESM-1b 650M + CaLM", "ESM-2 650M + ESM-1b 650M",
         "Triple model vs\nESM-2 (650M) + ESM-1b (650M)"),
    ]
    labels, fold_values, points, lows, highs = [], [], [], [], []
    for name, a, b, label in specs:
        row = boot.loc[name]
        labels.append(label)
        fold_values.append(fold_deltas(folds, a, b))
        points.append(row.pooled_oof_delta_auroc)
        lows.append(row.bootstrap_95ci_low)
        highs.append(row.bootstrap_95ci_high)
    fig, ax = plt.subplots(figsize=(12.2, 4.8))
    x = np.arange(len(labels))
    rng = np.random.default_rng(20260909)
    ax.axhline(0, color=ZERO, lw=0.8, ls="--", zorder=0)
    ax.bar(x, points, width=0.58, color=COLORS[:len(labels)],
           edgecolor=EDGE, linewidth=0.8, alpha=0.62, zorder=1)
    for xi, values, point, low, high, color in zip(
            x, fold_values, points, lows, highs, COLORS):
        jitter = rng.normal(0, 0.055, len(values))
        ax.scatter(np.full(len(values), xi) + jitter, values, s=27,
                   color=color, edgecolor=EDGE, linewidth=0.6,
                   alpha=0.95, zorder=3)
        ax.errorbar(xi, point, yerr=[[point - low], [high - point]],
                    fmt="none", color=DARK, ecolor=DARK,
                    elinewidth=1.2, capsize=4, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(-0.007, 0.070)
    ax.set_ylabel(r"$\Delta$AUROC")
    ax.tick_params(axis="x", length=0)
    style(ax)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.29)
    fig.savefig(OUT / "figS1.png", dpi=600, bbox_inches="tight")


def aligned_distribution(path, comparison):
    data = pd.read_csv(path)
    return data[data.comparison == comparison].sort_values("bootstrap_replicate").delta_auroc.to_numpy()


def fig_s3():
    folds = pd.read_csv(MODEL / "model_control_gene_heldout_fold_results.csv")
    boot = pd.read_csv(MODEL / "model_control_gene_bootstrap.csv").set_index("comparison")
    dist_path = MODEL / "model_control_gene_bootstrap_distributions.csv.gz"
    weak = "ESM-2 150M + ESM-2 650M vs ESM-2 650M"
    strong = "ESM-2 650M + ESM-1b 650M vs ESM-1b 650M"
    calm = "ESM-2 650M + CaLM vs ESM-2 650M"
    weak_boot = aligned_distribution(dist_path, weak)
    strong_boot = aligned_distribution(dist_path, strong)
    calm_boot = aligned_distribution(dist_path, calm)
    comparisons = [
        ("CaLM - ESM-2 (150M)", calm_boot - weak_boot,
         boot.loc[calm, "pooled_oof_delta_auroc"] - boot.loc[weak, "pooled_oof_delta_auroc"]),
        ("CaLM - ESM-1b (650M)", calm_boot - strong_boot,
         boot.loc[calm, "pooled_oof_delta_auroc"] - boot.loc[strong, "pooled_oof_delta_auroc"]),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), gridspec_kw={"width_ratios": [1.25, 1.0]})
    gains = [
        fold_deltas(folds, "ESM-2 150M + ESM-2 650M", "ESM-2 650M"),
        fold_deltas(folds, "ESM-2 650M + ESM-1b 650M", "ESM-1b 650M"),
        fold_deltas(folds, "ESM-2 650M + CaLM", "ESM-2 650M"),
    ]
    x = np.arange(3)
    rng = np.random.default_rng(17)
    panel_colors = ["#D7DEE2", "#F4D986", "#A7DAD2"]
    pooled = [
        boot.loc[weak, "pooled_oof_delta_auroc"],
        boot.loc[strong, "pooled_oof_delta_auroc"],
        boot.loc[calm, "pooled_oof_delta_auroc"],
    ]
    gain_lows = [boot.loc[weak, "bootstrap_95ci_low"],
                 boot.loc[strong, "bootstrap_95ci_low"],
                 boot.loc[calm, "bootstrap_95ci_low"]]
    gain_highs = [boot.loc[weak, "bootstrap_95ci_high"],
                  boot.loc[strong, "bootstrap_95ci_high"],
                  boot.loc[calm, "bootstrap_95ci_high"]]
    axes[0].bar(x, pooled, width=0.58, color=panel_colors,
                edgecolor=EDGE, linewidth=0.8, alpha=0.62, zorder=1)
    for xi, values, color in zip(x, gains, panel_colors):
        axes[0].scatter(np.full(len(values), xi) + rng.normal(0, 0.035, len(values)), values,
                        s=25, color=color, edgecolor=EDGE, linewidth=0.6,
                        alpha=0.95, zorder=3)
    axes[0].errorbar(x, pooled,
                     yerr=[np.asarray(pooled) - np.asarray(gain_lows),
                           np.asarray(gain_highs) - np.asarray(pooled)],
                     fmt="none", ecolor=DARK, elinewidth=1.2, capsize=4, zorder=4)
    axes[0].axhline(0, color=ZERO, lw=0.8, ls="--")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["ESM-2 (150M)\nprotein control",
                             "ESM-1b (650M)\nprotein control",
                             "CaLM\ncodon model"])
    axes[0].set_ylabel(r"$\Delta$AUROC")
    axes[0].text(-0.10, 1.02, "A", transform=axes[0].transAxes, fontweight="bold", fontsize=11)
    style(axes[0])

    fold_differences = [gains[2] - gains[0], gains[2] - gains[1]]
    labels, values, points, lows, highs = [], [], [], [], []
    for (label, delta, observed), fold_difference in zip(comparisons, fold_differences):
        labels.append(label)
        values.append(fold_difference)
        points.append(float(observed))
        low, high = np.quantile(delta, [0.025, 0.975])
        lows.append(low)
        highs.append(high)
    bx = np.arange(len(labels))
    bcolors = ["#A7DAD2", "#F4D986"]
    axes[1].axhline(0, color=ZERO, lw=0.8, ls="--", zorder=0)
    axes[1].bar(bx, points, width=0.58, color=bcolors, edgecolor=EDGE,
                linewidth=0.8, alpha=0.62, zorder=1)
    for xi, fold_values, color in zip(bx, values, bcolors):
        axes[1].scatter(np.full(len(fold_values), xi) + rng.normal(0, 0.035, len(fold_values)),
                        fold_values, s=25, color=color, edgecolor=EDGE,
                        linewidth=0.6, alpha=0.95, zorder=3)
    axes[1].errorbar(bx, points,
                     yerr=[np.asarray(points) - np.asarray(lows),
                           np.asarray(highs) - np.asarray(points)],
                     fmt="none", ecolor=DARK, elinewidth=1.2, capsize=4, zorder=4)
    axes[1].set_xticks(bx)
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel(r"Difference in gain ($\Delta$AUROC)")
    axes[1].tick_params(axis="x", length=0)
    style(axes[1])
    for ax in axes:
        ax.set_ylim(-0.010, 0.035)
        ax.set_yticks([-0.01, 0.00, 0.01, 0.02, 0.03])
    axes[1].text(-0.10, 1.02, "B", transform=axes[1].transAxes, fontweight="bold", fontsize=11)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.23, wspace=0.32)
    fig.savefig(OUT / "figS3.png", dpi=600, bbox_inches="tight")


def fig_s5():
    folds = pd.read_csv(CONTEXT / "context_control_fold_results.csv")
    boot = pd.read_csv(CONTEXT / "context_control_gene_bootstrap.csv").set_index("comparison")
    specs = [
        ("ESM-2 650M + mutational context + CaLM vs ESM-2 650M + mutational context",
         "ESM-2 650M + mutational context + CaLM", "ESM-2 650M + mutational context",
         "CaLM beyond mutational context"),
        ("ESM-2 650M + mutational context + CaLM vs ESM-2 650M + CaLM",
         "ESM-2 650M + mutational context + CaLM", "ESM-2 650M + CaLM",
         "Mutational context beyond CaLM"),
    ]
    labels, fold_values, points, lows, highs = [], [], [], [], []
    for name, a, b, label in specs:
        row = boot.loc[name]
        labels.append(label)
        fold_values.append(fold_deltas(folds.rename(columns={"test_auc": "test_auc"}), a, b))
        points.append(row.pooled_oof_delta_auroc)
        lows.append(row.bootstrap_95ci_low)
        highs.append(row.bootstrap_95ci_high)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    x = np.arange(len(labels))
    colors = ["#F4B7B5", "#C8DCA9"]
    rng = np.random.default_rng(20260909)
    ax.axhline(0, color=ZERO, lw=0.8, ls="--", zorder=0)
    ax.bar(x, points, width=0.58, color=colors, edgecolor=EDGE,
           linewidth=0.8, alpha=0.62, zorder=1)
    for xi, values, color in zip(x, fold_values, colors):
        ax.scatter(np.full(len(values), xi) + rng.normal(0, 0.04, len(values)),
                   values, s=27, color=color, edgecolor=EDGE,
                   linewidth=0.6, alpha=0.95, zorder=3)
    ax.errorbar(x, points,
                yerr=[np.asarray(points) - np.asarray(lows),
                      np.asarray(highs) - np.asarray(points)],
                fmt="none", ecolor=DARK, elinewidth=1.2, capsize=4, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(["CaLM beyond\nmutational context",
                        "Mutational context\nbeyond CaLM"])
    ax.set_ylim(-0.002, 0.022)
    ax.set_ylabel(r"Conditional $\Delta$AUROC")
    ax.tick_params(axis="x", length=0)
    style(ax)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.96, bottom=0.24)
    fig.savefig(OUT / "figS5.png", dpi=600, bbox_inches="tight")


def mave_fold_delta(folds, background, assay):
    sub = folds[(folds.assay == assay) & (folds.case_class == "gof")]
    column = "delta_combo_vs_esm2_650m" if background == "ESM-2 (650M)" else "delta_combo_vs_esm1b_650m"
    return sub[column].dropna().to_numpy()


def fig_s6():
    boot = pd.read_csv(MAVE / "clinmave_functional_class_gene_bootstrap.csv")
    fold_files = {
        "ESM-2 (650M)": MAVE / "clinmave_esm2_650m_calm_cv_fold_metrics.csv",
        "ESM-1b (650M)": MAVE / "clinmave_esm1b_650m_calm_cv_fold_metrics.csv",
    }
    fold_tables = {name: pd.read_csv(path) for name, path in fold_files.items()}
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharey=False)
    rng = np.random.default_rng(20260909)
    for ax, assay in zip(axes, ["DMS", "CBGE"]):
        labels, fold_values, points, lows, highs = [], [], [], [], []
        for background in ["ESM-2 (650M)", "ESM-1b (650M)"]:
            row = boot[(boot.plm_background == background) & (boot.assay == assay) & (boot.case_class == "gof")].iloc[0]
            labels.append(background)
            fold_values.append(mave_fold_delta(fold_tables[background], background, assay))
            points.append(row.pooled_oof_delta_auroc)
            lows.append(row.bootstrap_95ci_low)
            highs.append(row.bootstrap_95ci_high)
        x = np.arange(len(labels))
        colors = ["#9DD3E8", "#C5DDA9"]
        ax.axhline(0, color=ZERO, lw=0.8, ls="--", zorder=0)
        ax.bar(x, points, width=0.58, color=colors, edgecolor=EDGE,
               linewidth=0.8, alpha=0.62, zorder=1)
        for xi, values, color in zip(x, fold_values, colors):
            ax.scatter(np.full(len(values), xi) + rng.normal(0, 0.04, len(values)),
                       values, s=27, color=color, edgecolor=EDGE,
                       linewidth=0.6, alpha=0.95, zorder=3)
        ax.errorbar(x, points,
                    yerr=[np.asarray(points) - np.asarray(lows),
                          np.asarray(highs) - np.asarray(points)],
                    fmt="none", ecolor=DARK, elinewidth=1.2, capsize=4, zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylim((-0.13, 0.19) if assay == "DMS" else (-0.05, 0.03))
        ax.set_ylabel(r"$\Delta$AUROC")
        ax.tick_params(axis="x", length=0)
        style(ax)
        ax.set_title(assay, fontsize=10)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.22, wspace=0.28)
    fig.savefig(OUT / "figS6.png", dpi=600, bbox_inches="tight")


def write_table():
    tables = [
        pd.read_csv(MODEL / "model_control_gene_bootstrap.csv"),
        pd.read_csv(CONTEXT / "context_control_gene_bootstrap.csv"),
        pd.read_csv(MAVE / "clinmave_functional_class_gene_bootstrap.csv"),
    ]
    keep = ["analysis", "comparison", "plm_background", "assay", "case_class", "n_variants", "n_genes",
            "pooled_oof_auroc_a", "pooled_oof_auroc_b", "pooled_oof_delta_auroc",
            "bootstrap_95ci_low", "bootstrap_95ci_high", "n_bootstrap_valid", "resampling_unit", "ci_method"]
    table = pd.concat(tables, ignore_index=True, sort=False)
    table = table.reindex(columns=keep)
    table.to_csv(TABLE_OUT / "Supplementary_Table_gene_cluster_bootstrap.csv", index=False)


def main():
    plt.rcParams.update({"font.family": "Arial", "font.size": 8.5, "axes.linewidth": 0.8, "ps.fonttype": 42})
    fig_s1()
    fig_s3()
    fig_s5()
    fig_s6()
    write_table()
    for name in ["figS1.png", "figS3.png", "figS5.png", "figS6.png", "Supplementary_Table_gene_cluster_bootstrap.csv"]:
        print((TABLE_OUT if name.startswith("Supplementary_Table") else OUT) / name)


if __name__ == "__main__":
    main()
