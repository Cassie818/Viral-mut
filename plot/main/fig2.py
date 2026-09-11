#!/usr/bin/env python3
"""Composite Fig. 2: codon information versus PLMs and mutational context."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


MODEL_DIR = Path("Results/ClinVar/model_control")
CONTEXT_DIR = Path("Results/ClinVar/context_control")
FIG_DIR = Path("Figure")

MODEL_FOLDS = MODEL_DIR / "model_control_gene_heldout_fold_results.csv"
MODEL_AUROCS = MODEL_DIR / "model_control_gene_bootstrap_aurocs.csv"
MODEL_BOOTSTRAP = MODEL_DIR / "model_control_gene_bootstrap.csv"
CONTEXT_AUROCS = CONTEXT_DIR / "context_control_gene_bootstrap_aurocs.csv"
CONTEXT_BOOTSTRAP = CONTEXT_DIR / "context_control_gene_bootstrap.csv"


MODEL_ORDER = [
    "ESM-2 150M",
    "ESM-2 650M",
    "ESM-1b 650M",
    "CaLM",
    "ESM-2 150M + ESM-2 650M",
    "ESM-2 650M + ESM-1b 650M",
    "ESM-2 150M + CaLM",
    "ESM-2 650M + CaLM",
    "ESM-1b 650M + CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM",
]
MODEL_LABELS = {
    "ESM-2 150M": "ESM-2\n(150M)",
    "ESM-2 650M": "ESM-2\n(650M)",
    "ESM-1b 650M": "ESM-1b\n(650M)",
    "CaLM": "CaLM",
    "ESM-2 150M + ESM-2 650M": "ESM-2 (150M)\n+ ESM-2 (650M)",
    "ESM-2 650M + ESM-1b 650M": "ESM-2 (650M)\n+ ESM-1b (650M)",
    "ESM-2 150M + CaLM": "ESM-2 (150M)\n+ CaLM",
    "ESM-2 650M + CaLM": "ESM-2 (650M)\n+ CaLM",
    "ESM-1b 650M + CaLM": "ESM-1b (650M)\n+ CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM": "ESM-2 (650M)\n+ ESM-1b (650M)\n+ CaLM",
}
MODEL_COLORS = {
    "ESM-2 150M": "#D7DEE2",
    "ESM-2 650M": "#9BC9E7",
    "ESM-1b 650M": "#AFCF9F",
    "CaLM": "#EBCB70",
    "ESM-2 150M + ESM-2 650M": "#8ED1CC",
    "ESM-2 650M + ESM-1b 650M": "#C4B2E6",
    "ESM-2 150M + CaLM": "#E7A96B",
    "ESM-2 650M + CaLM": "#F0B7A4",
    "ESM-1b 650M + CaLM": "#E59AB2",
    "ESM-2 650M + ESM-1b 650M + CaLM": "#8FA3B3",
}
MODEL_SCORE_COLUMNS = {
    "ESM-2 150M": "score_esm2_150m",
    "ESM-2 650M": "score_esm2_650m",
    "ESM-1b 650M": "score_esm1b_650m",
    "CaLM": "score_calm",
    "ESM-2 150M + ESM-2 650M": "score_esm2_150m_esm2_650m",
    "ESM-2 650M + ESM-1b 650M": "score_esm2_650m_esm1b_650m",
    "ESM-2 150M + CaLM": "score_esm2_150m_calm",
    "ESM-2 650M + CaLM": "score_esm2_650m_calm",
    "ESM-1b 650M + CaLM": "score_esm1b_650m_calm",
    "ESM-2 650M + ESM-1b 650M + CaLM": "score_esm2_650m_esm1b_650m_calm",
}

GAIN_ORDER = [
    ("ESM-2 150M + ESM-2 650M", "ESM-2 650M", "ESM-2 (150M)\n+ ESM-2 (650M)"),
    ("ESM-2 650M + ESM-1b 650M", "ESM-2 650M", "ESM-2 (650M)\n+ ESM-1b (650M)"),
    ("ESM-2 650M + CaLM", "ESM-2 650M", "ESM-2 (650M)\n+ CaLM"),
    ("ESM-1b 650M + CaLM", "ESM-1b 650M", "ESM-1b (650M)\n+ CaLM"),
    (
        "ESM-2 650M + ESM-1b 650M + CaLM",
        "ESM-2 650M + ESM-1b 650M",
        "Triple model",
    ),
]
GAIN_COLORS = [
    MODEL_COLORS["ESM-2 150M + ESM-2 650M"],
    MODEL_COLORS["ESM-2 650M + ESM-1b 650M"],
    MODEL_COLORS["ESM-2 650M + CaLM"],
    MODEL_COLORS["ESM-1b 650M + CaLM"],
    MODEL_COLORS["ESM-2 650M + ESM-1b 650M + CaLM"],
]

CONTEXT_ORDER = [
    "ESM-2 650M",
    "ESM-2 650M + mutational context",
    "ESM-2 650M + CaLM",
    "ESM-2 650M + mutational context + CaLM",
]
CONTEXT_LABELS = [
    "ESM-2\n(650M)",
    "+ context",
    "+ CaLM",
    "+ context\n+ CaLM",
]
CONTEXT_COLORS = ["#A8D5F0", "#F6E6A3", "#F8D5C7", "#DCE2E7"]
CONTEXT_SCORE_COLUMNS = {
    "ESM-2 650M": "score_esm2_650m",
    "ESM-2 650M + mutational context": "score_esm2_650m_context",
    "ESM-2 650M + CaLM": "score_esm2_650m_calm",
    "ESM-2 650M + mutational context + CaLM": "score_esm2_650m_context_calm",
}
SOFT_EDGE = "#9A9A9A"
DARK_EDGE = "#4A4A4A"

WEIGHT_ORDER = [
    "ESM-2 150M + ESM-2 650M",
    "ESM-2 650M + ESM-1b 650M",
    "ESM-2 150M + CaLM",
    "ESM-2 650M + CaLM",
    "ESM-1b 650M + CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM",
]
WEIGHT_LABELS = {
    "ESM-2 150M + ESM-2 650M": "ESM-2 (150M)\n+ ESM-2 (650M)",
    "ESM-2 650M + ESM-1b 650M": "ESM-2 (650M)\n+ ESM-1b (650M)",
    "ESM-2 150M + CaLM": "ESM-2 (150M)\n+ CaLM",
    "ESM-2 650M + CaLM": "ESM-2 (650M)\n+ CaLM",
    "ESM-1b 650M + CaLM": "ESM-1b (650M)\n+ CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM": "ESM-2 (650M)\n+ ESM-1b (650M)\n+ CaLM",
}
WEIGHT_COMPONENTS = [
    ("w_esm2_150m", "ESM-2 (150M)", MODEL_COLORS["ESM-2 150M"]),
    ("w_esm2_650m", "ESM-2 (650M)", MODEL_COLORS["ESM-2 650M"]),
    ("w_esm1b_650m", "ESM-1b (650M)", MODEL_COLORS["ESM-1b 650M"]),
    ("w_calm", "CaLM", MODEL_COLORS["CaLM"]),
]

MODEL_COMPARISONS = [
    ("ESM-2 650M", "ESM-2 150M + ESM-2 650M", 0),
    ("ESM-2 650M", "ESM-2 650M + ESM-1b 650M", 1),
    ("ESM-2 150M", "ESM-2 150M + CaLM", 2),
    ("ESM-2 650M", "ESM-2 650M + CaLM", 3),
    ("ESM-1b 650M", "ESM-1b 650M + CaLM", 4),
    ("ESM-2 650M + ESM-1b 650M", "ESM-2 650M + ESM-1b 650M + CaLM", 5),
]


def add_panel_label(ax, label: str, title: str) -> None:
    ax.text(
        -0.095,
        1.035,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=10.5,
        clip_on=False,
    )


def add_bracket(ax, x1: float, x2: float, y: float, text: str, height: float = 0.0016) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color=DARK_EDGE, linewidth=0.75)
    ax.text((x1 + x2) / 2, y + height + 0.0008, text, ha="center", va="bottom", fontsize=7.2)


def bootstrap_delta(bootstrap: pd.DataFrame, combo: str, baseline: str) -> float:
    comparison = f"{combo} vs {baseline}"
    hit = bootstrap.loc[bootstrap["comparison"] == comparison, "pooled_oof_delta_auroc"]
    if len(hit) != 1:
        raise ValueError(f"Expected one bootstrap comparison for {comparison}, found {len(hit)}")
    return float(hit.iloc[0])


def add_model_comparison_bracket(
    ax,
    x_lookup: dict[str, int],
    y_lookup: dict[str, float],
    yerr_lookup: dict[str, float],
    baseline: str,
    combo: str,
    delta: float,
    level: int,
) -> None:
    x1 = x_lookup[baseline]
    x2 = x_lookup[combo]
    low, high = sorted([x1, x2])
    y = 0.952 + level * 0.010
    h = 0.0025
    ax.plot([low, low, high, high], [y, y + h, y + h, y], color=DARK_EDGE, linewidth=0.72, clip_on=False)
    ax.text((low + high) / 2, y + h + 0.0012, rf"$\Delta=+{delta:.4f}$", ha="center", va="bottom", fontsize=6.3)


def plot_auroc_panel(ax, summary: pd.DataFrame, bootstrap: pd.DataFrame) -> None:
    summary = summary.set_index("score_column")
    x = np.arange(len(MODEL_ORDER))
    score_columns = [MODEL_SCORE_COLUMNS[model] for model in MODEL_ORDER]
    means = summary.loc[score_columns, "pooled_oof_auroc"].to_numpy()
    lows = summary.loc[score_columns, "bootstrap_95ci_low"].to_numpy()
    highs = summary.loc[score_columns, "bootstrap_95ci_high"].to_numpy()
    yerr = np.vstack([means - lows, highs - means])
    ax.bar(
        x,
        means,
        yerr=yerr,
        capsize=2.3,
        width=0.72,
        color=[MODEL_COLORS[model] for model in MODEL_ORDER],
        edgecolor=SOFT_EDGE,
        alpha=0.76,
        linewidth=0.9,
        error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": DARK_EDGE},
        zorder=2,
    )
    for idx, model in enumerate(MODEL_ORDER):
        ax.text(idx, highs[idx] + 0.0025, f"{means[idx]:.3f}", ha="center", va="bottom", fontsize=7.0, color="#2F2F2F")

    x_lookup = {model: idx for idx, model in enumerate(MODEL_ORDER)}
    y_lookup = {model: mean for model, mean in zip(MODEL_ORDER, means)}
    yerr_lookup = {model: err for model, err in zip(MODEL_ORDER, yerr[1])}
    for baseline, combo, level in MODEL_COMPARISONS:
        delta = bootstrap_delta(bootstrap, combo, baseline)
        add_model_comparison_bracket(ax, x_lookup, y_lookup, yerr_lookup, baseline, combo, delta, level)

    add_panel_label(ax, "A", "Single-model and combined-model AUROC")
    ax.set_ylabel("Pooled out-of-fold AUROC")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_LABELS[model] for model in MODEL_ORDER], rotation=34, ha="right", fontsize=6.4)
    ax.set_ylim(0.80, 1.020)
    ax.set_yticks(np.arange(0.80, 1.011, 0.035))
    ax.tick_params(axis="x", length=0)


def plot_weight_panel(ax, folds: pd.DataFrame) -> None:
    standardized = pd.read_csv(MODEL_DIR / "standardized_equivalent_weights_by_fold.csv")
    standardized = standardized.pivot_table(
        index=["ensemble", "fold"], columns="component",
        values="standardized_equivalent_weight", fill_value=0.0,
    ).reset_index()
    rows = []
    for model in WEIGHT_ORDER:
        sub = standardized[standardized["ensemble"] == model]
        row = {"model": model}
        component_names = {"w_esm2_150m": "ESM-2 150M", "w_esm2_650m": "ESM-2 650M",
                           "w_esm1b_650m": "ESM-1b 650M", "w_calm": "CaLM"}
        for col, _, _ in WEIGHT_COMPONENTS:
            component = component_names[col]
            row[col] = float(sub[component].mean()) if component in sub.columns else 0.0
        rows.append(row)
    summary = pd.DataFrame(rows)

    y = np.arange(len(summary))
    left = np.zeros(len(summary))
    for col, label, color in WEIGHT_COMPONENTS:
        vals = summary[col].to_numpy()
        ax.barh(
            y,
            vals,
            left=left,
            height=0.62,
            color=color,
            edgecolor=SOFT_EDGE,
            linewidth=0.9,
            label=label,
            zorder=2,
        )
        for yi, val, start in zip(y, vals, left):
            if val >= 0.07:
                ax.text(start + val / 2, yi, f"{val:.2f}", ha="center", va="center", fontsize=7.0, color="#2F2F2F")
        left += vals

    add_panel_label(ax, "B", "Scale-adjusted ensemble weights")
    ax.set_xlabel("Scale-adjusted ensemble weight")
    ax.set_xlim(0, 1)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(y)
    ax.set_yticklabels([WEIGHT_LABELS[model] for model in WEIGHT_ORDER], fontsize=6.8, linespacing=0.9)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.legend(
        frameon=False,
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.19),
        fontsize=7.0,
        handlelength=1.15,
        columnspacing=1.0,
    )


def plot_gain_panel(ax, bootstrap: pd.DataFrame) -> None:
    bootstrap = bootstrap.set_index("comparison")
    rows = []
    for a, b, label in GAIN_ORDER:
        row = bootstrap.loc[f"{a} vs {b}"]
        mean = float(row["pooled_oof_delta_auroc"])
        low = float(row["bootstrap_95ci_low"])
        high = float(row["bootstrap_95ci_high"])
        rows.append((label, mean, low, high))
    y = np.arange(len(rows))[::-1]
    means = np.array([row[1] for row in rows])
    lows = np.array([row[2] for row in rows])
    highs = np.array([row[3] for row in rows])
    ax.axvline(0, color=DARK_EDGE, linewidth=0.8)
    ax.barh(
        y,
        means,
        xerr=np.vstack([means - lows, highs - means]),
        color=GAIN_COLORS,
        edgecolor=SOFT_EDGE,
        alpha=0.76,
        linewidth=0.9,
        capsize=2.5,
        error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": DARK_EDGE},
        zorder=2,
    )
    for yi, (_, mean, _, high) in zip(y, rows):
        ax.text(high + 0.002, yi, rf"$\Delta=+{mean:.3f}$", va="center", fontsize=7.4)
    add_panel_label(ax, "B", "Incremental gain over stronger constituent")
    ax.set_xlabel("ΔAUROC")
    ax.set_yticks(y)
    ax.set_yticklabels([row[0] for row in rows], fontsize=7.0, linespacing=0.9)
    ax.set_xlim(-0.004, 0.064)
    ax.set_xticks(np.arange(0, 0.065, 0.01))
    ax.tick_params(axis="y", length=0)


def plot_context_panel(ax, summary: pd.DataFrame, bootstrap: pd.DataFrame) -> None:
    summary = summary.set_index("score_column")
    x = np.arange(len(CONTEXT_ORDER))
    score_columns = [CONTEXT_SCORE_COLUMNS[model] for model in CONTEXT_ORDER]
    means = summary.loc[score_columns, "pooled_oof_auroc"].to_numpy()
    lows = summary.loc[score_columns, "bootstrap_95ci_low"].to_numpy()
    highs = summary.loc[score_columns, "bootstrap_95ci_high"].to_numpy()
    yerr = np.vstack([means - lows, highs - means])
    ax.bar(
        x,
        means,
        yerr=yerr,
        capsize=2.5,
        width=0.68,
        color=CONTEXT_COLORS,
        edgecolor=SOFT_EDGE,
        alpha=0.76,
        linewidth=0.9,
        error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": DARK_EDGE},
        zorder=2,
    )
    for idx, model in enumerate(CONTEXT_ORDER):
        ax.text(idx, means[idx] + yerr[1, idx] + 0.004, f"{means[idx]:.3f}", ha="center", fontsize=7.6)

    bootstrap = bootstrap.set_index("comparison")
    context_gain = bootstrap.loc["ESM-2 650M + mutational context vs ESM-2 650M"]
    calm_gain = bootstrap.loc["ESM-2 650M + CaLM vs ESM-2 650M"]
    calm_after_context = bootstrap.loc[
        "ESM-2 650M + mutational context + CaLM vs ESM-2 650M + mutational context"
    ]
    context_after_calm = bootstrap.loc[
        "ESM-2 650M + mutational context + CaLM vs ESM-2 650M + CaLM"
    ]
    add_bracket(ax, 0, 1, 0.927, rf"$\Delta=+{context_gain['pooled_oof_delta_auroc']:.3f}$")
    add_bracket(ax, 0, 2, 0.933, rf"$\Delta=+{calm_gain['pooled_oof_delta_auroc']:.3f}$")
    add_bracket(ax, 1, 3, 0.939, rf"$\Delta=+{calm_after_context['pooled_oof_delta_auroc']:.3f}$")
    add_bracket(ax, 2, 3, 0.945, rf"$\Delta=+{context_after_calm['pooled_oof_delta_auroc']:.3f}$")

    add_panel_label(ax, "C", "Fixed ESM-2 (650M) baseline with context and CaLM")
    ax.set_ylabel("Pooled out-of-fold AUROC")
    ax.set_xticks(x)
    ax.set_xticklabels(CONTEXT_LABELS, fontsize=8)
    ax.set_ylim(0.875, 0.953)
    ax.set_yticks(np.arange(0.88, 0.954, 0.01))
    ax.tick_params(axis="x", length=0)


def plot_conditional_panel(ax, bootstrap: pd.DataFrame) -> None:
    bootstrap = bootstrap.set_index("comparison")
    specs = [
        (
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + mutational context",
            "CaLM\ncontext-adjusted",
            "#F8D5C7",
        ),
        (
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + CaLM",
            "Context\nCaLM-adjusted",
            "#F6E6A3",
        ),
    ]
    positions = np.array([0.0, 0.72])
    labels = []
    colors = []
    rows = []
    for a, b, label, color in specs:
        row = bootstrap.loc[f"{a} vs {b}"]
        labels.append(label)
        colors.append(color)
        rows.append(row)

    means = np.array([float(row["pooled_oof_delta_auroc"]) for row in rows])
    lows = np.array([float(row["bootstrap_95ci_low"]) for row in rows])
    highs = np.array([float(row["bootstrap_95ci_high"]) for row in rows])
    bars = ax.bar(
        positions,
        means,
        width=0.38,
        yerr=np.vstack([means - lows, highs - means]),
        color=colors,
        edgecolor=SOFT_EDGE,
        linewidth=0.9,
        capsize=2.5,
        error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": DARK_EDGE},
        zorder=2,
    )
    for bar, mean, high in zip(bars, means, highs):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            high + 0.0007,
            rf"$\Delta=+{mean:.3f}$",
            ha="center",
            va="bottom",
            fontsize=7.6,
        )
    add_panel_label(ax, "D", "Conditional gain after adjustment")
    ax.axhline(0, color=DARK_EDGE, linewidth=0.8)
    ax.set_ylabel("Pooled out-of-fold ΔAUROC")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlim(-0.30, 1.05)
    ax.set_ylim(0.002, 0.023)
    ax.set_yticks(np.arange(0.004, 0.023, 0.004))
    ax.tick_params(axis="x", length=0)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    model_summary = pd.read_csv(MODEL_AUROCS)
    model_folds = pd.read_csv(MODEL_FOLDS)
    model_bootstrap = pd.read_csv(MODEL_BOOTSTRAP)
    context_summary = pd.read_csv(CONTEXT_AUROCS)
    context_bootstrap = pd.read_csv(CONTEXT_BOOTSTRAP)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.45), gridspec_kw={"width_ratios": [1.56, 1.0]})
    plot_auroc_panel(axes[0, 0], model_summary, model_bootstrap)
    plot_weight_panel(axes[0, 1], model_folds)
    plot_context_panel(axes[1, 0], context_summary, context_bootstrap)
    plot_conditional_panel(axes[1, 1], context_bootstrap)

    fig.subplots_adjust(left=0.062, right=0.985, top=0.94, bottom=0.08, hspace=0.56, wspace=0.25)

    path = FIG_DIR / "fig2.png"
    fig.savefig(path, dpi=600, bbox_inches="tight")
    print(path)


if __name__ == "__main__":
    main()
