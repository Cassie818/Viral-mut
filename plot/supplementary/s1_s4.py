#!/usr/bin/env python3
"""Supplementary materials for ClinVar model-control analyses.

Outputs:
- all pairwise fold-level statistics, including non-significant comparisons
- fold-level delta AUROC distributions for planned contrasts
- gene-wise versus variant-wise CV comparison
- same-modality placeholder control for generic ensembling
- per-fold ensemble-weight distributions
"""

from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/viral_mut_mpl_cache")

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


MODEL_DIR = Path("Results/ClinVar/model_control")
OUT_DIR = MODEL_DIR
FIG_DIR = Path("Figure")

FOLD_RESULTS = MODEL_DIR / "model_control_gene_heldout_fold_results.csv"
SCORE_TABLE = Path(
    "Results/ClinVar/current_missense_cohort/"
    "clinvar_missense_len1022_complete_cases_model_scores.csv"
)

SCORE_COLUMNS = {
    "ESM-2 150M": "esm2_150m_score",
    "ESM-2 650M": "esm2_650m_score",
    "ESM-1b 650M": "esm1b_650m_score",
    "CaLM": "calm_score",
}

MODEL_SPECS = [
    ("ESM-2 150M", ["ESM-2 150M"]),
    ("ESM-2 650M", ["ESM-2 650M"]),
    ("ESM-1b 650M", ["ESM-1b 650M"]),
    ("CaLM", ["CaLM"]),
    ("ESM-2 150M + ESM-2 650M", ["ESM-2 150M", "ESM-2 650M"]),
    ("ESM-2 650M + ESM-1b 650M", ["ESM-2 650M", "ESM-1b 650M"]),
    ("ESM-2 150M + ESM-1b 650M", ["ESM-2 150M", "ESM-1b 650M"]),
    ("ESM-2 150M + CaLM", ["ESM-2 150M", "CaLM"]),
    ("ESM-2 650M + CaLM", ["ESM-2 650M", "CaLM"]),
    ("ESM-1b 650M + CaLM", ["ESM-1b 650M", "CaLM"]),
    ("ESM-2 650M + ESM-1b 650M + CaLM", ["ESM-2 650M", "ESM-1b 650M", "CaLM"]),
]

PLANNED_CONTRASTS = [
    ("ESM-2 150M + ESM-2 650M", "ESM-2 650M", "ESM-2 150M + ESM-2 650M vs ESM-2 650M"),
    ("ESM-2 650M + ESM-1b 650M", "ESM-2 650M", "ESM-2 650M + ESM-1b 650M vs ESM-2 650M"),
    ("ESM-2 650M + CaLM", "ESM-2 650M", "ESM-2 650M + CaLM vs ESM-2 650M"),
    (
        "ESM-2 650M + CaLM",
        "ESM-2 150M + ESM-2 650M",
        "ESM-2 650M + CaLM vs same-modality placeholder",
    ),
    (
        "ESM-2 650M + ESM-1b 650M",
        "ESM-1b 650M",
        "Strong same-modality placeholder vs ESM-1b 650M",
    ),
    ("ESM-1b 650M + CaLM", "ESM-1b 650M", "ESM-1b 650M + CaLM vs ESM-1b 650M"),
    (
        "ESM-2 650M + ESM-1b 650M + CaLM",
        "ESM-2 650M + ESM-1b 650M",
        "Triple model vs ESM-2 650M + ESM-1b 650M",
    ),
]

MODEL_LABELS = {
    "ESM-2 150M": "ESM-2\n(150M)",
    "ESM-2 650M": "ESM-2\n(650M)",
    "ESM-1b 650M": "ESM-1b\n(650M)",
    "CaLM": "CaLM",
    "ESM-2 150M + ESM-2 650M": "ESM-2 (150M)\n+ ESM-2 (650M)",
    "ESM-2 650M + ESM-1b 650M": "ESM-2 (650M)\n+ ESM-1b (650M)",
    "ESM-2 150M + ESM-1b 650M": "ESM-2 (150M)\n+ ESM-1b (650M)",
    "ESM-2 150M + CaLM": "ESM-2 (150M)\n+ CaLM",
    "ESM-2 650M + CaLM": "ESM-2 (650M)\n+ CaLM",
    "ESM-1b 650M + CaLM": "ESM-1b (650M)\n+ CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM": "ESM-2 (650M)\n+ ESM-1b (650M)\n+ CaLM",
}

COLORS = {
    "grey": "#DDE2E5",
    "blue": "#B9DDE9",
    "green": "#BFE3D8",
    "yellow": "#F2E5A6",
    "rose": "#F0B9AD",
    "lilac": "#D8CAE8",
    "teal": "#A9D8D1",
    "text": "#2F2F2F",
    "edge": "#8E8E8A",
}


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
            "ps.fonttype": 42,
        }
    )


def format_ax(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_color(COLORS["edge"])
        spine.set_linewidth(0.9)
    ax.tick_params(axis="both", color=COLORS["edge"], labelcolor=COLORS["text"], width=0.8, length=3)
    ax.grid(False)


def pair_weights(step: float = 0.01) -> list[tuple[float, float]]:
    values = np.arange(0.0, 1.0 + step / 2, step)
    return [(float(1.0 - w), float(w)) for w in values]


def triple_weights(step: float = 0.05) -> list[tuple[float, float, float]]:
    values = np.arange(0.0, 1.0 + step / 2, step)
    weights = []
    for w1 in values:
        for w2 in values:
            w3 = 1.0 - w1 - w2
            if w3 >= -1e-9:
                weights.append((float(w1), float(w2), float(max(0.0, w3))))
    return weights


def best_weights(y_train: np.ndarray, train_scores: list[np.ndarray], candidates: list[tuple[float, ...]]) -> tuple[float, ...]:
    best_auc = -np.inf
    best_candidate = candidates[0]
    for weights in candidates:
        score = sum(w * s for w, s in zip(weights, train_scores))
        auc = roc_auc_score(y_train, score)
        if auc > best_auc:
            best_auc = auc
            best_candidate = weights
    return best_candidate


def all_pairwise_tests(fold_df: pd.DataFrame) -> pd.DataFrame:
    wide = fold_df.pivot(index="fold", columns="model", values="test_auc")
    summary = fold_df.groupby("model")["test_auc"].agg(["mean", "std"]).rename(columns={"mean": "mean_auc", "std": "sd_auc"})
    planned_unordered = {frozenset((a, b)) for a, b, _ in PLANNED_CONTRASTS}
    rows = []

    def add_row(a: str, b: str, planned_contrast: bool, contrast_label: str) -> None:
        delta = wide[a] - wide[b]
        try:
            wp = float(wilcoxon(delta).pvalue)
        except ValueError:
            wp = np.nan
        paired_p = float(ttest_rel(wide[a], wide[b]).pvalue)
        rows.append(
            {
                "contrast_label": contrast_label,
                "model_a": a,
                "model_b": b,
                "mean_auc_a": float(summary.loc[a, "mean_auc"]),
                "mean_auc_b": float(summary.loc[b, "mean_auc"]),
                "mean_delta_a_minus_b": float(delta.mean()),
                "sd_delta": float(delta.std(ddof=1)),
                "sem_delta": float(delta.sem()),
                "paired_t_p": paired_p,
                "wilcoxon_p": wp,
                "significance_paired_t": p_label(paired_p),
                "planned_contrast": planned_contrast,
                "fold_deltas": ";".join(f"{v:.5f}" for v in delta),
            }
        )

    for a, b, label in PLANNED_CONTRASTS:
        add_row(a, b, True, label)

    for a, b in combinations([name for name, _ in MODEL_SPECS], 2):
        if frozenset((a, b)) in planned_unordered:
            continue
        add_row(a, b, False, f"{a} vs {b}")

    out = pd.DataFrame(rows)
    return out


def p_label(p: float) -> str:
    if not np.isfinite(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def evaluate_variantwise(df: pd.DataFrame, n_splits: int = 10, seed: int = 16) -> pd.DataFrame:
    y = df["label"].to_numpy(dtype=int)
    splits = list(StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(df, y))
    pair_candidates = pair_weights(0.01)
    triple_candidates = triple_weights(0.05)
    rows = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        y_train = y[train_idx]
        y_test = y[test_idx]
        for model_name, components in MODEL_SPECS:
            if len(components) == 1:
                weights = (1.0,)
            else:
                candidates = pair_candidates if len(components) == 2 else triple_candidates
                train_scores = [df.iloc[train_idx][SCORE_COLUMNS[c]].to_numpy(float) for c in components]
                weights = best_weights(y_train, train_scores, candidates)
            test_scores = [df.iloc[test_idx][SCORE_COLUMNS[c]].to_numpy(float) for c in components]
            pred = sum(w * s for w, s in zip(weights, test_scores))
            weight_by_component = dict(zip(components, weights))
            rows.append(
                {
                    "fold": fold,
                    "model": model_name,
                    "test_auc": roc_auc_score(y_test, pred),
                    "w_esm2_150m": weight_by_component.get("ESM-2 150M", 0.0),
                    "w_esm2_650m": weight_by_component.get("ESM-2 650M", 0.0),
                    "w_esm1b_650m": weight_by_component.get("ESM-1b 650M", 0.0),
                    "w_calm": weight_by_component.get("CaLM", 0.0),
                    "n_test": len(test_idx),
                    "n_test_genes": df.iloc[test_idx]["Gene_prot"].nunique(),
                    "n_test_pathogenic": int(y_test.sum()),
                }
            )
    return pd.DataFrame(rows)


def plot_fold_delta_distribution(gene_folds: pd.DataFrame) -> None:
    wide = gene_folds.pivot(index="fold", columns="model", values="test_auc")
    rows = []
    full_labels = {
        "ESM-2 150M + ESM-2 650M vs ESM-2 650M": "ESM-2 (150M) + ESM-2 (650M) vs ESM-2 (650M)",
        "ESM-2 650M + ESM-1b 650M vs ESM-2 650M": "ESM-2 (650M) + ESM-1b (650M) vs ESM-2 (650M)",
        "ESM-2 650M + CaLM vs ESM-2 650M": "ESM-2 (650M) + CaLM vs ESM-2 (650M)",
        "ESM-2 650M + CaLM vs same-modality placeholder": (
            "ESM-2 (650M) + CaLM vs ESM-2 (150M) + ESM-2 (650M)"
        ),
        "Strong same-modality placeholder vs ESM-1b 650M": (
            "ESM-1b (650M) + ESM-2 (650M) vs ESM-1b (650M)"
        ),
        "ESM-1b 650M + CaLM vs ESM-1b 650M": "ESM-1b (650M) + CaLM vs ESM-1b (650M)",
        "Triple model vs ESM-2 650M + ESM-1b 650M": (
            "ESM-2 (650M) + ESM-1b (650M) + CaLM vs "
            "ESM-2 (650M) + ESM-1b (650M)"
        ),
    }
    for a, b, label in PLANNED_CONTRASTS:
        for fold, delta in (wide[a] - wide[b]).items():
            rows.append({"contrast": label, "fold": fold, "delta": delta})
    df = pd.DataFrame(rows)
    order = [label for _, _, label in PLANNED_CONTRASTS]
    fig, ax = plt.subplots(figsize=(6.6, 3.35))
    rng = np.random.default_rng(3)
    colors = [
        COLORS["grey"],
        COLORS["lilac"],
        COLORS["green"],
        COLORS["teal"],
        COLORS["yellow"],
        COLORS["rose"],
        COLORS["blue"],
    ]
    for idx, contrast in enumerate(order):
        vals = df[df["contrast"] == contrast]["delta"].to_numpy(float)
        y = np.full_like(vals, idx, dtype=float) + rng.normal(0, 0.055, size=len(vals))
        ax.scatter(vals, y, s=20, color=colors[idx], edgecolor=COLORS["edge"], linewidth=0.45, zorder=3)
        mean = vals.mean()
        sd = vals.std(ddof=1)
        ax.errorbar(mean, idx, xerr=sd, fmt="o", color=COLORS["text"], markersize=3.2, capsize=4, linewidth=0.9, zorder=4)
    ax.axvline(0, color="#9B9B96", linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([full_labels[s] for s in order])
    ax.invert_yaxis()
    ax.set_xlabel("Fold-level ΔAUROC")
    format_ax(ax)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figS1.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_score_correlation_audit(score_df: pd.DataFrame) -> pd.DataFrame:
    esm2_baseline = SCORE_COLUMNS["ESM-2 650M"]
    esm1b_baseline = SCORE_COLUMNS["ESM-1b 650M"]
    rows = []
    for baseline_name, baseline_col, scorer, col in [
        ("ESM-2 650M", esm2_baseline, "ESM-2 150M", SCORE_COLUMNS["ESM-2 150M"]),
        ("ESM-2 650M", esm2_baseline, "ESM-1b 650M", SCORE_COLUMNS["ESM-1b 650M"]),
        ("ESM-2 650M", esm2_baseline, "CaLM", SCORE_COLUMNS["CaLM"]),
        ("ESM-1b 650M", esm1b_baseline, "ESM-2 650M", SCORE_COLUMNS["ESM-2 650M"]),
    ]:
        rows.append(
            {
                "baseline_scorer": baseline_name,
                "added_scorer": scorer,
                "pearson_r": float(score_df[baseline_col].corr(score_df[col], method="pearson")),
                "spearman_rho": float(score_df[baseline_col].corr(score_df[col], method="spearman")),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "model_control_added_scorer_correlation_audit.csv", index=False)

    focused = (
        out[out["baseline_scorer"].eq("ESM-2 650M")]
        .assign(
            scorer_modality=lambda x: np.where(
                x["added_scorer"].eq("CaLM"), "codon-level LM", "protein LM"
            ),
            control_role=lambda x: x["added_scorer"].map(
                {
                    "ESM-2 150M": "weak PLM placeholder",
                    "ESM-1b 650M": "strong PLM placeholder",
                    "CaLM": "cross-modal codon scorer",
                }
            ),
            n_variants=int(len(score_df)),
        )
        .sort_values("spearman_rho", ascending=False)
        .reset_index(drop=True)
    )
    focused.insert(0, "redundancy_rank", np.arange(1, len(focused) + 1))
    focused.to_csv(OUT_DIR / "supp_table_model_control_scorer_correlations.csv", index=False)
    return out


def same_modality_placeholder_control(gene_folds: pd.DataFrame) -> pd.DataFrame:
    wide = gene_folds.pivot(index="fold", columns="model", values="test_auc")
    baseline = "ESM-2 650M"
    weak_placeholder = "ESM-2 150M + ESM-2 650M"
    strong_placeholder_baseline = "ESM-1b 650M"
    strong_placeholder = "ESM-2 650M + ESM-1b 650M"
    calm = "ESM-2 650M + CaLM"
    out = pd.DataFrame(
        {
            "fold": wide.index,
            "esm2_650m_auc": wide[baseline].to_numpy(float),
            "weak_same_modality_placeholder_auc": wide[weak_placeholder].to_numpy(float),
            "esm1b_650m_auc": wide[strong_placeholder_baseline].to_numpy(float),
            "strong_same_modality_placeholder_auc": wide[strong_placeholder].to_numpy(float),
            "esm2_650m_plus_calm_auc": wide[calm].to_numpy(float),
        }
    )
    out["weak_same_modality_placeholder_gain"] = (
        out["weak_same_modality_placeholder_auc"] - out["esm2_650m_auc"]
    )
    out["strong_same_modality_placeholder_gain"] = (
        out["strong_same_modality_placeholder_auc"] - out["esm1b_650m_auc"]
    )
    out["calm_gain"] = out["esm2_650m_plus_calm_auc"] - out["esm2_650m_auc"]
    out["calm_minus_weak_placeholder_gain"] = (
        out["calm_gain"] - out["weak_same_modality_placeholder_gain"]
    )
    out["calm_minus_strong_placeholder_gain"] = (
        out["calm_gain"] - out["strong_same_modality_placeholder_gain"]
    )
    out.to_csv(OUT_DIR / "model_control_same_modality_placeholder_fold_deltas.csv", index=False)

    try:
        wilcoxon_weak_p = float(wilcoxon(out["calm_minus_weak_placeholder_gain"]).pvalue)
    except ValueError:
        wilcoxon_weak_p = np.nan
    try:
        wilcoxon_strong_p = float(wilcoxon(out["calm_minus_strong_placeholder_gain"]).pvalue)
    except ValueError:
        wilcoxon_strong_p = np.nan
    summary = pd.DataFrame(
        [
            {
                "n_folds": int(len(out)),
                "mean_esm2_650m_auc": float(out["esm2_650m_auc"].mean()),
                "mean_weak_same_modality_placeholder_auc": float(
                    out["weak_same_modality_placeholder_auc"].mean()
                ),
                "mean_esm1b_650m_auc": float(out["esm1b_650m_auc"].mean()),
                "mean_strong_same_modality_placeholder_auc": float(
                    out["strong_same_modality_placeholder_auc"].mean()
                ),
                "mean_esm2_650m_plus_calm_auc": float(
                    out["esm2_650m_plus_calm_auc"].mean()
                ),
                "mean_weak_same_modality_placeholder_gain": float(
                    out["weak_same_modality_placeholder_gain"].mean()
                ),
                "sd_weak_same_modality_placeholder_gain": float(
                    out["weak_same_modality_placeholder_gain"].std(ddof=1)
                ),
                "mean_strong_same_modality_placeholder_gain": float(
                    out["strong_same_modality_placeholder_gain"].mean()
                ),
                "sd_strong_same_modality_placeholder_gain": float(
                    out["strong_same_modality_placeholder_gain"].std(ddof=1)
                ),
                "mean_calm_gain": float(out["calm_gain"].mean()),
                "sd_calm_gain": float(out["calm_gain"].std(ddof=1)),
                "mean_calm_minus_weak_placeholder_gain": float(
                    out["calm_minus_weak_placeholder_gain"].mean()
                ),
                "sd_calm_minus_weak_placeholder_gain": float(
                    out["calm_minus_weak_placeholder_gain"].std(ddof=1)
                ),
                "mean_calm_minus_strong_placeholder_gain": float(
                    out["calm_minus_strong_placeholder_gain"].mean()
                ),
                "sd_calm_minus_strong_placeholder_gain": float(
                    out["calm_minus_strong_placeholder_gain"].std(ddof=1)
                ),
                "paired_t_p_calm_gain_vs_weak_placeholder_gain": float(
                    ttest_rel(out["calm_gain"], out["weak_same_modality_placeholder_gain"]).pvalue
                ),
                "wilcoxon_p_calm_gain_vs_weak_placeholder_gain": wilcoxon_weak_p,
                "paired_t_p_calm_gain_vs_strong_placeholder_gain": float(
                    ttest_rel(out["calm_gain"], out["strong_same_modality_placeholder_gain"]).pvalue
                ),
                "wilcoxon_p_calm_gain_vs_strong_placeholder_gain": wilcoxon_strong_p,
                "fraction_folds_calm_gain_exceeds_weak_placeholder_gain": float(
                    (out["calm_gain"] > out["weak_same_modality_placeholder_gain"]).mean()
                ),
                "fraction_folds_calm_gain_exceeds_strong_placeholder_gain": float(
                    (out["calm_gain"] > out["strong_same_modality_placeholder_gain"]).mean()
                ),
            }
        ]
    )
    summary.to_csv(OUT_DIR / "model_control_same_modality_placeholder_summary.csv", index=False)
    return out


def plot_same_modality_placeholder_control(placeholder_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.7), gridspec_kw={"width_ratios": [1.2, 1.15]})
    ax = axes[0]
    x = np.array([0, 1, 2])
    rng = np.random.default_rng(11)
    for _, row in placeholder_df.iterrows():
        vals = [
            row["weak_same_modality_placeholder_gain"],
            row["strong_same_modality_placeholder_gain"],
            row["calm_gain"],
        ]
        jitter = rng.normal(0, 0.018, size=3)
        ax.plot(x + jitter, vals, color="#B8B8B2", linewidth=0.7, zorder=1)
        ax.scatter(
            x + jitter,
            vals,
            s=22,
            color=[COLORS["grey"], COLORS["yellow"], COLORS["green"]],
            edgecolor=COLORS["edge"],
            linewidth=0.45,
            zorder=3,
        )
    means = [
        placeholder_df["weak_same_modality_placeholder_gain"].mean(),
        placeholder_df["strong_same_modality_placeholder_gain"].mean(),
        placeholder_df["calm_gain"].mean(),
    ]
    ax.scatter(x, means, s=34, color=COLORS["text"], zorder=4)
    ax.axhline(0, color="#9B9B96", linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_xticks(x)
    ax.set_xticklabels(["Weak PLM\nplaceholder", "Strong PLM\nplaceholder", "CaLM"])
    ax.set_ylabel("Fold-level gain over matched PLM baseline")
    ax.set_xlim(-0.38, 2.38)
    ax.text(-0.13, 1.04, "A", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    format_ax(ax)

    ax = axes[1]
    diff_specs = [
        ("CaLM - weak", "calm_minus_weak_placeholder_gain", COLORS["teal"], 0.13),
        ("CaLM - strong", "calm_minus_strong_placeholder_gain", COLORS["yellow"], -0.13),
    ]
    for label, col, color, y_center in diff_specs:
        vals = placeholder_df[col].to_numpy(float)
        y = np.full_like(vals, y_center, dtype=float) + rng.normal(0, 0.024, size=len(vals))
        ax.scatter(vals, y, s=24, color=color, edgecolor=COLORS["edge"], linewidth=0.45, zorder=3, label=label)
        mean = vals.mean()
        sd = vals.std(ddof=1)
        ax.errorbar(mean, y_center, xerr=sd, fmt="o", color=COLORS["text"], markersize=3.8, capsize=4, linewidth=0.9, zorder=4)
    ax.axvline(0, color="#9B9B96", linewidth=0.8, linestyle=(0, (3, 2)))
    ax.set_yticks([0.13, -0.13])
    ax.set_yticklabels(["CaLM - weak", "CaLM - strong"])
    ax.set_xlabel("Gain difference (ΔAUROC)")
    ax.set_ylim(-0.36, 0.36)
    ax.text(-0.11, 1.04, "B", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    format_ax(ax)

    fig.tight_layout(w_pad=1.35)
    fig.savefig(FIG_DIR / "figS3.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_cv_comparison(gene_folds: pd.DataFrame, variant_folds: pd.DataFrame) -> pd.DataFrame:
    gene = gene_folds.groupby("model")["test_auc"].agg(gene_wise_mean_auc="mean", gene_wise_sd_auc="std").reset_index()
    var = variant_folds.groupby("model")["test_auc"].agg(variant_wise_mean_auc="mean", variant_wise_sd_auc="std").reset_index()
    out = gene.merge(var, on="model")
    out["variant_minus_gene_mean_auc"] = out["variant_wise_mean_auc"] - out["gene_wise_mean_auc"]
    out["model"] = pd.Categorical(out["model"], categories=[m for m, _ in MODEL_SPECS], ordered=True)
    out = out.sort_values("model")
    out.to_csv(OUT_DIR / "supp_table_model_control_gene_vs_variantwise_cv.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.0, 3.35))
    x = np.arange(len(out))
    offset = 0.14
    ax.errorbar(
        x - offset,
        out["gene_wise_mean_auc"],
        yerr=out["gene_wise_sd_auc"],
        fmt="o",
        markersize=5.0,
        color=COLORS["blue"],
        markeredgecolor=COLORS["edge"],
        markeredgewidth=0.75,
        ecolor=COLORS["edge"],
        elinewidth=0.85,
        capsize=2.8,
        label="Gene-wise CV",
        zorder=3,
    )
    ax.errorbar(
        x + offset,
        out["variant_wise_mean_auc"],
        yerr=out["variant_wise_sd_auc"],
        fmt="o",
        markersize=5.0,
        color=COLORS["yellow"],
        markeredgecolor=COLORS["edge"],
        markeredgewidth=0.75,
        ecolor=COLORS["edge"],
        elinewidth=0.85,
        capsize=2.8,
        label="Variant-wise CV",
        zorder=3,
    )
    for idx, row in out.reset_index(drop=True).iterrows():
        ax.plot(
            [idx - offset, idx + offset],
            [row["gene_wise_mean_auc"], row["variant_wise_mean_auc"]],
            color="#B9B9B4",
            linewidth=0.55,
            zorder=2,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_LABELS[str(m)] for m in out["model"]], rotation=35, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.78, 0.965)
    ax.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.31),
        fontsize=7.2,
        ncol=2,
        handletextpad=0.5,
        columnspacing=1.4,
    )
    format_ax(ax)
    fig.subplots_adjust(bottom=0.205, top=0.98, left=0.075, right=0.995)
    fig.savefig(FIG_DIR / "figS2.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_weight_distribution(gene_folds: pd.DataFrame) -> None:
    models = [
        "ESM-2 150M + ESM-2 650M",
        "ESM-2 650M + ESM-1b 650M",
        "ESM-2 150M + CaLM",
        "ESM-2 650M + CaLM",
        "ESM-1b 650M + CaLM",
        "ESM-2 650M + ESM-1b 650M + CaLM",
    ]
    component_cols = [
        ("w_esm2_150m", "ESM-2 (150M)", COLORS["grey"]),
        ("w_esm2_650m", "ESM-2 (650M)", COLORS["blue"]),
        ("w_esm1b_650m", "ESM-1b (650M)", COLORS["lilac"]),
        ("w_calm", "CaLM", COLORS["green"]),
    ]
    sub = gene_folds[gene_folds["model"].isin(models)].copy()
    summary = (
        sub.groupby("model")[[col for col, _, _ in component_cols]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary.columns = ["_".join([str(c) for c in col if c]) for col in summary.columns.to_flat_index()]
    summary.to_csv(OUT_DIR / "supp_table_model_control_fold_weight_summary.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.2), sharey=True)
    rng = np.random.default_rng(4)
    for ax, model in zip(axes.ravel(), models):
        m = sub[sub["model"] == model]
        xs = np.arange(len(component_cols))
        for idx, (col, label, color) in enumerate(component_cols):
            vals = m[col].to_numpy(float)
            jitter = rng.normal(0, 0.045, len(vals))
            ax.scatter(np.full(len(vals), idx) + jitter, vals, s=16, color=color, edgecolor=COLORS["edge"], linewidth=0.4, zorder=3)
            if vals.max() > 0:
                ax.plot([idx - 0.18, idx + 0.18], [vals.mean(), vals.mean()], color=COLORS["text"], linewidth=1.1, zorder=4)
        ax.set_title(MODEL_LABELS[model].replace("\n", " "), fontsize=7.1, pad=3)
        ax.set_xticks(xs)
        ax.set_xticklabels([label for _, label, _ in component_cols], rotation=38, ha="right", fontsize=6.3)
        ax.set_ylim(-0.03, 1.03)
        format_ax(ax)
    axes[0, 0].set_ylabel("Optimized weight")
    axes[1, 0].set_ylabel("Optimized weight")
    fig.tight_layout(w_pad=0.65, h_pad=0.85)
    fig.savefig(FIG_DIR / "figS4.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    setup()

    gene_folds = pd.read_csv(FOLD_RESULTS)
    score_df = pd.read_csv(SCORE_TABLE)

    all_pairwise = all_pairwise_tests(gene_folds)
    all_pairwise.to_csv(OUT_DIR / "supp_table_model_control_all_pairwise_tests.csv", index=False)
    score_correlation = write_score_correlation_audit(score_df)
    placeholder_df = same_modality_placeholder_control(gene_folds)

    variant_fold_path = OUT_DIR / "model_control_variantwise_fold_results.csv"
    if variant_fold_path.exists():
        variant_folds = pd.read_csv(variant_fold_path)
    else:
        variant_folds = evaluate_variantwise(score_df)
        variant_folds.to_csv(variant_fold_path, index=False)

    cv_comparison = plot_cv_comparison(gene_folds, variant_folds)
    plot_fold_delta_distribution(gene_folds)
    plot_same_modality_placeholder_control(placeholder_df)
    plot_weight_distribution(gene_folds)

    print(OUT_DIR / "supp_table_model_control_all_pairwise_tests.csv")
    print(OUT_DIR / "model_control_added_scorer_correlation_audit.csv")
    print(OUT_DIR / "supp_table_model_control_scorer_correlations.csv")
    print(OUT_DIR / "model_control_same_modality_placeholder_fold_deltas.csv")
    print(OUT_DIR / "model_control_same_modality_placeholder_summary.csv")
    print(OUT_DIR / "supp_table_model_control_gene_vs_variantwise_cv.csv")
    print(OUT_DIR / "supp_table_model_control_fold_weight_summary.csv")
    print(FIG_DIR / "figS1.png")
    print(FIG_DIR / "figS2.png")
    print(FIG_DIR / "figS4.png")
    print(FIG_DIR / "figS3.png")
    print(score_correlation.to_string(index=False))
    print(cv_comparison[["model", "gene_wise_mean_auc", "variant_wise_mean_auc", "variant_minus_gene_mean_auc"]].to_string(index=False))


if __name__ == "__main__":
    main()
