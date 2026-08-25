#!/usr/bin/env python3
"""Re-optimise frozen ClinMAVE functional-class analyses with Bayesian search."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ensemble_optim import bayes_optimize_weights


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--plm-column", required=True)
    parser.add_argument("--plm-name", required=True)
    parser.add_argument("--old-folds", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument("--bayes-init", type=int, default=10)
    parser.add_argument("--bayes-iter", type=int, default=20)
    return parser.parse_args()


def standardize(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = train.mean(axis=0)
    sd = train.std(axis=0)
    sd[sd == 0] = 1.0
    return (train - mean) / sd, (test - mean) / sd


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    return np.nan if len(np.unique(y)) < 2 else float(roc_auc_score(y, score))


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input)
    rows = []

    for comparison, df in data.groupby("comparison", sort=False):
        df = df.dropna(subset=["label", "Gene", args.plm_column, "calm_score"]).copy()
        y = df["label"].to_numpy(int)
        groups = df["Gene"].astype(str).to_numpy()
        splits = list(
            StratifiedGroupKFold(
                n_splits=args.folds,
                shuffle=True,
                random_state=args.seed,
            ).split(df, y, groups=groups)
        )
        for fold, (train_idx, test_idx) in enumerate(splits, start=1):
            train_raw = df.iloc[train_idx][[args.plm_column, "calm_score"]].to_numpy(float)
            test_raw = df.iloc[test_idx][[args.plm_column, "calm_score"]].to_numpy(float)
            train, test = standardize(train_raw, test_raw)
            weights, train_auc, n_evals = bayes_optimize_weights(
                y[train_idx],
                [train[:, 0], train[:, 1]],
                seed=args.seed + fold + 1000 * len(rows),
                n_init=args.bayes_init,
                n_iter=args.bayes_iter,
                candidate_step=0.001,
            )
            plm_auc = safe_auc(y[test_idx], test[:, 0])
            calm_auc = safe_auc(y[test_idx], test[:, 1])
            combo_auc = safe_auc(y[test_idx], weights[0] * test[:, 0] + weights[1] * test[:, 1])
            assay, case_class, _, _ = comparison.split("_")
            rows.append(
                {
                    "assay": assay,
                    "case_class": case_class,
                    "comparison": comparison,
                    "fold": fold,
                    "optimizer": "bayes",
                    "calm_weight": weights[1],
                    "plm_weight": weights[0],
                    "train_auc_at_selected_weights": train_auc,
                    "n_objective_evals": n_evals,
                    "auroc_plm": plm_auc,
                    "auroc_calm": calm_auc,
                    "auroc_combo": combo_auc,
                    "delta_combo_vs_plm": combo_auc - plm_auc,
                    "n_test": len(test_idx),
                    "n_test_genes": len(np.unique(groups[test_idx])),
                }
            )

    folds = pd.DataFrame(rows)
    folds.to_csv(args.outdir / "bayesian_fold_metrics.csv", index=False)

    summary_rows = []
    for (assay, case_class), sub in folds.groupby(["assay", "case_class"], sort=False):
        valid = sub.dropna(subset=["auroc_combo", "auroc_plm"])
        p_value = ttest_rel(valid["auroc_combo"], valid["auroc_plm"]).pvalue if len(valid) >= 2 else np.nan
        summary_rows.append(
            {
                "plm": args.plm_name,
                "assay": assay,
                "case_class": case_class,
                "mean_calm_weight": sub["calm_weight"].mean(),
                "sd_calm_weight": sub["calm_weight"].std(ddof=1),
                "mean_auroc_plm": sub["auroc_plm"].mean(),
                "mean_auroc_combo": sub["auroc_combo"].mean(),
                "mean_delta_combo_vs_plm": sub["delta_combo_vs_plm"].mean(),
                "n_evaluable_folds": len(valid),
                "paired_t_p": p_value,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.outdir / "bayesian_summary.csv", index=False)

    old = pd.read_csv(args.old_folds)
    # Old files have model-specific names, so only the common weight is required for comparison.
    comparison = old.merge(
        folds,
        on=["assay", "case_class", "fold"],
        suffixes=("_grid", "_bayes"),
    )
    comparison["bayes_minus_grid_weight"] = comparison["calm_weight_bayes"] - comparison["calm_weight_grid"]
    comparison.to_csv(args.outdir / "bayesian_vs_grid_fold_comparison.csv", index=False)

    print(summary.to_string(index=False, float_format=lambda value: f"{value:.6g}"))
    print(f"\nWrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
