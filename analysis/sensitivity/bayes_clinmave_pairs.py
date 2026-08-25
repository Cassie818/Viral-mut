#!/usr/bin/env python3
"""Re-optimise frozen matched ClinMAVE DMS-CBGE pairs with Bayesian search."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ensemble_optim import bayes_optimize_weights


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--plm-llr-column", required=True)
    parser.add_argument("--old-folds", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=16)
    parser.add_argument("--bayes-init", type=int, default=10)
    parser.add_argument("--bayes-iter", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input)
    rows = []
    for pair_index, (pair_id, table) in enumerate(data.groupby("pair_id", sort=False)):
        split_y = table["DMS_label"].astype(str) + "_" + table["CBGE_label"].astype(str)
        n_splits = min(10, int(split_y.value_counts().min()))
        if n_splits < 2:
            continue
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.seed)
        for fold, (train_idx, test_idx) in enumerate(cv.split(table, split_y), start=1):
            train = table.iloc[train_idx]
            test = table.iloc[test_idx]
            for platform in ["DMS", "CBGE"]:
                label_col = f"{platform}_label"
                y_train = train[label_col].to_numpy(int)
                y_test = test[label_col].to_numpy(int)
                plm_train = -train[args.plm_llr_column].to_numpy(float)
                calm_train = -train["calm_llr"].to_numpy(float)
                weights, train_auc, n_evals = bayes_optimize_weights(
                    y_train,
                    [plm_train, calm_train],
                    seed=args.seed + pair_index * 1000 + fold * 10 + int(platform == "CBGE"),
                    n_init=args.bayes_init,
                    n_iter=args.bayes_iter,
                    candidate_step=0.001,
                )
                evaluable = len(np.unique(y_test)) == 2
                if evaluable:
                    plm_test = -test[args.plm_llr_column].to_numpy(float)
                    calm_test = -test["calm_llr"].to_numpy(float)
                    combo = weights[0] * plm_test + weights[1] * calm_test
                    plm_auc = roc_auc_score(y_test, plm_test)
                    calm_auc = roc_auc_score(y_test, calm_test)
                    combo_auc = roc_auc_score(y_test, combo)
                else:
                    plm_auc = calm_auc = combo_auc = np.nan
                rows.append(
                    {
                        "pair_id": pair_id,
                        "Gene": table["Gene"].iloc[0],
                        "DMS_dataset": table["DMS_dataset"].iloc[0],
                        "CBGE_dataset": table["CBGE_dataset"].iloc[0],
                        "fold": fold,
                        "platform": platform,
                        "optimizer": "bayes",
                        "calm_weight": weights[1],
                        "plm_weight": weights[0],
                        "train_auc_at_selected_weights": train_auc,
                        "n_objective_evals": n_evals,
                        "auroc_plm": plm_auc,
                        "auroc_calm": calm_auc,
                        "auroc_combo": combo_auc,
                        "n_test": len(test),
                    }
                )

    folds = pd.DataFrame(rows)
    folds.to_csv(args.outdir / "bayesian_pair_fold_metrics.csv", index=False)
    summary = (
        folds.groupby(["pair_id", "Gene", "DMS_dataset", "CBGE_dataset", "platform"], as_index=False)
        .agg(
            mean_calm_weight=("calm_weight", "mean"),
            sd_calm_weight=("calm_weight", "std"),
            mean_auroc_plm=("auroc_plm", "mean"),
            mean_auroc_combo=("auroc_combo", "mean"),
            n_folds=("fold", "nunique"),
        )
    )
    summary.to_csv(args.outdir / "bayesian_pair_summary.csv", index=False)
    wide = summary.pivot(index=["pair_id", "Gene"], columns="platform", values="mean_calm_weight").reset_index()
    wide["delta_CBGE_minus_DMS"] = wide["CBGE"] - wide["DMS"]
    wide.to_csv(args.outdir / "bayesian_pair_weight_shifts.csv", index=False)

    old = pd.read_csv(args.old_folds)
    comparison = old.merge(
        folds,
        on=["pair_id", "Gene", "fold", "platform"],
        suffixes=("_grid", "_bayes"),
    )
    comparison["bayes_minus_grid_weight"] = comparison["calm_weight_bayes"] - comparison["calm_weight_grid"]
    comparison.to_csv(args.outdir / "bayesian_vs_grid_pair_folds.csv", index=False)

    print(
        pd.DataFrame(
            [{
                "n_pairs": wide["pair_id"].nunique(),
                "positive_shifts": int((wide["delta_CBGE_minus_DMS"] > 0).sum()),
                "median_shift": wide["delta_CBGE_minus_DMS"].median(),
                "mean_abs_fold_weight_change_vs_grid": comparison["bayes_minus_grid_weight"].abs().mean(),
            }]
        ).to_string(index=False, float_format=lambda value: f"{value:.6g}")
    )
    print(f"\nWrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
