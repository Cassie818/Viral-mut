#!/usr/bin/env python3
"""Re-optimise ClinVar gene-level ensembles with Bayesian search."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ensemble_optim import bayes_optimize_weights


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "Results/Revision/len1022_model_control/"
            "len1022_model_control_score_table_complete_cases.csv"
        ),
    )
    parser.add_argument(
        "--grid-summary",
        type=Path,
        default=Path("Results/Revision/len1022_core_clinvar/per_gene_ensemble_summary.csv"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("Results/ClinVar/gene_level_bayesian"),
    )
    parser.add_argument("--seed", type=int, default=16)
    parser.add_argument("--bayes-init", type=int, default=10)
    parser.add_argument("--bayes-iter", type=int, default=20)
    return parser.parse_args()


def optimise(y: np.ndarray, first: np.ndarray, second: np.ndarray, seed: int, args: argparse.Namespace):
    weights, auc, n_evals = bayes_optimize_weights(
        y,
        [first, second],
        seed=seed,
        n_init=args.bayes_init,
        n_iter=args.bayes_iter,
        candidate_step=0.001,
    )
    return float(weights[1]), float(auc), int(n_evals)


def regression(df: pd.DataFrame, x_col: str, y_col: str, label: str) -> dict[str, float | str | int]:
    clean = df.dropna(subset=[x_col, y_col]).copy()
    model = sm.OLS(
        clean[y_col].to_numpy(float),
        sm.add_constant(clean[x_col].to_numpy(float)),
    ).fit(cov_type="HC3")
    low, high = model.conf_int(alpha=0.05)[1]
    rho, rho_p = spearmanr(clean[x_col], clean[y_col])
    return {
        "analysis": label,
        "n_genes": len(clean),
        "slope": model.params[1],
        "slope_95ci_low": low,
        "slope_95ci_high": high,
        "ols_p": model.pvalues[1],
        "r_squared": model.rsquared,
        "spearman_rho": rho,
        "spearman_p": rho_p,
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input).drop_duplicates(["variant_id", "Gene_gene"])
    rows = []
    for gene_index, (gene, sub) in enumerate(df.groupby("Gene_gene", sort=True)):
        y = sub["label"].to_numpy(int)
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        if n_pos < 5 or n_neg < 5:
            continue
        s150 = sub["esm2_150m_score"].to_numpy(float)
        s650 = sub["esm2_650m_score"].to_numpy(float)
        s1b = sub["esm1b_650m_score"].to_numpy(float)
        calm = sub["calm_score"].to_numpy(float)
        w_150_650, auc_150_650, n_eval_1 = optimise(y, s150, s650, args.seed + gene_index * 10, args)
        w_650_calm, auc_650_calm, n_eval_2 = optimise(y, s650, calm, args.seed + gene_index * 10 + 1, args)
        w_1b_calm, auc_1b_calm, n_eval_3 = optimise(y, s1b, calm, args.seed + gene_index * 10 + 2, args)
        auc_650 = float(roc_auc_score(y, s650))
        auc_1b = float(roc_auc_score(y, s1b))
        auc_calm = float(roc_auc_score(y, calm))
        rows.append(
            {
                "gene": gene,
                "n": len(sub),
                "n_pathogenic": n_pos,
                "n_benign": n_neg,
                "optimizer": "bayes",
                "esm2_calm_weight": w_650_calm,
                "esm2_calm_auc": auc_650_calm,
                "esm2_150m_650m_weight_second": w_150_650,
                "esm2_150m_650m_auc": auc_150_650,
                "esm1b_calm_weight": w_1b_calm,
                "esm1b_calm_auc": auc_1b_calm,
                "esm2_650m_auc": auc_650,
                "esm1b_650m_auc": auc_1b,
                "calm_auc": auc_calm,
                "esm2_cross_modal_gain": auc_650_calm - auc_650,
                "esm2_cross_modal_advantage": auc_650_calm - auc_150_650,
                "esm2_gain_over_protein": auc_650_calm - auc_650,
                "esm1b_gain_over_protein": auc_1b_calm - auc_1b,
                "esm2_independent_calm_signal": auc_calm - auc_650,
                "esm1b_independent_calm_signal": auc_calm - auc_1b,
                "n_objective_evals": n_eval_1 + n_eval_2 + n_eval_3,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "bayesian_gene_level_summary.csv", index=False)

    regressions = pd.DataFrame(
        [
            regression(result, "esm2_calm_weight", "esm2_cross_modal_gain", "esm2_weight_predicts_cross_modal_gain"),
            regression(result, "esm2_independent_calm_signal", "esm2_cross_modal_gain", "esm2_independent_signal_predicts_cross_modal_gain"),
            regression(result, "esm1b_calm_weight", "esm1b_gain_over_protein", "esm1b_weight_predicts_gain"),
            regression(result, "esm1b_independent_calm_signal", "esm1b_gain_over_protein", "esm1b_independent_signal_predicts_gain"),
        ]
    )
    regressions.to_csv(args.outdir / "bayesian_gene_level_regressions.csv", index=False)

    grid = pd.read_csv(args.grid_summary)
    compare = grid.merge(result, on="gene", suffixes=("_grid", "_bayes"))
    compare["esm2_calm_weight_bayes_minus_grid"] = (
        compare["esm2_calm_weight"] - compare["ESM-2 650M + CaLM weight_second"]
    )
    compare["cross_modal_gain_bayes_minus_grid"] = (
        compare["esm2_cross_modal_gain"] - compare["650MCaLM_minus_650M"]
    )
    compare["cross_modal_advantage_bayes_minus_grid"] = (
        compare["esm2_cross_modal_advantage"] - compare["650MCaLM_minus_150M650M"]
    )
    compare.to_csv(args.outdir / "bayesian_vs_grid_gene_comparison.csv", index=False)

    print(regressions.to_string(index=False, float_format=lambda value: f"{value:.6g}"))
    print("\nBayesian versus grid")
    print(
        compare[[
            "esm2_calm_weight_bayes_minus_grid",
            "cross_modal_gain_bayes_minus_grid",
            "cross_modal_advantage_bayes_minus_grid",
        ]]
        .describe()
        .to_string(float_format=lambda value: f"{value:.6g}")
    )
    print(f"\nWrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
