#!/usr/bin/env python3
"""Gene-level codon contribution summaries for the ClinVar cohort."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm


IN_DIR = Path("Results/Revision/len1022_core_clinvar")
OUT_DIR = IN_DIR / "gene_level_dosage"


def fit_regression(df: pd.DataFrame, x_col: str, y_col: str, label: str, weighted: bool) -> dict[str, object]:
    reg_df = df.dropna(subset=[x_col, y_col, "n_pathogenic", "n_benign"]).copy()
    x = sm.add_constant(reg_df[x_col].to_numpy(float), has_constant="add")
    y = reg_df[y_col].to_numpy(float)
    if weighted:
        weights = (
            reg_df["n_pathogenic"].to_numpy(float)
            * reg_df["n_benign"].to_numpy(float)
            / (reg_df["n_pathogenic"].to_numpy(float) + reg_df["n_benign"].to_numpy(float))
        )
        model = sm.WLS(y, x, weights=weights).fit(cov_type="HC3")
        method = "WLS_HC3"
        weight_definition = "n_pathogenic*n_benign/(n_pathogenic+n_benign)"
    else:
        model = sm.OLS(y, x).fit(cov_type="HC3")
        method = "OLS_HC3"
        weight_definition = ""
    ci_low, ci_high = model.conf_int(alpha=0.05)[1]
    return {
        "analysis": label,
        "method": method,
        "weighted": int(weighted),
        "n_genes": int(len(reg_df)),
        "x_column": x_col,
        "y_column": y_col,
        "weight_definition": weight_definition,
        "intercept": float(model.params[0]),
        "slope": float(model.params[1]),
        "slope_95ci_low": float(ci_low),
        "slope_95ci_high": float(ci_high),
        "p_value": float(model.pvalues[1]),
        "r_squared": float(model.rsquared),
    }


def fit_spearman_and_top5_sensitivity(df: pd.DataFrame, x_col: str, y_col: str, label: str) -> dict[str, object]:
    reg_df = df.dropna(subset=[x_col, y_col, "n_pathogenic", "n_benign"]).copy()
    rho, rho_p = spearmanr(reg_df[x_col], reg_df[y_col])
    cutoff = float(reg_df[x_col].quantile(0.95))
    trimmed = reg_df[reg_df[x_col] <= cutoff].copy()
    x = sm.add_constant(trimmed[x_col].to_numpy(float), has_constant="add")
    y = trimmed[y_col].to_numpy(float)
    model = sm.OLS(y, x).fit(cov_type="HC3")
    ci_low, ci_high = model.conf_int(alpha=0.05)[1]
    return {
        "analysis": label,
        "n_genes": int(len(reg_df)),
        "x_column": x_col,
        "y_column": y_col,
        "spearman_rho": float(rho),
        "spearman_p_value": float(rho_p),
        "top5_percent_cutoff": cutoff,
        "n_excluded_top5_percent": int(len(reg_df) - len(trimmed)),
        "trimmed_n_genes": int(len(trimmed)),
        "trimmed_ols_slope": float(model.params[1]),
        "trimmed_ols_slope_95ci_low": float(ci_low),
        "trimmed_ols_slope_95ci_high": float(ci_high),
        "trimmed_ols_p_value": float(model.pvalues[1]),
        "trimmed_ols_r_squared": float(model.rsquared),
    }


def add_rankings(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["esm2_calm_weight"] = out["ESM-2 650M + CaLM weight_second"]
    out["esm2_cross_modal_gain"] = out["650MCaLM_minus_150M650M"]
    out["esm2_gain_over_protein"] = out["650MCaLM_minus_650M"]
    out["esm2_independent_calm_signal"] = out["CaLM AUROC"] - out["ESM-2 650M AUROC"]

    out["esm1b_calm_weight"] = out["ESM-1b 650M + CaLM weight_second"]
    out["esm1b_gain_over_protein"] = out["1bCaLM_minus_1b"]
    out["esm1b_independent_calm_signal"] = out["CaLM AUROC"] - out["ESM-1b 650M AUROC"]

    out["esm2_codon_rank"] = out["esm2_calm_weight"].rank(method="first", ascending=False).astype(int)
    out["esm2_gain_rank"] = out["esm2_cross_modal_gain"].rank(method="first", ascending=False).astype(int)
    out["esm1b_codon_rank"] = out["esm1b_calm_weight"].rank(method="first", ascending=False).astype(int)
    out["esm1b_gain_rank"] = out["esm1b_gain_over_protein"].rank(method="first", ascending=False).astype(int)

    n_top = max(1, int(np.ceil(len(out) * 0.10)))
    out["esm2_top_codon_weight_decile"] = out["esm2_codon_rank"] <= n_top
    out["esm2_top_cross_modal_gain_decile"] = out["esm2_gain_rank"] <= n_top
    out["esm1b_top_codon_weight_decile"] = out["esm1b_codon_rank"] <= n_top
    out["esm1b_top_gain_decile"] = out["esm1b_gain_rank"] <= n_top
    return out


def write_gene_list(df: pd.DataFrame, mask_col: str, metric_col: str, filename: str) -> None:
    cols = [
        "gene",
        "n",
        "n_pathogenic",
        "n_benign",
        "ESM-2 650M AUROC",
        "ESM-1b 650M AUROC",
        "CaLM AUROC",
        "esm2_calm_weight",
        "esm2_cross_modal_gain",
        "esm2_independent_calm_signal",
        "esm1b_calm_weight",
        "esm1b_gain_over_protein",
        "esm1b_independent_calm_signal",
    ]
    df.loc[df[mask_col]].sort_values(metric_col, ascending=False)[cols].to_csv(OUT_DIR / filename, index=False)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    per_gene = pd.read_csv(IN_DIR / "per_gene_ensemble_summary.csv")
    ranked = add_rankings(per_gene)
    ranked.to_csv(OUT_DIR / "gene_level_codon_contribution_summary.csv", index=False)

    regressions = pd.DataFrame(
        [
            fit_regression(ranked, "esm2_calm_weight", "esm2_cross_modal_gain", "esm2_weight_predicts_cross_modal_gain", False),
            fit_regression(ranked, "esm2_calm_weight", "esm2_cross_modal_gain", "weighted_esm2_weight_predicts_cross_modal_gain", True),
            fit_regression(
                ranked,
                "esm2_independent_calm_signal",
                "esm2_cross_modal_gain",
                "esm2_independent_calm_signal_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                ranked,
                "esm2_independent_calm_signal",
                "esm2_cross_modal_gain",
                "weighted_esm2_independent_calm_signal_predicts_cross_modal_gain",
                True,
            ),
            fit_regression(ranked, "esm1b_calm_weight", "esm1b_gain_over_protein", "esm1b_weight_predicts_gain_over_esm1b", False),
            fit_regression(
                ranked,
                "esm1b_calm_weight",
                "esm1b_gain_over_protein",
                "weighted_esm1b_weight_predicts_gain_over_esm1b",
                True,
            ),
            fit_regression(
                ranked,
                "esm1b_independent_calm_signal",
                "esm1b_gain_over_protein",
                "esm1b_independent_calm_signal_predicts_gain_over_esm1b",
                False,
            ),
            fit_regression(
                ranked,
                "esm1b_independent_calm_signal",
                "esm1b_gain_over_protein",
                "weighted_esm1b_independent_calm_signal_predicts_gain_over_esm1b",
                True,
            ),
        ]
    )
    regressions.to_csv(OUT_DIR / "gene_level_codon_contribution_regressions.csv", index=False)

    sensitivity = pd.DataFrame(
        [
            fit_spearman_and_top5_sensitivity(
                ranked,
                "esm2_independent_calm_signal",
                "esm2_cross_modal_gain",
                "esm2_independent_calm_signal_spearman_and_top5_sensitivity",
            ),
            fit_spearman_and_top5_sensitivity(
                ranked,
                "esm1b_independent_calm_signal",
                "esm1b_gain_over_protein",
                "esm1b_independent_calm_signal_spearman_and_top5_sensitivity",
            ),
        ]
    )
    sensitivity.to_csv(OUT_DIR / "gene_level_independent_signal_sensitivity.csv", index=False)

    write_gene_list(ranked, "esm2_top_codon_weight_decile", "esm2_calm_weight", "esm2_top_codon_weight_decile_genes.csv")
    write_gene_list(ranked, "esm2_top_cross_modal_gain_decile", "esm2_cross_modal_gain", "esm2_top_cross_modal_gain_decile_genes.csv")
    write_gene_list(ranked, "esm1b_top_codon_weight_decile", "esm1b_calm_weight", "esm1b_top_codon_weight_decile_genes.csv")
    write_gene_list(ranked, "esm1b_top_gain_decile", "esm1b_gain_over_protein", "esm1b_top_gain_decile_genes.csv")

    summary = pd.DataFrame(
        [
            {
                "analysis": "gene_level_codon_contribution",
                "n_genes": int(len(ranked)),
                "n_top_decile": int(np.ceil(len(ranked) * 0.10)),
                "mean_esm2_calm_weight": float(ranked["esm2_calm_weight"].mean()),
                "median_esm2_calm_weight": float(ranked["esm2_calm_weight"].median()),
                "mean_esm2_cross_modal_gain": float(ranked["esm2_cross_modal_gain"].mean()),
                "median_esm2_cross_modal_gain": float(ranked["esm2_cross_modal_gain"].median()),
                "mean_esm1b_calm_weight": float(ranked["esm1b_calm_weight"].mean()),
                "median_esm1b_calm_weight": float(ranked["esm1b_calm_weight"].median()),
                "mean_esm1b_gain_over_protein": float(ranked["esm1b_gain_over_protein"].mean()),
                "median_esm1b_gain_over_protein": float(ranked["esm1b_gain_over_protein"].median()),
            }
        ]
    )
    summary.to_csv(OUT_DIR / "gene_level_codon_contribution_run_summary.csv", index=False)

    print(summary.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nRegressions")
    print(regressions.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nIndependent signal sensitivity")
    print(sensitivity.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nTop ESM-2 codon-weight genes")
    print(
        ranked.sort_values("esm2_calm_weight", ascending=False)[
            ["gene", "n", "n_pathogenic", "n_benign", "esm2_calm_weight", "esm2_cross_modal_gain", "esm2_independent_calm_signal"]
        ]
        .head(15)
        .to_string(index=False)
    )
    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
