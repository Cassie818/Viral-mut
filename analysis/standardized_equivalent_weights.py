#!/usr/bin/env python3
"""Convert raw ClinVar ensemble weights to standardized-equivalent weights.

For each cross-validation fold, component standard deviations are estimated on
the training variants only.  If alpha_j is a raw-score ensemble weight and
sigma_j is the corresponding training-fold standard deviation, the equivalent
weight after z-standardization is

    beta_j = alpha_j * sigma_j / sum_l(alpha_l * sigma_l).

The transformed coefficient is used only for reporting. Original ensemble
predictions are retained, so pooled out-of-fold rankings and AUROC are unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


SCORES = Path(
    "Results/ClinVar/current_missense_cohort/"
    "clinvar_missense_len1022_complete_cases_model_scores.csv"
)
OOF = Path("Results/ClinVar/model_control/model_control_oof_predictions.csv.gz")
FOLD_WEIGHTS = Path(
    "Results/ClinVar/model_control/model_control_gene_heldout_fold_results.csv"
)
OUTPUT_CSV = Path(
    "Results/ClinVar/model_control/standardized_equivalent_weights.csv"
)
OUTPUT_FOLD_CSV = Path(
    "Results/ClinVar/model_control/standardized_equivalent_weights_by_fold.csv"
)
OUTPUT_TEX = Path("manuscript/Supplementary_Table_standardized_equivalent_weights.tex")

COMPONENTS = {
    "ESM-2 150M": ("esm2_150m_score", "w_esm2_150m"),
    "ESM-2 650M": ("esm2_650m_score", "w_esm2_650m"),
    "ESM-1b 650M": ("esm1b_650m_score", "w_esm1b_650m"),
    "CaLM": ("calm_score", "w_calm"),
}

ENSEMBLES = [
    "ESM-2 150M + ESM-2 650M",
    "ESM-2 650M + ESM-1b 650M",
    "ESM-2 150M + CaLM",
    "ESM-2 650M + CaLM",
    "ESM-1b 650M + CaLM",
    "ESM-2 650M + ESM-1b 650M + CaLM",
]


def latex_escape(value: str) -> str:
    return value.replace("&", r"\&").replace("%", r"\%")


def main() -> None:
    scores = pd.read_csv(SCORES)
    oof = pd.read_csv(OOF, usecols=["variant_id", "fold"])
    weights = pd.read_csv(FOLD_WEIGHTS)

    if scores["variant_id"].duplicated().any() or oof["variant_id"].duplicated().any():
        raise ValueError("Variant identifiers must be unique in score and OOF tables")

    data = scores.merge(oof, on="variant_id", how="inner", validate="one_to_one")
    if len(data) != len(scores):
        raise ValueError("OOF fold assignments do not cover the complete score table")

    weights = weights[weights["model"].isin(ENSEMBLES)].copy()
    expected = len(ENSEMBLES) * data["fold"].nunique()
    if len(weights) != expected:
        raise ValueError(f"Expected {expected} ensemble-fold rows, found {len(weights)}")

    fold_rows: list[dict[str, object]] = []
    for _, row in weights.iterrows():
        fold = int(row["fold"])
        ensemble = str(row["model"])
        component_names = [name.strip() for name in str(row["components"]).split("+")]
        train = data[data["fold"] != fold]

        numerators: dict[str, float] = {}
        sds: dict[str, float] = {}
        raw: dict[str, float] = {}
        for component in component_names:
            score_col, weight_col = COMPONENTS[component]
            sigma = float(train[score_col].to_numpy(float).std(ddof=0))
            alpha = float(row[weight_col])
            sds[component] = sigma
            raw[component] = alpha
            numerators[component] = alpha * sigma

        denominator = sum(numerators.values())
        if denominator <= 0:
            raise ValueError(f"Invalid standardized-weight denominator for {ensemble}, fold {fold}")

        std_weights = {
            component: numerators[component] / denominator
            for component in component_names
        }

        # Confirm algebraic equivalence on held-out variants.
        test = data[data["fold"] == fold]
        raw_score = np.zeros(len(test), dtype=float)
        standardized_score = np.zeros(len(test), dtype=float)
        for component in component_names:
            score_col, _ = COMPONENTS[component]
            values = test[score_col].to_numpy(float)
            mean = float(train[score_col].mean())
            raw_score += raw[component] * values
            standardized_score += std_weights[component] * ((values - mean) / sds[component])
        expected_standardized = raw_score / denominator - sum(
            std_weights[component]
            * float(train[COMPONENTS[component][0]].mean())
            / sds[component]
            for component in component_names
        )
        max_affine_error = float(
            np.max(np.abs(standardized_score - expected_standardized))
        )

        for component in component_names:
            fold_rows.append(
                {
                    "ensemble": ensemble,
                    "fold": fold,
                    "component": component,
                    "training_score_sd": sds[component],
                    "raw_weight": raw[component],
                    "standardized_equivalent_weight": std_weights[component],
                    "weight_difference": std_weights[component] - raw[component],
                    "max_heldout_affine_equivalence_error": max_affine_error,
                }
            )

    fold_df = pd.DataFrame(fold_rows)
    summary = (
        fold_df.groupby(["ensemble", "component"], sort=False)
        .agg(
            n_folds=("fold", "nunique"),
            mean_training_score_sd=("training_score_sd", "mean"),
            sd_training_score_sd=("training_score_sd", "std"),
            mean_raw_weight=("raw_weight", "mean"),
            sd_raw_weight=("raw_weight", "std"),
            mean_standardized_equivalent_weight=("standardized_equivalent_weight", "mean"),
            sd_standardized_equivalent_weight=("standardized_equivalent_weight", "std"),
            mean_weight_difference=("weight_difference", "mean"),
            max_heldout_affine_equivalence_error=(
                "max_heldout_affine_equivalence_error",
                "max",
            ),
        )
        .reset_index()
    )
    summary["raw_weight_mean_sd"] = summary.apply(
        lambda r: f"{r['mean_raw_weight']:.3f} ({r['sd_raw_weight']:.3f})", axis=1
    )
    summary["standardized_weight_mean_sd"] = summary.apply(
        lambda r: (
            f"{r['mean_standardized_equivalent_weight']:.3f} "
            f"({r['sd_standardized_equivalent_weight']:.3f})"
        ),
        axis=1,
    )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_TEX.parent.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(OUTPUT_FOLD_CSV, index=False)
    summary.to_csv(OUTPUT_CSV, index=False)

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\caption{Scale sensitivity of optimised ensemble weights in the ClinVar analysis.}",
        r"\label{tab:standardized_equivalent_weights}",
        r"\begin{tabular}{llccc}",
        r"\toprule",
        r"Ensemble & Component & \shortstack{Training-score\\SD} & \shortstack{Raw weight\\mean (SD)} & \shortstack{Standardised-equivalent weight\\mean (SD)} \\",
        r"\midrule",
    ]
    for ensemble, group in summary.groupby("ensemble", sort=False):
        first = True
        for _, r in group.iterrows():
            ensemble_label = latex_escape(ensemble) if first else ""
            lines.append(
                f"{ensemble_label} & {latex_escape(str(r['component']))} & "
                f"{r['mean_training_score_sd']:.3f} & "
                f"{r['raw_weight_mean_sd']} & "
                f"{r['standardized_weight_mean_sd']} \\\\"
            )
            first = False
        lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\begin{minipage}{0.98\textwidth}",
            r"\footnotesize Values are means across the 10 gene-held-out cross-validation folds; values in parentheses are standard deviations across folds. Training-score standard deviations were calculated using the training variants in each fold. Standardised-equivalent weights were calculated as $w_j\sigma_j/\sum_l w_l\sigma_l$ and were used only to report coefficients; all out-of-fold predictions and AUROCs retain the original raw-score ensembles.",
            r"\end{minipage}",
            r"\end{table}",
            "",
        ]
    )
    OUTPUT_TEX.write_text("\n".join(lines))

    if not np.allclose(
        fold_df.groupby(["ensemble", "fold"])["standardized_equivalent_weight"].sum().to_numpy(),
        1.0,
        atol=1e-12,
    ):
        raise AssertionError("Standardized-equivalent weights do not sum to one")
    if fold_df["max_heldout_affine_equivalence_error"].max() > 1e-12:
        raise AssertionError("Weight conversion did not preserve held-out scores up to an affine transformation")

    print(summary.to_string(index=False))
    print(f"\nWrote {OUTPUT_CSV}")
    print(f"Wrote {OUTPUT_FOLD_CSV}")
    print(f"Wrote {OUTPUT_TEX}")


if __name__ == "__main__":
    main()
