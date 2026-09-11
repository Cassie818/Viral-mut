#!/usr/bin/env python3
"""Create scale-adjusted, standardized-equivalent weights for Figs. 3, 6, and 7.

For a two-model raw-score ensemble with CaLM weight w, the equivalent CaLM
weight after z-standardizing both scores on the training data is

    w_std = w * sd_calm / ((1 - w) * sd_plm + w * sd_calm).

The transformed coefficient is used only for reporting. Original ensemble
predictions are retained, so pooled out-of-fold rankings and AUROC are unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold


ROOT = Path(__file__).resolve().parents[1]


def standardized_equivalent(w_calm: float, sd_plm: float, sd_calm: float) -> float:
    denominator = (1.0 - w_calm) * sd_plm + w_calm * sd_calm
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("Score standard deviations must give a positive finite denominator")
    return float(w_calm * sd_calm / denominator)


def fig3_weights() -> None:
    base = ROOT / "Results/ClinMAVE/functional_effects"
    specs = {
        "ESM-2 (650M)": (
            base / "clinmave_esm2_650m_calm_cv_fold_metrics.csv",
            base / "clinmave_esm2_650m_calm_oof_predictions.csv.gz",
        ),
        "ESM-1b (650M)": (
            base / "clinmave_esm1b_650m_calm_cv_fold_metrics.csv",
            base / "clinmave_esm1b_650m_calm_oof_predictions.csv.gz",
        ),
    }
    rows: list[dict[str, object]] = []
    for model, (fold_path, oof_path) in specs.items():
        folds = pd.read_csv(fold_path)
        oof = pd.read_csv(oof_path)
        for _, fold_row in folds.iterrows():
            subset = oof[
                (oof["assay"] == fold_row["assay"])
                & (oof["case_class"] == fold_row["case_class"])
            ]
            train = subset[subset["fold"] != int(fold_row["fold"])]
            sd_plm = float(train["score_plm"].std(ddof=1))
            sd_calm = float(train["score_calm"].std(ddof=1))
            raw = float(fold_row["calm_weight"])
            std = standardized_equivalent(raw, sd_plm, sd_calm)
            rows.append(
                {
                    "model": model,
                    "assay": fold_row["assay"],
                    "case_class": fold_row["case_class"],
                    "fold": int(fold_row["fold"]),
                    "raw_calm_weight": raw,
                    "training_sd_plm": sd_plm,
                    "training_sd_calm": sd_calm,
                    "standardized_equivalent_calm_weight": std,
                }
            )
    fold_out = pd.DataFrame(rows)
    fold_out.to_csv(base / "fig3_standardized_equivalent_fold_weights.csv", index=False)
    summary = (
        fold_out.groupby(["model", "assay", "case_class"], as_index=False)
        .agg(
            mean_raw_calm_weight=("raw_calm_weight", "mean"),
            sd_raw_calm_weight=("raw_calm_weight", "std"),
            mean_standardized_equivalent_calm_weight=("standardized_equivalent_calm_weight", "mean"),
            sd_standardized_equivalent_calm_weight=("standardized_equivalent_calm_weight", "std"),
        )
    )
    summary.to_csv(base / "fig3_standardized_equivalent_weight_summary.csv", index=False)


def fig6_weights() -> None:
    source_path = ROOT / "Results/Revision/len1022_model_control/len1022_model_control_score_table_complete_cases.csv"
    fold_path = ROOT / "Results/ClinVar/gene_level/gene_level_nested_cv_fold_metrics.csv"
    out_base = ROOT / "Results/ClinVar/gene_level"
    source = pd.read_csv(source_path).drop_duplicates(subset=["variant_id", "Gene_gene"])
    folds = pd.read_csv(fold_path)
    score_pairs = {
        "ESM-2 650M + CaLM": ("esm2_650m_score", "calm_score", 2),
        "ESM-1b 650M + CaLM": ("esm1b_650m_score", "calm_score", 4),
    }
    eligible_genes = sorted(folds["gene"].unique())
    rows: list[dict[str, object]] = []
    for gene_index, gene in enumerate(eligible_genes):
        gene_df = source[source["Gene_gene"] == gene]
        y = gene_df["label"].to_numpy(int)
        for model, (plm_col, calm_col, offset) in score_pairs.items():
            model_folds = folds[(folds["gene"] == gene) & (folds["model"] == model)].sort_values("fold")
            if model_folds.empty:
                raise ValueError(f"Missing fold weights for {gene}, {model}")
            splitter = StratifiedKFold(
                n_splits=len(model_folds),
                shuffle=True,
                random_state=16 + 1009 * offset + gene_index,
            )
            splits = list(splitter.split(np.zeros(len(y)), y))
            for (_, fold_row), (train_idx, test_idx) in zip(model_folds.iterrows(), splits):
                if len(train_idx) != int(fold_row["n_train"]) or len(test_idx) != int(fold_row["n_test"]):
                    raise AssertionError(f"Fold reconstruction mismatch for {gene}, {model}")
                sd_plm = float(gene_df.iloc[train_idx][plm_col].std(ddof=1))
                sd_calm = float(gene_df.iloc[train_idx][calm_col].std(ddof=1))
                raw = float(fold_row["weight_second"])
                rows.append(
                    {
                        "gene": gene,
                        "model": model,
                        "fold": int(fold_row["fold"]),
                        "raw_calm_weight": raw,
                        "training_sd_plm": sd_plm,
                        "training_sd_calm": sd_calm,
                        "standardized_equivalent_calm_weight": standardized_equivalent(raw, sd_plm, sd_calm),
                    }
                )
    fold_out = pd.DataFrame(rows)
    fold_out.to_csv(out_base / "fig6_standardized_equivalent_fold_weights.csv", index=False)
    summary = (
        fold_out.groupby(["gene", "model"], as_index=False)
        .agg(
            mean_raw_calm_weight=("raw_calm_weight", "mean"),
            sd_raw_calm_weight=("raw_calm_weight", "std"),
            mean_standardized_equivalent_calm_weight=("standardized_equivalent_calm_weight", "mean"),
            sd_standardized_equivalent_calm_weight=("standardized_equivalent_calm_weight", "std"),
        )
    )
    summary.to_csv(out_base / "fig6_standardized_equivalent_weight_summary.csv", index=False)


def fig7_weights() -> None:
    revision = ROOT / "Results/Revision/ClinMAVE_api_cross_platform"
    out_base = ROOT / "Results/ClinMAVE/cross_platform_context"
    specs = {
        "150M": revision / "dataset_pair_modality_weights_esm2_150m",
        "650M": revision / "dataset_pair_modality_weights",
    }
    all_fold_rows: list[dict[str, object]] = []
    all_summary_rows: list[dict[str, object]] = []
    for model_short, base in specs.items():
        matched = pd.read_csv(base / "dataset_pair_matched_variants_with_scores.csv")
        folds = pd.read_csv(base / "dataset_pair_fold_metrics.csv")
        plm_col = f"esm2_{model_short.lower()}_llr"
        for pair_id, pair_df in matched.groupby("pair_id", sort=False):
            pair_folds = folds[folds["pair_id"] == pair_id]
            if pair_folds.empty:
                continue
            split_y = pair_df["DMS_label"].astype(str) + "_" + pair_df["CBGE_label"].astype(str)
            n_splits = min(10, int(split_y.value_counts().min()))
            if n_splits < 2:
                raise AssertionError(f"Stored fold rows exist for unevaluable pair {pair_id}")
            splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=16)
            splits = list(splitter.split(pair_df, split_y))
            for platform in ["DMS", "CBGE"]:
                platform_rows: list[dict[str, object]] = []
                platform_folds = pair_folds[pair_folds["platform"] == platform].sort_values("fold")
                for (_, fold_row), (train_idx, test_idx) in zip(platform_folds.iterrows(), splits):
                    if len(test_idx) != int(fold_row["n_test"]):
                        raise AssertionError(f"Fold reconstruction mismatch for {pair_id}, {model_short}")
                    train = pair_df.iloc[train_idx]
                    sd_plm = float(train[plm_col].std(ddof=1))
                    sd_calm = float(train["calm_llr"].std(ddof=1))
                    raw = float(fold_row["calm_weight"])
                    row = {
                        "pair_id": pair_id,
                        "Gene": fold_row["Gene"],
                        "DMS_dataset": fold_row["DMS_dataset"],
                        "CBGE_dataset": fold_row["CBGE_dataset"],
                        "model_short": model_short,
                        "platform": platform,
                        "fold": int(fold_row["fold"]),
                        "raw_calm_weight": raw,
                        "training_sd_plm": sd_plm,
                        "training_sd_calm": sd_calm,
                        "standardized_equivalent_calm_weight": standardized_equivalent(raw, sd_plm, sd_calm),
                    }
                    platform_rows.append(row)
                    all_fold_rows.append(row)
                platform_df = pd.DataFrame(platform_rows)
                all_summary_rows.append(
                    {
                        "pair_id": pair_id,
                        "Gene": platform_df["Gene"].iloc[0],
                        "DMS_dataset": platform_df["DMS_dataset"].iloc[0],
                        "CBGE_dataset": platform_df["CBGE_dataset"].iloc[0],
                        "model_short": model_short,
                        "platform": platform,
                        "mean_raw_calm_weight": platform_df["raw_calm_weight"].mean(),
                        "mean_standardized_equivalent_calm_weight": platform_df["standardized_equivalent_calm_weight"].mean(),
                    }
                )
    pd.DataFrame(all_fold_rows).to_csv(out_base / "fig7_standardized_equivalent_fold_weights.csv", index=False)
    summary = pd.DataFrame(all_summary_rows)
    wide = summary.pivot(
        index=["pair_id", "Gene", "DMS_dataset", "CBGE_dataset", "model_short"],
        columns="platform",
        values="mean_standardized_equivalent_calm_weight",
    ).reset_index()
    wide = wide.rename(columns={"DMS": "DMS_weight_standardized", "CBGE": "CBGE_weight_standardized"})
    wide["delta_weight_standardized"] = wide["CBGE_weight_standardized"] - wide["DMS_weight_standardized"]
    raw_table = pd.read_csv(out_base / "dataset_pair_modality_weights_11_gene_pairs_table.csv")
    merge_keys = ["Gene", "DMS_dataset", "CBGE_dataset"]
    for model_short in ["150M", "650M"]:
        add = wide[wide["model_short"] == model_short][
            merge_keys
            + ["DMS_weight_standardized", "CBGE_weight_standardized", "delta_weight_standardized"]
        ].rename(
            columns={
                "DMS_weight_standardized": f"DMS_weight_standardized_{model_short}",
                "CBGE_weight_standardized": f"CBGE_weight_standardized_{model_short}",
                "delta_weight_standardized": f"delta_weight_standardized_{model_short}",
            }
        )
        raw_table = raw_table.merge(add, on=merge_keys, how="left", validate="one_to_one")
    raw_table.to_csv(out_base / "dataset_pair_modality_weights_11_gene_pairs_table_standardized.csv", index=False)


def main() -> None:
    fig3_weights()
    fig6_weights()
    fig7_weights()
    print("Wrote standardized-equivalent weight tables for Figs. 3, 6, and 7")


if __name__ == "__main__":
    main()
