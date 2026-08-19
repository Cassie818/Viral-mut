#!/usr/bin/env python3
"""Evaluate ClinMAVE functional-effect prediction with ESM-2 650M and CaLM."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


BASE = Path("Results/ClinMAVE")
ESM650 = Path("Results/Revision/ClinMAVE_ESM2_650M/clinmave_missense_all_with_esm2_650m.csv")
OUTDIR = Path("Results/Revision/ClinMAVE_ESM2_650M/functional_effects")

MERGE_KEYS = ["Identifier", "Gene", "Transcriptid", "Site"]
ESM_COLS = MERGE_KEYS + ["Ref", "Mut"]
CALM_COLS = MERGE_KEYS + ["Ref", "Mut"]
WEIGHTS = np.linspace(0.0, 1.0, 101)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--esm650", type=Path, default=ESM650)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument(
        "--drop-global-conflicts",
        action="store_true",
        help="Drop identifiers assigned to more than one functional class across all ClinMAVE missense files.",
    )
    return parser.parse_args()


def read_calm(base: Path, assay: str, effect_class: str) -> pd.DataFrame:
    path = base / assay / "missense" / f"{effect_class}_{assay.lower()}_LLR_CaLM_results.csv"
    df = pd.read_csv(path)
    keep = CALM_COLS + ["LLR"]
    df = df[keep].rename(columns={"Ref": "calm_ref_codon", "Mut": "calm_mut_codon", "LLR": "calm_llr"})
    df["assay"] = assay
    df["effect_class"] = effect_class
    return df


def drop_ambiguous_and_duplicate_variants(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    before = len(df)
    ambiguous = df.groupby("Identifier")["label"].nunique()
    ambiguous_ids = set(ambiguous[ambiguous > 1].index)
    if ambiguous_ids:
        df = df[~df["Identifier"].isin(ambiguous_ids)].copy()
    after_ambiguous = len(df)
    df = df.sort_values(["Identifier", "label"]).drop_duplicates("Identifier", keep="first")
    return df, {
        "input_rows": before,
        "ambiguous_identifier_rows": before - after_ambiguous,
        "duplicate_identifier_rows": after_ambiguous - len(df),
        "analysis_rows": len(df),
    }


def standardize_train_test(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    train_x = train[cols].to_numpy(float)
    test_x = test[cols].to_numpy(float)
    mean = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd[sd == 0] = 1.0
    return (train_x - mean) / sd, (test_x - mean) / sd


def best_weight(y: np.ndarray, esm: np.ndarray, calm: np.ndarray) -> float:
    best_w = 0.0
    best_auc = -np.inf
    for w in WEIGHTS:
        score = w * calm + (1.0 - w) * esm
        auc = roc_auc_score(y, score)
        if auc > best_auc:
            best_auc = auc
            best_w = float(w)
    return best_w


def safe_roc_auc(y: np.ndarray, score: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return np.nan
    return float(roc_auc_score(y, score))


def paired_t_p(a: pd.Series, b: pd.Series) -> float:
    pairs = pd.concat([a, b], axis=1).dropna()
    if len(pairs) < 2:
        return np.nan
    diff = pairs.iloc[:, 0] - pairs.iloc[:, 1]
    if np.allclose(diff, diff.iloc[0]):
        return np.nan
    return float(ttest_rel(pairs.iloc[:, 0], pairs.iloc[:, 1]).pvalue)


def make_splits(df: pd.DataFrame, folds: int):
    y = df["label"].to_numpy()
    groups = df["Gene"].astype(str).to_numpy()
    n_splits = min(folds, int(np.bincount(y).min()), int(pd.Series(groups).nunique()))
    if n_splits < 2:
        raise ValueError("Not enough data for cross-validation")
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=7)
    return list(splitter.split(df, y, groups=groups)), "stratified_gene"


def evaluate_dataset(df: pd.DataFrame, assay: str, case_class: str, folds: int) -> tuple[pd.DataFrame, dict]:
    splits, split_type = make_splits(df, folds)
    rows = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train = df.iloc[train_idx].copy()
        test = df.iloc[test_idx].copy()
        y_train = train["label"].to_numpy()
        y_test = test["label"].to_numpy()

        train_scaled, test_scaled = standardize_train_test(train, test, ["esm2_650m_score", "calm_score"])
        esm_train, calm_train = train_scaled[:, 0], train_scaled[:, 1]
        esm_test, calm_test = test_scaled[:, 0], test_scaled[:, 1]
        w = best_weight(y_train, esm_train, calm_train)

        esm_auc = safe_roc_auc(y_test, esm_test)
        calm_auc = safe_roc_auc(y_test, calm_test)
        combo_auc = safe_roc_auc(y_test, w * calm_test + (1.0 - w) * esm_test)

        fold_row = {
            "assay": assay,
            "case_class": case_class,
            "fold": fold,
            "split_type": split_type,
            "n_train": len(train),
            "n_test": len(test),
            "n_test_genes": test["Gene"].nunique(),
            "n_test_cases": int(y_test.sum()),
            "n_test_controls": int((1 - y_test).sum()),
            "calm_weight": w,
            "esm2_650m_weight": 1.0 - w,
            "auroc_esm2_650m": esm_auc,
            "auroc_calm": calm_auc,
            "auroc_esm2_650m_calm": combo_auc,
        }
        fold_row["delta_combo_vs_esm2_650m"] = fold_row["auroc_esm2_650m_calm"] - fold_row["auroc_esm2_650m"]
        fold_row["delta_combo_vs_calm"] = fold_row["auroc_esm2_650m_calm"] - fold_row["auroc_calm"]
        rows.append(fold_row)

    fold_df = pd.DataFrame(rows)
    summary = {
        "assay": assay,
        "case_class": case_class,
        "split_type": split_type,
        "n_variants": len(df),
        "n_cases": int(df["label"].sum()),
        "n_controls": int((1 - df["label"]).sum()),
        "n_genes": df["Gene"].nunique(),
        "n_folds": fold_df["fold"].nunique(),
        "n_evaluable_folds": int(fold_df["auroc_esm2_650m_calm"].notna().sum()),
        "mean_calm_weight": fold_df["calm_weight"].mean(),
        "sd_calm_weight": fold_df["calm_weight"].std(ddof=1),
    }
    for col in ["auroc_esm2_650m", "auroc_calm", "auroc_esm2_650m_calm", "delta_combo_vs_esm2_650m", "delta_combo_vs_calm"]:
        summary[f"{col}_mean"] = fold_df[col].mean()
        summary[f"{col}_sd"] = fold_df[col].std(ddof=1)
    summary["p_delta_combo_vs_esm2_650m_paired_t"] = paired_t_p(
        fold_df["auroc_esm2_650m_calm"], fold_df["auroc_esm2_650m"]
    )
    summary["p_delta_combo_vs_calm_paired_t"] = paired_t_p(
        fold_df["auroc_esm2_650m_calm"], fold_df["auroc_calm"]
    )
    return fold_df, summary


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    esm = pd.read_csv(args.esm650)
    esm = esm[ESM_COLS + ["assay", "effect_class", "LLR", "esm2_650m_llr", "esm2_650m_score"]]
    esm = esm.rename(columns={"Ref": "esm_ref_aa", "Mut": "esm_mut_aa", "LLR": "esm2_150m_llr"})
    global_conflict_ids: set[str] = set()
    if args.drop_global_conflicts:
        class_counts = esm.dropna(subset=["esm2_650m_score"]).groupby("Identifier")["effect_class"].nunique()
        global_conflict_ids = set(class_counts[class_counts > 1].index)
        esm = esm[~esm["Identifier"].isin(global_conflict_ids)].copy()

    all_folds = []
    summaries = []
    audits = []
    audits.append({"global_conflicting_identifiers_removed": len(global_conflict_ids)})
    merged_outputs = []

    for assay in ["DMS", "CBGE"]:
        for case_class in ["lof", "gof"]:
            parts = []
            for effect_class, label in [("normal", 0), (case_class, 1)]:
                e = esm[(esm["assay"] == assay) & (esm["effect_class"] == effect_class)].copy()
                c = read_calm(args.base, assay, effect_class)
                merged = e.merge(c, on=MERGE_KEYS + ["assay", "effect_class"], how="inner")
                merged["label"] = label
                parts.append(merged)

                audits.append({
                    "assay": assay,
                    "case_class": case_class,
                    "effect_class": effect_class,
                    "esm2_650m_rows": len(e),
                    "calm_rows": len(c),
                    "merged_rows": len(merged),
                })

            df = pd.concat(parts, ignore_index=True)
            df, dup_audit = drop_ambiguous_and_duplicate_variants(df)
            df["calm_score"] = -df["calm_llr"]
            df["esm2_150m_score"] = -df["esm2_150m_llr"]
            before_complete = len(df)
            df = df.dropna(subset=["esm2_650m_score", "calm_score", "label", "Gene"]).copy()
            dup_audit["incomplete_score_rows"] = before_complete - len(df)
            df["comparison"] = f"{assay}_{case_class}_vs_normal"
            dup_audit.update({"assay": assay, "case_class": case_class})
            audits.append(dup_audit)
            merged_outputs.append(df)

            fold_df, summary = evaluate_dataset(df, assay, case_class, args.folds)
            all_folds.append(fold_df)
            summaries.append(summary)

    pd.concat(all_folds, ignore_index=True).to_csv(args.outdir / "clinmave_esm2_650m_calm_cv_fold_metrics.csv", index=False)
    pd.DataFrame(summaries).to_csv(args.outdir / "clinmave_esm2_650m_calm_cv_summary.csv", index=False)
    pd.DataFrame(audits).to_csv(args.outdir / "clinmave_esm2_650m_calm_merge_audit.csv", index=False)
    pd.concat(merged_outputs, ignore_index=True).to_csv(args.outdir / "clinmave_esm2_650m_calm_analysis_table.csv", index=False)

    print(pd.DataFrame(summaries).to_string(index=False))
    print(f"\nWrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
