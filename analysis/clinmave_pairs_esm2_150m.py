#!/usr/bin/env python3
"""Evaluate CaLM contribution across ClinMAVE DMS/CBGE pairs with ESM-2 150M."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


RAW_JSONL = Path("Results/Revision/ClinMAVE_api_cross_platform/clinmave_dms_cbge_gene_variants.jsonl")
PAIRS = Path("Results/Revision/ClinMAVE_api_cross_platform/clinmave_dataset_pair_eligible_old_style.csv")
SCORE_BASE = Path("Results/ClinMAVE")
OUTDIR = Path("Results/Revision/ClinMAVE_api_cross_platform/dataset_pair_modality_weights_esm2_150m")
WEIGHTS = np.linspace(0.0, 1.0, 201)


def platform(value: str) -> str | None:
    if value == "Deep Mutational Scanning":
        return "DMS"
    if value == "CRISPR-Based Genome Editing":
        return "CBGE"
    return None


def label(value: str) -> int | None:
    if value == "Functionally normal":
        return 0
    if value == "Loss-of-function":
        return 1
    return None


def load_records() -> pd.DataFrame:
    rows = []
    with RAW_JSONL.open() as handle:
        for line in handle:
            record = json.loads(line)
            y = label(record.get("consequenceClass"))
            assay = platform(record.get("maveTechnique"))
            if y is None or assay is None or record.get("molecularConsequence") != "Missense":
                continue
            rows.append(
                {
                    "Gene": record.get("geneName"),
                    "Identifier": record.get("identifier"),
                    "datasetId": record.get("datasetId"),
                    "platform": assay,
                    "label": y,
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["Gene", "Identifier", "datasetId", "platform", "label"])


def load_model_scores(prefix: str, score_name: str) -> pd.DataFrame:
    frames = []
    for assay in ["DMS", "CBGE"]:
        for effect in ["normal", "lof", "gof"]:
            suffix = "_LLR_CaLM_results.csv" if prefix == "calm" else "_LLR_results.csv"
            path = SCORE_BASE / assay / "missense" / f"{effect}_{assay.lower()}{suffix}"
            if not path.exists():
                continue
            df = pd.read_csv(path, usecols=["Identifier", "LLR"]).rename(columns={"LLR": score_name})
            frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    return out.dropna(subset=["Identifier", score_name]).drop_duplicates("Identifier")


def load_scores() -> pd.DataFrame:
    esm = load_model_scores("esm2_150m", "esm2_150m_llr")
    calm = load_model_scores("calm", "calm_llr")
    return esm.merge(calm, on="Identifier", how="inner")


def pair_table(records: pd.DataFrame, pair: pd.Series, scores: pd.DataFrame) -> pd.DataFrame:
    gene = pair["Gene"]
    dms = records[
        (records["Gene"] == gene)
        & (records["platform"] == "DMS")
        & (records["datasetId"] == pair["DMS_dataset"])
    ]
    cbge = records[
        (records["Gene"] == gene)
        & (records["platform"] == "CBGE")
        & (records["datasetId"] == pair["CBGE_dataset"])
    ]

    def strict_labels(df: pd.DataFrame, col: str) -> pd.DataFrame:
        grouped = (
            df.groupby("Identifier")["label"]
            .agg(lambda x: "|".join(map(str, sorted(set(x)))))
            .reset_index()
        )
        grouped = grouped[grouped["label"].isin(["0", "1"])].copy()
        grouped[col] = grouped["label"].astype(int)
        return grouped[["Identifier", col]]

    wide = strict_labels(dms, "DMS_label").merge(strict_labels(cbge, "CBGE_label"), on="Identifier", how="inner")
    return wide.merge(scores, on="Identifier", how="inner")


def best_weight(train: pd.DataFrame, label_col: str) -> float:
    y = train[label_col].to_numpy(int)
    if len(np.unique(y)) < 2:
        return np.nan
    prot = train["esm2_150m_llr"].to_numpy(float)
    calm = train["calm_llr"].to_numpy(float)
    best_w, best_auc = np.nan, -np.inf
    for w in WEIGHTS:
        combo = w * calm + (1.0 - w) * prot
        auc = roc_auc_score(y, -combo)
        if auc > best_auc:
            best_w = float(w)
            best_auc = auc
    return best_w


def evaluate_pair(table: pd.DataFrame, random_state: int = 16) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame()
    split_y = table["DMS_label"].astype(str) + "_" + table["CBGE_label"].astype(str)
    min_stratum = split_y.value_counts().min()
    if pd.isna(min_stratum):
        return pd.DataFrame()
    n_splits = min(10, int(min_stratum))
    if n_splits < 2:
        return pd.DataFrame()
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    rows = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(table, split_y), start=1):
        train = table.iloc[train_idx]
        test = table.iloc[test_idx]
        for assay, label_col in [("DMS", "DMS_label"), ("CBGE", "CBGE_label")]:
            y_test = test[label_col].to_numpy(int)
            w = best_weight(train, label_col)
            evaluable = len(np.unique(y_test)) == 2 and np.isfinite(w)
            if evaluable:
                prot = test["esm2_150m_llr"].to_numpy(float)
                calm = test["calm_llr"].to_numpy(float)
                combo = w * calm + (1.0 - w) * prot
                auroc_prot = roc_auc_score(y_test, -prot)
                auroc_calm = roc_auc_score(y_test, -calm)
                auroc_combo = roc_auc_score(y_test, -combo)
            else:
                auroc_prot = auroc_calm = auroc_combo = np.nan
            rows.append(
                {
                    "fold": fold,
                    "platform": assay,
                    "calm_weight": w,
                    "auroc_esm2_150m": auroc_prot,
                    "auroc_calm": auroc_calm,
                    "auroc_combo": auroc_combo,
                    "n_test": len(test),
                    "n_test_cases": int(y_test.sum()),
                    "n_test_controls": int((1 - y_test).sum()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = load_records()
    pairs = pd.read_csv(PAIRS)
    scores = load_scores()

    matched_rows = []
    fold_rows = []
    for _, pair in pairs.iterrows():
        table = pair_table(records, pair, scores)
        pair_id = f"{pair['Gene']}|{pair['DMS_dataset']}|{pair['CBGE_dataset']}"
        table["pair_id"] = pair_id
        table["Gene"] = pair["Gene"]
        table["DMS_dataset"] = pair["DMS_dataset"]
        table["CBGE_dataset"] = pair["CBGE_dataset"]
        matched_rows.append(table)
        folds = evaluate_pair(table)
        if folds.empty:
            continue
        folds["pair_id"] = pair_id
        folds["Gene"] = pair["Gene"]
        folds["DMS_dataset"] = pair["DMS_dataset"]
        folds["CBGE_dataset"] = pair["CBGE_dataset"]
        fold_rows.append(folds)

    matched = pd.concat(matched_rows, ignore_index=True) if matched_rows else pd.DataFrame()
    folds = pd.concat(fold_rows, ignore_index=True) if fold_rows else pd.DataFrame()

    summary_rows = []
    if not folds.empty:
        for (pair_id, assay), sub in folds.groupby(["pair_id", "platform"]):
            meta = sub.iloc[0][["Gene", "DMS_dataset", "CBGE_dataset"]].to_dict()
            mt = matched[matched["pair_id"] == pair_id]
            label_col = f"{assay}_label"
            summary_rows.append(
                {
                    "pair_id": pair_id,
                    **meta,
                    "platform": assay,
                    "n_variants": len(mt),
                    "n_cases": int(mt[label_col].sum()),
                    "n_controls": int((1 - mt[label_col]).sum()),
                    "label_concordance": float((mt["DMS_label"] == mt["CBGE_label"]).mean()),
                    "n_folds": sub["fold"].nunique(),
                    "mean_calm_weight": sub["calm_weight"].mean(),
                    "sd_calm_weight": sub["calm_weight"].std(ddof=1),
                    "mean_auroc_esm2_150m": sub["auroc_esm2_150m"].mean(),
                    "mean_auroc_calm": sub["auroc_calm"].mean(),
                    "mean_auroc_combo": sub["auroc_combo"].mean(),
                    "mean_delta_combo_vs_esm2_150m": (sub["auroc_combo"] - sub["auroc_esm2_150m"]).mean(),
                }
            )
    summary = pd.DataFrame(summary_rows)
    tests = []
    for pair_id, sub in folds.groupby("pair_id") if not folds.empty else []:
        wide = sub.pivot(index="fold", columns="platform", values="calm_weight").dropna()
        row = {"pair_id": pair_id, "paired_folds": len(wide)}
        if {"DMS", "CBGE"}.issubset(wide.columns) and len(wide) >= 2:
            diff = wide["CBGE"] - wide["DMS"]
            row["delta_CBGE_minus_DMS"] = float(diff.mean())
            row["p_wilcoxon"] = 1.0 if np.allclose(diff, 0) else wilcoxon(wide["CBGE"], wide["DMS"]).pvalue
        tests.append(row)
    tests = pd.DataFrame(tests)
    if not summary.empty and not tests.empty:
        summary = summary.merge(tests, on="pair_id", how="left")

    matched.to_csv(OUTDIR / "dataset_pair_matched_variants_with_scores.csv", index=False)
    folds.to_csv(OUTDIR / "dataset_pair_fold_metrics.csv", index=False)
    summary.to_csv(OUTDIR / "dataset_pair_summary.csv", index=False)

    compact = summary.pivot_table(
        index=["pair_id", "Gene", "DMS_dataset", "CBGE_dataset", "n_variants", "label_concordance"],
        columns="platform",
        values=["mean_calm_weight", "mean_auroc_esm2_150m", "mean_auroc_calm", "mean_auroc_combo"],
        aggfunc="first",
    )
    compact.columns = [f"{metric}_{assay}" for metric, assay in compact.columns]
    compact = compact.reset_index().merge(tests, on="pair_id", how="left")
    compact = compact.sort_values(["n_variants", "Gene"], ascending=[False, True])
    compact.to_csv(OUTDIR / "dataset_pair_summary_compact.csv", index=False)

    show_cols = [
        "Gene",
        "DMS_dataset",
        "CBGE_dataset",
        "n_variants",
        "label_concordance",
        "mean_calm_weight_DMS",
        "mean_calm_weight_CBGE",
        "delta_CBGE_minus_DMS",
        "p_wilcoxon",
    ]
    print(compact[show_cols].to_string(index=False))
    print(f"\nWrote outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
