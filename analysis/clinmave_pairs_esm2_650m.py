#!/usr/bin/env python3
"""Evaluate CaLM contribution across ClinMAVE DMS/CBGE dataset pairs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


RAW_JSONL = Path("Results/Revision/ClinMAVE_api_cross_platform/clinmave_dms_cbge_gene_variants.jsonl")
PAIRS = Path("Results/Revision/ClinMAVE_api_cross_platform/clinmave_dataset_pair_eligible_old_style.csv")
ESM2_650M = Path("Results/Revision/ClinMAVE_ESM2_650M/clinmave_missense_all_with_esm2_650m.csv")
CALM_BASE = Path("Results/ClinMAVE")
OUTDIR = Path("Results/Revision/ClinMAVE_api_cross_platform/dataset_pair_modality_weights")
WEIGHTS = np.linspace(0.0, 1.0, 201)
MAX_FOLDS = 10
RANDOM_SEED = 16
MIN_PREMATCH_CLASS_COUNT = 5
ELIGIBILITY_COUNT_COLUMNS = [
    "DMS_cases",
    "DMS_controls",
    "CBGE_cases",
    "CBGE_controls",
]


def load_pairs() -> pd.DataFrame:
    pairs = pd.read_csv(PAIRS)
    missing = set(ELIGIBILITY_COUNT_COLUMNS) - set(pairs.columns)
    if missing:
        raise ValueError(f"Pair eligibility table is missing columns: {sorted(missing)}")
    eligible = pairs[ELIGIBILITY_COUNT_COLUMNS].ge(MIN_PREMATCH_CLASS_COUNT).all(axis=1)
    if not eligible.all():
        raise ValueError(
            "Every candidate pair must have at least five normal and five LoF "
            "variants on each platform before score matching"
        )
    return pairs


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
            if y is None or assay is None:
                continue
            if record.get("molecularConsequence") != "Missense":
                continue
            rows.append(
                {
                    "Gene": record.get("geneName"),
                    "Identifier": record.get("identifier"),
                    "datasetId": record.get("datasetId"),
                    "platform": assay,
                    "label": y,
                    "phenotype": record.get("phenotype"),
                    "pmid": record.get("pmid"),
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["Gene", "Identifier", "datasetId", "platform", "label"])


def load_scores() -> pd.DataFrame:
    esm = pd.read_csv(ESM2_650M, usecols=["Identifier", "esm2_650m_llr"])
    esm = esm.dropna(subset=["Identifier", "esm2_650m_llr"]).drop_duplicates("Identifier")

    calm_frames = []
    for assay in ["DMS", "CBGE"]:
        for effect in ["normal", "lof", "gof"]:
            path = CALM_BASE / assay / "missense" / f"{effect}_{assay.lower()}_LLR_CaLM_results.csv"
            if not path.exists():
                continue
            df = pd.read_csv(path, usecols=["Identifier", "LLR"]).rename(columns={"LLR": "calm_llr"})
            calm_frames.append(df)
    calm = pd.concat(calm_frames, ignore_index=True)
    calm = calm.dropna(subset=["Identifier", "calm_llr"]).drop_duplicates("Identifier")
    scores = esm.merge(calm, on="Identifier", how="inner")
    return scores


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
    prot = train["esm2_650m_llr"].to_numpy(float)
    calm = train["calm_llr"].to_numpy(float)
    best_w, best_auc = np.nan, -np.inf
    for w in WEIGHTS:
        combo = w * calm + (1.0 - w) * prot
        auc = roc_auc_score(y, -combo)
        if auc > best_auc:
            best_w = float(w)
            best_auc = auc
    return best_w


def evaluate_pair(table: pd.DataFrame, random_state: int = RANDOM_SEED) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame()
    split_y = table["DMS_label"].astype(str) + "_" + table["CBGE_label"].astype(str)
    min_stratum = split_y.value_counts().min()
    if pd.isna(min_stratum):
        return pd.DataFrame()
    n_splits = min(MAX_FOLDS, int(min_stratum))
    if n_splits < 2:
        return pd.DataFrame()
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    rows = []
    for fold, (train_idx, test_idx) in enumerate(cv.split(table, split_y), start=1):
        train = table.iloc[train_idx]
        test = table.iloc[test_idx]
        for platform, label_col in [("DMS", "DMS_label"), ("CBGE", "CBGE_label")]:
            w = best_weight(train, label_col)
            y_test = test[label_col].to_numpy(int)
            rows.append(
                {
                    "fold": fold,
                    "platform": platform,
                    "calm_weight": w,
                    "n_test": len(test),
                    "n_test_cases": int(y_test.sum()),
                    "n_test_controls": int((1 - y_test).sum()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = load_records()
    pairs = load_pairs()
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

    matched.to_csv(OUTDIR / "dataset_pair_matched_variants_with_scores.csv", index=False)
    folds.to_csv(OUTDIR / "dataset_pair_fold_weights.csv", index=False)
    print(f"\nWrote outputs to {OUTDIR}")


if __name__ == "__main__":
    main()
