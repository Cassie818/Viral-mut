#!/usr/bin/env python3
"""ClinVar model-control ensembles.

This compares same-modality PLM+PLM ensembles against cross-modal PLM+CaLM
ensembles using gene-held-out folds. Ensemble weights are optimized only on
training genes and evaluated on held-out genes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


PATHOGENIC_LABELS = {"pathogenic", "likely_pathogenic"}
KEY_COLS = [
    "Gene_prot",
    "Site_prot",
    "Ref_prot",
    "Mut_prot",
    "Gene_gene",
    "Site_gene",
    "Ref_gene",
    "Mut_gene",
    "Label_prot",
    "Label_gene",
]
SCORE_COLUMNS = {
    "ESM-2 150M": "esm2_150m_score",
    "ESM-2 650M": "esm2_650m_score",
    "ESM-1b 650M": "esm1b_650m_score",
    "CaLM": "calm_score",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        default="Results/Revision/length_filtered_clinvar/clinvar_missense_len1022_scores.csv",
    )
    parser.add_argument(
        "--esm2-650m",
        default=(
            "Results/Revision/esm2_650m_len1022_full_clinvar/"
            "clinvar_missense_len1022_esm2_650m_scores_clean.csv"
        ),
    )
    parser.add_argument(
        "--esm1b-650m",
        default=(
            "Results/Revision/esm1b_650m_len1022_full_clinvar/"
            "clinvar_missense_len1022_esm1b_650m_scores_clean.csv"
        ),
    )
    parser.add_argument(
        "--out-dir",
        default="Results/Revision/len1022_model_control",
    )
    parser.add_argument("--n-splits", type=int, default=10)
    parser.add_argument("--seed", type=int, default=16)
    parser.add_argument("--pair-grid-step", type=float, default=0.01)
    parser.add_argument("--triple-grid-step", type=float, default=0.05)
    return parser.parse_args()


def pair_weights(step: float) -> list[tuple[float, float]]:
    values = np.arange(0.0, 1.0 + step / 2.0, step)
    return [(float(1.0 - w), float(w)) for w in values]


def triple_weights(step: float) -> list[tuple[float, float, float]]:
    values = np.arange(0.0, 1.0 + step / 2.0, step)
    weights = []
    for w1 in values:
        for w2 in values:
            w3 = 1.0 - w1 - w2
            if w3 >= -1e-9:
                weights.append((float(w1), float(w2), float(max(0.0, w3))))
    return weights


def best_weights(
    y_train: np.ndarray,
    train_scores: list[np.ndarray],
    candidates: list[tuple[float, ...]],
) -> tuple[float, ...]:
    best_auc = -np.inf
    best_candidate = candidates[0]
    for weights in candidates:
        train_mix = sum(w * score for w, score in zip(weights, train_scores))
        auc = roc_auc_score(y_train, train_mix)
        if auc > best_auc:
            best_auc = auc
            best_candidate = weights
    return best_candidate


def mix_scores(weights: tuple[float, ...], scores: list[np.ndarray]) -> np.ndarray:
    return sum(w * score for w, score in zip(weights, scores))


def load_scores(args: argparse.Namespace) -> pd.DataFrame:
    base = pd.read_csv(args.base)
    esm2 = pd.read_csv(args.esm2_650m)[KEY_COLS + ["esm2_650m_llr", "esm2_650m_score"]]
    esm1b = pd.read_csv(args.esm1b_650m)[KEY_COLS + ["esm1b_650m_llr", "esm1b_650m_score"]]

    df = base.copy()
    df["label"] = df["Label_prot"].astype(str).isin(PATHOGENIC_LABELS).astype(int)
    df["esm2_150m_llr"] = df["LLR_prot"].astype(float)
    df["esm2_150m_score"] = -df["esm2_150m_llr"]
    df["calm_llr"] = df["LLR_gene"].astype(float)
    df["calm_score"] = -df["calm_llr"]
    df["variant_id"] = (
        df["Gene_prot"].astype(str)
        + ":"
        + df["Site_prot"].astype(str)
        + ":"
        + df["Ref_prot"].astype(str)
        + ">"
        + df["Mut_prot"].astype(str)
        + ":"
        + df["Ref_gene"].astype(str)
        + ">"
        + df["Mut_gene"].astype(str)
    )
    df = df.merge(esm2, on=KEY_COLS, how="left", validate="one_to_one")
    df = df.merge(esm1b, on=KEY_COLS, how="left", validate="one_to_one")
    return df


def paired_tests(fold_df: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("ESM-2 150M + ESM-2 650M", "ESM-2 650M"),
        ("ESM-2 650M + ESM-1b 650M", "ESM-2 650M"),
        ("ESM-2 650M + CaLM", "ESM-2 650M"),
        ("ESM-1b 650M + CaLM", "ESM-1b 650M"),
        ("ESM-2 650M + CaLM", "ESM-2 650M + ESM-1b 650M"),
        (
            "ESM-2 650M + ESM-1b 650M + CaLM",
            "ESM-2 650M + ESM-1b 650M",
        ),
        (
            "ESM-2 650M + ESM-1b 650M + CaLM",
            "ESM-2 650M + CaLM",
        ),
    ]
    wide = fold_df.pivot(index="fold", columns="model", values="test_auc")
    rows = []
    for a, b in comparisons:
        delta = wide[a] - wide[b]
        try:
            wilcoxon_p = float(wilcoxon(delta).pvalue)
        except ValueError:
            wilcoxon_p = np.nan
        rows.append(
            {
                "a": a,
                "b": b,
                "mean_delta": float(delta.mean()),
                "sd_delta": float(delta.std(ddof=1)),
                "paired_t_p": float(ttest_rel(wide[a], wide[b]).pvalue),
                "wilcoxon_p": wilcoxon_p,
                "fold_deltas": ";".join(f"{value:.5f}" for value in delta),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_all = load_scores(args)
    score_cols = list(SCORE_COLUMNS.values())
    df = df_all.dropna(subset=["label", "Gene_prot", *score_cols]).copy()
    df_all.to_csv(out_dir / "len1022_model_control_score_table_all_variants.csv", index=False)
    df.to_csv(out_dir / "len1022_model_control_score_table_complete_cases.csv", index=False)

    y = df["label"].to_numpy(dtype=int)
    groups = df["Gene_prot"].astype(str).to_numpy()
    splits = list(
        StratifiedGroupKFold(
            n_splits=args.n_splits, shuffle=True, random_state=args.seed
        ).split(df, y, groups=groups)
    )

    pair_candidates = pair_weights(args.pair_grid_step)
    triple_candidates = triple_weights(args.triple_grid_step)
    model_specs = [
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
        (
            "ESM-2 650M + ESM-1b 650M + CaLM",
            ["ESM-2 650M", "ESM-1b 650M", "CaLM"],
        ),
    ]

    rows = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        y_train = y[train_idx]
        y_test = y[test_idx]
        for model_name, components in model_specs:
            if len(components) == 1:
                weights = (1.0,)
            else:
                candidates = pair_candidates if len(components) == 2 else triple_candidates
                train_scores = [
                    df.iloc[train_idx][SCORE_COLUMNS[component]].to_numpy()
                    for component in components
                ]
                weights = best_weights(y_train, train_scores, candidates)

            test_scores = [
                df.iloc[test_idx][SCORE_COLUMNS[component]].to_numpy()
                for component in components
            ]
            test_score = mix_scores(weights, test_scores)
            weight_by_component = dict(zip(components, weights))
            rows.append(
                {
                    "fold": fold,
                    "model": model_name,
                    "components": " + ".join(components),
                    "test_auc": roc_auc_score(y_test, test_score),
                    "w_esm2_150m": weight_by_component.get("ESM-2 150M", 0.0),
                    "w_esm2_650m": weight_by_component.get("ESM-2 650M", 0.0),
                    "w_esm1b_650m": weight_by_component.get("ESM-1b 650M", 0.0),
                    "w_calm": weight_by_component.get("CaLM", 0.0),
                    "n_test": len(test_idx),
                    "n_test_genes": len(np.unique(groups[test_idx])),
                    "n_test_pathogenic": int(y_test.sum()),
                }
            )

    fold_df = pd.DataFrame(rows)
    summary = (
        fold_df.groupby("model", sort=False)
        .agg(
            mean_auc=("test_auc", "mean"),
            sd_auc=("test_auc", "std"),
            mean_w_esm2_150m=("w_esm2_150m", "mean"),
            mean_w_esm2_650m=("w_esm2_650m", "mean"),
            mean_w_esm1b_650m=("w_esm1b_650m", "mean"),
            mean_w_calm=("w_calm", "mean"),
            n_folds=("fold", "count"),
        )
        .reset_index()
    )
    audit = pd.DataFrame(
        [
            {
                "all_variants": len(df_all),
                "all_genes": df_all["Gene_prot"].nunique(),
                "complete_case_variants": len(df),
                "complete_case_genes": df["Gene_prot"].nunique(),
                "dropped_variants_missing_scores": len(df_all) - len(df),
                "dropped_genes_missing_scores": df_all.loc[
                    df_all[score_cols].isna().any(axis=1), "Gene_prot"
                ].nunique(),
                "n_pathogenic_complete": int(df["label"].sum()),
                "n_benign_complete": int(len(df) - df["label"].sum()),
            }
        ]
    )

    fold_df.to_csv(out_dir / "model_control_gene_heldout_fold_results.csv", index=False)
    summary.to_csv(out_dir / "model_control_gene_heldout_summary.csv", index=False)
    paired_tests(fold_df).to_csv(out_dir / "model_control_paired_tests.csv", index=False)
    audit.to_csv(out_dir / "model_control_input_audit.csv", index=False)

    print(audit.to_string(index=False))
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(paired_tests(fold_df).to_string(index=False, float_format=lambda value: f"{value:.4g}"))


if __name__ == "__main__":
    main()
