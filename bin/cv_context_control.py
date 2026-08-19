#!/usr/bin/env python3
"""Fixed ESM-2 650M baseline with mutational-context and CaLM increments."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, wilcoxon
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PATHOGENIC_LABELS = {"pathogenic", "likely_pathogenic"}
BASES = ("A", "C", "G", "T")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scores",
        default=(
            "Results/Revision/len1022_model_control/"
            "len1022_model_control_score_table_all_variants.csv"
        ),
    )
    parser.add_argument("--out-dir", default="Results/Revision/fig3_fixed_esm2_context_calm")
    parser.add_argument("--n-splits", type=int, default=10)
    parser.add_argument("--seed", type=int, default=16)
    return parser.parse_args()


def mutation_features(ref: str, mut: str) -> dict[str, object]:
    ref = str(ref).upper().replace("U", "T")
    mut = str(mut).upper().replace("U", "T")
    out: dict[str, object] = {
        "valid_snv_codon": 0,
        "ref_codon": ref,
        "mut_codon": mut,
        "ref_base": "NA",
        "mut_base": "NA",
        "substitution": "NA",
        "mutated_codon_position": "NA",
        "is_transition": np.nan,
        "is_transversion": np.nan,
        "ref_codon_gc": np.nan,
        "mut_codon_gc": np.nan,
        "delta_codon_gc": np.nan,
        "ref_cpg": np.nan,
        "mut_cpg": np.nan,
        "delta_cpg": np.nan,
        "local_degeneracy": "NA",
    }
    if len(ref) != 3 or len(mut) != 3 or any(base not in BASES for base in ref + mut):
        return out
    diffs = [idx for idx, (a, b) in enumerate(zip(ref, mut)) if a != b]
    if len(diffs) != 1:
        return out

    pos = diffs[0]
    ref_aa = CODON_TABLE.get(ref)
    syn_count = 0
    if ref_aa is not None:
        for base in BASES:
            candidate = ref[:pos] + base + ref[pos + 1 :]
            if CODON_TABLE.get(candidate) == ref_aa:
                syn_count += 1
    ref_gc = sum(base in {"G", "C"} for base in ref) / 3
    mut_gc = sum(base in {"G", "C"} for base in mut) / 3
    is_transition = int((ref[pos], mut[pos]) in TRANSITIONS)
    ref_cpg = int("CG" in ref)
    mut_cpg = int("CG" in mut)
    out.update(
        {
            "valid_snv_codon": 1,
            "ref_base": ref[pos],
            "mut_base": mut[pos],
            "substitution": f"{ref[pos]}>{mut[pos]}",
            "mutated_codon_position": str(pos + 1),
            "is_transition": is_transition,
            "is_transversion": 1 - is_transition,
            "ref_codon_gc": ref_gc,
            "mut_codon_gc": mut_gc,
            "delta_codon_gc": mut_gc - ref_gc,
            "ref_cpg": ref_cpg,
            "mut_cpg": mut_cpg,
            "delta_cpg": mut_cpg - ref_cpg,
            "local_degeneracy": str(syn_count),
        }
    )
    return out


def add_context_features(df: pd.DataFrame) -> pd.DataFrame:
    features = pd.DataFrame(
        [mutation_features(ref, mut) for ref, mut in zip(df["Ref_gene"], df["Mut_gene"])],
        index=df.index,
    )
    return pd.concat([df, features], axis=1)


def make_pipeline(numeric_cols: list[str], categorical_cols: list[str]) -> Pipeline:
    preprocessor = ColumnTransformer(
        transformers=[
            ("num", StandardScaler(), numeric_cols),
            ("cat", OneHotEncoder(handle_unknown="ignore", sparse=False), categorical_cols),
        ],
        remainder="drop",
    )
    return Pipeline(
        steps=[
            ("preprocess", preprocessor),
            (
                "model",
                LogisticRegression(
                    solver="lbfgs",
                    max_iter=1000,
                    penalty="l2",
                ),
            ),
        ]
    )


def paired_tests(fold_df: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("ESM-2 650M + mutational context", "ESM-2 650M"),
        ("ESM-2 650M + CaLM", "ESM-2 650M"),
        (
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + mutational context",
        ),
        (
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + CaLM",
        ),
        (
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M",
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

    df = pd.read_csv(args.scores)
    if "label" not in df.columns:
        df["label"] = df["Label_prot"].astype(str).isin(PATHOGENIC_LABELS).astype(int)
    df = add_context_features(df)

    numeric_context = [
        "valid_snv_codon",
        "is_transition",
        "is_transversion",
        "ref_codon_gc",
        "mut_codon_gc",
        "delta_codon_gc",
        "ref_cpg",
        "mut_cpg",
        "delta_cpg",
    ]
    categorical_context = [
        "ref_codon",
        "mut_codon",
        "ref_base",
        "mut_base",
        "substitution",
        "mutated_codon_position",
        "local_degeneracy",
    ]
    model_specs = [
        ("ESM-2 650M", ["esm2_650m_score"], []),
        ("ESM-2 650M + mutational context", ["esm2_650m_score", *numeric_context], categorical_context),
        ("ESM-2 650M + CaLM", ["esm2_650m_score", "calm_score"], []),
        (
            "ESM-2 650M + mutational context + CaLM",
            ["esm2_650m_score", "calm_score", *numeric_context],
            categorical_context,
        ),
    ]
    required = sorted(
        set(["label", "Gene_prot"] + [col for _, nums, cats in model_specs for col in nums + cats])
    )
    complete = df.dropna(subset=required).copy()
    complete.to_csv(out_dir / "fig3_fixed_esm2_context_calm_score_table.csv", index=False)

    y = complete["label"].to_numpy(dtype=int)
    groups = complete["Gene_prot"].astype(str).to_numpy()
    splits = list(
        StratifiedGroupKFold(n_splits=args.n_splits, shuffle=True, random_state=args.seed).split(
            complete, y, groups
        )
    )

    rows = []
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        train = complete.iloc[train_idx]
        test = complete.iloc[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]
        for model_name, numeric_cols, categorical_cols in model_specs:
            pipeline = make_pipeline(numeric_cols, categorical_cols)
            pipeline.fit(train[numeric_cols + categorical_cols], y_train)
            pred = pipeline.predict_proba(test[numeric_cols + categorical_cols])[:, 1]
            rows.append(
                {
                    "fold": fold,
                    "model": model_name,
                    "test_auc": roc_auc_score(y_test, pred),
                    "n_test": len(test),
                    "n_test_genes": test["Gene_prot"].nunique(),
                    "n_test_pathogenic": int(y_test.sum()),
                }
            )

    fold_df = pd.DataFrame(rows)
    summary = (
        fold_df.groupby("model", sort=False)
        .agg(
            mean_auc=("test_auc", "mean"),
            sd_auc=("test_auc", "std"),
            n_folds=("fold", "count"),
        )
        .reset_index()
    )
    audit = pd.DataFrame(
        [
            {
                "all_variants": len(df),
                "all_genes": df["Gene_prot"].nunique(),
                "complete_case_variants": len(complete),
                "complete_case_genes": complete["Gene_prot"].nunique(),
                "dropped_variants": len(df) - len(complete),
                "n_pathogenic_complete": int(complete["label"].sum()),
                "n_benign_complete": int(len(complete) - complete["label"].sum()),
            }
        ]
    )
    fold_df.to_csv(out_dir / "fig3_fixed_esm2_context_calm_fold_results.csv", index=False)
    summary.to_csv(out_dir / "fig3_fixed_esm2_context_calm_summary.csv", index=False)
    paired_tests(fold_df).to_csv(out_dir / "fig3_fixed_esm2_context_calm_paired_tests.csv", index=False)
    audit.to_csv(out_dir / "fig3_fixed_esm2_context_calm_input_audit.csv", index=False)

    print(audit.to_string(index=False))
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(paired_tests(fold_df).to_string(index=False, float_format=lambda value: f"{value:.4g}"))


if __name__ == "__main__":
    main()
