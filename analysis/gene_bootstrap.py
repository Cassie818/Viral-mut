#!/usr/bin/env python3
"""Gene-cluster bootstrap inference for pooled out-of-fold predictions."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


CLINVAR_MODEL_OOF = Path("Results/ClinVar/model_control/model_control_oof_predictions.csv.gz")
CLINVAR_CONTEXT_OOF = Path("Results/ClinVar/context_control/context_control_oof_predictions.csv.gz")
CLINMAVE_ESM2_OOF = Path(
    "Results/ClinMAVE/functional_effects/clinmave_esm2_650m_calm_oof_predictions.csv.gz"
)
CLINMAVE_ESM1B_OOF = Path(
    "Results/ClinMAVE/functional_effects/clinmave_esm1b_650m_calm_oof_predictions.csv.gz"
)

MODEL_COMPARISONS = [
    (
        "CaLM vs ESM-2 650M",
        "score_calm",
        "score_esm2_650m",
    ),
    (
        "ESM-2 150M + ESM-2 650M vs ESM-2 650M",
        "score_esm2_150m_esm2_650m",
        "score_esm2_650m",
    ),
    (
        "ESM-2 650M + ESM-1b 650M vs ESM-2 650M",
        "score_esm2_650m_esm1b_650m",
        "score_esm2_650m",
    ),
    (
        "ESM-2 650M + ESM-1b 650M vs ESM-1b 650M",
        "score_esm2_650m_esm1b_650m",
        "score_esm1b_650m",
    ),
    (
        "ESM-2 650M + CaLM vs ESM-2 650M",
        "score_esm2_650m_calm",
        "score_esm2_650m",
    ),
    (
        "ESM-2 150M + CaLM vs ESM-2 150M",
        "score_esm2_150m_calm",
        "score_esm2_150m",
    ),
    (
        "ESM-1b 650M + CaLM vs ESM-1b 650M",
        "score_esm1b_650m_calm",
        "score_esm1b_650m",
    ),
    (
        "ESM-2 650M + ESM-1b 650M + CaLM vs ESM-2 650M + ESM-1b 650M",
        "score_esm2_650m_esm1b_650m_calm",
        "score_esm2_650m_esm1b_650m",
    ),
    (
        "ESM-2 650M + ESM-1b 650M + CaLM vs ESM-2 650M + CaLM",
        "score_esm2_650m_esm1b_650m_calm",
        "score_esm2_650m_calm",
    ),
]

CONTEXT_COMPARISONS = [
    (
        "ESM-2 650M + mutational context vs ESM-2 650M",
        "score_esm2_650m_context",
        "score_esm2_650m",
    ),
    (
        "ESM-2 650M + CaLM vs ESM-2 650M",
        "score_esm2_650m_calm",
        "score_esm2_650m",
    ),
    (
        "ESM-2 650M + mutational context + CaLM vs ESM-2 650M + mutational context",
        "score_esm2_650m_context_calm",
        "score_esm2_650m_context",
    ),
    (
        "ESM-2 650M + mutational context + CaLM vs ESM-2 650M + CaLM",
        "score_esm2_650m_context_calm",
        "score_esm2_650m_calm",
    ),
    (
        "ESM-2 650M + mutational context + CaLM vs ESM-2 650M",
        "score_esm2_650m_context_calm",
        "score_esm2_650m",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument(
        "--scope", choices=["all", "clinvar", "clinmave"], default="all"
    )
    return parser.parse_args()


def prepare_auc(y: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pre-sort one score vector and retain tie-group starts for weighted AUC."""
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_score) != 0) + 1]
    return order, starts, y[order].astype(float)


def weighted_auc_batch(
    cluster_counts: np.ndarray,
    variant_cluster: np.ndarray,
    prepared: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Calculate tie-aware weighted AUROC for a batch of cluster resamples."""
    order, starts, y_sorted = prepared
    weights = cluster_counts[:, variant_cluster[order]].astype(float, copy=False)
    positive = np.add.reduceat(weights * y_sorted, starts, axis=1)
    negative = np.add.reduceat(weights * (1.0 - y_sorted), starts, axis=1)
    negative_before = np.cumsum(negative, axis=1) - negative
    numerator = np.sum(positive * (negative_before + 0.5 * negative), axis=1)
    total_positive = positive.sum(axis=1)
    total_negative = negative.sum(axis=1)
    denominator = total_positive * total_negative
    return np.divide(
        numerator,
        denominator,
        out=np.full(len(cluster_counts), np.nan, dtype=float),
        where=denominator > 0,
    )


def bootstrap_auc_columns(
    df: pd.DataFrame,
    score_columns: list[str],
    *,
    n_bootstrap: int,
    seed: int,
    batch_size: int,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    required = ["gene", "label", *score_columns]
    data = df.dropna(subset=required).copy()
    y = data["label"].to_numpy(int)
    if len(np.unique(y)) != 2:
        raise ValueError("Both outcome classes are required for AUROC")

    genes, variant_cluster = np.unique(data["gene"].astype(str), return_inverse=True)
    n_clusters = len(genes)
    probability = np.full(n_clusters, 1.0 / n_clusters)
    rng = np.random.default_rng(seed)

    point = {
        column: float(roc_auc_score(y, data[column].to_numpy(float)))
        for column in score_columns
    }
    prepared = {
        column: prepare_auc(y, data[column].to_numpy(float))
        for column in score_columns
    }
    bootstrap = {
        column: np.full(n_bootstrap, np.nan, dtype=float)
        for column in score_columns
    }

    for start in range(0, n_bootstrap, batch_size):
        stop = min(start + batch_size, n_bootstrap)
        counts = rng.multinomial(n_clusters, probability, size=stop - start)
        for column in score_columns:
            bootstrap[column][start:stop] = weighted_auc_batch(
                counts, variant_cluster, prepared[column]
            )

    return point, bootstrap


def summarise_comparisons(
    df: pd.DataFrame,
    comparisons: list[tuple[str, str, str]],
    *,
    analysis: str,
    n_bootstrap: int,
    seed: int,
    batch_size: int,
    metadata: dict[str, object] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_columns = list(dict.fromkeys([column for _, a, b in comparisons for column in (a, b)]))
    point, bootstrap = bootstrap_auc_columns(
        df,
        score_columns,
        n_bootstrap=n_bootstrap,
        seed=seed,
        batch_size=batch_size,
    )
    score_rows = []
    for column in score_columns:
        valid = bootstrap[column][np.isfinite(bootstrap[column])]
        lower, upper = np.quantile(valid, [0.025, 0.975])
        score_row = {
            "analysis": analysis,
            "score_column": column,
            "n_variants": int(len(df)),
            "n_genes": int(df["gene"].nunique()),
            "n_bootstrap_requested": int(n_bootstrap),
            "n_bootstrap_valid": int(len(valid)),
            "pooled_oof_auroc": point[column],
            "bootstrap_mean_auroc": float(np.mean(valid)),
            "bootstrap_sd_auroc": float(np.std(valid, ddof=1)),
            "bootstrap_95ci_low": float(lower),
            "bootstrap_95ci_high": float(upper),
            "resampling_unit": "gene",
            "ci_method": "percentile cluster bootstrap",
        }
        if metadata:
            score_row.update(metadata)
        score_rows.append(score_row)

    rows = []
    distributions = []
    for comparison, model_a, model_b in comparisons:
        delta = bootstrap[model_a] - bootstrap[model_b]
        valid = delta[np.isfinite(delta)]
        lower, upper = np.quantile(valid, [0.025, 0.975])
        lower_tail = (np.count_nonzero(valid <= 0) + 1) / (len(valid) + 1)
        upper_tail = (np.count_nonzero(valid >= 0) + 1) / (len(valid) + 1)
        row = {
            "analysis": analysis,
            "comparison": comparison,
            "model_a_score": model_a,
            "model_b_score": model_b,
            "n_variants": int(len(df)),
            "n_genes": int(df["gene"].nunique()),
            "n_bootstrap_requested": int(n_bootstrap),
            "n_bootstrap_valid": int(len(valid)),
            "pooled_oof_auroc_a": point[model_a],
            "pooled_oof_auroc_b": point[model_b],
            "pooled_oof_delta_auroc": point[model_a] - point[model_b],
            "bootstrap_mean_delta_auroc": float(np.mean(valid)),
            "bootstrap_sd_delta_auroc": float(np.std(valid, ddof=1)),
            "bootstrap_95ci_low": float(lower),
            "bootstrap_95ci_high": float(upper),
            "bootstrap_two_sided_tail_p": float(min(1.0, 2.0 * min(lower_tail, upper_tail))),
            "ci_excludes_zero": bool(lower > 0 or upper < 0),
            "resampling_unit": "gene",
            "ci_method": "percentile cluster bootstrap",
        }
        if metadata:
            row.update(metadata)
        rows.append(row)
        dist = pd.DataFrame(
            {
                "analysis": analysis,
                "comparison": comparison,
                "bootstrap_replicate": np.arange(1, n_bootstrap + 1),
                "delta_auroc": delta,
            }
        )
        if metadata:
            for key, value in metadata.items():
                dist[key] = value
        distributions.append(dist)
    return (
        pd.DataFrame(rows),
        pd.concat(distributions, ignore_index=True),
        pd.DataFrame(score_rows),
    )


def run_clinvar(
    *, n_bootstrap: int, seed: int, batch_size: int
) -> tuple[list[pd.DataFrame], list[pd.DataFrame], list[pd.DataFrame]]:
    summaries, distributions, score_summaries = [], [], []
    model = pd.read_csv(CLINVAR_MODEL_OOF)
    summary, dist, score_summary = summarise_comparisons(
        model,
        MODEL_COMPARISONS,
        analysis="ClinVar model combinations",
        n_bootstrap=n_bootstrap,
        seed=seed,
        batch_size=batch_size,
    )
    summaries.append(summary)
    distributions.append(dist)
    score_summaries.append(score_summary)

    context = pd.read_csv(CLINVAR_CONTEXT_OOF)
    summary, dist, score_summary = summarise_comparisons(
        context,
        CONTEXT_COMPARISONS,
        analysis="ClinVar mutational-context control",
        n_bootstrap=n_bootstrap,
        seed=seed + 1,
        batch_size=batch_size,
    )
    summaries.append(summary)
    distributions.append(dist)
    score_summaries.append(score_summary)
    return summaries, distributions, score_summaries


def run_clinmave(
    *, n_bootstrap: int, seed: int, batch_size: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries, distributions = [], []
    for offset, (background, path) in enumerate(
        [("ESM-2 (650M)", CLINMAVE_ESM2_OOF), ("ESM-1b (650M)", CLINMAVE_ESM1B_OOF)]
    ):
        data = pd.read_csv(path)
        for task_index, ((assay, case_class), task) in enumerate(
            data.groupby(["assay", "case_class"], sort=True)
        ):
            comparison = f"{background} + CaLM vs {background}"
            summary, dist, _ = summarise_comparisons(
                task,
                [(comparison, "score_ensemble", "score_plm")],
                analysis="ClinMAVE functional-class comparison",
                n_bootstrap=n_bootstrap,
                seed=seed + 10 + 10 * offset + task_index,
                batch_size=batch_size,
                metadata={
                    "plm_background": background,
                    "assay": assay,
                    "case_class": case_class,
                },
            )
            summaries.append(summary)
            distributions.append(dist)
    return pd.concat(summaries, ignore_index=True), pd.concat(distributions, ignore_index=True)


def main() -> None:
    args = parse_args()
    if args.scope in {"all", "clinvar"}:
        clinvar_summaries, clinvar_distributions, clinvar_score_summaries = run_clinvar(
            n_bootstrap=args.n_bootstrap,
            seed=args.seed,
            batch_size=args.batch_size,
        )
        model_summary, context_summary = clinvar_summaries
        model_dist, context_dist = clinvar_distributions
        model_score_summary, context_score_summary = clinvar_score_summaries
        model_summary.to_csv(
            "Results/ClinVar/model_control/model_control_gene_bootstrap.csv", index=False
        )
        model_dist.to_csv(
            "Results/ClinVar/model_control/model_control_gene_bootstrap_distributions.csv.gz",
            index=False,
            compression="gzip",
        )
        model_score_summary.to_csv(
            "Results/ClinVar/model_control/model_control_gene_bootstrap_aurocs.csv",
            index=False,
        )
        context_summary.to_csv(
            "Results/ClinVar/context_control/context_control_gene_bootstrap.csv", index=False
        )
        context_dist.to_csv(
            "Results/ClinVar/context_control/context_control_gene_bootstrap_distributions.csv.gz",
            index=False,
            compression="gzip",
        )
        context_score_summary.to_csv(
            "Results/ClinVar/context_control/context_control_gene_bootstrap_aurocs.csv",
            index=False,
        )
        print("\nClinVar model combinations")
        print(model_summary.to_string(index=False, float_format=lambda value: f"{value:.6g}"))
        print("\nClinVar mutational-context control")
        print(context_summary.to_string(index=False, float_format=lambda value: f"{value:.6g}"))

    if args.scope in {"all", "clinmave"}:
        clinmave_summary, clinmave_dist = run_clinmave(
            n_bootstrap=args.n_bootstrap,
            seed=args.seed,
            batch_size=args.batch_size,
        )
        clinmave_summary.to_csv(
            "Results/ClinMAVE/functional_effects/clinmave_functional_class_gene_bootstrap.csv",
            index=False,
        )
        clinmave_dist.to_csv(
            "Results/ClinMAVE/functional_effects/clinmave_functional_class_gene_bootstrap_distributions.csv.gz",
            index=False,
            compression="gzip",
        )
        print("\nClinMAVE functional classes")
        print(clinmave_summary.to_string(index=False, float_format=lambda value: f"{value:.6g}"))


if __name__ == "__main__":
    main()
