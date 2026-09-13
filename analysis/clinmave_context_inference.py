#!/usr/bin/env python3
"""Gene-block inference for matched ClinMAVE DMS-CBGE weight shifts."""

from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path("Results/ClinMAVE/cross_platform_context")
PAIR_TABLE = BASE / "dataset_pair_modality_weights_11_gene_pairs_table_standardized.csv"
OUTPUT = BASE / "dataset_pair_modality_weight_gene_block_inference.csv"
BOOTSTRAP_REPLICATES = 10_000
RANDOM_SEED = 16


def gene_blocks(data: pd.DataFrame, value_column: str) -> list[np.ndarray]:
    return [
        group[value_column].dropna().to_numpy(float)
        for _, group in data.groupby("Gene", sort=True)
    ]


def cluster_bootstrap_ci(
    blocks: list[np.ndarray],
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = RANDOM_SEED,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates, dtype=float)
    n_genes = len(blocks)
    for index in range(replicates):
        sampled = rng.integers(0, n_genes, size=n_genes)
        estimates[index] = np.median(np.concatenate([blocks[i] for i in sampled]))
    low, high = np.percentile(estimates, [2.5, 97.5])
    return float(low), float(high)


def exact_gene_sign_flip_p(blocks: list[np.ndarray]) -> tuple[float, int]:
    observed = float(np.median(np.concatenate(blocks)))
    null_statistics = np.fromiter(
        (
            np.median(
                np.concatenate([sign * block for sign, block in zip(signs, blocks)])
            )
            for signs in product((-1.0, 1.0), repeat=len(blocks))
        ),
        dtype=float,
        count=2 ** len(blocks),
    )
    tolerance = np.finfo(float).eps * max(1.0, abs(observed)) * 8
    p_value = np.mean(np.abs(null_statistics) >= abs(observed) - tolerance)
    return float(p_value), int(len(null_statistics))


def summarize(
    data: pd.DataFrame,
    value_column: str,
    analysis: str,
    role: str,
) -> dict[str, object]:
    subset = data[["Gene", "DMS_dataset", "CBGE_dataset", value_column]].dropna()
    blocks = gene_blocks(subset, value_column)
    values = np.concatenate(blocks)
    ci_low, ci_high = cluster_bootstrap_ci(blocks)
    p_value, n_patterns = exact_gene_sign_flip_p(blocks)
    return {
        "analysis": analysis,
        "analysis_role": role,
        "estimand": "median pair-level CBGE minus DMS scale-adjusted CaLM weight",
        "n_weight_differences": int(len(values)),
        "n_dataset_pairs": int(
            subset[["Gene", "DMS_dataset", "CBGE_dataset"]].drop_duplicates().shape[0]
        ),
        "n_genes": int(len(blocks)),
        "n_positive_differences": int((values > 0).sum()),
        "median_delta_weight_standardized": float(np.median(values)),
        "gene_cluster_bootstrap_95ci_low": ci_low,
        "gene_cluster_bootstrap_95ci_high": ci_high,
        "exact_gene_block_sign_flip_p_two_sided": p_value,
        "sign_flip_patterns": n_patterns,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "random_seed": RANDOM_SEED,
    }


def main() -> None:
    data = pd.read_csv(PAIR_TABLE)
    rows = [
        summarize(
            data,
            "delta_weight_standardized_150M",
            "ESM-2 (150M)",
            "background-specific sensitivity analysis",
        ),
        summarize(
            data,
            "delta_weight_standardized_650M",
            "ESM-2 (650M)",
            "background-specific sensitivity analysis",
        ),
    ]
    combined = data[["Gene", "DMS_dataset", "CBGE_dataset"]].copy()
    combined["delta_weight_standardized_combined"] = (
        data["delta_weight_standardized_150M"]
        + data["delta_weight_standardized_650M"]
    ) / 2.0
    rows.insert(
        0,
        summarize(
            combined,
            "delta_weight_standardized_combined",
            "Combined PLM backgrounds",
            "combined-background analysis",
        ),
    )
    result = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False)
    print(result.to_string(index=False))
    print(f"\nWrote {OUTPUT}")


if __name__ == "__main__":
    main()
