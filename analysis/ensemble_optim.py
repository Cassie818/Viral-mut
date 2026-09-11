#!/usr/bin/env python3
"""Weight optimisers for convex score ensembles."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score


def _objective(
    y: np.ndarray,
    scores: list[np.ndarray],
    weights: np.ndarray,
) -> float:
    mixed = sum(float(w) * score for w, score in zip(weights, scores))
    return float(roc_auc_score(y, mixed))


def _candidate_pool(n_components: int, step: float) -> np.ndarray:
    if n_components == 2:
        second = np.arange(0.0, 1.0 + step / 2.0, step)
        return np.column_stack([1.0 - second, second])
    if n_components == 3:
        values = np.arange(0.0, 1.0 + step / 2.0, step)
        candidates = []
        for first in values:
            for second in values:
                third = 1.0 - first - second
                if third >= -1e-9:
                    candidates.append((first, second, max(0.0, third)))
        return np.asarray(candidates, dtype=float)
    raise ValueError(f"Only two- and three-model ensembles are supported, got {n_components}")


def grid_optimize_weights(
    y: np.ndarray,
    scores: list[np.ndarray],
    step: float,
) -> tuple[tuple[float, ...], float, int]:
    candidates = _candidate_pool(len(scores), step)
    values = np.asarray([_objective(y, scores, w) for w in candidates])
    best = int(np.argmax(values))
    return tuple(candidates[best]), float(values[best]), int(len(candidates))

