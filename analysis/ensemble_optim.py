#!/usr/bin/env python3
"""Weight optimisers for convex score ensembles."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.stats import norm
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern
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


def bayes_optimize_weights(
    y: np.ndarray,
    scores: list[np.ndarray],
    *,
    seed: int,
    n_init: int = 10,
    n_iter: int = 20,
    candidate_step: float = 0.001,
    xi: float = 0.001,
) -> tuple[tuple[float, ...], float, int]:
    """Optimise AUROC with a Gaussian-process expected-improvement search.

    The search is performed over a dense finite pool on the convex simplex.
    Simplex vertices and the equal-weight point are always evaluated before
    random initialisation so boundary optima are not missed solely by chance.
    """

    candidates = _candidate_pool(len(scores), candidate_step)
    rng = np.random.default_rng(seed)

    anchors = [int(np.argmin(np.linalg.norm(candidates - target, axis=1))) for target in np.eye(len(scores))]
    anchors.append(
        int(
            np.argmin(
                np.linalg.norm(candidates - np.repeat(1.0 / len(scores), len(scores)), axis=1)
            )
        )
    )
    evaluated = list(dict.fromkeys(anchors))
    remaining = np.setdiff1d(np.arange(len(candidates)), evaluated, assume_unique=False)
    n_random = max(0, min(n_init, len(candidates)) - len(evaluated))
    if n_random:
        evaluated.extend(rng.choice(remaining, size=n_random, replace=False).tolist())

    values = [_objective(y, scores, candidates[idx]) for idx in evaluated]
    kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(
        length_scale=np.repeat(0.2, len(scores)),
        length_scale_bounds=(1e-3, 10.0),
        nu=2.5,
    )

    for iteration in range(n_iter):
        available = np.setdiff1d(np.arange(len(candidates)), evaluated, assume_unique=False)
        if len(available) == 0:
            break
        gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-7,
            normalize_y=True,
            random_state=seed + iteration,
            optimizer=None,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            gp.fit(candidates[evaluated], np.asarray(values))
        mean, std = gp.predict(candidates[available], return_std=True)
        improvement = mean - max(values) - xi
        z = np.divide(improvement, std, out=np.zeros_like(improvement), where=std > 0)
        expected_improvement = improvement * norm.cdf(z) + std * norm.pdf(z)
        expected_improvement[std <= 1e-12] = 0.0
        chosen = int(available[int(np.argmax(expected_improvement))])
        evaluated.append(chosen)
        values.append(_objective(y, scores, candidates[chosen]))

    best = int(np.argmax(values))
    return tuple(candidates[evaluated[best]]), float(values[best]), int(len(evaluated))
