#!/usr/bin/env python3
"""Pair-level concordance of platform shifts across ESM-2 backgrounds."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


DATA = Path(
    "Results/ClinMAVE/cross_platform_context/"
    "dataset_pair_modality_weights_150m_vs_650m_comparison.csv"
)
OUT = Path("Figure/figS_platform_baseline_concordance.png")

POINT = "#B9DDE9"
EDGE = "#6F8F99"
TEXT = "#30302E"
AXIS = "#8D8E8A"


def p_text(p: float) -> str:
    if p >= 0.001:
        return rf"$p={p:.3f}$"
    exponent = int(np.floor(np.log10(p)))
    mantissa = p / 10**exponent
    return rf"$p={mantissa:.1f}\times10^{{{exponent}}}$"


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    df = pd.read_csv(DATA)
    df = df[
        (df["Gene"] != "BRCA1")
        & df["delta_150m"].notna()
        & df["delta_650m"].notna()
    ].copy()
    rho, pval = spearmanr(df["delta_150m"], df["delta_650m"])

    fig, ax = plt.subplots(figsize=(4.35, 4.05))
    sizes = 22 + np.clip(df["n_variants"].to_numpy(float), 0, 250) * 0.16
    ax.scatter(
        df["delta_150m"],
        df["delta_650m"],
        s=sizes,
        facecolor=POINT,
        edgecolor=EDGE,
        linewidth=0.8,
        alpha=0.88,
        zorder=3,
    )

    lim = (-0.36, 0.76)
    ax.plot(lim, lim, color="#767773", linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="#C6C6C1", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#C6C6C1", linewidth=0.8, zorder=1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"CBGE$-$DMS CaLM-weight shift with ESM-2 (150M)")
    ax.set_ylabel(r"CBGE$-$DMS CaLM-weight shift with ESM-2 (650M)")

    ax.text(
        0.04,
        0.96,
        rf"Spearman $\rho={rho:.2f}$" + "\n" + p_text(float(pval)) + f"\n$n={len(df)}$ pairs",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color=TEXT,
        bbox={"facecolor": (1, 1, 1, 0.86), "edgecolor": "none", "pad": 2.2},
        zorder=5,
    )

    for spine in ax.spines.values():
        spine.set_color(AXIS)
        spine.set_linewidth(1.0)
    ax.tick_params(color=AXIS, labelcolor=TEXT, width=0.9, length=3)
    fig.subplots_adjust(left=0.19, right=0.97, top=0.97, bottom=0.16)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=600, bbox_inches="tight")
    print(OUT)


if __name__ == "__main__":
    main()
