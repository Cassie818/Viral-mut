#!/usr/bin/env python3
"""Supplementary full-view probability-space alignment control for Fig. 5."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import matplotlib.pyplot as plt


INPUT = Path("Results/Revision/len1022_aa_aggregation/len1022_calm_aa_aggregation_variant_scores.csv")
OUTPUT_DIR = Path("Results/ClinVar/substitution_discordance")
FIGURE = Path("Figure/figS7_discordance_probability_space_alignment.png")


def main() -> None:
    data = pd.read_csv(
        INPUT,
        usecols=["esm2_650m_llr", "calm_codon_llr_recomputed", "calm_aa_agg_llr"],
    ).dropna()
    codon = data["esm2_650m_llr"] - data["calm_codon_llr_recomputed"]
    aa_aggregated = data["esm2_650m_llr"] - data["calm_aa_agg_llr"]
    pearson = pearsonr(codon, aa_aggregated)
    spearman = spearmanr(codon, aa_aggregated)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "n_variants": len(data),
                "pearson_r": pearson.statistic,
                "pearson_p": pearson.pvalue,
                "spearman_rho": spearman.correlation,
                "spearman_p": spearman.pvalue,
            }
        ]
    ).to_csv(OUTPUT_DIR / "discordance_probability_space_correlations.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    figure, axis = plt.subplots(figsize=(5.1, 5.0))
    limit = float(np.ceil(np.nanpercentile(np.abs(np.r_[codon, aa_aggregated]), 99.7)))

    # Pale points retain the full cohort while making the below-diagonal shift visible.
    axis.scatter(codon, aa_aggregated, s=2.4, color="#8FB8C9", alpha=0.075, edgecolors="none", rasterized=True)
    axis.fill_between([-limit, limit], [-limit, limit], [-limit, -limit], color="#F7E6E2", alpha=0.65, zorder=0)
    axis.plot([-limit, limit], [-limit, limit], color="#575757", linewidth=0.85, linestyle=(0, (3, 2)), zorder=2)
    axis.axhline(0, color="#D7D7D7", linewidth=0.6, zorder=1)
    axis.axvline(0, color="#D7D7D7", linewidth=0.6, zorder=1)

    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Codon-space discordance")
    axis.set_ylabel("AA-space discordance")
    axis.text(
        0.04,
        0.96,
        f"Pearson $r$ = {pearson.statistic:.2f}\nSpearman $\\rho$ = {spearman.correlation:.2f}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        color="#3E3E3E",
    )
    axis.text(
        0.96,
        0.08,
        "Below diagonal:\naggregation reduces\ndiscordance magnitude",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#87645E",
    )
    axis.tick_params(width=0.8, length=3.2, color="#A0A0A0", labelcolor="#3E3E3E")
    for spine in axis.spines.values():
        spine.set_color("#A0A0A0")
        spine.set_linewidth(0.85)

    figure.tight_layout(pad=0.5)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE, dpi=600, bbox_inches="tight")
    print(FIGURE)


if __name__ == "__main__":
    main()
