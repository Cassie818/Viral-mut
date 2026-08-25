#!/usr/bin/env python3
"""Plot Fig. 7 experimental-context modality contribution."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


BASE = Path("Results/ClinMAVE/cross_platform_context")
PATHS = {
    "ESM-2 (150M)": BASE / "dataset_pair_summary_compact_esm2_150m.csv",
    "ESM-2 (650M)": BASE / "dataset_pair_summary_compact_esm2_650m.csv",
}
OUTS = [Path("Figure/fig7.png")]


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 9,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "axes.linewidth": 0.9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

COLORS = {
    "150": "#9CCDDD",
    "650": "#E7D27A",
    "pos": "#79B3A6",
    "neg": "#D99C91",
    "neutral": "#8E9390",
    "line": "#C9CECB",
    "grid": "#ECEBE7",
    "text": "#333333",
}


def load_data() -> pd.DataFrame:
    frames = []
    for model, path in PATHS.items():
        df = pd.read_csv(path)
        df = df[df["Gene"] != "BRCA1"].copy()
        df["model"] = model
        df["model_short"] = "150M" if "150M" in model else "650M"
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def style_ax(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color("#8F918E")
        ax.spines[spine].set_linewidth(0.9)
    ax.grid(False)
    ax.tick_params(width=0.85, colors=COLORS["text"], length=3)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    y = 1.11 if label in {"C", "D"} else 1.12
    x = -0.13 if label in {"B", "D"} else -0.20
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=11.5,
        fontweight="bold",
        ha="left",
        va="top",
        color="black",
    )


def paired_panel(ax: plt.Axes, df: pd.DataFrame, model: str, color: str) -> None:
    sub = df[df["model"] == model].copy()
    x = np.array([0, 1])
    for _, row in sub.iterrows():
        dms = row["mean_calm_weight_DMS"]
        cbge = row["mean_calm_weight_CBGE"]
        delta = row["delta_CBGE_minus_DMS"]
        line_color = COLORS["pos"] if delta >= 0 else COLORS["neg"]
        ax.plot(x, [dms, cbge], color=line_color, linewidth=0.8, alpha=0.52, zorder=1)
        ax.scatter(0, dms, s=24, facecolor=line_color, edgecolor="white", linewidth=0.45, alpha=0.9, zorder=2)
        ax.scatter(1, cbge, s=24, facecolor=line_color, edgecolor="white", linewidth=0.45, alpha=0.9, zorder=3)

    ax.set_xlim(-0.28, 1.28)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, 1], ["DMS", "CBGE"])
    ax.set_ylabel("Optimised CaLM weight")
    ax.set_title(model, fontsize=9.2, pad=5)
    style_ax(ax)


def gene_delta_panel(ax: plt.Axes, df: pd.DataFrame) -> None:
    gene = (
        df.groupby(["Gene", "model_short"], as_index=False)
        .agg(
            mean_delta=("delta_CBGE_minus_DMS", "mean"),
            n_pairs=("pair_id", "count"),
            mean_n=("n_variants", "mean"),
        )
        .pivot(index="Gene", columns="model_short", values=["mean_delta", "n_pairs", "mean_n"])
    )
    gene.columns = [f"{metric}_{model}" for metric, model in gene.columns]
    gene = gene.reset_index()
    gene["sort_delta"] = gene[["mean_delta_150M", "mean_delta_650M"]].mean(axis=1, skipna=True)
    gene = gene.sort_values("sort_delta", ascending=True).reset_index(drop=True)
    y = np.arange(len(gene))

    ax.axvline(0, color="#9C9C98", linewidth=1.1, linestyle=(0, (4, 3)), zorder=1)
    for idx, row in gene.iterrows():
        values = []
        if np.isfinite(row.get("mean_delta_150M", np.nan)):
            values.append(row["mean_delta_150M"])
        if np.isfinite(row.get("mean_delta_650M", np.nan)):
            values.append(row["mean_delta_650M"])
        if len(values) == 2:
            ax.plot(values, [idx, idx], color=COLORS["line"], linewidth=1.0, zorder=2)

    for model_short, color, offset in [("150M", COLORS["150"], -0.13), ("650M", COLORS["650"], 0.13)]:
        col = f"mean_delta_{model_short}"
        ncol = f"n_pairs_{model_short}"
        subset = gene[np.isfinite(gene[col])].copy()
        ax.scatter(
            subset[col],
            subset.index + offset,
            s=34,
            facecolor=color,
            edgecolor="#7E8784",
            linewidth=0.7,
            alpha=0.95,
            zorder=3,
            label=f"ESM-2 ({model_short})",
        )

    ax.set_yticks(y, gene["Gene"])
    ax.set_xlim(-0.34, 0.74)
    ax.set_ylim(-0.6, len(gene) - 0.4)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(frameon=False, loc="lower right", fontsize=7.3, handletextpad=0.25)
    style_ax(ax)


def baseline_sensitivity(ax: plt.Axes, df: pd.DataFrame) -> None:
    shared = (
        df.pivot_table(
            index=["pair_id", "Gene", "n_variants"],
            columns="model_short",
            values="delta_CBGE_minus_DMS",
            aggfunc="first",
        )
        .reset_index()
        .dropna(subset=["150M", "650M"])
    )
    ax.scatter(
        shared["150M"],
        shared["650M"],
        s=28,
        facecolor="#B9D8E1",
        edgecolor="#6F8F99",
        linewidth=0.8,
        alpha=0.86,
        zorder=3,
    )
    lim = (-0.36, 0.76)
    ax.plot(lim, lim, color="#8B8E8C", linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
    ax.axhline(0, color="#C9C9C4", linewidth=1.0, zorder=1)
    ax.axvline(0, color="#C9C9C4", linewidth=1.0, zorder=1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"$\Delta w$ with ESM-2 (150M)")
    ax.set_ylabel(r"$\Delta w$ with ESM-2 (650M)")
    style_ax(ax)
    pearson_r = shared["150M"].corr(shared["650M"], method="pearson")
    ax.text(0.05, 0.94, rf"$r={pearson_r:.2f}$", transform=ax.transAxes, ha="left", va="top", fontsize=8)



def main() -> None:
    df = load_data()
    fig = plt.figure(figsize=(6.10, 5.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.03], wspace=0.27, hspace=0.24)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    paired_panel(ax_a, df, "ESM-2 (150M)", COLORS["150"])
    paired_panel(ax_b, df, "ESM-2 (650M)", COLORS["650"])
    gene_delta_panel(ax_c, df)
    baseline_sensitivity(ax_d, df)

    for label, ax in zip("ABCD", [ax_a, ax_b, ax_c, ax_d]):
        add_panel_label(ax, label)

    fig.subplots_adjust(left=0.13, right=0.98, top=0.93, bottom=0.11)
    for out in OUTS:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=450)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
