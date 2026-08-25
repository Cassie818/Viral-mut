#!/usr/bin/env python3
"""Create supplementary tables and a fold-level plot for the context control."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {
    "context": "#BFD8A8",
    "calm": "#F1B9BE",
    "text": "#363636",
    "edge": "#8F8F8A",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        default="Results/ClinVar/context_control",
    )
    parser.add_argument(
        "--table-dir",
        default="Results/ClinVar/context_control",
    )
    parser.add_argument("--figure-dir", default="Figure")
    return parser.parse_args()


def context_covariates() -> pd.DataFrame:
    rows = [
        ("Reference codon", "Categorical", "Wild-type DNA codon (U converted to T)."),
        ("Mutant codon", "Categorical", "Mutant DNA codon (U converted to T)."),
        ("Mutated codon position", "Categorical", "Position 1, 2, or 3 within the codon."),
        ("Reference nucleotide", "Categorical", "Wild-type nucleotide at the mutated position."),
        ("Mutant nucleotide", "Categorical", "Mutant nucleotide at the mutated position."),
        ("Nucleotide substitution", "Categorical", "Directed substitution, for example C>T."),
        ("Transition status", "Numeric", "1 for A>G, G>A, C>T, or T>C; otherwise 0."),
        ("Transversion status", "Numeric", "1 minus transition status."),
        ("Reference codon GC content", "Numeric", "Fraction of G or C nucleotides in the reference codon."),
        ("Mutant codon GC content", "Numeric", "Fraction of G or C nucleotides in the mutant codon."),
        ("Change in codon GC content", "Numeric", "Mutant minus reference codon GC content."),
        ("Reference within-codon CpG status", "Numeric", "1 if the reference codon contains CG; otherwise 0."),
        ("Mutant within-codon CpG status", "Numeric", "1 if the mutant codon contains CG; otherwise 0."),
        ("Change in within-codon CpG status", "Numeric", "Mutant minus reference within-codon CpG status."),
        (
            "Local codon degeneracy",
            "Categorical",
            "Number of possible nucleotides at the mutated position that preserve the reference amino acid.",
        ),
    ]
    return pd.DataFrame(rows, columns=["covariate", "encoding", "definition"])


def write_latex_covariate_table(table: pd.DataFrame, output_path: Path) -> None:
    lines = [
        "% Requires \\usepackage{booktabs,longtable}",
        "\\begin{longtable}{p{0.24\\textwidth} p{0.13\\textwidth} p{0.53\\textwidth}}",
        "\\caption{Mutational-context covariates used in the context-control analysis.}\\label{tab:context_covariates}\\\\",
        "\\toprule",
        "Covariate & Encoding & Definition \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        "Covariate & Encoding & Definition \\\\",
        "\\midrule",
        "\\endhead",
    ]
    for record in table.itertuples(index=False):
        lines.append(f"{record.covariate} & {record.encoding} & {record.definition} \\\\")
    lines.extend(["\\bottomrule", "\\end{longtable}", ""])
    output_path.write_text("\n".join(lines), encoding="utf-8")


def plot_conditional_deltas(paired: pd.DataFrame, figure_path: Path) -> None:
    contrasts = [
        (
            "CaLM beyond mutational context",
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + mutational context",
            COLORS["calm"],
        ),
        (
            "Mutational context beyond CaLM",
            "ESM-2 650M + mutational context + CaLM",
            "ESM-2 650M + CaLM",
            COLORS["context"],
        ),
    ]
    rows = []
    for label, a, b, color in contrasts:
        record = paired[(paired["a"] == a) & (paired["b"] == b)].iloc[0]
        for fold, value in enumerate(str(record["fold_deltas"]).split(";"), start=1):
            rows.append({"contrast": label, "fold": fold, "delta": float(value), "color": color})
    values = pd.DataFrame(rows)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, ax = plt.subplots(figsize=(5.2, 2.45))
    rng = np.random.default_rng(16)
    for index, (label, _, _, color) in enumerate(contrasts):
        subset = values[values["contrast"] == label]
        y = np.full(len(subset), index, dtype=float) + rng.normal(0, 0.055, len(subset))
        x = subset["delta"].to_numpy(float)
        ax.scatter(x, y, s=24, color=color, edgecolor=COLORS["edge"], linewidth=0.45, zorder=3)
        ax.errorbar(
            x.mean(),
            index,
            xerr=x.std(ddof=1),
            fmt="o",
            color=COLORS["text"],
            markersize=3.5,
            capsize=3,
            linewidth=0.9,
            zorder=4,
        )
    ax.axvline(0, color="#A8A8A3", linestyle=(0, (3, 2)), linewidth=0.9)
    ax.set_yticks(range(len(contrasts)))
    ax.set_yticklabels([item[0] for item in contrasts])
    ax.invert_yaxis()
    ax.set_xlabel("Fold-level conditional $\\Delta$AUROC")
    ax.tick_params(color=COLORS["edge"], labelcolor=COLORS["text"], width=0.8, length=3)
    for spine in ax.spines.values():
        spine.set_color(COLORS["edge"])
        spine.set_linewidth(0.9)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    table_dir = Path(args.table_dir)
    figure_dir = Path(args.figure_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    folds = pd.read_csv(input_dir / "context_control_fold_results.csv")
    paired = pd.read_csv(input_dir / "context_control_paired_tests.csv")

    write_latex_covariate_table(
        context_covariates(),
        table_dir / "Supplementary_Table_S4_mutational_context_covariates.tex",
    )
    plot_conditional_deltas(paired, figure_dir / "figS4.png")


if __name__ == "__main__":
    main()
