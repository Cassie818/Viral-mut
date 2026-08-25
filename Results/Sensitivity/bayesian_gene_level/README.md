# Bayesian gene-level optimisation sensitivity analysis

This directory preserves the Bayesian-optimisation sensitivity analysis for
the ClinVar gene-level CaLM--PLM ensembles. It is separate from the primary
deterministic grid-search workflow and does not overwrite its outputs.

## Files

- `bayesian_gene_level_summary.csv`: Bayesian weights, AUROCs, and gene-level
  gain measures for 466 genes.
- `bayesian_gene_level_regressions.csv`: OLS regression and Spearman results.
- `bayesian_vs_grid_gene_comparison.csv`: gene-by-gene Bayesian-versus-grid
  comparison.

## Optimisation

The Bayesian search uses a Gaussian-process surrogate with a Matern kernel and
expected improvement over the convex weight simplex. Each optimisation uses
10 initial evaluations and 20 optimisation iterations. The random seed is 16
by default and is varied deterministically across genes and model combinations.

## Reproduction

Run from the repository root:

```bash
python analysis/sensitivity/bayes_clinvar_gene.py
```

The default input is the matched ClinVar model-score table used by the
gene-level analysis. The primary grid-search summary is used only for the
Bayesian-versus-grid comparison. Command-line options can override the input,
grid-summary, output directory, seed, and optimisation budget.

## Interpretation

Bayesian and grid search produce nearly identical cross-modal gains, whereas
exact optimised weights and weight-based gene rankings are less stable. These
outputs therefore support optimiser robustness of performance-based conclusions
but should not be used to treat an exact weight as an optimiser-independent
biological quantity.
