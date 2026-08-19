# Current ClinVar analysis files

This directory contains the current ClinVar files.

## Subdirectories

- `current_missense_cohort/`: current length-compatible ClinVar missense cohort,
  cohort counts, and model-score complete cases used for the main model-control
  analyses.
- `model_control/`: PLM/CaLM model-combination summaries, fold metrics, paired
  tests, ensemble weights, and supplementary model-control statistics.
- `context_control/`: fixed ESM-2 (650M) mutational-context control results,
  including fold-level AUROCs and all prespecified paired comparisons.
- `substitution_discordance/`: CaLM amino-acid aggregation, degeneracy, and
  substitution enrichment results used for Fig. 5.
- `gene_level/`: gene-level codon contribution and top-gene analyses used for
  Fig. 6.
