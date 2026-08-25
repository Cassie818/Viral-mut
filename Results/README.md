# Results

Final manuscript result tables are organised by dataset:

- `ClinVar/`: the current matched missense cohort and results for model control,
  mutational context, substitution discordance, and gene-level analyses.
- `ClinMAVE/`: results for functional-class, nonsense/synonymous, and matched
  cross-platform analyses.
- `Sensitivity/`: results from prespecified optimisation-sensitivity analyses.

Each dataset directory contains a short README describing its analysis-level
subdirectories. All manuscript figure scripts read inputs from the public
`Results/ClinVar/` and `Results/ClinMAVE/` directories.

`Revision/` is an ignored local workspace containing large model caches,
variant-level score tables, and upstream intermediate files. It is retained so
that expensive model-scoring steps do not need to be repeated, but it is not
part of the public GitHub result set.
