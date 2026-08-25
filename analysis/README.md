# Analysis scripts

Top-level ClinVar analyses:

- `clinvar_ensemble.py`: primary PLM/CaLM ensemble comparison
- `clinvar_context.py`: explicit mutational-context control
- `clinvar_discordance.py`: amino-acid aggregation and substitution discordance
- `clinvar_gene_prepare.py`: preparation of gene-level source tables
- `clinvar_gene.py`: gene-level codon-contribution analysis

Top-level ClinMAVE analyses:

- `clinmave_classes_esm2.py`: functional-class analysis with ESM-2
- `clinmave_classes_esm1b.py`: functional-class analysis with ESM-1b
- `clinmave_pairs_fetch.py`: retrieval of matched DMS--CBGE pairs
- `clinmave_pairs_esm2_150m.py`: matched-pair analysis with ESM-2 (150M)
- `clinmave_pairs_esm2_650m.py`: matched-pair analysis with ESM-2 (650M)

Subdirectories:

- `scoring/`: sequence preparation and model scoring
- `sensitivity/`: Bayesian-optimisation sensitivity analyses

Run scripts from the repository root so relative `Results/` paths resolve
consistently.
