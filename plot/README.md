# Figure scripts

- `main/`: one script per current main figure (`fig2`, `fig3`, `fig5`, `fig6`, and `fig7`)
- `supplementary/`: explicitly numbered scripts for Supplementary Figs. S1--S9

All scripts read final tables from `Results/` and write PNG files to `Figure/`.
They do not generate PDF duplicates.

Supplementary script mapping:

- `s1_s3.py`: model-control fold gains, cross-validation comparison, and ensemble weights
- `s4.py`: mutational-context conditional gains
- `s5_s6.py`: ClinMAVE GoF fold gains and optimised weights
- `s7.py`: amino-acid probability-alignment control
- `s8.py`: gene-level PLM-background controls
- `s9.py`: matched-platform PLM-background concordance
