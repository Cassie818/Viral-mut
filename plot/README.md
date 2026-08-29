# Figure scripts

- `main/`: one script per current main figure (`fig2`, `fig3`, `fig5`, `fig6`, and `fig7`)
- `supplementary/`: explicitly numbered scripts for Supplementary Figs. S1--S10

All scripts read final tables from `Results/` and write PNG files to `Figure/`.
They do not generate PDF duplicates.

Supplementary script mapping:

- `s1_s4.py`: model-control fold gains, cross-validation comparison, same-modality placeholder control, and ensemble weights
- `s5.py`: mutational-context conditional gains
- `s6_s7.py`: ClinMAVE GoF fold gains and optimised weights
- `s8.py`: amino-acid probability-alignment control
- `s9.py`: gene-level PLM-background controls
- `s10.py`: matched-platform PLM-background concordance
