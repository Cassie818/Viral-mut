# Codon-level information complements protein language models for missense variant effect prediction

The project compares codon-level language-model scores from CaLM with residue-level protein language-model scores from ESM-2 and ESM-1b to test whether codon-resolution information contributes to missense-variant interpretation beyond protein-level constraints.

<br>
<img src="Figure/fig1.jpeg" width="850">
<br>

## Overview

The main analyses use ClinVar and ClinMAVE variant datasets to evaluate:

- whether CaLM contributes signal beyond PLM baselines and PLM-PLM ensembles
- whether the CaLM contribution persists after explicit mutational-context controls
- how modality contribution varies across functional classes in ClinMAVE
- whether codon-level discordance remains after aggregating CaLM probabilities into amino-acid space
- whether gene-level CaLM contribution varies across genes and PLM backgrounds
- whether inferred modality contribution differs across matched DMS and CBGE assay contexts

## Repository Structure

```text
Figure/                 Main and supplementary manuscript figures
Results/Main_Analysis/  Result tables used in the manuscript
bin/                    Analysis and scoring scripts
plot/                   Figure generation scripts
```
## Key Scripts

ClinVar model and context analyses:

- `bin/cv_model_control.py`: core ClinVar PLM/CaLM model-combination analysis for Fig. 2A-B
- `bin/cv_context_control.py`: ESM-2 (650M), mutational-context, and CaLM comparison for Fig. 2C-D
- `bin/cv_aa_agg.py`: CaLM amino-acid aggregation and substitution-discordance analyses for Fig. 5
- `bin/cv_gene_codon.py`: gene-level codon contribution analyses for Fig. 6
- `bin/cv_gene_controls.py`: ClinVar gene-level source tables and sequence-confounder controls

ClinMAVE analyses:

- `bin/cm_function_esm2_650m.py`: ClinMAVE functional-class analyses using ESM-2 (650M) and CaLM
- `bin/cm_function_esm1b_650m.py`: ClinMAVE functional-class analyses using ESM-1b (650M) and CaLM
- `bin/cm_fetch_pairs.py`: ClinMAVE DMS-CBGE dataset-pair retrieval
- `bin/cm_pair_context_650m.py`: matched DMS-CBGE dataset-pair analysis using ESM-2 (650M) and CaLM
- `bin/cm_pair_context_150m.py`: matched DMS-CBGE sensitivity analysis using ESM-2 (150M) and CaLM

Scoring utilities:

- `bin/score_calm_codon_logits.py` and `bin/score_calm_codon_llr.py`: CaLM codon-level scoring
- `bin/score_plm_residue_logits.py` and `bin/score_plm_residue_llr.py`: generic PLM residue-level scoring
- `bin/score_cv_esm2_650m.py`: ClinVar ESM-2 (650M) scoring
- `bin/score_cm_esm2_650m.py`: ClinMAVE ESM-2 (650M) scoring
- `bin/score_cm_esm1b_650m.py`: ClinMAVE ESM-1b (650M) scoring

## Dependencies

The analysis scripts are Python-based and use common scientific packages, including:

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `statsmodels`
- `matplotlib`
- `seaborn`
- `torch`
- `fair-esm`

CaLM should be installed from the original CaLM repository:

```bash
pip install torch torchvision torchaudio
pip install einops rotary_embedding_torch biopython
python setup.py install
```

ESM models can be installed using:

```bash
pip install fair-esm
```

## Data Sources

Variant annotations were obtained from ClinVar and ClinMAVE:

- ClinVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar
- ClinMAVE: https://ngdc.cncb.ac.cn/clinmave/

Language-model resources:

- CaLM: https://github.com/oxpig/CaLM
- ESM: https://github.com/facebookresearch/esm

## References

1. Rives, Alexander, et al. "Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences." Proceedings of the national academy of sciences 118.15 (2021): e2016239118.
2. Lin, Zeming, et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." Science 379.6637 (2023): 1123-1130.
3. Outeiral, Carlos, and Charlotte M. Deane. "Codon language embeddings provide strong signals for use in protein engineering." Nature Machine Intelligence 6.2 (2024): 170-179. 
4. Meier, Joshua, et al. "Language models enable zero-shot prediction of the effects of mutations on protein function." Advances in neural information processing systems 34 (2021): 29287-29303.
