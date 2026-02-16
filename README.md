# Coupling codon and protein constraints decouples drivers of variant pathogenicity

## Background
Predicting the functional impact of genetic variants remains a fundamental challenge in genomics. While existing models excel at identifying protein-intrinsic defects, they often overlook the regulatory syntax embedded within coding sequences. Here, we employ single-sequence codon (CaLM) and protein (ESM-2) language models as orthogonal probes to dissect the drivers of variant pathogenicity. Benchmarking against ClinVar variants demonstrated the nearly equal contributions of the two modalities. Crucially, dissecting the ClinMAVE dataset elucidated a mechanistic dichotomy: loss-of-function variants are predominantly governed by residue-level features, whereas gain-of-function variants and dosage-sensitive genes rely heavily on codon-level information that regulates protein biogenesis. These results highlight that pathogenicity stems from both the "product" and the "process," offering a biological framework to decode the mechanisms underlying disease.
<br> <img src="https://github.com/Cassie818/Viral-mut/blob/main/Figure/fig1.jpeg" width=850> <br>


## Contents

- [Notebooks/Scripts](#notebooksscripts)
- [Results](#results)
- [Dataset](#datasets)
- [Installation](#installation)
  - [CaLM](#1-install-calm)
  - [ESM-2](#2-install-esm-2)
- [References](#references)


## Notebooks/Scripts
- `cal_codon_logits.py`: code for computing logits from CaLM
- `cal_codon_scores.py`: code for computing effect scores at codon-level
- `cal_residue_logits.py`: code for computing logits from ESM-2
- `cal_residue_scores.py`: code for computing effect scores at residue-level
- `figx.ipynb`: code for all relevant analysis involving figx
- `get_fasta.py`: code for downloading cDNA and protein sequences from NCBI
- `go.R` code for GO enrichment analysis


## Results
- `missense`: contain all the prediction results for missense variants
- `nonsense`: contain all the prediction results for nonsense variants
- `synonymous`: contain all the prediction results for synonymous variants


## Datasets
All single-point mutations were downloaded from ClinVar (https://ftp.ncbi.nlm.nih.gov/pub/clinvar) and ClinMAVE (https://ngdc.cncb.ac.cn/clinmave/).


## Installation
To replicate our experiments, following the ESM-2 and CaLM first:
### 1. Install CaLM
```
pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu 
pip3 install einops
pip3 install rotary_embedding_torch
pip3 install biopython
python setup.py install 
```
### 2. Install ESM-2
```
pip install fair-esm  # latest release, OR:
pip install git+https://github.com/facebookresearch/esm.git  # bleeding edge, current repo main branch
```

## References
1. Meier, J., Rao, R., Verkuil, R., Liu, J., Sercu, T., Rives, A.: Language models enable zero-shot prediction of the effects of mutations on protein function. Advances in neural information processing systems 34, 29287–29303 (2021)
2. Outeiral, C., Deane, C.M.: Codon language embeddings provide strong signals for use in protein engineering. Nature Machine Intelligence 6(2), 170–179 (2024)
3. Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., et al.: Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379(6637), 1123–1130 (2023)




