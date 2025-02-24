# Multilingual model improves zero-shot prediction of disease effects on proteins

## Background
Models for mutation effect prediction in coding sequences rely on sequence-, structure-, or homology-based features. Here, we introduce a novel method that combines a codon language model with a protein language model, providing a dual representation of disease effects. By capturing contextual dependencies at both the genetic and protein level, our approach achieves a 3% increase in ROC-AUC classifying disease effects for 137,350 ClinVar missense variants across 13,791 genes, outperforming two single-sequence-based language models. Obviously the codon language model can uniquely differentiate synonymous from nonsense mutations. Our strategy of using information at complementary biological scales (akin to human multilingual models) extends to protein fitness landscape modeling and evolutionary studies, with potential applications in precision medicine, protein engineering, and genomics.
<br> <img src="https://github.com/Cassie818/Viral-mut/blob/main/Figure/fig1.png" width=900> <br>

## Notebooks/Scripts
- `cal_codon_logits.py`: code for computing logits from CaLM
- `cal_codon_scores.py`: code for computing effect scores at codon-level
- `cal_residue_logits.py`: code for computing logits from ESM-2
- `cal_residue_scores.py`: code for computing effect scores at residue-level
- `codon_preference.ipynb`: code for identifying codon preference
- `get_fasta.py`: code for downloading sequences from NCBI


## Results
- `missense`: contain all the prediction results for missense variants
- `nonsense`: contain all the prediction results for nonsense variants
- `synonymous`: contain all the prediction results for synonymous variants
- `codon preference.csv`: codon preference results

## Installation
### 1. Install CaLM
```
! pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu \
! pip3 install einops \
! pip3 install rotary_embedding_torch \
! pip3 install biopython \
! python setup.py install 
```
### 2. Install ESM-2
```
pip install fair-esm  # latest release, OR:
pip install git+https://github.com/facebookresearch/esm.git  # bleeding edge, current repo main branch
```

## Reference
1. Meier, J., Rao, R., Verkuil, R., Liu, J., Sercu, T., Rives, A.: Language models enable zero-shot prediction of the effects of mutations on protein function. Advances in neural information processing systems 34, 29287–29303 (2021)
2. Outeiral, C., Deane, C.M.: Codon language embeddings provide strong signals for use in protein engineering. Nature Machine Intelligence 6(2), 170–179 (2024)
3. Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., et al.: Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379(6637), 1123–1130 (2023)




