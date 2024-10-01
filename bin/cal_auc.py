import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
import pickle


# Define a function to load data
def load_data(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

# Load the data
LLR_benign_gene = load_data('./data/CaLM_LLR_benign.txt')
LLR_Pathogenic_gene = load_data('./data/CaLM_LLR_Pathogenic.txt')
LLR_benign_prot = load_data('./data/ESM2_LLR_benign.txt')
LLR_Pathogenic_prot = load_data('./data/ESM2_LLR_Pathogenic.txt')

# **Assign labels: 1 for pathogenic samples, 0 for benign samples**
labels_benign = np.zeros(len(LLR_benign_gene))          # Label for benign samples is 0
labels_pathogenic = np.ones(len(LLR_Pathogenic_gene))   # Label for pathogenic samples is 1

# Combine labels
labels_gene = np.concatenate([labels_benign, labels_pathogenic])
labels_prot = np.concatenate([labels_benign, labels_pathogenic])
labels_combined = labels_gene

# Combine scores
scores_gene = np.concatenate([LLR_benign_gene, LLR_Pathogenic_gene])
scores_prot = np.concatenate([LLR_benign_prot, LLR_Pathogenic_prot])
scores_combined = 0.265 * scores_gene + 0.006 * (scores_prot ** 3) + 0.04 * (scores_gene * scores_prot)

# **If higher scores indicate benign samples, invert scores**
# This is to ensure that higher scores correspond to pathogenic samples (label 1)
scores_gene_inverted = -scores_gene
scores_prot_inverted = -scores_prot
scores_combined_inverted = -scores_combined

# Calculate ROC curve and AUC
fpr_gene, tpr_gene, _ = roc_curve(labels_gene, scores_gene_inverted)
auc_gene = roc_auc_score(labels_gene, scores_gene_inverted)

fpr_prot, tpr_prot, _ = roc_curve(labels_prot, scores_prot_inverted)
auc_prot = roc_auc_score(labels_prot, scores_prot_inverted)

fpr_combined, tpr_combined, _ = roc_curve(labels_combined, scores_combined_inverted)
auc_combined = roc_auc_score(labels_combined, scores_combined_inverted)

# Calculate Precision-Recall (PR) curve and AUC
precision_gene, recall_gene, _ = precision_recall_curve(labels_gene, scores_gene_inverted)
pr_auc_gene = auc(recall_gene, precision_gene)

precision_prot, recall_prot, _ = precision_recall_curve(labels_prot, scores_prot_inverted)
pr_auc_prot = auc(recall_prot, precision_prot)

precision_combined, recall_combined, _ = precision_recall_curve(labels_combined, scores_combined_inverted)
pr_auc_combined = auc(recall_combined, precision_combined)

# Plot ROC curve
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(fpr_gene,
         tpr_gene,
         label=f'Gene data ROC (AUC = {auc_gene:.4f})',
         color='blue')
plt.plot(fpr_prot,
         tpr_prot,
         label=f'Protein data ROC (AUC = {auc_prot:.4f})',
         color='green')
plt.plot(fpr_combined,
         tpr_combined,
         label=f'Combined data ROC (AUC = {auc_combined:.4f})',
         color='red',
         linestyle='--')
plt.plot([0, 1], [0, 1], 'k--', label='Random Guess')

plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.grid(True)

# Plot Precision-Recall (PR) curve
plt.subplot(1, 2, 2)
plt.plot(recall_gene,
         precision_gene,
         label=f'Gene data PR (AUC = {pr_auc_gene:.4f})',
         color='blue')
plt.plot(recall_prot,
         precision_prot,
         label=f'Protein data PR (AUC = {pr_auc_prot:.4f})',
         color='green')
plt.plot(recall_combined,
         precision_combined,
         label=f'Combined data PR (AUC = {pr_auc_combined:.4f})',
         color='red',
         linestyle='--')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc='lower left')
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()
