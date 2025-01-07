import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

np.random.seed(56)

# Protein sequence and amino acids for the axes
gene_sequence = ['TAC', 'TTC', 'GAC', 'GAA', 'CTG', 'GCC', 'TGC', 'CGG', 'GAC', 'ACA']
codon_vertical = ['GAT', 'GAA', 'CGA', 'CAC', 'AAA']
log_scores = np.random.rand(len(codon_vertical), len(gene_sequence))

cell_width = 0.8
cell_height = 0.8
fig_width = len(gene_sequence) * cell_width
fig_height = len(codon_vertical) * cell_height

# Configure plot
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['arial']
plt.figure(figsize=(fig_width, fig_height))

# Create a light green color palette
colors = ['#e8f5e9', '#c8e6c9', '#a5d6a7', '#81c784', '#66bb6a']
cmap = ListedColormap(sns.color_palette(colors, n_colors=256))

# Plot heatmap
ax = sns.heatmap(
    log_scores,
    annot=False,
    cmap=cmap,
    cbar=False,
    xticklabels=gene_sequence,
    yticklabels=codon_vertical,
    linewidth=0.5,
    linecolor='grey'
)

# Highlight a specific cell
highlight_position = (3, 2)
ax.add_patch(
    plt.Rectangle(
        (highlight_position[1], highlight_position[0]), 1, 1, fill=False, edgecolor='red', lw=3
    )
)
plt.text(
    highlight_position[1] + 0.5,
    highlight_position[0] + 0.5,
    "GAC73CAC",
    color="black",
    ha="center",
    va="center",
    fontsize=10,
)

ax.set_xticklabels(gene_sequence, fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)

# Add position numbers at the start and end of the sequence
plt.text(0, 5.7, f"{210}", ha="center", va="center", fontsize=14, transform=ax.transData)
plt.text(len(gene_sequence), 5.7, f"{240}", ha="center", va="center", fontsize=14, transform=ax.transData)

plt.tight_layout()
plt.show()
