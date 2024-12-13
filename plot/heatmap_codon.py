import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap


np.random.seed(56)

# Original protein sequence
protein_sequence = ['Y', 'F', 'D', 'E', 'L', 'A', 'C', 'R', 'D', 'T']

# Generate corresponding gene sequence
gene_sequence = ['TAC', 'TTC',
                 'GAC', 'GAA',
                 'CTG', 'GCC',
                 'TGC', 'CGG',
                 'GAC','ACA']

positions = list(range(1, len(gene_sequence) + 1))

amino_acids_vertical = ['GAT', 'GAA', 'CGA', 'CAC', 'AAA']
log_scores = np.random.rand(len(amino_acids_vertical), len(gene_sequence))

cell_width = 0.9
cell_height = 0.8
fig_width = len(gene_sequence) * cell_width
fig_height = len(amino_acids_vertical) * cell_height

# Configure plot
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['arial']
plt.figure(figsize=(fig_width, fig_height))

# Create a colormap with transparency
cmap = sns.color_palette('coolwarm', as_cmap=True)
transparent_cmap = cmap(np.linspace(0, 1, 256))
transparent_cmap[:, -1] = 0.8
transparent_cmap = ListedColormap(transparent_cmap)

# Plot heatmap
ax = sns.heatmap(
    log_scores,
    annot=False,
    cmap=transparent_cmap,
    cbar=False,
    xticklabels=gene_sequence,
    yticklabels=amino_acids_vertical,
)

# Highlight a specific cell
highlight_position = (3, 2)
ax.add_patch(
    plt.Rectangle(
        (highlight_position[1], highlight_position[0]), 1, 1, fill=False, edgecolor='yellow', lw=2
    )
)
plt.text(
    highlight_position[1] + 0.5,
    highlight_position[0] + 0.5,
    "GAC73CAC",
    color="red",
    ha="center",
    va="center",
    fontsize=10,
)


ax.set_xticklabels(gene_sequence, fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)

# Add position numbers at the start and end of the sequence
plt.text(-0.8, 5.7, f"{210}", ha="center", va="center", fontsize=18, transform=ax.transData)
plt.text(len(gene_sequence), 5.7, f"{240}", ha="center", va="center", fontsize=18, transform=ax.transData)

plt.tight_layout()
plt.show()