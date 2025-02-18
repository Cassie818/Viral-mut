import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

np.random.seed(52)

# Protein sequence and amino acids for the axes
protein_sequence = ['Y', 'F', 'D', 'E', 'L', 'A', 'C', 'R', 'D', 'T']
amino_acids_vertical = ['D', 'E', 'R', 'H', 'K']

# Generate random scores for the heatmap
log_scores = np.random.rand(len(amino_acids_vertical), len(protein_sequence))

# Figure dimensions
fig_width = len(protein_sequence) * 0.8
fig_height = len(amino_acids_vertical) * 0.8
plt.figure(figsize=(fig_width, fig_height))

# Aesthetic Morandi-inspired color palette with modern twist
color_palette = ['#e6f7ff', '#b3e0ff', '#80ccff', '#4db8ff', '#1aa3ff']
cmap = ListedColormap(sns.color_palette(color_palette, n_colors=256))

# Create the heatmap
ax = sns.heatmap(log_scores,
                 annot=False,
                 cmap=cmap,
                 cbar=False,
                 xticklabels=protein_sequence,
                 yticklabels=amino_acids_vertical,
                 linewidth=0.5,
                 linecolor='grey')

# Customize the font size for xticks and yticks
plt.xticks(fontsize=20, rotation=0)
plt.yticks(fontsize=20, rotation=0)

# Optional: Highlight a specific mutation
highlight_position = (3, 2)
plt.gca().add_patch(plt.Rectangle((highlight_position[1],
                                   highlight_position[0]),
                                  1,
                                  1,
                                  fill=False,
                                  edgecolor='red',
                                  lw=3))
plt.text(highlight_position[1] + 0.5,
         highlight_position[0] + 0.5,
         'D73H',
         color='black',
         ha='center',
         va='center',
         fontsize=20)
plt.text(0.0,
         5.2,
         f"{70}",
         ha="center",
         va="center",
         fontsize=14,
         transform=ax.transData)
plt.text(9.85,
         5.2,
         f"{80}",
         ha="center",
         va="center",
         fontsize=14,
         transform=ax.transData)

plt.subplots_adjust(right=0.9)

# Apply tight layout to optimize space
plt.tight_layout()
plt.show()
