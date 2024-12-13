import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap


np.random.seed(56)


# Protein sequences: YFDELACRDT
protein_sequence = ['Y', 'F', 'D', 'E', 'L', 'A', 'C', 'R', 'D', 'T']
amino_acids_vertical = ['D', 'E', 'R', 'H', 'K']
positions = list(range(1, len(protein_sequence) + 1))

# Generate random log scores for heatmap
log_scores = np.random.rand(len(amino_acids_vertical), len(protein_sequence))


cell_width = 0.8
cell_height = 0.8
fig_width = len(protein_sequence) * cell_width
fig_height = len(amino_acids_vertical) * cell_height

# Configure plot
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['arial']
plt.figure(figsize=(fig_width, fig_height))

# Create a colormap with transparency
cmap = sns.color_palette("viridis", as_cmap=True)
transparent_cmap = cmap(np.linspace(0, 1, 256))
transparent_cmap[:, -1] = 0.8
transparent_cmap = ListedColormap(transparent_cmap)


ax = sns.heatmap(
    log_scores,
    annot=False,
    cmap=transparent_cmap,
    cbar=False,  # Disable color bar
    xticklabels=protein_sequence,
    yticklabels=amino_acids_vertical,
)


highlight_position = (3, 2)
ax.add_patch(
    plt.Rectangle(
        (highlight_position[1], highlight_position[0]), 1, 1, fill=False, edgecolor='yellow', lw=2
    )
)
plt.text(
    highlight_position[1] + 0.5,
    highlight_position[0] + 0.5,
    "D73H",
    color="red",
    ha="center",
    va="center",
    fontsize=18,
)


ax.set_xticklabels(protein_sequence, fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)


plt.text(0, 5.5, f"{70}", ha="center", va="center", fontsize=18, transform=ax.transData)
plt.text(len(protein_sequence), 5.5, f"{80}", ha="center", va="center", fontsize=18, transform=ax.transData)


plt.tight_layout()
plt.show()