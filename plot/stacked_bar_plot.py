import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from typing import List


def plot_distribution_by_mutation_types(
        missense_pathogenic: List,
        missense_likely_pathogenic: List,
        missense_benign: List,
        missense_likely_benign: List,
        nonsense_pathogenic: List,
        nonsense_likely_pathogenic: List,
        nonsense_benign: List,
        nonsense_likely_benign: List,
        synonymous_pathogenic: List,
        synonymous_likely_pathogenic: List,
        synonymous_benign: List,
        synonymous_likely_benign: List):
    """
    Plots stacked bar charts showing the distribution of mutation effects
    for each mutation type (Missense, Nonsense, Synonymous, SNP Combined).
    """
    # Calculate counts for each mutation type and effect
    missense_counts = [
        len(missense_pathogenic),
        len(missense_likely_pathogenic),
        len(missense_benign),
        len(missense_likely_benign)
    ]
    nonsense_counts = [
        len(nonsense_pathogenic),
        len(nonsense_likely_pathogenic),
        len(nonsense_benign),
        len(nonsense_likely_benign)
    ]
    synonymous_counts = [
        len(synonymous_pathogenic),
        len(synonymous_likely_pathogenic),
        len(synonymous_benign),
        len(synonymous_likely_benign)
    ]
    combined_counts = [
        sum(x) for x in zip(missense_counts, nonsense_counts, synonymous_counts)
    ]

    # Data preparation
    mutation_types = ['Missense', 'Nonsense', 'Synonymous', 'All']
    effects = ['Pathogenic', 'Likely Pathogenic', 'Benign', 'Likely Benign']
    data = [missense_counts, nonsense_counts, synonymous_counts, combined_counts]
    colors = ['#c1121f', '#ff9896', '#023e8a', '#aec7e8']  # Colors for each effect

    # Transpose data for stacked bar chart
    data_transposed = np.array(data).T

    # Plotting stacked bar chart
    x = np.arange(len(mutation_types))  # Position of bars for each mutation type
    bar_width = 0.5

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['arial']

    bottom_values = np.zeros(len(mutation_types))
    for i, effect in enumerate(effects):
        ax.bar(x, data_transposed[i], bar_width, label=effect, alpha=0.7, color=colors[i], bottom=bottom_values)
        bottom_values += data_transposed[i]

    # Format y-axis to include commas for thousands
    def thousand_separator_formatter(x, _):
        if x == 0:
            return "0"
        return f"{int(x):,}"

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_linewidth(1.5)

    ax.yaxis.set_major_formatter(FuncFormatter(thousand_separator_formatter))

    # Adding labels, title, and legend
    ax.set_xticks(x)
    ax.set_xticklabels(mutation_types, fontsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # Adjusting the legend
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=2, fontsize=16, frameon=False)

    ax.grid(False)

    plt.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()
