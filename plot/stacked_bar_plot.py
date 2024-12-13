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
    colors = ['#d62728', '#ff9896', '#1f77b4', '#aec7e8']  # Colors for each effect

    # Transpose data for stacked bar chart
    data_transposed = np.array(data).T

    # Plotting stacked bar chart
    x = np.arange(len(mutation_types))  # Position of bars for each mutation type
    bar_width = 0.6

    fig, ax = plt.subplots(figsize=(10, 6))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['arial']

    bottom_values = np.zeros(len(mutation_types))
    for i, effect in enumerate(effects):
        ax.bar(x, data_transposed[i], bar_width, label=effect, color=colors[i], bottom=bottom_values)
        bottom_values += data_transposed[i]

    # Format y-axis to include commas for thousands
    def thousand_separator_formatter(x, _):
        if x == 0:
            return "0"
        return f"{int(x):,}"

    ax.yaxis.set_major_formatter(FuncFormatter(thousand_separator_formatter))

    # Adding labels, title, and legend
    ax.set_xticks(x)
    ax.set_xticklabels(mutation_types, fontsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Adjusting the legend
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=14, frameon=False)

    ax.grid(False)

    plt.tight_layout(rect=[0, 0.01, 1, 0.95])
    plt.show()

