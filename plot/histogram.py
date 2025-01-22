import matplotlib.pyplot as plt
import seaborn as sns
from typing import List
import matplotlib.ticker as ticker


def plot_histogram(pathogenic_lists: List[List[float]],
                   likely_pathogenic_lists: List[List[float]],
                   benign_lists: List[List[float]],
                   likely_benign_lists: List[List[float]],
                   model_names: List[str],
                   bins=50,
                   x_range=(-15, 2.5)):
    """
    Plots histograms for multiple models

    Parameters:
    - pathogenic_lists: list of lists, each containing pathogenic values for a model.
    - likely_pathogenic_lists: list of lists, each containing likely pathogenic values for a model.
    - benign_lists: list of lists, each containing benign values for a model.
    - likely_benign_lists: list of lists, each containing likely benign values for a model.
    - model_names: list of model names corresponding to the data.
    - bins: number of bins to use in the histograms (default is 50).
    - x_range: tuple specifying the fixed x-axis range (default is (-15, 10)).
    """
    sns.set(style="whitegrid", context="talk")

    # Create subplots
    num_models = len(model_names)
    fig, axes = plt.subplots(num_models, 1, figsize=(12, 9), sharex=True)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['arial']

    for i, ax in enumerate(axes):
        ax.hist(pathogenic_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.5,
                color='#d62728',
                label='Pathogenic',
                edgecolor='white',
                linewidth=1)

        ax.hist(likely_pathogenic_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.5,
                color='#ff9896',
                label='Likely Pathogenic',
                edgecolor='white',
                linewidth=1)

        ax.hist(benign_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.5,
                color='#1f77b4',
                label='Benign',
                edgecolor='white',
                linewidth=1)

        ax.hist(likely_benign_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.5,
                color='#aec7e8',
                label='Likely Benign',
                edgecolor='white',
                linewidth=1)

        # Add titles and labels
        ax.set_title(model_names[i])
        ax.set_ylim(0, 12000)
        ax.set_xlim(-15, 2.5)

        # Format y-axis with comma-separated labels every 2000
        ax.set_yticks(range(0, 12001, 4000))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

        ax.grid(True, linestyle='--', alpha=0.5)

    # Add legend to the last plot only
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=14)

    # Set a common x-axis label and range
    # Add a common x-axis label
    fig.text(0.5, 0.08, 'LLR', ha='center', va='center', fontsize=16)
    axes[-1].set_xlim(x_range)

    # Adjust layout for better spacing
    plt.tight_layout(rect=[0, 0.07, 1, 1])

    # Show the plot
    plt.show()



if __name__ == "__main__":
    pathogenic_lists = [LLR_path_prot,
                        LLR_path_gene,
                        LLR_path_combined]
    likely_pathogenic_lists = [LLR_likely_path_prot,
                               LLR_likely_path_gene,
                               LLR_likely_path_combined]
    benign_lists = [LLR_benign_prot,
                    LLR_benign_gene,
                    LLR_benign_combined]

    likely_benign_lists = [LLR_likely_benign_prot,
                           LLR_likely_benign_gene,
                           LLR_likely_benign_combined]

    model_names = ['ESM-2', 'CaLM', 'Hybrid']

    plot_histogram(pathogenic_lists,
                   likely_pathogenic_lists,
                   benign_lists,
                   likely_benign_lists,
                   model_names)
