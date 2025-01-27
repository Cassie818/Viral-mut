import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from typing import List


def plot_histogram(pathogenic_lists: List[List[float]],
                   likely_pathogenic_lists: List[List[float]],
                   benign_lists: List[List[float]],
                   likely_benign_lists: List[List[float]],
                   all_pathogenic_lists: List[List[float]],
                   all_benign_lists: List[List[float]],
                   model_names: List[str],
                   bins=50,
                   x_range=(-15, 2.5)):
    """
    Plots histograms and density plots for multiple models with shared x-axis and fixed range.
    Highlights the overlapping areas between Pathogenic and Benign density plots.
    """
    sns.set(style="whitegrid", context="talk")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

    # Create subplots
    num_models = len(model_names)
    fig, axes = plt.subplots(num_models, 1, figsize=(12, 12), sharex=True)
    if num_models == 1:
        axes = [axes]

    all_handles = []
    all_labels = []

    for i, ax in enumerate(axes):
        # Plot histograms on the left y-axis
        ax.hist(pathogenic_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.6,
                color='#d62728',
                label='Pathogenic',
                edgecolor='white',
                linewidth=1.5)
        ax.hist(likely_pathogenic_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.6,
                color='#ff9896',
                label='Likely Pathogenic',
                edgecolor='white',
                linewidth=1.5)
        ax.hist(benign_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.6,
                color='#1f77b4',
                label='Benign',
                edgecolor='white',
                linewidth=1.5)
        ax.hist(likely_benign_lists[i],
                bins=bins,
                range=x_range,
                alpha=0.6,
                color='#aec7e8',
                label='Likely Benign',
                edgecolor='white',
                linewidth=1.5)

        # Create a secondary y-axis for density plots
        ax_density = ax.twinx()

        kde_pathogenic = sns.kdeplot(all_pathogenic_lists[i], ax=ax_density, color='#d62728', linewidth=2.5,
                                     linestyle="--", label='Pathogenic')
        kde_benign = sns.kdeplot(all_benign_lists[i], ax=ax_density, color='#1f77b4', linewidth=2.5, linestyle="--",
                                 label='Benign')

        kde_pathogenic_x, kde_pathogenic_y = kde_pathogenic.get_lines()[-2].get_data()  # Pathogenic line data
        kde_benign_x, kde_benign_y = kde_benign.get_lines()[-1].get_data()  # Benign line data

        common_x = np.linspace(x_range[0], x_range[1], 5000)
        kde_pathogenic_interp = np.interp(common_x, kde_pathogenic_x, kde_pathogenic_y)
        kde_benign_interp = np.interp(common_x, kde_benign_x, kde_benign_y)

        overlap_density = np.minimum(kde_pathogenic_interp, kde_benign_interp)
        overlap_area = np.trapz(overlap_density, common_x)

        total_area_pathogenic = np.trapz(kde_pathogenic_interp, common_x)
        total_area_benign = np.trapz(kde_benign_interp, common_x)

        ax_density.fill_between(common_x, 0, overlap_density, color='#FFD700', alpha=0.3, label='Overlap')

        ax.set_ylim(0, 12000)
        ax.set_yticks(range(0, 12001, 4000))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax.grid(True, linestyle='--', alpha=0.5)

        ax_density.set_ylabel('')
        ax_density.set_ylim(0, 0.4)
        ax_density.grid(False)

        # Add title and annotate overlap area
        ax.set_title(f"{model_names[i]}", fontsize=16, pad=20)
        ax_density.text(0.95, 0.95,
                        f"Overlap: {overlap_area:.4f}",
                        transform=ax_density.transAxes, ha='right', va='top', fontsize=14)

        # Collect legend handles and labels from the first subplot
        if i == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles_density, labels_density = ax_density.get_legend_handles_labels()
            all_handles.extend(handles + handles_density)
            all_labels.extend(labels + labels_density)

    fig.text(0.5, 0.05, 'LLR', ha='center', va='center', fontsize=18)
    axes[-1].set_xlim(x_range)

    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.95)

    fig.legend(all_handles, all_labels, loc='lower center', ncol=4, fontsize=14, bbox_to_anchor=(0.5, -0.05))

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

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
