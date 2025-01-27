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
    palette = {
        'Pathogenic': '#c1121f',
        'Likely Pathogenic': '#ff9896',
        'Benign': '#023e8a',
        'Likely Benign': '#aec7e8',
        'Overlap': '#ff9f1c'
    }

    sns.set(style="white", context="paper")
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'legend.fontsize': 10
    })

    num_models = len(model_names)
    fig, axes = plt.subplots(num_models, 1, figsize=(9, 3 * num_models), sharex=True)
    plt.subplots_adjust(hspace=0.4)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['arial']

    for i, ax in enumerate(axes):
        hist_params = {
            'bins': bins,
            'range': x_range,
            'alpha': 0.8,
            'edgecolor': 'gray',
            'linewidth': 0.5,
            'density': False
        }

        ax.hist(pathogenic_lists[i], **hist_params, color=palette['Pathogenic'])
        ax.hist(likely_pathogenic_lists[i], **hist_params, color=palette['Likely Pathogenic'])
        ax.hist(benign_lists[i], **hist_params, color=palette['Benign'])
        ax.hist(likely_benign_lists[i], **hist_params, color=palette['Likely Benign'])
        ax.tick_params(axis='y',
                       labelsize=12,
                       length=4,
                       width=1,
                       color='#404040')

        ax_density = ax.twinx()
        kde_params = {
            'linewidth': 2,
            'linestyle': '-',
            'alpha': 0.8
        }

        kde_pathogenic = sns.kdeplot(all_pathogenic_lists[i], ax=ax_density, color='#d62728', linewidth=1.5,
                                     linestyle="--", label='Pathogenic')
        kde_benign = sns.kdeplot(all_benign_lists[i], ax=ax_density, color='#1f77b4', linewidth=1.5, linestyle="--",
                                 label='Benign')

        kde_pathogenic_x, kde_pathogenic_y = kde_pathogenic.get_lines()[-2].get_data()  # Pathogenic line data
        kde_benign_x, kde_benign_y = kde_benign.get_lines()[-1].get_data()  # Benign line data

        common_x = np.linspace(x_range[0], x_range[1], 5000)
        kde_pathogenic_interp = np.interp(common_x, kde_pathogenic_x, kde_pathogenic_y)
        kde_benign_interp = np.interp(common_x, kde_benign_x, kde_benign_y)

        overlap_density = np.minimum(kde_pathogenic_interp, kde_benign_interp)
        overlap_area = np.trapz(overlap_density, common_x)

        ax_density.fill_between(common_x, 0, overlap_density,
                                color=palette['Overlap'], alpha=0.2, lw=0)

        ax.set_ylim(0, 12000)
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
        ax.tick_params(axis='x', labelsize=12)

        ax_density.set_ylim(0, 0.4)
        ax_density.tick_params(axis='y', labelsize=12, right=False, color='black')
        ax_density.spines['right'].set_visible(False)
        ax_density.set_ylabel('')

        for ax in [ax, ax_density]:
            for spine in ax.spines.values():
                spine.set_color('#808080')
                spine.set_linewidth(0.8)

        ax.set_title(f"{model_names[i]}", fontweight='semibold', pad=10)
        ax_density.text(0.99, 0.95,
                        f"overlap: {overlap_area:.4f}",
                        transform=ax_density.transAxes,
                        ha='right', va='top',
                        fontsize=10)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=palette['Pathogenic'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Likely Pathogenic'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Benign'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Likely Benign'], alpha=0.8),
        plt.Line2D([], [], color=palette['Pathogenic'], linewidth=2),
        plt.Line2D([], [], color=palette['Benign'], linewidth=2),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Overlap'], alpha=0.3)
    ]

    labels = ['Pathogenic', 'Likely Pathogenic', 'Benign', 'Likely Benign',
              'Pathogenic density', 'Benign density', 'Overlap']

    # 优化布局
    fig.supxlabel('LLR', y=0.08, fontsize=12, fontweight='normal')
    axes[-1].set_xlim(x_range)
    ax.tick_params(axis='x', labelsize=12)
    fig.align_ylabels()

    fig.legend(handles, labels,
               loc='lower center',
               ncol=4,
               bbox_to_anchor=(0.5, -0.01),
               frameon=True,
               framealpha=0.9)

    plt.tight_layout(rect=[0, 0.05, 1, 0.99])


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
