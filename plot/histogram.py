import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from typing import List


def plot_histogram(
        pathogenic_lists: List[List[float]],
        likely_pathogenic_lists: List[List[float]],
        benign_lists: List[List[float]],
        likely_benign_lists: List[List[float]],
        all_pathogenic_lists: List[List[float]],
        all_benign_lists: List[List[float]],
        model_names: List[str],
        bins: int = 50,
        x_range: tuple = (-15, 2.5)
) -> None:
    """
    Plots histograms and density curves (KDE) for multiple models.
    Each model includes LLR values for Pathogenic, Likely Pathogenic,
    Benign, and Likely Benign. The overlapping area of the KDE curves
    between pathogenic-related and benign-related Figure is also visualized.

    Parameters
    ----------
    pathogenic_lists : List[List[float]]
        LLR values for the 'Pathogenic' category, grouped by models.
    likely_pathogenic_lists : List[List[float]]
        LLR values for the 'Likely Pathogenic' category, grouped by models.
    benign_lists : List[List[float]]
        LLR values for the 'Benign' category, grouped by models.
    likely_benign_lists : List[List[float]]
        LLR values for the 'Likely Benign' category, grouped by models.
    all_pathogenic_lists : List[List[float]]
        Combined LLR values for Pathogenic and Likely Pathogenic,
        used for density estimation for each model.
    all_benign_lists : List[List[float]]
        Combined LLR values for Benign and Likely Benign,
        used for density estimation for each model.
    model_names : List[str]
        Names of the models, used for subplot titles.
    bins : int, optional
        Number of bins for the histogram, default is 50.
    x_range : tuple, optional
        The x-axis range for both histogram and density plots,
        default is (-15, 2.5).

    """

    palette = {
        'Pathogenic': '#c1121f',
        'Likely Pathogenic': '#ff9896',
        'Benign': '#023e8a',
        'Likely Benign': '#aec7e8',
        'Overlap': '#FFD700'
    }

    sns.set(style="white", context="paper")
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'legend.fontsize': 10,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']
    })

    num_models = len(model_names)
    fig, axes = plt.subplots(num_models,
                             1,
                             figsize=(9, 3 * num_models),
                             sharex=True)
    plt.subplots_adjust(hspace=0.4)

    # Common histogram parameters
    hist_params = {
        'bins': bins,
        'range': x_range,
        'alpha': 0.8,
        'edgecolor': 'gray',
        'linewidth': 0.5,
        'density': False
    }

    # Iterate through each model
    for i, ax in enumerate(axes):
        # Plot histograms for each category
        ax.hist(pathogenic_lists[i], **hist_params, color=palette['Pathogenic'])
        ax.hist(likely_pathogenic_lists[i], **hist_params, color=palette['Likely Pathogenic'])
        ax.hist(benign_lists[i], **hist_params, color=palette['Benign'])
        ax.hist(likely_benign_lists[i], **hist_params, color=palette['Likely Benign'])

        # Format y-axis for histograms
        ax.tick_params(axis='y', labelsize=12, length=4, width=1, color='#404040')

        # Create a secondary axis for density plots
        ax_density = ax.twinx()

        # KDE plots for Pathogenic and Benign
        kde_pathogenic = sns.kdeplot(
            all_pathogenic_lists[i],
            ax=ax_density,
            color='#d62728',
            linewidth=1.5,
            linestyle="--",
            label='Pathogenic'
        )
        kde_benign = sns.kdeplot(
            all_benign_lists[i],
            ax=ax_density,
            color='#1f77b4',
            linewidth=1.5,
            linestyle="--",
            label='Benign'
        )

        # Extract (x, y) Figure from the last two lines: pathogenic and benign
        pathogenic_x, pathogenic_y = kde_pathogenic.get_lines()[-2].get_data()
        benign_x, benign_y = kde_benign.get_lines()[-1].get_data()

        # Interpolate on a common x range for overlap calculation
        common_x = np.linspace(x_range[0], x_range[1], 5000)
        pathogenic_interp = np.interp(common_x, pathogenic_x, pathogenic_y)
        benign_interp = np.interp(common_x, benign_x, benign_y)

        # Calculate overlap area
        overlap_density = np.minimum(pathogenic_interp, benign_interp)
        overlap_area = np.trapz(overlap_density, common_x)

        # Fill the overlapping area
        ax_density.fill_between(
            common_x,
            0,
            overlap_density,
            color=palette['Overlap'],
            alpha=0.2,
            lw=0
        )

        # Set y-axis limits for histograms and format
        ax.set_ylim(0, 12000)
        ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
        ax.tick_params(axis='x', labelsize=12)

        # Set y-axis limits for the density plot
        ax_density.set_ylim(0, 0.4)

        ax_density.tick_params(
            axis='y',
            which='both',
            direction='in',
            pad=-25,
            labelright=True,
            labelleft=False,
            color='black',
            labelsize=12
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color('#808080')
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_color('#808080')
        ax.spines["bottom"].set_linewidth(0.8)


        ax_density.spines["top"].set_visible(False)
        ax_density.spines["left"].set_visible(False)
        ax_density.spines["bottom"].set_visible(False)
        ax_density.spines["right"].set_visible(True)
        ax_density.spines["right"].set_color('#808080')
        ax_density.spines["right"].set_linewidth(0.8)
        ax_density.set_ylabel('')

        # Title and text box for overlap area
        ax.set_title(f"{model_names[i]}", fontweight='semibold', pad=10)
        ax_density.text(
            0.99,
            0.14,
            f"overlap area: {overlap_area:.4f}",
            transform=ax_density.transAxes,
            color='gray',
            ha='right',
            va='top',
            fontsize=10
        )

    # Create legend handles
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=palette['Pathogenic'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Likely Pathogenic'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Benign'], alpha=0.8),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Likely Benign'], alpha=0.8),
        plt.Line2D([], [], color='#d62728', linewidth=2, linestyle='--'),
        plt.Line2D([], [], color='#1f77b4', linewidth=2, linestyle='--'),
        plt.Rectangle((0, 0), 1, 1, fc=palette['Overlap'], alpha=0.3)
    ]
    labels = [
        'Pathogenic', 'Likely Pathogenic',
        'Benign', 'Likely Benign',
        'Combined pathogenic density', 'Combined benign density',
        'Overlap area'
    ]

    # X-axis range for the last subplot
    axes[-1].set_xlim(x_range)

    # Overall labeling and legend
    fig.supxlabel('LLR', y=0.08, fontsize=12, fontweight='normal')
    fig.align_ylabels()

    fig.legend(
        handles,
        labels,
        loc='lower center',
        ncol=4,
        bbox_to_anchor=(0.5, 0),
        frameon=False
    )

    # Tight layout adjustment
    plt.tight_layout(rect=[0, 0.05, 1, 0.99])
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


    all_pathogenic_lists = [np.append(LLR_path_prot, LLR_likely_path_prot),
                            np.append(LLR_path_gene, LLR_likely_path_gene),
                            np.append(LLR_path_combined, LLR_likely_path_combined)]

    all_benign_lists = [np.append(LLR_benign_prot, LLR_likely_benign_prot),
                        np.append(LLR_benign_gene, LLR_likely_benign_gene),
                        np.append(LLR_benign_combined, LLR_likely_benign_combined)]

    model_names = ['ESM-2', 'CaLM', 'Hybrid']

    plot_histogram(pathogenic_lists,
                   likely_pathogenic_lists,
                   benign_lists,
                   likely_benign_lists,
                   all_pathogenic_lists,
                   all_benign_lists,
                   model_names)
