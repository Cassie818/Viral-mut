import matplotlib.pyplot as plt
import seaborn as sns
from typing import List


def plot_histogram(pathogenic_list: List,
                   likely_pathogenic_list: List,
                   benign_list: List,
                   likely_benign_list: List,
                   model_name,
                   bins=50):
    """
    Plots a histogram comparing pathogenic and benign ClinVar data with enhanced aesthetics.

    Parameters:
    - pathogenic_list: list or array of values representing pathogenic data.
    - benign_list: list or array of values representing benign data.
    - bins: number of bins to use in the histograms (default is 50).
    """
    sns.set(style="whitegrid", context="talk")

    # Create the plot
    plt.figure(figsize=(7, 6))

    # Plot the histograms for both lists with enhanced transparency and modern colors
    plt.hist(pathogenic_list,
             bins=bins,
             alpha=0.5,
             color='#FF0000',  # red
             label='Pathogenic',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(likely_pathogenic_list,
             bins=bins,
             alpha=0.5,
             color='#FF6666',  # a lighter red
             label='Likely Pathogenic',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(benign_list,
             bins=bins,
             alpha=0.5,
             color='#0000FF',  # a soft blue color
             label='Benign',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(likely_benign_list,
             bins=bins,
             alpha=0.5,
             color='#99CCFF',  # a lighter blue color
             label='Likely Benign',
             edgecolor='black',
             linewidth=1.2)

    # Add titles and labels with larger font sizes
    plt.title(f'{model_name}')
    plt.xlabel('LLR')

    plt.ylim(0, 12000)  # Set y-axis range
    plt.yticks(range(0, 12001, 2000))

    # Display grid for better readability
    plt.grid(True, linestyle='--', alpha=0.6)

    # Tight layout for better spacing
    plt.tight_layout()

    # Show the plot
    plt.show()


if __name__ == "__main__":
    plot_histogram(LLR_path_prot,
                   LLR_likely_path_prot,
                   LLR_benign_prot,
                   LLR_likely_benign_prot,
                   'ESM-2 (150M)')
