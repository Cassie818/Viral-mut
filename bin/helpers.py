import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List


# Function to plot the distribution of grammaticality at a specific site
def plot_distribution(df: pd.DataFrame,
                      site: int) -> None:
    """
    Plots the distribution of grammaticality values for a given site.

    Args:
        df (pd.DataFrame): Dataframe containing grammaticality values.
        site (int): The site (1-based index) for which the distribution is plotted.
    """
    row = df.iloc[site - 1]

    plt.figure(figsize=(12, 6))
    sns.barplot(x=row.index, y=row.values)
    plt.title(f"Distribution of grammaticality in Site {site}", fontsize=16)
    plt.xlabel("Tokens", fontsize=12)
    plt.ylabel("Values", fontsize=12)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plot_histogram(pathogenic_list: List,
                   benign_list: List,
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
    plt.figure(figsize=(8, 6))

    # Plot the histograms for both lists with enhanced transparency and modern colors
    plt.hist(pathogenic_list,
             bins=bins,
             alpha=0.5,
             color='#ff7f0e',  # a soft orange color
             label='ClinVar: Pathogenic',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(benign_list,
             bins=bins,
             alpha=0.5,
             color='#1f77b4',  # a soft blue color
             label='ClinVar: Benign',
             edgecolor='black',
             linewidth=1.2)

    # Add titles and labels with larger font sizes
    plt.title('Combined')
    plt.xlabel('LLR')
    plt.ylabel('Frequency')
    plt.legend(fontsize=12)

    # Display grid for better readability
    plt.grid(True, linestyle='--', alpha=0.6)

    # Tight layout for better spacing
    plt.tight_layout()

    # Show the plot
    plt.show()


def flatten(nested_list: List) -> List:
    return [item for sublist in nested_list for item in flatten(sublist)] if isinstance(nested_list, list) else [
        nested_list]
