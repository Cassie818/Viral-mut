import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List


# Function to plot the distribution of grammaticality at a specific site
def plot_distribution(df: pd.DataFrame, site: int) -> None:
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


def plot_histogram(Pathogenic_list: List, benign_list: List, bins=50):
    """
    Plots a histogram comparing pathogenic and benign ClinVar data.

    Parameters:
    - Pathogenic_list: list or array of values representing pathogenic data.
    - benign_list: list or array of values representing benign data.
    - bins: number of bins to use in the histograms (default is 50).
    """
    plt.hist(Pathogenic_list, bins=bins, alpha=0.6, color='orange', label='ClinVar: pathogenic', edgecolor='black')
    plt.hist(benign_list, bins=bins, alpha=0.5, color='blue', label='ClinVar: benign', edgecolor='black')

    plt.title('ESM2-150M')
    plt.xlabel('LLR')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()
