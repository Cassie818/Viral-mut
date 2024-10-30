import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List


def load_data(gene: str):
    protein_data = pd.read_csv(f"Results/Protein/{gene}_ESM2_grammaticality.csv")
    gene_data = pd.read_csv(f"Results/Gene/{gene}_CaLM_grammaticality.csv")
    return protein_data, gene_data


def plot_distribution(df: pd.DataFrame, site: int, typ: str) -> None:
    """
    Plots the distribution of grammaticality values for a given site.

    Args:
        df (pd.DataFrame): Dataframe containing grammaticality values.
        site (int): The site (1-based index) for which the distribution is plotted.
        typ (str): The type of data, either "protein" or "gene".
    """
    if typ == "protein":
        df = df[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']]
    else:
        df = df[['AAA', 'AAU', 'AAC', 'AAG', 'AUA', 'AUU', 'AUC', 'AUG',
                 'ACA', 'ACU', 'ACC', 'ACG', 'AGA', 'AGU', 'AGC', 'AGG',
                 'UAA', 'UAU', 'UAC', 'UAG', 'UUA', 'UUU', 'UUC', 'UUG',
                 'UCA', 'UCU', 'UCC', 'UCG', 'UGA', 'UGU', 'UGC', 'UGG',
                 'GAA', 'GAU', 'CAC', 'CAG', 'CUA', 'CUU', 'CUC', 'CUG',
                 'CCA', 'CCU', 'CCC', 'CCG', 'CGA', 'CGU', 'CGC', 'CGG',
                 'GAA', 'GAU', 'GAC', 'GAG', 'GUA', 'GUU', 'GUC', 'GUG',
                 'GCA', 'GCU', 'GCC', 'GCG', 'GGA', 'GGU', 'GGC', 'GGG']]

    row = df.iloc[site - 1]

    # Use a more refined theme and customize it further
    sns.set_theme(style="whitegrid")

    # Increase figure size for better readability
    plt.figure(figsize=(12, 6))

    # Use a color palette with a gradient for better contrast
    colors = sns.color_palette("mako", n_colors=len(row))

    # Create the bar plot
    sns.barplot(x=row.index, y=row.values, palette=colors)

    # Add a title with increased font size and padding, without bold or italics
    plt.title(f"Grammaticality Distribution at Site {site}", fontsize=12, pad=15)

    # Set x and y labels with increased font size and padding, no bold or italics
    plt.xlabel("Tokens", fontsize=10, labelpad=15)
    plt.ylabel("Grammaticality Values", fontsize=10, labelpad=15)

    # Rotate the x-axis labels for better readability
    plt.xticks(rotation=45, ha="right", fontsize=10)

    # Add grid lines for y-axis only with enhanced visibility
    plt.grid(axis='y', linestyle='--', alpha=0.6)

    # Add a custom background color to the plot
    plt.gca().set_facecolor('#f7f7f7')

    # Tight layout for avoiding overlap
    plt.tight_layout()

    # Show the plot
    plt.show()


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
    plt.figure(figsize=(8, 6))

    # Plot the histograms for both lists with enhanced transparency and modern colors
    plt.hist(pathogenic_list,
             bins=bins,
             alpha=0.5,
             color='#ff7f0e',  # a soft orange color
             label='ClinVar: Pathogenic',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(likely_pathogenic_list,
             bins=bins,
             alpha=0.5,
             color='#ffbb78',  # a lighter orange color
             label='ClinVar: Likely Pathogenic',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(benign_list,
             bins=bins,
             alpha=0.5,
             color='#1f77b4',  # a soft blue color
             label='ClinVar: Benign',
             edgecolor='black',
             linewidth=1.2)

    plt.hist(likely_benign_list,
             bins=bins,
             alpha=0.5,
             color='#aec7e8',  # a lighter blue color
             label='ClinVar: Likely Benign',
             edgecolor='black',
             linewidth=1.2)

    # Add titles and labels with larger font sizes
    plt.title(f'{model_name}')
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
