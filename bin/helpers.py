import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List


def load_data(gene: str):
    protein_data = pd.read_csv(f"Results/Protein/{gene}_ESM2_grammaticality.csv")
    gene_data = pd.read_csv(f"Results/Gene/{gene}_CaLM_grammaticality.csv")
    return protein_data, gene_data


def plot_distribution(df: pd.DataFrame,
                      site: int,
                      typ: str) -> None:
    """
    Plots the distribution of grammaticality values for a given site.

    Args:
        df (pd.DataFrame): Dataframe containing grammaticality values.
        site (int): The site (1-based index) for which the distribution is plotted.
        typ (str): The type of data, either "protein" or "gene".
    """
    if typ == "protein":
        df = df[['K', 'N', 'I', 'M', 'T',
                 'R', 'S', 'Y', 'L', 'F',
                 'C', 'W', 'Q', 'H', 'P',
                 'E', 'D', 'V', 'A', 'G']]
    else:
        df = df[['AAA', 'AAG',
                 'AAU', 'AAC',
                 'AUA', 'AUU', 'AUC',
                 'AUG',
                 'ACA', 'ACC', 'ACG', 'ACU',
                 'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU',
                 'AGU', 'AGC', 'UCA', 'UCC', 'UCG', 'UCU',
                 'UAU', 'UAC',
                 'UUA', 'UUG', 'CUA', 'CUC', 'CUG', 'CUU',
                 'UUU', 'UUC',
                 'UGU', 'UGC',
                 'UGG',
                 'CAA', 'CAG',
                 'CAU', 'CAC',
                 'CCA', 'CCC', 'CCG', 'CCU',
                 'GAA', 'GAG',
                 'GAU', 'GAC',
                 'GUA', 'GUG', 'GUC', 'GUU',
                 'GCA', 'GCG', 'GCU', 'GCC',
                 'GGA', 'GGG', 'GGU', 'GGC',
                 'UAA', 'UAG', 'UGA']]

    df = df.rename(columns={'AAA': 'AAA (K)', 'AAG': 'AAG (K)',
                            'AAU': 'AAU (N)', 'AAC': 'AAC (N)',
                            'AUA': 'AUA (I)', 'AUU': 'AUU (I)', 'AUC': 'AUC (I)',
                            'AUG': 'AUG (M)',
                            'ACA': 'ACA (T)', 'ACC': 'ACC (T)', 'ACG': 'ACG (T)', 'ACU': 'ACU (T)',
                            'AGA': 'AGA (R)', 'AGG': 'AGG (R)', 'CGA': 'CGA (R)', 'CGC': 'CGC (R)', 'CGG': 'CGG (R)',
                            'CGU': 'CGU (R)',
                            'AGU': 'AGU (S)', 'AGC': 'AGC (S)', 'UCA': 'UCA (S)', 'UCC': 'UCC (S)', 'UCG': 'UCG (S)',
                            'UCU': 'UCU (S)',
                            'UAU': 'UAU (Y)', 'UAC': 'UAC (Y)',
                            'UUA': 'UUA (L)', 'UUG': 'UUG (L)', 'CUA': 'CUA (L)', 'CUC': 'CUC (L)', 'CUG': 'CUG (L)',
                            'CUU': 'CUU (L)',
                            'UUU': 'UUU (F)', 'UUC': 'UUC (F)',
                            'UGU': 'UGU (C)', 'UGC': 'UGC (C)',
                            'UGG': 'UGG (W)',
                            'CAA': 'CAA (Q)', 'CAG': 'CAG (Q)',
                            'CAU': 'CAU (H)', 'CAC': 'CAC (H)',
                            'CCA': 'CCA (P)', 'CCC': 'CCC (P)', 'CCG': 'CCG (P)', 'CCU': 'CCU (P)',
                            'GAA': 'GAA (E)', 'GAG': 'GAG (E)',
                            'GAU': 'GAU (D)', 'GAC': 'GAC (D)',
                            'GUA': 'GUA (V)', 'GUG': 'GUG (V)', 'GUC': 'GUC (V)', 'GUU': 'GUU (V)',
                            'GCA': 'GCA (A)', 'GCG': 'GCG (A)', 'GCU': 'GCU (A)', 'GCC': 'GCC (A)',
                            'GGA': 'GGA (G)', 'GGG': 'GGG (G)', 'GGU': 'GGU (G)', 'GGC': 'GGC (G)',
                            'UAA': 'UAA (Stop)', 'UAG': 'UAG (Stop)', 'UGA': 'UGA (Stop)',
                            })

    row = df.iloc[site - 1]

    # Use a more refined theme and customize it further
    sns.set_theme(style="whitegrid")

    # Increase figure size for better readability
    plt.figure(figsize=(12, 6))

    # Use a color palette with a gradient for better contrast
    colors = sns.color_palette("mako", n_colors=len(row))

    # Create the barplot
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
