import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List
import matplotlib.ticker as ticker


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


def plot_histogram(
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
        synonymous_likely_benign: List,
        model_name: str,
        bins=50):
    """
    Plots separate histograms as individual charts. One chart each for Missense, Nonsense, Synonymous mutations,
    and one for SNPs Combined.

    Parameters:
    - Each parameter corresponds to a specific type and pathogenicity level.
    - model_name: title of the plot.
    - bins: number of bins to use in the histograms (default is 50).
    """
    sns.set(style="whitegrid", context="talk")

    # Define colors for each mutation type
    colors = {
        'pathogenic': '#d62728',  # red
        'likely_pathogenic': '#ff9896',  # light red
        'benign': '#1f77b4',  # blue
        'likely_benign': '#aec7e8'  # light blue
    }

    # Combine data for SNPs
    snp_pathogenic = missense_pathogenic + nonsense_pathogenic + synonymous_pathogenic
    snp_likely_pathogenic = missense_likely_pathogenic + nonsense_likely_pathogenic + synonymous_likely_pathogenic
    snp_benign = missense_benign + nonsense_benign + synonymous_benign
    snp_likely_benign = missense_likely_benign + nonsense_likely_benign + synonymous_likely_benign

    # Helper function to format the y-axis
    def format_axis(ax):
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
        ax.set_xlabel('LLR', fontsize=14)
        ax.set_ylabel('Frequency', fontsize=14)

    # Plot Missense mutations
    plt.figure(figsize=(8, 6))
    plt.hist(missense_pathogenic, bins=bins, alpha=0.6, color=colors['pathogenic'],
             label='Pathogenic', edgecolor='black', linewidth=1)
    plt.hist(missense_likely_pathogenic, bins=bins, alpha=0.6, color=colors['likely_pathogenic'],
             label='Likely Pathogenic', edgecolor='black', linewidth=1)
    plt.hist(missense_benign, bins=bins, alpha=0.6, color=colors['benign'],
             label='Benign', edgecolor='black', linewidth=1)
    plt.hist(missense_likely_benign, bins=bins, alpha=0.6, color=colors['likely_benign'],
             label='Likely Benign', edgecolor='black', linewidth=1)
    plt.title('Missense mutations', fontsize=14)
    plt.legend(fontsize=10)
    format_axis(plt.gca())
    plt.tight_layout()
    plt.show()

    # Plot Nonsense mutations
    plt.figure(figsize=(8, 6))
    plt.hist(nonsense_pathogenic, bins=bins, alpha=0.6, color=colors['pathogenic'], edgecolor='black', linewidth=1,
             label='Pathogenic')
    plt.hist(nonsense_likely_pathogenic, bins=bins, alpha=0.6, color=colors['likely_pathogenic'], edgecolor='black',
             linewidth=1, label='Likely Pathogenic')
    plt.hist(nonsense_benign, bins=bins, alpha=0.6, color=colors['benign'], edgecolor='black', linewidth=1,
             label='Benign')
    plt.hist(nonsense_likely_benign, bins=bins, alpha=0.6, color=colors['likely_benign'], edgecolor='black',
             linewidth=1, label='Likely Benign')
    plt.title('Nonsense mutations', fontsize=14)
    plt.legend(fontsize=10)
    format_axis(plt.gca())
    plt.tight_layout()
    plt.show()

    # Plot Synonymous mutations
    plt.figure(figsize=(8, 6))
    plt.hist(synonymous_pathogenic, bins=bins, alpha=0.6, color=colors['pathogenic'],
             edgecolor='black', linewidth=1, label='Pathogenic')
    plt.hist(synonymous_likely_pathogenic, bins=bins, alpha=0.6, color=colors['likely_pathogenic'],
             edgecolor='black', linewidth=1, label='Likely Pathogenic')
    plt.hist(synonymous_benign, bins=bins, alpha=0.6, color=colors['benign'],
             edgecolor='black', linewidth=1, label='Benign')
    plt.hist(synonymous_likely_benign, bins=bins, alpha=0.6, color=colors['likely_benign'],
             edgecolor='black', linewidth=1, label='Likely Benign')
    plt.title('Synonymous mutations', fontsize=14)
    plt.legend(fontsize=10)
    format_axis(plt.gca())
    plt.tight_layout()
    plt.show()

    # Plot SNPs Combined
    plt.figure(figsize=(8, 6))
    plt.hist(snp_pathogenic, bins=bins, alpha=0.4, color=colors['pathogenic'],
             edgecolor='black', linewidth=1,label='Pathogenic')
    plt.hist(snp_likely_pathogenic, bins=bins, alpha=0.4, color=colors['likely_pathogenic'],
             edgecolor='black', linewidth=1, label='Likely Pathogenic')
    plt.hist(snp_benign, bins=bins, alpha=0.4, color=colors['benign'],
             edgecolor='black', linewidth=1, label='Benign')
    plt.hist(snp_likely_benign, bins=bins, alpha=0.4, color=colors['likely_benign'],
             edgecolor='black', linewidth=1, label='Likely Benign')
    plt.title('SNPs', fontsize=14)
    plt.legend(fontsize=10)
    format_axis(plt.gca())
    plt.tight_layout()
    plt.show()


def flatten(nested_list: List) -> List:
    return [item for sublist in nested_list for item in flatten(sublist)] if isinstance(nested_list, list) else [
        nested_list]
