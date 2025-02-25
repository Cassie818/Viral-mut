import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List
import matplotlib.ticker as ticker

def load_data(gene: str):
    protein_data = pd.read_csv(f"Results/Protein/{gene}_ESM2_grammaticality.csv",
                               sep = ',')
    gene_data = pd.read_csv(f"Results/Gene/{gene}_CaLM_grammaticality.csv",
                            sep = ',')
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
                            'AGA': 'AGA (R)', 'AGG': 'AGG (R)', 'CGA': 'CGA (R)', 'CGC': 'CGC (R)', 'CGG': 'CGG (R)', 'CGU': 'CGU (R)',
                            'AGU': 'AGU (S)', 'AGC': 'AGC (S)', 'UCA': 'UCA (S)', 'UCC': 'UCC (S)', 'UCG': 'UCG (S)', 'UCU': 'UCU (S)',
                            'UAU': 'UAU (Y)', 'UAC': 'UAC (Y)',
                            'UUA': 'UUA (L)', 'UUG': 'UUG (L)', 'CUA': 'CUA (L)', 'CUC': 'CUC (L)', 'CUG': 'CUG (L)', 'CUU': 'CUU (L)',
                            'UUU': 'UUU (F)', 'UUC': 'UUC (F)',
                            'UGU': 'UGU (C)', 'UGC': 'UGC (C)',
                            'UGG': 'UGG (W)',
                            'CAA': 'CAA (Q)', 'CAG': 'CAG (Q)',
                            'CAU': 'CAU (H)', 'CAC': 'CAC (H)',
                            'CCA': 'CCA (P)', 'CCC': 'CCC (P)', 'CCG': 'CCG (P)', 'CCU': 'CCU (P)',
                            'GAA': 'GAA (E)', 'GAG': 'GAG (E)',
                            'GAU': 'GAU (D)', 'GAC': 'GAC (D)',
                            'GUA': 'GUA (V)', 'GUG': 'GUG (V)',  'GUC': 'GUC (V)', 'GUU': 'GUU (V)',
                            'GCA': 'GCA (A)', 'GCG': 'GCG (A)', 'GCU': 'GCU (A)', 'GCC': 'GCC (A)',
                            'GGA': 'GGA (G)', 'GGG': 'GGG (G)', 'GGU': 'GGU (G)', 'GGC': 'GGC (G)',
                            'UAA': 'UAA (Stop)', 'UAG': 'UAG (Stop)', 'UGA': 'UGA (Stop)',
                            })

    row = df.iloc[site - 1]

    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(12, 6))

    colors = sns.color_palette("mako", n_colors=len(row))

    sns.barplot(x=row.index, y=row.values, palette=colors)

    plt.title(f"Grammaticality at site {site}", fontsize=12, pad=15)

    plt.xlabel("Tokens", fontsize=10, labelpad=15)
    plt.ylabel("Grammaticality values", fontsize=10, labelpad=15)

    plt.xticks(rotation=45, ha="right", fontsize=10)

    plt.grid(axis='y', linestyle='--', alpha=0.6)

    plt.gca().set_facecolor('#f7f7f7')

    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    protein_data, gene_data = load_data("PAH")

    plot_distribution(protein_data, 261, 'protein')
    plot_distribution(gene_data, 261, 'cDNA')