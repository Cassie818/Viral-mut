import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from bin.params import amino_acid_dict
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


# Function to extract mutation information from the dataframe
def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts mutation information including reference and mutant amino acids, and site positions.

    Args:
        df (pd.DataFrame): Dataframe containing mutation data.

    Returns:
        pd.DataFrame: Updated dataframe with additional columns for mutation details.
    """
    # Extract amino acid mutation information
    df['aaMut'] = df['Name'].str.extract(r'(p\.\w+\d+\w+)')
    df['aaSite'] = df['aaMut'].str.extract(r'(\d+)').astype(int)  # Numeric site positions
    df['Ref'] = df['aaMut'].str.extract(r'p\.([A-Za-z]{3})\d+').map(amino_acid_dict)  # Reference amino acid
    df['Mut'] = df['aaMut'].str.extract(r'p\.[A-Za-z]{3}\d+([A-Za-z]{3})').map(amino_acid_dict)  # Mutant amino acid
    return df


# Function to calculate log-likelihood ratios
def calculate_llr(df: pd.DataFrame, grammaticality: pd.DataFrame) -> List[float]:
    """
    Calculates the log-likelihood ratio (LLR) for each mutation.

    Args:
        df (pd.DataFrame): Dataframe containing mutation information.
        grammaticality (pd.DataFrame): Dataframe with grammaticality values for each site and amino acid.

    Returns:
        List[float]: List of LLR values for the given mutations.
    """
    LLR = []
    for _, row in df.iterrows():
        aasite = row['aaSite']
        ref = row['Ref']
        mut = row['Mut']

        # Retrieve the grammaticality values for the reference and mutant codons
        wt = grammaticality.iloc[aasite - 1][ref]
        mt = grammaticality.iloc[aasite - 1][mut]

        # Calculate log-likelihood ratio (LLR)
        llr = np.log(mt) - np.log(wt)
        LLR.append(llr)

    return LLR


if __name__ == "__main__":
    # List of genes to be processed
    Gene_list = ['IRF6']

    # Lists to store LLR values for benign and pathogenic mutations
    LLR_Pathogenic_list = []
    LLR_benign_list = []

    for gene in Gene_list:
        # Load grammaticality data
        grammaticality_file = f"./data/{gene}_ESM2_grammaticality.csv"
        grammaticality = pd.read_csv(grammaticality_file, sep=',')

        # Load benign and pathogenic mutation data
        benign_info_file = f"./data/{gene}_benign.txt"
        pathogenic_info_file = f"./data/{gene}.txt"

        benign_info = pd.read_csv(benign_info_file, sep='\t')
        pathogenic_info = pd.read_csv(pathogenic_info_file, sep='\t')

        # Extract mutation information
        benign_mutations = extract_mutation_info(benign_info)
        pathogenic_mutations = extract_mutation_info(pathogenic_info)

        # Calculate LLR for benign and pathogenic mutations
        LLR_benign = calculate_llr(benign_mutations, grammaticality)
        LLR_pathogenic = calculate_llr(pathogenic_mutations, grammaticality)

        # Append results to the lists
        LLR_benign_list.append(LLR_benign)
        LLR_Pathogenic_list.append(LLR_pathogenic)
