import pandas as pd
import numpy as np
from bin.helpers import flatten
import pickle
from bin.cal_gene_grammaticality import read_fasta_nuc

from typing import List, Tuple


def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts detailed mutation information from the 'Name' column of a DataFrame.

    Parameters:
    - df: pandas DataFrame containing a column 'Name' with mutation information.

    Returns:
    - df: pandas DataFrame with additional columns for nucleotide mutation,
          amino acid mutation, numeric sites, and reference/mutant nucleotides.
    """
    # Extract nucleotide mutation in format 'c.<number><nucleotide1>><nucleotide2>'
    df['ncMut'] = df['Name'].str.extract(r'(c\.\d+[A-Z]>[A-Z])')

    # Extract amino acid mutation in format 'p.<amino acid><number><amino acid>'
    df['aaMut'] = df['Name'].str.extract(r'(p\.\w+\d+\w+)')

    # Extract the numeric part of nucleotide and amino acid mutation sites
    df['ncSite'] = df['ncMut'].str.extract(r'(\d+)')
    df['aaSite'] = df['aaMut'].str.extract(r'(\d+)')

    # Extract the reference nucleotide and mutated nucleotide from ncMut
    df['Ref'] = df['ncMut'].str.extract(r'([A-Z])>')
    df['Mut'] = df['ncMut'].str.extract(r'>([A-Z])')

    return df


def calculate_llr(df: pd.DataFrame,
                  sequence: str,
                  grammaticality: pd.DataFrame) -> List[float]:
    """
    Calculates the log-likelihood ratio (LLR) for each mutation in the given DataFrame.

    Parameters:
    - df: pandas DataFrame containing mutation information.
    - sequence: str, the nucleotide sequence of the gene.
    - grammaticality: pandas DataFrame containing grammaticality scores for each codon.

    Returns:
    - LLR: List of log-likelihood ratios (LLRs) for each mutation.
    """
    LLR = []

    # Define the split_into_codons function to divide the sequence into triplets (codons)
    def split_into_codons(seq: str) -> List[str]:
        """
        Splits the nucleotide sequence into codons (3 bases each).

        Parameters:
        - seq: str, nucleotide sequence.

        Returns:
        - List of 3-base strings (codons).
        """
        return [seq[i:i + 3] for i in range(0, len(seq), 3)]

    # Iterate over each row of the mutation DataFrame
    for index, row in df.iterrows():
        # Extract numeric part of 'ncSite' and 'aaSite' columns (mutation locations)
        ncs = ''.join(filter(str.isdigit, row['ncSite']))
        aas = ''.join(filter(str.isdigit, row['aaSite']))

        # Replace thymine (T) with uracil (U) to represent RNA mutation
        mtaa = row['Mut'].replace('T', 'U')

        # Calculate mutation site position within the codon
        mtsite = int(ncs) % 3 - 1  # Find the position within the codon (0, 1, or 2)
        aasite = int(aas)  # Amino acid site index

        # Get reference codon and replace 'T' with 'U' to represent RNA
        ref_codon = split_into_codons(sequence)[aasite - 1].replace('T', 'U')

        # Create the mutant codon by modifying the appropriate base
        if mtsite < 0:
            # If mtsite is negative, modify the third base of the codon
            mut_codon = ref_codon[:2] + mtaa + ref_codon[3:]
        else:
            # Otherwise, modify the corresponding base position in the codon
            mut_codon = ref_codon[:mtsite] + mtaa + ref_codon[mtsite + 1:]

        # Retrieve grammaticality scores for the reference and mutant codons
        wt = grammaticality.iloc[aasite - 1][ref_codon]
        mt = grammaticality.iloc[aasite - 1][mut_codon]

        # Calculate log-likelihood ratio (LLR)
        llr = np.log(mt) - np.log(wt)
        LLR.append(llr)

    return LLR


if __name__ == "__main__":
    # List of genes to be processed
    Gene_list = ['IRF6']
    LLR_Pathogenic_list = []
    LLR_benign_list = []

    # Process each gene in the Gene_list
    for gene in Gene_list:
        # Read grammaticality scores for the given gene from a CSV file
        file = f"./data/{gene}_CaLM_grammaticality.csv"
        grammaticality = pd.read_csv(file, sep=',')

        # Read benign mutation information from text file
        benign_info = pd.read_csv(f"./data/{gene}_benign.txt", sep='\t')

        # Read pathogenic mutation information from text file
        Pathogenic_info = pd.read_csv(f"./data/{gene}.txt", sep='\t')

        # Extract detailed mutation information for benign and pathogenic mutations
        benign = extract_mutation_info(benign_info)
        Pathogenic = extract_mutation_info(Pathogenic_info)

        # Read the nucleotide sequence from a FASTA file
        seq_path = f"./data/{gene}.faa"
        sequence = read_fasta_nuc(seq_path)[0][1]  # Extract the sequence part

        # Calculate LLR for benign and pathogenic mutations
        LLR_benign = calculate_llr(benign, sequence, grammaticality)
        LLR_Pathogenic = calculate_llr(Pathogenic, sequence, grammaticality)

        # Append LLR values to respective lists
        LLR_benign_list.append(LLR_benign)
        LLR_Pathogenic_list.append(LLR_Pathogenic)

    benign_list = flatten(LLR_benign_list)
    Pathogenic_list = flatten(LLR_Pathogenic_list)

    # Save files
    with open('./data/CaLM_LLR_benign.txt', 'wb') as f:
        pickle.dump(benign_list, f)

    with open('./data/CaLM_LLR_Pathogenic.txt', 'wb') as f:
        pickle.dump(Pathogenic_list, f)
