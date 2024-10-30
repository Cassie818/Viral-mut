from bin.cal_gene_grammaticality import read_fasta_nuc
from bin.params import codon_list
import pandas as pd
import numpy as np
from typing import List
import re
import os
import csv
import logging


def flatten(nested_list: List) -> List:
    """
    Flattens a nested list into a single-level list.

    Args:
        nested_list (List): The list to flatten.

    Returns:
        List: A flattened list.
    """
    return [item for sublist in nested_list for item in flatten(sublist)] if isinstance(nested_list, list) else [
        nested_list]


def split_into_codons(seq: str) -> List[str]:
    """
    Splits the nucleotide sequence into codons (3 bases each).

    Parameters:
    - seq: str, nucleotide sequence.

    Returns:
    - List of 3-base strings (codons).
    """
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts mutation information from the 'Name' column of a DataFrame,
    including reference and mutant codons, and site positions.

    Parameters:
    - df (pd.DataFrame): DataFrame containing a 'Name' column with mutation data.

    Returns:
    - pd.DataFrame: Updated DataFrame with additional columns for mutation details:
        - 'ncMut': Nucleotide mutation (e.g., 'c.123A>T')
        - 'ncSite': Nucleotide mutation position (integer)
        - 'Ref': Reference nucleotide (single letter)
        - 'Mut': Mutant nucleotide (single letter)
        - 'aaMut': Amino acid mutation (e.g., 'p.Arg259Pro')
        - 'aaSite': Amino acid mutation position (integer)
        - 'Ref_AA': Reference amino acid (single letter)
        - 'Mut_AA': Mutant amino acid (single letter)
    """

    # Precompile regular expressions for efficiency
    nc_mut_pattern = re.compile(r'(c\.\d+[A-Z]>[A-Z])')
    nc_site_pattern = re.compile(r'c\.(\d+)')
    ref_nc_pattern = re.compile(r'c\.\d+([A-Z])>')
    mut_nc_pattern = re.compile(r'>([A-Z])')
    aa_mut_pattern = re.compile(r'(p\.[A-Za-z]{3}\d+[A-Za-z]{3})')
    aa_site_pattern = re.compile(r'^p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}$')  # Refined pattern
    ref_aa_pattern = re.compile(r'p\.([A-Za-z]{3})\d+[A-Za-z]{3}')
    mut_aa_pattern = re.compile(r'p\.[A-Za-z]{3}\d+([A-Za-z]{3})')

    def process_row(row: pd.Series) -> pd.Series:
        """
        Processes a single row to extract mutation details.

        Parameters:
        - row (pd.Series): A row from the DataFrame.

        Returns:
        - pd.Series: The updated row with extracted mutation information.
        """
        name = row.get('Name', '')

        # Initialize new columns with NaN
        row['ncMut'] = np.nan
        row['ncSite'] = np.nan
        row['Ref'] = np.nan
        row['Mut'] = np.nan
        row['aaMut'] = np.nan
        row['aaSite'] = np.nan
        row['Ref_AA'] = np.nan
        row['Mut_AA'] = np.nan

        # Extract nucleotide mutation
        ncMut_match = nc_mut_pattern.search(name)
        if ncMut_match:
            row['ncMut'] = ncMut_match.group(1)

            # Extract nucleotide mutation position
            ncSite_match = nc_site_pattern.search(row['ncMut'])
            if ncSite_match:
                row['ncSite'] = int(ncSite_match.group(1))
            else:
                print(f"Warning: Nucleotide site not found in '{row['ncMut']}'")

            # Extract reference and mutant nucleotides
            ref_nc_match = ref_nc_pattern.search(row['ncMut'])
            mut_nc_match = mut_nc_pattern.search(row['ncMut'])
            if ref_nc_match:
                row['Ref'] = ref_nc_match.group(1)
            else:
                print(f"Warning: Reference nucleotide not found in '{row['ncMut']}'")

            if mut_nc_match:
                row['Mut'] = mut_nc_match.group(1)
            else:
                print(f"Warning: Mutant nucleotide not found in '{row['ncMut']}'")

        # Extract amino acid mutation
        aaMut_match = aa_mut_pattern.search(name)
        if aaMut_match:
            row['aaMut'] = aaMut_match.group(1)

            # Extract amino acid mutation position using refined pattern
            aaSite_match = aa_site_pattern.search(row['aaMut'])
            if aaSite_match:
                row['aaSite'] = int(aaSite_match.group(1))
            else:
                print(f"Warning: Amino acid site not found or pattern mismatch in '{row['aaMut']}'")

        return row

    # Apply the row-wise processing function
    df = df.apply(process_row, axis=1)

    return df


def calculate_llr(row: pd.Series,
                  grammaticality: pd.DataFrame,
                  gene: str) -> float:
    """
    Calculates the log-likelihood ratio (LLR) for a single mutation.

    Args:
        row (pd.Series): A row from the DataFrame containing mutation information.
        grammaticality (pd.DataFrame): DataFrame with grammaticality values for each site and amino acid.
        gene (str): The gene associated with the mutation.

    Returns:
        float: The LLR value for the given mutation.
    """
    ncsite = row['ncSite']  # The nucleotide site where the mutation occurs
    aasite = row['aaSite']  # The amino acid site in the sequence
    mtsite = (int(ncsite) - 1) % 3

    ref = row['Ref']  # Reference nucleotide
    mut = row['Mut']  # Mutant nucleotide

    # Replace thymine (T) with uracil (U) to represent RNA mutation
    mtnc = mut.replace('T', 'U')

    # Read the nucleotide sequence from a FASTA file
    seq_path = f"./data/Gene/{gene}.fasta"
    sequence = read_fasta_nuc(seq_path)[0][1]  # Extract the sequence part

    # Get reference codon by splitting the sequence into codons
    ref_codon = split_into_codons(sequence)[aasite - 1].replace('T', 'U')  # Ensure 0-based indexing for amino acid site

    # Create the mutant codon by modifying the appropriate base
    mut_codon = list(ref_codon)  # Convert to list for mutating a single nucleotide

    if 0 <= mtsite < len(mut_codon):  # Ensure the mutation is within the codon range (0, 1, 2)
        mut_codon[mtsite] = mtnc  # Modify the nucleotide at the appropriate site
    else:
        raise ValueError(f"Invalid mutation site index {mtsite} for codon '{ref_codon}'")

    mut_codon = ''.join(mut_codon)  # Convert back to string

    # Store the reference and mutant codons in the row
    row['Ref_codon'] = ref_codon
    row['Mut_codon'] = mut_codon

    # Check if codons are valid
    if ref_codon not in codon_list or mut_codon not in codon_list:
        print(
            f"Warning: Reference codon '{ref_codon}' or mutant codon '{mut_codon}' not found in codon list. Skipping!!!")
        return None
    else:
        # Retrieve grammaticality scores for the reference and mutant codons
        wt = grammaticality.iloc[aasite - 1][ref_codon]
        mt = grammaticality.iloc[aasite - 1][mut_codon]

        # Calculate log-likelihood ratio (LLR)
        llr = np.log(mt) - np.log(wt)
        print(f"Calculated LLR for Gene={gene}, Site={aasite}, Ref={ref_codon}, Mut={mut_codon}, LLR={llr}")

        return llr


def initialize_output_file(output_path: str):
    """
    Initializes the output CSV file with headers.

    Args:
        output_path (str): Path to the output CSV file.
    """
    headers = ['Label', 'Gene', 'Site', 'Ref', 'Mut', 'LLR']
    with open(output_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)


def append_batch_to_output_file(output_path: str,
                                batch_data: List):
    """
    Appends a batch of data to the output CSV file.

    Args:
        output_path (str): Path to the output CSV file.
        batch_data (List): List of rows to append, each row is a list.
    """
    with open(output_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(batch_data)


def setup_logging():
    """
    Configures the logging settings.
    """
    logging.basicConfig(
        filename='processing.log',
        filemode='a',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def process_data(data: pd.DataFrame,
                 label: str,
                 output_dir: str,
                 batch_size: int = 100):
    """
    Processes the mutation data by extracting mutation information and calculating LLR.
    Immediately saves gene site and LLR to an output file in batches.

    Args:
        data (pd.DataFrame): The input mutation data.
        label (str): Label to identify the dataset (e.g., 'benign', 'pathogenic').
        output_dir (str): Directory where the output file will be saved.
        batch_size (int): Number of rows to write in each batch.

    Returns:
        None
    """

    processed_data = extract_mutation_info(data)

    # Define output file path
    output_path = os.path.join(output_dir, f'{label}_LLR_CaLM_results.csv')

    # Initialize the output file with headers
    initialize_output_file(output_path)

    batch_data = []
    total_processed = 0
    total_skipped = 0

    for idx, row in processed_data.iterrows():
        # Use a regular expression to extract the content inside parentheses
        gene_match = re.search(r'\(([^)]+)\)', row['Name'])
        gene = gene_match.group(1) if gene_match else None

        if gene:
            grammaticality_file = f"./Results/Gene/{gene}_CaLM_grammaticality.csv"
            if not os.path.isfile(grammaticality_file):
                logging.error(f"File not found: {grammaticality_file}")
                total_skipped += 1
                continue

            try:
                grammaticality = pd.read_csv(grammaticality_file, sep=',')
            except Exception as e:
                logging.error(f"Error reading {grammaticality_file}: {e}")
                total_skipped += 1
                continue

            # Calculate LLR for the current row, passing the gene name
            llr = calculate_llr(row, grammaticality, gene)

            # Prepare data to append
            output_row = [
                label,
                gene,
                row['aaSite'],
                row['Ref_codon'],
                row['Mut_codon'],
                llr
            ]

            batch_data.append(output_row)
            total_processed += 1

            # Write batch to file if batch size is reached
            if len(batch_data) >= batch_size:
                append_batch_to_output_file(output_path, batch_data)
                logging.info(f"Appended batch of {len(batch_data)} rows to {output_path}")
                batch_data = []

    # Write any remaining data in the batch
    if batch_data:
        append_batch_to_output_file(output_path, batch_data)
        logging.info(f"Appended final batch of {len(batch_data)} rows to {output_path}")

    print(f"Completed processing for '{label}' dataset.")
    print(f"Total LLR calculations performed: {total_processed}")
    print(f"Total mutations skipped: {total_skipped}")


def main():
    # Setup logging
    setup_logging()

    # Define the directory where output files will be saved
    output_dir = "./LLR"
    os.makedirs(output_dir, exist_ok=True)

    try:
        benign_data = pd.read_csv('./data/benign_data.csv')
        logging.info("Loaded benign_data.csv successfully.")
    except Exception as e:
        print(f"Error loading 'benign_data.csv': {e}")
        logging.error(f"Error loading benign_data.csv: {e}")
        benign_data = pd.DataFrame()

    try:
        likely_benign_data = pd.read_csv('./data/likely_benign_data.csv')
        logging.info("Loaded benign_data.csv successfully.")
    except Exception as e:
        print(f"Error loading 'benign_data.csv': {e}")
        logging.error(f"Error loading benign_data.csv: {e}")
        likely_benign_data = pd.DataFrame()

    try:
        pathogenic_data = pd.read_csv('./data/pathogenic_data.csv')
        logging.info("Loaded pathogenic_data.csv successfully.")
    except Exception as e:
        print(f"Error loading 'pathogenic_data.csv': {e}")
        logging.error(f"Error loading pathogenic_data.csv: {e}")
        pathogenic_data = pd.DataFrame()

    try:
        likely_pathogenic_data = pd.read_csv('./data/likely_pathogenic_data.csv')
        logging.info("Loaded likely_pathogenic_data.csv successfully.")
    except Exception as e:
        print(f"Error loading 'likely_pathogenic_data.csv': {e}")
        logging.error(f"Error loading likely_pathogenic_data.csv: {e}")
        likely_pathogenic_data = pd.DataFrame()

    # Process both benign and pathogenic datasets
    if not benign_data.empty:
        process_data(benign_data, 'benign', output_dir)
    else:
        print("Benign data is empty. Skipping processing for benign dataset.")
        logging.warning("Benign data is empty. Skipping processing for benign dataset.")

    if not likely_benign_data.empty:
        process_data(likely_benign_data, 'likely_benign', output_dir)
    else:
        print("Likely benign data is empty. Skipping processing for likely benign dataset.")
        logging.warning("Likely benign data is empty. Skipping processing for likely benign dataset.")

    if not pathogenic_data.empty:
        process_data(pathogenic_data, 'pathogenic', output_dir)
    else:
        print("Pathogenic data is empty. Skipping processing for pathogenic dataset.")
        logging.warning("Pathogenic data is empty. Skipping processing for pathogenic dataset.")

    if not likely_pathogenic_data.empty:
        process_data(likely_pathogenic_data, 'likely_pathogenic', output_dir)
    else:
        print("Likely pathogenic data is empty. Skipping processing for likely pathogenic dataset.")
        logging.warning("Likely pathogenic data is empty. Skipping processing for likely pathogenic dataset.")

    print("All datasets have been processed.")
    logging.info("Completed processing of mutation data.")


if __name__ == "__main__":
    main()
