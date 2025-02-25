import pandas as pd
import numpy as np
from typing import List
import re
import os
import csv
import logging
from params import amino_acid_dict


def flatten(nested_list: List) -> List:
    """
    Flattens a nested list into a single-level list.
    """
    return [item for sublist in nested_list for item in flatten(sublist)] if isinstance(nested_list, list) else [
        nested_list]


def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts mutation information including reference and mutant amino acids, and site positions.

    Args: DataFrame containing mutation data.

    Returns: Updated DataFrame with additional columns for mutation details.
    """

    def process_row(row):
        # Extract aaMut
        aaMut_match = re.search(r'(p\.\w+\d+\w+)', row['Name'])
        row['aaMut'] = aaMut_match.group(1) if aaMut_match else np.nan

        # Extract Site
        aaSite_match = re.search(r'(\d+)', row['aaMut']) if pd.notna(row['aaMut']) else None
        row['aaSite'] = int(aaSite_match.group(1)) if aaSite_match else np.nan

        # Extract Ref
        ref_match = re.search(r'p\.([A-Za-z]{3})\d+', row['aaMut']) if pd.notna(row['aaMut']) else None
        row['Ref'] = amino_acid_dict.get(ref_match.group(1), np.nan) if ref_match else np.nan

        # Extract Mut
        mut_match = re.search(r'p\.[A-Za-z]{3}\d+([A-Za-z]{3})', row['aaMut']) if pd.notna(row['aaMut']) else None
        row['Mut'] = amino_acid_dict.get(mut_match.group(1), np.nan) if mut_match else np.nan

        return row

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

    Returns: Float: The LLR value for the given mutation.
    """
    aasite = row['aaSite']
    ref = row['Ref']
    mut = row['Mut']

    # Skip nonsense mutation
    if pd.isna(mut) or pd.isna(ref):
        print(f"Skipping termination mutation for Gene={gene}, Site={aasite}.")

    else:
        # Access grammaticality values using row index and column label
        wt = grammaticality.loc[aasite - 1, ref]
        mt = grammaticality.loc[aasite - 1, mut]

        # Calculate log-likelihood ratio (LLR)
        llr = np.log(mt) - np.log(wt)
        print(f"Calculated LLR for Gene={gene}, Site={aasite}, Ref={ref}, Mut={mut}: LLR={llr}")
        return llr


def initialize_output_file(output_path: str):
    """
    Initializes the output CSV file with headers.

    Args: output_path (str): Path to the output CSV file.
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
    Extracting mutation information and calculating LLR.

    Args:
        data (pd.DataFrame): The input mutation data.
        label (str): Label to identify the dataset (e.g., 'benign', 'pathogenic').
        output_dir (str): Directory where the output file will be saved.
        batch_size (int): Number of rows to write in each batch.
    """

    processed_data = extract_mutation_info(data)

    output_path = os.path.join(output_dir, f'{label}_LLR_results.csv')

    initialize_output_file(output_path)

    batch_data = []
    total_processed = 0
    total_skipped = 0

    for idx, row in processed_data.iterrows():
        # Use a regular expression to extract the content inside parentheses
        gene_match = re.search(r'\(([^)]+)\)', row['Name'])
        gene = gene_match.group(1) if gene_match else None

        if gene:
            grammaticality_file = f"./Results/Protein/{gene}_ESM2_grammaticality.csv"
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

            # Prepare Figure to append
            output_row = [
                label,
                gene,
                row['aaSite'],
                row['Ref'],
                row['Mut'],
                llr
            ]

            batch_data.append(output_row)
            total_processed += 1

            # Write batch to file if batch size is reached
            if len(batch_data) >= batch_size:
                append_batch_to_output_file(output_path, batch_data)
                logging.info(f"Appended batch of {len(batch_data)} rows to {output_path}")
                batch_data = []

    if batch_data:
        append_batch_to_output_file(output_path, batch_data)
        logging.info(f"Appended final batch of {len(batch_data)} rows to {output_path}")

    print(f"Completed processing for '{label}' dataset.")
    print(f"Total LLR calculations performed: {total_processed}")
    print(f"Total mutations skipped: {total_skipped}")


def main():

    setup_logging()

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
