import csv
from typing import List, Tuple
import pandas as pd
import torch
import torch.nn.functional as F
import numpy as np
from itertools import groupby


def load_esm_model(model_name: str):
    repr_layer = int(model_name.split('_')[1][1:])
    model, alphabet = torch.hub.load("facebookresearch/esm:main", model_name)
    batch_converter = alphabet.get_batch_converter()
    return model.eval(), alphabet, batch_converter, repr_layer


def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Reads a FASTA file and returns a list of tuples containing protein names and sequences.

    Args: file_path (str): Path to the FASTA file.

    Returns: List[Tuple[str, str]]: A list of tuples where each tuple contains a protein name and its corresponding sequence.
    """
    with open(file_path, 'r') as fasta_file:
        fasta_entries = groupby(fasta_file, lambda line: line.startswith(">"))

        return [
            (header[1:].strip(), ''.join(seq.strip() for seq in sequence_group).upper())
            for is_header, header_lines in fasta_entries if is_header
            for header in header_lines
            for _, sequence_group in fasta_entries
        ]


def prepare_grammaticality_data(model_name: str,
                                seq_path: str,
                                output_csv_path: str):
    """
    Load an ESM-2 model, reads protein sequence, calculates grammaticality probabilities, and saves to CSV.

    Args:
        model_name (str): The name of the ESM model to load.
        seq_path (str): Path to the input protein sequence file in FASTA format.
        output_csv_path (str): Path to the output CSV file to save grammaticality probabilities.

    Returns: pd.DataFrame: A DataFrame with the calculated grammaticality for selected amino acids.
    """
    # Load ESM model and prepare data
    model, alphabet, batch_converter, repr_layer = load_esm_model(model_name)

    # Read protein sequence and get batch tokens
    # List of tuples (protein_name, sequence)
    data: List[Tuple[str, str]] = read_fasta(seq_path)

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    # Length of each sequence in the batch
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[repr_layer], return_contacts=False)

    # Calculate grammaticality probabilities and save to CSV
    with open(output_csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Header for amino acids
        header = [alphabet.get_tok(i) for i in range(len(alphabet))]
        csv_writer.writerow(header)

        for i, tokens_len in enumerate(batch_lens):
            # Shape: (batch_size, seq_length, alphabet_size)
            logits = results["logits"]
            # Shape: (sequence_length, alphabet_size)
            grammaticality = F.softmax(logits[i, 1:tokens_len - 1], dim=-1).cpu().numpy()

            csv_writer.writerows(grammaticality)


if __name__ == "__main__":
    Gene_list = pd.read_csv("./bin/gene_info.txt", sep="\t", header=None)[0].tolist()
    model_name = "esm2_t30_150M_UR50D"

    for gene in Gene_list:
        seq_path = f"./data/Protein/{gene}_protein.fasta"
        output_csv_path = f"./Results/{gene}_ESM2_grammaticality.csv"
        # Prepare grammaticality data
        prepare_grammaticality_data(model_name, seq_path, output_csv_path)

        # Load grammaticality data from CSV
        grammaticality: pd.DataFrame = pd.read_csv(output_csv_path)
        # Save grammaticality data for 20 common amino acids
        grammaticality = grammaticality[['A', 'C', 'D', 'E', 'F',
                                         'G', 'H', 'I', 'K', 'L',
                                         'M', 'N', 'P', 'Q', 'R',
                                         'S', 'T', 'V', 'W', 'Y']]
