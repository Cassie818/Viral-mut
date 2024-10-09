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

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        List[Tuple[str, str]]: A list of tuples where each tuple contains a protein name and its corresponding sequence.
    """
    with open(file_path, 'r') as fasta_file:
        # Group the lines by whether they start with '>', which indicates a header line.
        fasta_entries = groupby(fasta_file, lambda line: line.startswith(">"))

        # Process each group to create the final list of protein names and sequences.
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
    Loads an ESM model, reads protein sequence, extracts per-residue representations,
    calculates grammaticality probabilities, and saves to CSV.

    Args:
        model_name (str): The name of the ESM model to load.
        seq_path (str): Path to the input protein sequence file in FASTA format.
        output_csv_path (str): Path to the output CSV file to save grammaticality probabilities.

    Returns:
        pd.DataFrame: A DataFrame with the calculated grammaticality for selected amino acids.
    """
    # Load ESM model and prepare data
    model, alphabet, batch_converter, repr_layer = load_esm_model(model_name)

    # Read protein sequence and get batch tokens
    data: List[Tuple[str, str]] = read_fasta(seq_path)  # List of tuples (protein_name, sequence)
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    # batch_labels: List[str] -> Protein labels
    # batch_strs: List[str] -> Protein sequences as strings
    # batch_tokens: torch.Tensor -> Tokenized protein sequences, shape [batch_size, max_seq_length]
    batch_lens: torch.Tensor = (batch_tokens != alphabet.padding_idx).sum(1)  # Length of each sequence in the batch

    # Extract per-residue representations
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[repr_layer], return_contacts=False)
        # results["logits"]: torch.Tensor -> Shape [batch_size, max_seq_length, alphabet_size]

    # Calculate grammaticality probabilities and save to CSV
    with open(output_csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        header: List[str] = [alphabet.get_tok(i) for i in range(len(alphabet))]  # Header for amino acids
        csv_writer.writerow(header)

        for i, tokens_len in enumerate(batch_lens):
            logits: torch.Tensor = results["logits"]  # Shape [batch_size, seq_length, alphabet_size]
            grammaticality: np.ndarray = F.softmax(logits[i, 1:tokens_len - 1], dim=-1).cpu().numpy()
            # grammaticality: numpy.ndarray -> Shape [sequence_length, alphabet_size]

            csv_writer.writerows(grammaticality)


if __name__ == "__main__":
    Gene_list = pd.read_csv("./data/gene_info.txt", sep="\t", header=None)[0].tolist()
    model_name = "esm2_t30_150M_UR50D"

    for gene in Gene_list:
        seq_path = f"./data/Protein/{gene}_protein.fasta"
        output_csv_path = f"./Results/{gene}_ESM2_grammaticality.csv"
        # Prepare grammaticality data
        prepare_grammaticality_data(model_name, seq_path, output_csv_path)

        # Load grammaticality data from CSV
        grammaticality: pd.DataFrame = pd.read_csv(output_csv_path)
        grammaticality = grammaticality[['A', 'C', 'D', 'E', 'F',
                                         'G', 'H', 'I', 'K', 'L',
                                         'M', 'N', 'P', 'Q', 'R',
                                         'S', 'T', 'V', 'W', 'Y']]
