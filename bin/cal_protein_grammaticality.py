import csv
from typing import List, Tuple
import pandas as pd
import torch
import torch.nn.functional as F
import numpy as np


def load_esm_model(model_name: str):
    repr_layer = int(model_name.split('_')[1][1:])
    model, alphabet = torch.hub.load("facebookresearch/esm:main", model_name)
    batch_converter = alphabet.get_batch_converter()
    return model.eval(), alphabet, batch_converter, repr_layer


def read_fasta_prot(file_path: str):
    data = []
    with open(file_path, 'r') as fasta_file:
        protein_name = ""
        sequence = []
        for line in fasta_file:
            if line.startswith(">"):
                if sequence:
                    data.append((protein_name, ''.join(sequence)))
                protein_name = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
        if sequence:
            data.append((protein_name, ''.join(sequence)))
    return data


def prepare_grammaticality_data(model_name: str,
                                seq_path: str,
                                output_csv_path: str) -> pd.DataFrame:
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
    data: List[Tuple[str, str]] = read_fasta_prot(seq_path)  # List of tuples (protein_name, sequence)
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

        for pos, tokens_len in enumerate(batch_lens):
            logits: torch.Tensor = results["logits"]  # Shape [batch_size, max_seq_length, alphabet_size]
            grammaticality: np.ndarray = F.softmax(logits[pos, 1:tokens_len - 1], dim=-1).cpu().numpy()
            # grammaticality: numpy.ndarray -> Shape [sequence_length, alphabet_size]

            csv_writer.writerows(grammaticality)

    # Load grammaticality data from CSV
    grammaticality: pd.DataFrame = pd.read_csv(output_csv_path)
    grammaticality = grammaticality[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']]
    return grammaticality


# Example usage in the main script
if __name__ == "__main__":
    Gene_list = ['IRF6']
    model_name = "esm2_t33_650M_UR50D"

    for gene in Gene_list:
        seq_path = f"./data/{gene}_protein.faa"
        output_csv_path = f"./data/{gene}_ESM2_grammaticality.csv"
        # Prepare grammaticality data
        grammaticality: pd.DataFrame = prepare_grammaticality_data(model_name, seq_path, output_csv_path)

