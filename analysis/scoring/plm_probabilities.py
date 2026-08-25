import argparse
import csv
from pathlib import Path
from typing import List, Tuple
import torch
import torch.nn.functional as F
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export residue probabilities from an ESM model.")
    parser.add_argument("--gene-list", required=True, type=Path, help="One gene identifier per line.")
    parser.add_argument("--fasta-dir", required=True, type=Path, help="Directory containing <gene>_protein.fasta files.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--model", default="esm2_t30_150M_UR50D")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    genes = [line.strip().split("\t")[0] for line in args.gene_list.read_text().splitlines() if line.strip()]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for gene in genes:
        prepare_grammaticality_data(
            args.model,
            str(args.fasta_dir / f"{gene}_protein.fasta"),
            str(args.output_dir / f"{gene}_ESM2_grammaticality.csv"),
        )
