from typing import List, Tuple, Union
from calm.sequence import CodonSequence
import torch.nn.functional as F
from calm import CaLM
import numpy as np
import pandas as pd
import csv
import torch


class CaLMPluS(CaLM):

    def get_logits(self,
                   sequence: Union[str, 'CodonSequence']) -> np.ndarray:
        """
        Calculate and return the logits for a given input sequence.

        Args:
        - sequence: Union[str, CodonSequence]
            The input sequence for which logits are to be calculated. It can be either
            a string (representing nucleotide bases) or an instance of the CodonSequence class.

        Returns:
        - logits: np.ndarray
            The predicted logits after applying softmax, indicating probabilities across possible classes.
        """
        # Check if the input sequence is a string or CodonSequence instance.
        if isinstance(sequence, str):
            seq: CodonSequence = CodonSequence(sequence)  # Convert string to CodonSequence type.
        elif isinstance(sequence, CodonSequence):
            seq: CodonSequence = sequence
        else:
            # Raise an error if the input is neither a string nor a CodonSequence.
            raise ValueError('The input sequence must be a string or a CodonSequence instance.')

        # Tokenize the codon sequence into a format that the model can understand.
        tokens = self.tokenize(seq)

        with torch.no_grad():

            output = self.model(tokens)
            logits = output['logits']
            logits = F.softmax(logits, dim=-1)

            return logits.detach().cpu().numpy()


def read_fasta_nuc(file_path: str) -> List[Tuple[str, str]]:
    """
    Reads nucleotide sequences from a FASTA file and returns a list of tuples containing IDs and sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        List[Tuple[str, str]]: A list of tuples where each tuple contains an ID and its corresponding nucleotide sequence.
    """
    sequences: List[Tuple[str, str]] = []
    current_id: str = ""
    current_sequence: List[str] = []

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences.append((current_id, "".join(current_sequence)))
                current_id = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_id:
            sequences.append((current_id, "".join(current_sequence)))

    return sequences


if __name__ == "__main__":
    calm = CaLMPluS()
    gene_list = pd.read_csv("../data/gene_info.txt", sep="\t", header=None)[0].tolist()

    for gene in gene_list:

        seq_path = f"../data/Gene/{gene}.fasta"
        sequence = read_fasta_nuc(seq_path)[0][1]
        # remove the start token and end token
        logits = calm.get_logits(sequence)[:, 1:-1, ]

        tok_to_idx = calm.alphabet.tok_to_idx
        codons = [i for i in tok_to_idx]

        csv_fname = f"../Results/{gene}_CaLM_grammaticality.csv"

        with open(csv_fname, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            header = codons
            csv_writer.writerow(header)
            for probs in logits:
                for row in probs:
                    csv_writer.writerow(row)
