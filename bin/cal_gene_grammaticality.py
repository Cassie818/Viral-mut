from typing import List, Tuple, Union
from calm.sequence import CodonSequence
import torch.nn.functional as F
from calm import CaLM
import numpy as np
import pandas as pd
import csv


def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts nucleotide and amino acid mutation information from the 'Name' column.

    Parameters:
    - df: pandas DataFrame containing a column 'Name' with mutation information.

    Returns:
    - df: pandas DataFrame with additional columns for extracted mutation details.
    """
    # Extract the nucleotide mutation information
    df['ncMut'] = df['Name'].str.extract(r'(c\.\d+[A-Z]>[A-Z])')

    # Extract the amino acid mutation information
    df['aaMut'] = df['Name'].str.extract(r'(p\.\w+\d+\w+)')

    # Extract the numeric part (site) from ncMut and aaMut
    df['ncSite'] = df['ncMut'].str.extract(r'(\d+)')
    df['aaSite'] = df['aaMut'].str.extract(r'(\d+)')

    # Extract the reference nucleotide and the mutated nucleotide from ncMut
    df['Ref'] = df['ncMut'].str.extract(r'([A-Z])>')
    df['Mut'] = df['ncMut'].str.extract(r'>([A-Z])')

    return df


def read_fasta_nuc(file_path: str) -> List[Tuple[str, str]]:
    """
    Reads a nucleotide sequence from a FASTA file.

    Parameters:
    - file_path: string, path to the FASTA file.

    Returns:
    - sequences: list of tuples containing the ID and sequence from the FASTA file.
    """
    sequences: List[Tuple[str, str]] = []
    current_id: str = ""
    current_sequence: List[str] = []

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()  # Remove whitespace
            if line.startswith(">"):
                # Save the previous sequence before starting a new one
                if current_id and current_sequence:
                    sequences.append((current_id, "".join(current_sequence)))
                current_id = line[1:]  # Extract ID from the FASTA header
                current_sequence = []  # Reset sequence
            else:
                current_sequence.append(line)

        # Add the last sequence after the loop ends
        if current_id and current_sequence:
            sequences.append((current_id, "".join(current_sequence)))

    return sequences


def _split_into_codons(seq: str) -> List[str]:
    for i in range(0, len(seq), 3):
        yield seq[i:i + 3]


def split_into_codons(seq: str) -> List[str]:
    """Yield successive 3-letter chunks of a string/sequence."""
    return list(_split_into_codons(seq))


class CaLMPluS(CaLM):

    def get_logits(self, sequence: Union[str, 'CodonSequence']) -> np.ndarray:
        """
        Calculate and return the logits for a given input sequence.

        Parameters:
        - sequence: Union[str, CodonSequence]
            The input sequence for which logits are to be calculated. It can be either
            a string (representing nucleotide bases) or an instance of the CodonSequence class.

        Returns:
        - logits: np.ndarray
            The predicted logits after applying softmax, indicating probabilities across
            possible classes.

        Raises:
        - ValueError: If the input sequence is not a string or CodonSequence.
        """
        # Check if the input sequence is a string or CodonSequence instance.
        # Convert string to a CodonSequence if necessary.
        if isinstance(sequence, str):
            seq: CodonSequence = CodonSequence(sequence)  # Convert string to CodonSequence type.
        elif isinstance(sequence, CodonSequence):
            seq: CodonSequence = sequence
        else:
            # Raise an error if the input is neither a string nor a CodonSequence.
            raise ValueError('The input sequence must be a string or a CodonSequence instance.')

        # Tokenize the codon sequence into a format that the model can understand.
        tokens = self.tokenize(seq)  # Expected to be a tensor or similar type that can be used by the model.

        # Pass the tokens to the model to generate output.
        output = self.model(tokens)

        # Extract logits (raw prediction scores) from the model output.
        logits = output['logits']

        # Apply softmax to convert logits into probabilities
        logits = F.softmax(logits, dim=-1).detach().numpy()

        return logits


if __name__ == "__main__":
    calm_model = CaLMPluS()
    gene_list = ['IRF6']

    for gene in gene_list:
        seq_path = f"./data/{gene}.faa"
        sequence = read_fasta_nuc(seq_path)[0][1]
        # remove the start token and end token
        logits = calm_model.get_logits(sequence)[:, 1:-1, ]
        print("Logits shape:", logits.shape)

        tok_to_idx = calm_model.alphabet.tok_to_idx
        codons = [i for i in tok_to_idx]

        csv_fname = f"./data/{gene}_CaLM_grammaticality.csv"

        with open(csv_fname, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            header = codons
            csv_writer.writerow(header)
            for probs in logits:
                for row in probs:
                    csv_writer.writerow(row)
