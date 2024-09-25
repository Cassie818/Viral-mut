import pandas as pd
import torch
import esm
import torch.nn.functional as F
import csv
from params import amino_acid_dict


def extract_mutation_info(df: pd.DataFrame) -> pd.DataFrame:
    # Extract the amino acid mutation information
    df['aaMut'] = df['Name'].str.extract(r'(p\.\w+\d+\w+)')

    # Extract the numeric part site from aaMut
    df['aaSite'] = df['aaMut'].str.extract(r'(\d+)')

    # Extract the reference nucleotide and the mutated nucleotide from ncMut
    df['Ref'] = df['aaMut'].str.extract(r'p\.([A-Za-z]{3})\d+')
    df['Mut'] = df['aaMut'].str.extract(r'p\.[A-Za-z]{3}\d+([A-Za-z]{3})')
    df['Ref'] = df['Ref'].map(amino_acid_dict)
    df['Mut'] = df['Mut'].map(amino_acid_dict)

    return df


def read_fasta_prot(file_path):
    data = []
    protein_name = ""
    sequence = []

    with open(file_path, 'r') as fasta_file:

        for line in fasta_file:
            if line.startswith(">"):
                protein_name = line[1:-1]
                sequence = []
            else:
                sequence.append(line)

        data.append((protein_name, ''.join(sequence)))

    return data


def load_esm_model(model_name):
    repr_layer = int(model_name.split('_')[1][1:])
    model, alphabet = torch.hub.load("facebookresearch/esm:main", model_name)
    batch_converter = alphabet.get_batch_converter()
    return model.eval(), alphabet, batch_converter, repr_layer


model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t30_150M_UR50D")
# model,alphabet,batch_converter,repr_layer = load_esm_model(model_name='esm1b_t33_650M_UR50S')
batch_converter = alphabet.get_batch_converter()
model.eval()

batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

probabilities = []
for pos, tokens_len in enumerate(batch_lens):
    logit = results["logits"]
    # remove start token
    grammaticality = F.softmax(logit[pos, 1: tokens_len - 1], dim=-1).cpu().numpy()
    probabilities.append(grammaticality)

csv_fname = "./data/BRCA1_ESM2_grammaticality.csv"

with open(csv_fname, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    header = amino_acids
    csv_writer.writerow(header)
    for probs in probabilities:

        for row in probs:
            csv_writer.writerow(row)
