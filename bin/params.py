amino_acid_dict = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
}


codon_to_amino_acid = {
    'AAA': 'K',   # Lysine
    'AAU': 'N',   # Asparagine
    'AAC': 'N',   # Asparagine
    'AAG': 'K',   # Lysine
    'AUA': 'I',   # Isoleucine
    'AUU': 'I',   # Isoleucine
    'AUC': 'I',   # Isoleucine
    'AUG': 'M',   # Methionine (Start codon)
    'ACA': 'T',   # Threonine
    'ACU': 'T',   # Threonine
    'ACC': 'T',   # Threonine
    'ACG': 'T',   # Threonine
    'AGA': 'R',   # Arginine
    'AGU': 'S',   # Serine
    'AGC': 'S',   # Serine
    'AGG': 'R',   # Arginine
    'UAA': '*',   # Stop codon
    'UAU': 'Y',   # Tyrosine
    'UAC': 'Y',   # Tyrosine
    'UAG': '*',   # Stop codon
    'UUA': 'L',   # Leucine
    'UUU': 'F',   # Phenylalanine
    'UUC': 'F',   # Phenylalanine
    'UUG': 'L',   # Leucine
    'UCA': 'S',   # Serine
    'UCU': 'S',   # Serine
    'UCC': 'S',   # Serine
    'UCG': 'S',   # Serine
    'UGA': '*',   # Stop codon
    'UGU': 'C',   # Cysteine
    'UGC': 'C',   # Cysteine
    'UGG': 'W',   # Tryptophan
    'CAA': 'Q',   # Glutamine
    'CAU': 'H',   # Histidine
    'CAC': 'H',   # Histidine
    'CAG': 'Q',   # Glutamine
    'CUA': 'L',   # Leucine
    'CUU': 'L',   # Leucine
    'CUC': 'L',   # Leucine
    'CUG': 'L',   # Leucine
    'CCA': 'P',   # Proline
    'CCU': 'P',   # Proline
    'CCC': 'P',   # Proline
    'CCG': 'P',   # Proline
    'CGA': 'R',   # Arginine
    'CGU': 'R',   # Arginine
    'CGC': 'R',   # Arginine
    'CGG': 'R',   # Arginine
    'GAA': 'E',   # Glutamic
    'GAU': 'D',   # Aspartic
    'GAC': 'D',   # Aspartic
    'GAG': 'E',   # Glutamic
    'GUA': 'V',   # Valine
    'GUU': 'V',   # Valine
    'GUC': 'V',   # Valine
    'GUG': 'V',   # Valine
    'GCA': 'A',   # Alanine
    'GCU': 'A',   # Alanine
    'GCC': 'A',   # Alanine
    'GCG': 'A',   # Alanine
    'GGA': 'G',   # Glycine
    'GGU': 'G',   # Glycine
    'GGC': 'G',   # Glycine
    'GGG': 'G'    # Glycine
}

amino_acid_to_codons = {
    'K': ['AAA', 'AAG'],                                 # Lysine
    'N': ['AAU', 'AAC'],                                 # Asparagine
    'I': ['AUA', 'AUU', 'AUC'],                          # Isoleucine
    'M': ['AUG'],                                        # Methionine
    'T': ['ACA', 'ACU', 'ACC', 'ACG'],                   # Threonine
    'R': ['AGA', 'AGG', 'CGA', 'CGU', 'CGC', 'CGG'],     # Arginine
    'S': ['AGU', 'AGC', 'UCA', 'UCU', 'UCC', 'UCG'],     # Serine
    '*': ['UAA', 'UAG', 'UGA'],                          # Stop codons
    'Y': ['UAU', 'UAC'],                                 # Tyrosine
    'L': ['UUA', 'UUG', 'CUA', 'CUU', 'CUC', 'CUG'],     # Leucine
    'F': ['UUU', 'UUC'],                                 # Phenylalanine
    'C': ['UGU', 'UGC'],                                 # Cysteine
    'W': ['UGG'],                                        # Tryptophan
    'Q': ['CAA', 'CAG'],                                 # Glutamine
    'H': ['CAU', 'CAC'],                                 # Histidine
    'P': ['CCA', 'CCU', 'CCC', 'CCG'],                   # Proline
    'E': ['GAA', 'GAG'],                                 # Glutamic Acid
    'D': ['GAU', 'GAC'],                                 # Aspartic Acid
    'V': ['GUA', 'GUU', 'GUC', 'GUG'],                   # Valine
    'A': ['GCA', 'GCU', 'GCC', 'GCG'],                   # Alanine
    'G': ['GGA', 'GGU', 'GGC', 'GGG']                    # Glycine
}


# 64 codons
codon_list = ['AAA', 'AAU', 'AAC', 'AAG',
              'AUA', 'AUU', 'AUC', 'AUG',
              'ACA', 'ACU', 'ACC', 'ACG',
              'AGA', 'AGU', 'AGC', 'AGG',
              'UAA', 'UAU', 'UAC', 'UAG',
              'UUA', 'UUU', 'UUC', 'UUG',
              'UCA', 'UCU', 'UCC', 'UCG',
              'UGA', 'UGU', 'UGC', 'UGG',
              'CAA', 'CAU', 'CAC', 'CAG',
              'CUA', 'CUU', 'CUC', 'CUG',
              'CCA', 'CCU', 'CCC', 'CCG',
              'CGA', 'CGU', 'CGC', 'CGG',
              'GAA', 'GAU', 'GAC', 'GAG',
              'GUA', 'GUU', 'GUC', 'GUG',
              'GCA', 'GCU', 'GCC', 'GCG',
              'GGA', 'GGU', 'GGC', 'GGG']