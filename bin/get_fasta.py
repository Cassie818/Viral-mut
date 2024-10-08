from Bio import Entrez, SeqIO, Seq
from Bio.SeqRecord import SeqRecord

# Set your email. NCBI requires this.
Entrez.email = "ruyic818@gmail.com"  # Replace with your actual email


def write_fasta_file(records, filename):
    """
    Writes a list of SeqRecord objects to a FASTA file with sequences on a single line.

    Args:
        records (list of SeqRecord): The sequences to write.
        filename (str): The path to the output FASTA file.
    """
    try:
        with open(filename, 'w') as f:
            for record in records:
                # Write the header
                f.write(f">{record.id} {record.description}\n")
                # Write the sequence as a single line
                f.write(f"{record.seq}\n")
        print(f"Sequences have been saved to {filename}")
    except Exception as e:
        print(f"An error occurred while writing to {filename}: {e}")


def fetch_cds_and_protein_sequences(accession):
    try:
        # Fetch the GenBank record
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")

        # Extract CDS features
        cds_features = [feature for feature in record.features if feature.type == "CDS"]
        if not cds_features:
            print("No CDS features found in the record.")
            return

        # Prepare lists to hold CDS and protein SeqRecords
        cds_records = []
        protein_records = []

        # Iterate over CDS features and extract information
        for idx, cds in enumerate(cds_features, 1):
            qualifiers = cds.qualifiers
            protein_id = qualifiers.get("protein_id", [f"CDS_{idx}"])[0]
            product = qualifiers.get("product", ["N/A"])[0]
            gene = qualifiers.get("gene", ["N/A"])[0]
            location = str(cds.location)

            # Extract CDS nucleotide sequence
            cds_seq = cds.extract(record.seq)

            # Create a SeqRecord for CDS
            cds_record = SeqRecord(
                cds_seq,
                id=protein_id,
                description=f"Gene: {gene}; Product: {product}; Location: {location}"
            )
            cds_records.append(cds_record)

            # Attempt to extract the protein translation from qualifiers
            if "translation" in qualifiers:
                protein_seq = qualifiers["translation"][0]
                protein_record = SeqRecord(
                    Seq.Seq(protein_seq),
                    id=protein_id,
                    description=f"Gene: {gene}; Product: {product}; Location: {location}"
                )
            else:
                # If translation is not available, perform translation
                protein_seq = cds_seq.translate(to_stop=True)
                protein_record = SeqRecord(
                    protein_seq,
                    id=protein_id,
                    description=f"Gene: {gene}; Product: {product}; Location: {location}; Note: Translation performed manually."
                )

            protein_records.append(protein_record)

        # Determine the gene name for file naming
        # If multiple genes are present, you might need to adjust this logic
        # Here, we'll use the first gene name found; otherwise, default to 'UnknownGene'
        gene_names = [record.qualifiers.get("gene", ["UnknownGene"])[0] for record in cds_features]
        gene = gene_names[0] if gene_names else "UnknownGene"

        # Define output file paths
        cds_output_file = f"../data/ClinVar/gene/{gene}.fasta"
        protein_output_file = f"../data/ClinVar/protein/{gene}_protein.fasta"

        # Write CDS sequences to FASTA without wrapping
        write_fasta_file(cds_records, cds_output_file)

        # Write protein sequences to FASTA without wrapping
        write_fasta_file(protein_records, protein_output_file)

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    # Open the text file and read all the accession numbers
    with open("accession_numbers.txt", "r") as file:
        accession_numbers = file.read().splitlines()  # Read each line, remove newline characters

    # Loop through each accession number and call the function
    for accession_number in accession_numbers:
        fetch_cds_and_protein_sequences(accession_number)
