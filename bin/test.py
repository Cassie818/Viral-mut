import argparse
import sys


def clean_fasta(input_file, output_file):
    """
    Process a FASTA file by removing spaces, newline characters, and numbers,
    and converting nucleotide characters to uppercase.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output cleaned FASTA file.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            sequence = []
            header = None
            line_number = 0
            for line in infile:
                line_number += 1
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header:
                        # Write the previous sequence
                        outfile.write(header + '\n')
                        outfile.write(''.join(sequence) + '\n')
                        print(f"Processed sequence: {header}, Length: {len(''.join(sequence))}")
                        sequence = []
                    header = line
                else:
                    # Remove all spaces, digits, and convert to uppercase
                    cleaned_seq = ''.join(filter(lambda x: not x.isdigit(), line.replace(" ", "").upper()))
                    sequence.append(cleaned_seq)
            # Write the last sequence
            if header and sequence:
                outfile.write(header + '\n')
                outfile.write(''.join(sequence) + '\n')
                print(f"Processed sequence: {header}, Length: {len(''.join(sequence))}")
        print(f"\nProcessing complete! Output file generated: {output_file}")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found. Please check the file path.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Process a FASTA file.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input FASTA file")
    parser.add_argument('-o', '--output', required=True, help="Path to the output cleaned FASTA file")
    args = parser.parse_args()
    clean_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
