import pandas as pd


def load_clinvar_data(file_path):
    """Load ClinVar Figure from a compressed TSV file."""
    return pd.read_csv(file_path, sep='\t')


def filter_clinvar_data(df):
    """Filter ClinVar Figure based on specific criteria."""

    germline_classifications = ['Benign', 'Likely benign', 'Likely pathogenic', 'Pathogenic']
    variation_type = 'single nucleotide variant'

    filtered_df = df[
        (df['ClinicalSignificance'].isin(germline_classifications)) &  # Filter classification
        (df['Type'].str.contains(variation_type, case=False, na=False)) &  # Single nucleotide variation
        (
                (df['ReviewStatus'].str.contains('criteria provided, multiple submitters, no conflicts',
                                                 case=False, na=False)) |
                (df['ReviewStatus'].str.contains('criteria provided, single submitter',
                                                 case=False, na=False)) |
                (df['ReviewStatus'].str.contains('reviewed by expert panel',
                                                 case=False, na=False)) |
                (df['ReviewStatus'].str.contains('practice guideline',
                                                 case=False, na=False))
        )
        ]

    return filtered_df


# Extract transcript and gene info
def extract_transcript_gene(filtered_df, output_filtered_file):
    """Extract transcript and gene symbol from the Name column."""
    # Define the regular expression pattern for the type of annotation you're looking for
    pattern = r'NM_\d+\.\d+\([A-Za-z0-9]+\):c\.\d+[A-Z]?>[A-Z]?\s\(p\.[A-Za-z0-9]+\)'

    # Filter the dataframe based on the regex pattern in the "Name" column
    filtered_df = filtered_df[filtered_df['Name'].str.contains(pattern, regex=True, na=False)]
    save_filtered_data(filtered_df, output_filtered_file)
    print(f"Filtered Figure has been saved to {output_filtered_file}")

    # Define the regex pattern to extract transcript and gene symbol
    transcript_gene_pattern = r'(NM_\d+\.\d+)\(([A-Za-z0-9]+)\)'

    # Extract transcript and gene symbol
    df_transcript_info = filtered_df['Name'].str.extract(transcript_gene_pattern, expand=True)

    # Rename columns for clarity
    df_transcript_info.columns = ['Transcript', 'Gene']

    # Merge with original dataframe
    df_combined = pd.concat([filtered_df, df_transcript_info], axis=1)

    # Remove duplicates based on the Transcript and Gene columns
    df_deduplicated = df_combined.drop_duplicates(subset=['Transcript', 'Gene'])

    return df_deduplicated


# Save filtered Figure to a new CSV file
def save_filtered_data(filtered_df, output_file):
    """Save the filtered Figure to a new CSV file."""
    filtered_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    file_path = "./Figure/variant_summary.txt"
    output_filtered_file = 'filtered_clinvar_data.csv'
    output_grouped_sorted_file = "clinvar_gene_info.csv"

    clinvar_df = load_clinvar_data(file_path)

    filtered_clinvar_df = filter_clinvar_data(clinvar_df)

    df_deduplicated = extract_transcript_gene(filtered_clinvar_df, output_filtered_file)

    df_grouped_sorted = df_deduplicated.sort_values(by=['Gene'], ascending=True)[['Gene', 'Transcript']]

    df_grouped_sorted.to_csv(output_grouped_sorted_file, index=False)
    print(f"Grouped and sorted Figure has been saved to {output_grouped_sorted_file}")

    print(df_grouped_sorted)
