from Bio import AlignIO
import pandas as pd
import numpy as np
from datetime import datetime
import os
import argparse

NUM_SITES = 29903  # Total number of nucleotides in the original Wuhan COVID 2019 reference genome

def calculate_evolutionary_metrics(variantName):
    # Define the path to the alignment file
    alignment_file = os.path.join("..", "data", "aligned-sequences", f"{variantName}-aligned.txt")
    
    # Print the absolute path for debugging
    print("Alignment file path:", os.path.abspath(alignment_file))

    # Check if the file exists
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(f"The file {alignment_file} does not exist.")

    # Read the alignment file
    alignment = AlignIO.read(alignment_file, "clustal")
    reference_seq = alignment[0].seq  # Use the first sequence as reference
    data = []

    for record in alignment[1:]:
        sequence_id = record.id
        print("Processing Sequence ID:", record.id)  # Check each ID

        # Extract the date from the sequence ID
        sequence_date = record.id.split("###")[-1]  
        sequence_date = datetime.strptime(sequence_date, "%Y-%m-%d")

        # Calculate identity, substitutions, insertions, deletions
        substitutions = sum(1 for ref, base in zip(reference_seq, record.seq) if ref != base and base != "-")
        insertions = record.seq.count("-")  # Count gaps in the aligned sequence
        deletions = reference_seq.count("-")  # Count gaps in the reference sequence

        # Calculate percentage identity
        matches = sum(1 for ref, base in zip(reference_seq, record.seq) if ref == base)
        percent_identity = (matches / len(reference_seq)) * 100

        # Calculate mutation density (mutations per kb)
        mutation_density = (substitutions + insertions) / (len(reference_seq) / 1000)

        # Append the metrics
        data.append({
            "Date": sequence_date,
            "Sequence_ID": sequence_id,
            "Percent_Identity": percent_identity,
            "Substitutions": substitutions,
            "Insertions": insertions,
            "Deletions": deletions,
            "Total_Mutations": substitutions + insertions,
            "Mutation_Density_per_kb": mutation_density
        })

    # Sort by date and reset the index
    df = pd.DataFrame(data).sort_values("Date").reset_index(drop=True) 
    return df

def calculate_jukes_cantor_distance(df, total_sites):
    # Calculate the Jukes-Cantor distance for each sequence
    df['Jukes_Cantor_Distance'] = df['Substitutions'].apply(lambda x: jukes_cantor(x, total_sites))
    return df

def jukes_cantor(substitutions, total_sites):
    # Proportion of substitutions (p)
    p = substitutions / total_sites
    
    # Check if p exceeds theoretical limit for Jukes-Cantor
    if p >= 0.75:
        return np.inf  # Beyond model's accuracy range, assume infinite distance
    
    # Apply Jukes-Cantor formula
    try:
        d = -(3/4) * np.log(1 - (4/3) * p)
    except ValueError:
        d = np.nan  # Handle cases where log input might be out of range
    
    return d

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate evolutionary metrics for a given variant.")
    parser.add_argument("variantName", type=str, help="Name of the variant (e.g., Omicron)")

    # Parse arguments
    args = parser.parse_args()

    # Usage
    df = calculate_evolutionary_metrics(args.variantName)
    total_sites = NUM_SITES 
    df = calculate_jukes_cantor_distance(df, total_sites)

    # Display the DataFrame with the new Jukes-Cantor distance column
    pd.set_option('display.max_columns', None)

    print("DataFrame with Jukes-Cantor distance:")
    print("-" * 100)
    print(df)

    # Save the DataFrame in JSON format
    output_dir = os.path.join("..", "data", "analyzed-sequences")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"analyzed-{args.variantName}.json")
    df.to_json(output_file, orient='records', date_format='iso')

    print(f"DataFrame saved to {output_file}")