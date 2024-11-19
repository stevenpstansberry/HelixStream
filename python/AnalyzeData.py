"""
This script calculates evolutionary metrics for viral genome sequences, including percentage identity,
substitutions, insertions, deletions, total mutations, mutation density per kilobase, and applies the
Jukes-Cantor model to estimate evolutionary distances. It can analyze a specific variant or all variants
in the aligned sequences directory.

Dependencies:
    - Biopython
    - pandas
    - numpy

Usage:
    python AnalyzeData.py [variantName]

    If 'variantName' is provided, it analyzes the specified variant. Otherwise, it analyzes all variants.
"""

from Bio import AlignIO
import pandas as pd
import numpy as np
from datetime import datetime
import os
import argparse

NUM_SITES = 29903  # Total number of nucleotides in the original Wuhan COVID 2019 reference genome

def calculate_evolutionary_metrics(variantName):
    """
    Calculate evolutionary metrics for a given variant.

    Reads the alignment file for the specified variant and computes metrics such as percent identity,
    substitutions, insertions, deletions, total mutations, and mutation density per kilobase for each sequence.

    Parameters:
        variantName (str): The name of the variant to analyze (e.g., 'Omicron').

    Returns:
        pandas.DataFrame: A DataFrame containing the calculated metrics for each sequence.
    """

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
    """
    Apply the Jukes-Cantor model to calculate evolutionary distances.

    Adds a new column 'Jukes_Cantor_Distance' to the DataFrame based on the number of substitutions.

    Parameters:
        df (pandas.DataFrame): The DataFrame containing evolutionary metrics.
        total_sites (int): The total number of nucleotide sites considered in the analysis.

    Returns:
        pandas.DataFrame: The DataFrame with the Jukes-Cantor distances added.
    """    
    # Calculate the Jukes-Cantor distance for each sequence
    df['Jukes_Cantor_Distance'] = df['Substitutions'].apply(lambda x: jukes_cantor(x, total_sites))
    return df

def jukes_cantor(substitutions, total_sites):
    """
    Calculate the Jukes-Cantor distance for a given number of substitutions.

    Parameters:
        substitutions (int): The number of substitutions between sequences.
        total_sites (int): The total number of nucleotide sites.

    Returns:
        float: The Jukes-Cantor distance. Returns infinity if the proportion of substitutions exceeds 0.75.
    """
        
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

def analyze_all_variants():
    """
    Analyze all variants in the aligned sequences directory.

    Iterates over each alignment file in the directory, calculates evolutionary metrics,
    applies the Jukes-Cantor model, and saves the results to JSON files.
    """
    aligned_dir = os.path.join("..", "data", "aligned-sequences")
    for filename in os.listdir(aligned_dir):
        if filename.endswith("-aligned.txt"):
            variant_name = filename.replace("-aligned.txt", "")
            print(f"Analyzing variant: {variant_name}")
            df = calculate_evolutionary_metrics(variant_name)
            total_sites = NUM_SITES 
            df = calculate_jukes_cantor_distance(df, total_sites)

            # Save the DataFrame in JSON format
            output_dir = os.path.join("..", "data", "analyzed-sequences", variant_name)
            os.makedirs(output_dir, exist_ok=True)
            today_date = datetime.now().strftime("%Y-%m-%d")
            output_file = os.path.join(output_dir, f"analyzed-{variant_name}-{today_date}.json")
            df.to_json(output_file, orient='records', date_format='iso')

            print(f"DataFrame for {variant_name} saved to {output_file}")

if __name__ == "__main__":
    """
    Main execution block.

    This block sets up command-line argument parsing to accept an optional variant name.
    If a variant name is provided, it calculates evolutionary metrics for that variant,
    applies the Jukes-Cantor model, displays the results, and saves them to a JSON file.
    If no variant name is provided, it analyzes all available variants in the aligned sequences directory.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate evolutionary metrics for a given variant.")
    parser.add_argument("variantName", type=str, nargs='?', help="Name of the variant (e.g., Omicron)")

    # Parse arguments
    args = parser.parse_args()

    if args.variantName:
        # Analyze the specified variant
        df = calculate_evolutionary_metrics(args.variantName)
        total_sites = NUM_SITES 
        df = calculate_jukes_cantor_distance(df, total_sites)

        # Display the DataFrame with the new Jukes-Cantor distance column
        pd.set_option('display.max_columns', None)

        print("DataFrame with Jukes-Cantor distance:")
        print("-" * 100)
        print(df)

        # Save the DataFrame in JSON format
        output_dir = os.path.join("..", "data", "analyzed-sequences", args.variantName)
        os.makedirs(output_dir, exist_ok=True)
        today_date = datetime.now().strftime("%Y-%m-%d")
        output_file = os.path.join(output_dir, f"analyzed-{args.variantName}-{today_date}.json")
        df.to_json(output_file, orient='records', date_format='iso')

        print(f"DataFrame saved to {output_file}")
    else:
        # Analyze all variants
        analyze_all_variants()