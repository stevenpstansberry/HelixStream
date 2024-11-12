from Bio import AlignIO
import pandas as pd
from datetime import datetime
import os

def calculate_evolutionary_metrics():
    # Define the path to the alignment file
    alignment_file = os.path.join("..", "data", "aligned-sequences", "aligned.txt")
    
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

        # Extract the date from the sequence ID, assuming the format is ID###YYYY-MM-DD
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

# Usage
df = calculate_evolutionary_metrics()

# Display the data
print(df.head())
