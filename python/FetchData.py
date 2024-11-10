from Bio import Entrez
import os

Entrez.email = "stevenpstansberry@gmail.com"

def fetch_variant_sequences(variant, retmax=1):
    # Map variant names to Pango lineage codes
    variant_lineage_mapping = {
        'Alpha': 'B.1.1.7',
        'Beta': 'B.1.351',
        'Gamma': 'P.1',
        'Delta': 'B.1.617.2',
        'Omicron': 'B.1.1.529',
    }

    # Get the lineage code for the variant
    lineage = variant_lineage_mapping.get(variant, variant)

    # Construct the search query to fetch complete genomes
    search_query = f"SARS-CoV-2[Organism] AND ({lineage}[All Fields] OR {variant}[All Fields]) AND complete genome[Title]"

    # Fetch the IDs of matching records
    handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]

    print(f"Fetching {len(ids)} sequences for variant '{variant}'...")

    # Set up output directory one level above the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(script_dir)  # Go up one level
    output_dir = os.path.join(root_dir, 'data', 'fasta-sequences')
    os.makedirs(output_dir, exist_ok=True)

    # Initialize a list to store sequences for this variant
    variant_sequences = []

    for seq_id in ids:
        try:
            # Fetch the complete genome sequence
            seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            sequence_data = seq_handle.read()
            if sequence_data:
                variant_sequences.append(sequence_data)
            else:
                print(f"Warning: Sequence data for ID {seq_id} is empty.")
        except Exception as e:
            print(f"Error fetching sequence ID {seq_id}: {e}")
            continue

    # Save all sequences into a single FASTA file named after the variant
    if variant_sequences:
        variant_filename = f"{variant}.fasta"
        file_path = os.path.join(output_dir, variant_filename)
        with open(file_path, 'w') as f:
            f.write('\n'.join(variant_sequences))
        print(f"Saved {len(variant_sequences)} sequences to {file_path}")
    else:
        print(f"No sequences were fetched for variant '{variant}'.")

def fetch_wuhan_sequence():
    seq_id = 'NC_045512.2'  # Reference sequence ID for the original Wuhan strain
    try:
        seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        sequence_data = seq_handle.read()
        if sequence_data:
            # Save the sequence as 'Wuhan.fasta'
            script_dir = os.path.dirname(os.path.abspath(__file__))
            root_dir = os.path.dirname(script_dir)
            output_dir = os.path.join(root_dir, 'data', 'fasta-sequences')
            os.makedirs(output_dir, exist_ok=True)
            file_path = os.path.join(output_dir, 'Wuhan.fasta')
            with open(file_path, 'w') as f:
                f.write(sequence_data)
            print(f"Saved Wuhan reference sequence to {file_path}")
        else:
            print("Warning: Wuhan sequence data is empty.")
    except Exception as e:
        print(f"Error fetching Wuhan sequence: {e}")


# Fetch sequences for specific variants
print("Fetching sequences for selected SARS-CoV-2 variants...")
print("-" * 50)
print("Variant: Alpha")
fetch_variant_sequences('Alpha', retmax=1)
print("\nvariant: Betta")
fetch_variant_sequences('Beta', retmax=1)
print("\nvariant: Gamma")
fetch_variant_sequences('Gamma', retmax=1)
print("\nvariant: Delta")
fetch_variant_sequences('Delta', retmax=1)
print("\nvariant: Omicron")
fetch_variant_sequences('Omicron', retmax=1)
print ("\n Fetching base Wuhan sequence")
fetch_wuhan_sequence()
print("-" * 50)
print("Done fetching sequences.")





