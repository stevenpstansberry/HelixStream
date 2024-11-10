from Bio import Entrez, SeqIO
import os
import time

Entrez.email = "stevenpstansberry@gmail.com"

def fetch_variant_sequences(variant, retmax=20):
    from datetime import datetime

    # Map variant names to Pango lineage codes
    variant_lineage_mapping = {
        'Alpha': 'B.1.1.7',
        'Beta': 'B.1.351',
        'Gamma': 'P.1',
        'Delta': 'B.1.617.2',
        'Omicron': 'B.1.1.529',
        # Add more variants and lineages as needed
    }

    # Get the lineage code for the variant
    lineage = variant_lineage_mapping.get(variant, variant)

    # Construct the search query to fetch recent complete genomes
    search_query = (
        f"SARS-CoV-2[Organism] AND ({lineage}[All Fields] OR {variant}[All Fields]) "
        "AND complete genome[Title] AND 10000:50000[Sequence Length]"
    )

    # Fetch the IDs of matching records
    handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=retmax, sort='most recent')
    record = Entrez.read(handle)
    ids = record["IdList"]

    print(f"Fetching {len(ids)} sequences for variant '{variant}'...")

    # Set up output directory one level above the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(script_dir)  # Go up one level

    # Create a directory for the variant
    variant_dir = os.path.join(root_dir, 'data', 'fasta-sequences', variant)
    os.makedirs(variant_dir, exist_ok=True)

    for seq_id in ids:
        try:
            # Fetch the sequence in GenBank format
            seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            seq_record = SeqIO.read(seq_handle, "genbank")
            seq_handle.close()

            # Initialize collection_date as None
            collection_date = None

            # Extract collection date from features
            for feature in seq_record.features:
                if feature.type == 'source':
                    # Check for 'collection_date' key
                    if 'collection_date' in feature.qualifiers:
                        collection_date = feature.qualifiers['collection_date'][0]
                        break  # Exit loop after finding the date

            if not collection_date:
                # Fallback to publication date from Entrez summary
                summary_handle = Entrez.esummary(db="nucleotide", id=seq_id)
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()
                collection_date = summary_record[0].get('PubDate', 'Unknown')

            if collection_date and collection_date != 'Unknown':
                # Prepare the FASTA-formatted sequence with updated header
                new_header = f">{seq_record.description} | Date: {collection_date}"
                sequence_data_with_date = f"{new_header}\n{str(seq_record.seq)}\n"

                # Create a filename using the accession number and collection date
                sanitized_date = collection_date.replace(' ', '_').replace('/', '-').replace(':', '-')
                filename = f"{seq_id}_{sanitized_date}.fasta"
                file_path = os.path.join(variant_dir, filename)

                # Save the sequence to a file
                with open(file_path, 'w') as f:
                    f.write(sequence_data_with_date)
                print(f"Saved sequence ID {seq_id} to {file_path}")
            else:
                print(f"Skipped sequence ID {seq_id} due to missing date.")

        except Exception as e:
            print(f"Error fetching sequence ID {seq_id}: {e}")
            continue
        finally:
            time.sleep(0.5)  # Be nice to the NCBI servers

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
fetch_variant_sequences('Alpha', retmax=20)
print("\nvariant: Beta")
fetch_variant_sequences('Beta', retmax=20)
print("\nvariant: Gamma")
fetch_variant_sequences('Gamma', retmax=20)
print("\nvariant: Delta")
fetch_variant_sequences('Delta', retmax=20)
print("\nvariant: Omicron")
fetch_variant_sequences('Omicron', retmax=20)
print ("\n Fetching base Wuhan sequence")
fetch_wuhan_sequence()
print("-" * 50)
print("Done fetching sequences.")





