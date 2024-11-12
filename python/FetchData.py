from Bio import Entrez, SeqIO
import os
import time
from datetime import datetime


Entrez.email = "stevenpstansberry@gmail.com"

#Todo addd logic to record when last fetched

def fetch_variant_sequences(variant, retmax=50):
    # Map variant names to Pango lineage codes, can easily add more
    variant_lineage_mapping = {
        'Alpha': 'B.1.1.7',
        'Beta': 'B.1.351',
        'Gamma': 'P.1',
        'Delta': 'B.1.617.2',
        'Omicron': 'B.1.1.529',
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

    # Create a directory for the variants if it doesn't exist
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

            # Handle no collection date found in features
            if not collection_date:
                # Fallback to publication date from Entrez summary
                summary_handle = Entrez.esummary(db="nucleotide", id=seq_id)
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()
                collection_date = summary_record[0].get('PubDate', 'Unknown')

            # Function to parse various date formats
            def parse_date(date_str):
                for fmt in ("%d-%m-%Y", "%Y-%m-%d", "%d-%b-%Y", "%Y-%m"):  # Try multiple formats
                    try:
                        return datetime.strptime(date_str, fmt).strftime("%Y-%m-%d")
                    except ValueError:
                        continue
                return date_str  # If no format matches, return original

            # Format collection_date to "year-month-day"
            formatted_date = parse_date(collection_date)

            if formatted_date and formatted_date != 'Unknown':
                # Prepare the FASTA-formatted sequence with updated header
                new_header = f">{seq_record.description} | Date: {formatted_date}"
                sequence_data_with_date = f"{new_header}\n{str(seq_record.seq)}\n"

                # Create a filename using the accession number and formatted collection date
                sanitized_date = formatted_date.replace(' ', '_').replace('/', '-').replace(':', '-')
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
        # Fetch the sequence in GenBank format
        seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(seq_handle, "genbank")
        seq_handle.close()

        # Initialize collection_date as None
        collection_date = None

        # Extract collection date from features
        for feature in seq_record.features:
            if feature.type == 'source':
                if 'collection_date' in feature.qualifiers:
                    collection_date = feature.qualifiers['collection_date'][0]
                    print(f"Found collection_date: {collection_date}")
                    break
                else:
                    print(f"No 'collection_date' in qualifiers for ID {seq_id}")

        if not collection_date:
            # Fallback to publication date from Entrez summary
            summary_handle = Entrez.esummary(db="nucleotide", id=seq_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            collection_date = summary_record[0].get('PubDate', 'Unknown')
            print(f"Using PubDate as collection_date: {collection_date}")

        if collection_date and collection_date != 'Unknown':
            # Prepare the FASTA-formatted sequence with updated header
            new_header = f">{seq_record.description} | Date: {collection_date}"
            sequence_data_with_date = f"{new_header}\n{str(seq_record.seq)}\n"

            # Create a file for the Wuhan reference sequence
            filename = f"original_wuhan_strain.fasta"

            # Save the sequence to the 'Wuhan' directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            root_dir = os.path.dirname(script_dir)
            variant_dir = os.path.join(root_dir, 'data', 'fasta-sequences', 'Wuhan')
            os.makedirs(variant_dir, exist_ok=True)
            file_path = os.path.join(variant_dir, filename)

            with open(file_path, 'w') as f:
                f.write(sequence_data_with_date)
            print(f"Saved Wuhan reference sequence to {file_path}")
        else:
            print("Warning: Wuhan sequence data is missing date information.")
    except Exception as e:
        print(f"Error fetching Wuhan sequence: {e}")


# Fetch sequences for specific variants
print("Fetching sequences for selected SARS-CoV-2 variants...")
print("-" * 50)
# variants = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron']
variants = ['Omicron']


for variant in variants:
    print(f"\nVariant: {variant}")
    fetch_variant_sequences(variant, retmax=5)

print("\nFetching base Wuhan sequence")
fetch_wuhan_sequence()
print("-" * 50)
print("Done fetching sequences.")





