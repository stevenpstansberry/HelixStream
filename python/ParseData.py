"""
This script reads and processes SARS-CoV-2 genome sequences, including the Wuhan reference sequence and variant sequences.

It performs the following tasks:
- Reads the Wuhan reference sequence from the specified directory.
- Reads variant sequences for specified variants.
- Extracts metadata such as collection dates and custom identifiers from sequence headers.
- Aggregates sequences and saves them to a FASTA file for further analysis.

Dependencies:
    - Biopython

Usage:
    python ParseData.py [variantName]

    If 'variantName' is provided, it processes the specified variant.
    Otherwise, it processes all predefined variants.

Note:
    Ensure that the required directories and files are correctly set up before running the script.
"""
import os
from datetime import datetime
from Bio import SeqIO

def read_wuhan_sequence():
    """
    Read the Wuhan reference sequence from the designated directory.

    The function searches for a FASTA file containing the Wuhan reference sequence,
    reads it, and extracts a custom identifier and description for consistency.

    Returns:
        Bio.SeqRecord.SeqRecord: The Wuhan reference sequence record.

    Raises:
        FileNotFoundError: If the Wuhan sequence file is not found in the directory.
    """    
    print("Reading Wuhan reference sequence...")
    wuhan_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'fasta-sequences', 'Wuhan')
    for filename in os.listdir(wuhan_dir):
        if filename.endswith('.fasta'):
            wuhan_path = os.path.join(wuhan_dir, filename)
            print(f"Found Wuhan sequence file: {filename}")
            wuhan_record = SeqIO.read(wuhan_path, 'fasta')
            # Extract custom_id and assign to record.id
            wuhan_record.id = "Severe/acute/respiratory/syndrome/coronavirus/2/isolate/Wuhan/Hu/1,/complete/genome"
            wuhan_record.description = "Severe/acute/respiratory/syndrome/coronavirus/2/isolate/Wuhan/Hu/1,/complete/genome| Date: 2019-12-29"
            print(f"Wuhan sequence ID: {wuhan_record.id}")
            return wuhan_record
    raise FileNotFoundError("Wuhan sequence not found in directory.")

def read_variant_sequences(variant):
    """
    Read and process sequences for a specified SARS-CoV-2 variant.

    The function reads all FASTA files for the given variant, extracts the collection date
    and custom identifier from the sequence headers, and stores the sequences along with their dates.

    Parameters:
        variant (str): The name of the variant to process (e.g., 'Alpha').

    Returns:
        list: A list of dictionaries containing sequence records and their collection dates.
    """    
    print(f"Reading sequences for variant: {variant}")
    variant_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'fasta-sequences', variant)
    fasta_files = [filename for filename in os.listdir(variant_dir) if filename.endswith('.fasta')]
    total_sequences = len(fasta_files)
    print(f"{total_sequences} sequences found for {variant}.")

    variant_sequences = []

    for idx, filename in enumerate(fasta_files, start=1):
        print(f"Reading sequence {idx}/{total_sequences}: {filename}")
        file_path = os.path.join(variant_dir, filename)
        try:
            record = SeqIO.read(file_path, 'fasta')
            date_str = extract_date_from_header(record.description)
            custom_id = extract_custom_id(record.description)
            record.id = custom_id  # Assign custom_id to record.id
            if date_str:
                collection_date = parse_date(date_str)
                if collection_date:
                    variant_sequences.append({
                        'record': record,
                        'date': collection_date,
                    })
        except ValueError:
            print(f"Warning: No valid sequence found in file {filename}. Skipping.")
        except Exception as e:
            print(f"Error reading file {filename}: {e}")

    # Order sequences chronologically
    variant_sequences.sort(key=lambda x: x['date'])
    print(f"Finished reading {len(variant_sequences)} sequences for {variant}.\n")
    return variant_sequences

def extract_date_from_header(header):
    """
    Extract the collection date from a sequence header.

    Parameters:
        header (str): The header string from which to extract the date.

    Returns:
        str or None: The extracted date string, or None if not found.
    """    
    if "| Date: " in header:
        return header.split("| Date: ")[1].strip()
    return None

def parse_date(date_str):
    """
    Parse a date string into a datetime object.

    Tries multiple date formats to handle various possible date representations.

    Parameters:
        date_str (str): The date string to parse.

    Returns:
        datetime.datetime or None: The parsed datetime object, or None if parsing fails.
    """    
    for fmt in ('%Y-%m-%d', '%Y-%m', '%Y'):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None

def extract_custom_id(header):
    """
    Extract a custom identifier from a sequence header.

    The function looks for the pattern "SARS-CoV-2/" and extracts the subsequent text up to the first comma.

    Parameters:
        header (str): The header string from which to extract the custom ID.

    Returns:
        str: The extracted custom ID, or the original header if the pattern is not found.
    """    
    # Extracts the portion after "SARS-CoV-2/" up to the first comma
    if "SARS-CoV-2/" in header:
        start = header.find("SARS-CoV-2/") + len("SARS-CoV-2/")
        end = header.find(",", start)
        return header[start:end].strip()
    return header  # Fallback to entire header if pattern not found

def save_aggregated_sequences(variant, records):
    """
    Save aggregated sequences to a FASTA file.

    The function writes the given sequence records to a FASTA file, formatting the headers to include
    the custom IDs and collection dates.

    Parameters:
        variant (str): The name of the variant being processed.
        records (list): A list of SeqRecord objects to save.
    """    
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'aggregated-sequences')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{variant}-aggregated-sequences.fasta")

    with open(output_file, 'w') as f:
        for record in records:
            date_str = record.description.split("| Date: ")[1].strip() if "| Date: " in record.description else "Unknown"
            record_id = record.id.replace(" ", "/")
            f.write(f">{record_id}###{date_str}\n")
            f.write(f"{record.seq}\n")

    print(f"Aggregated sequences saved to {output_file}")

def main(variant=None):
    """
    Main function to process sequences.

    If a variant name is provided, it processes that variant. Otherwise, it processes all predefined variants.

    Parameters:
        variant (str, optional): The name of the variant to process.

    Returns:
        list: A list of all SeqRecord objects processed.
    """    
    if variant:
        print(f"Starting analysis for variant: {variant}")
        wuhan_record = read_wuhan_sequence()
        variant_sequences = read_variant_sequences(variant)
        all_records = [wuhan_record] + [entry['record'] for entry in variant_sequences]
        save_aggregated_sequences(variant, all_records)
    else:
        print("Starting analysis for all variants...")
        wuhan_record = read_wuhan_sequence()
        variants = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Omicron']
        for variant in variants:
            variant_sequences = read_variant_sequences(variant)
            all_records = [wuhan_record] + [entry['record'] for entry in variant_sequences]
            save_aggregated_sequences(variant, all_records)
    return all_records

if __name__ == "__main__":
    """
    Main execution block.

    Parses command-line arguments to optionally accept a variant name.
    Calls the main function with the specified variant or processes all variants if none is provided.
    """    
    import sys
    if len(sys.argv) > 1:
        variant_name = sys.argv[1]
        sequences = main(variant_name)
    else:
        sequences = main()