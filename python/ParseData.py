# ParseData.py

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import os
from datetime import datetime

def read_wuhan_sequence():
    print("Reading Wuhan reference sequence...")
    wuhan_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'fasta-sequences', 'Wuhan')
    for filename in os.listdir(wuhan_dir):
        if filename.endswith('.fasta'):
            wuhan_path = os.path.join(wuhan_dir, filename)
            print(f"Found Wuhan sequence file: {filename}")
            wuhan_record = SeqIO.read(wuhan_path, 'fasta')
            return wuhan_record
    raise FileNotFoundError("Wuhan sequence not found in directory.")

def read_variant_sequences(variant):
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
            if date_str:
                collection_date = parse_date(date_str)
                if collection_date:
                    variant_sequences.append({
                        'record': record,
                        'date': collection_date
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
    if "| Date: " in header:
        return header.split("| Date: ")[1].strip()
    return None

def parse_date(date_str):
    for fmt in ('%Y-%m-%d', '%Y-%m', '%Y'):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None

def align_sequences(wuhan_seq, variant_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
    alignments = aligner.align(wuhan_seq, variant_seq)
    
    # Get the best alignment
    aligned_wuhan_seq, aligned_variant_seq = alignments[0].aligned[0], alignments[0].aligned[1]
    aligned_wuhan = []
    aligned_variant = []
    wuhan_index = 0
    variant_index = 0
    
    for (start_wuhan, end_wuhan), (start_variant, end_variant) in zip(aligned_wuhan_seq, aligned_variant_seq):
        # Add gaps where needed
        while wuhan_index < start_wuhan:
            aligned_wuhan.append(wuhan_seq[wuhan_index])
            aligned_variant.append('-')
            wuhan_index += 1
        while variant_index < start_variant:
            aligned_wuhan.append('-')
            aligned_variant.append(variant_seq[variant_index])
            variant_index += 1
        # Add aligned segment
        while wuhan_index < end_wuhan:
            aligned_wuhan.append(wuhan_seq[wuhan_index])
            aligned_variant.append(variant_seq[variant_index])
            wuhan_index += 1
            variant_index += 1
    
    aligned_wuhan.extend(wuhan_seq[wuhan_index:])
    aligned_variant.extend(variant_seq[variant_index:])
    
    print("Alignment complete.")
    return ''.join(aligned_wuhan), ''.join(aligned_variant)

def main(variant):
    print(f"Starting analysis for variant: {variant}")
    
    wuhan_record = read_wuhan_sequence()
    wuhan_seq = str(wuhan_record.seq).upper()
    variant_sequences = read_variant_sequences(variant)
    
    print("Beginning alignment of each variant sequence with respect to the Wuhan reference...")
    print("This may take a while depending on the number of sequences.")
    print("-" * 100)
    
    aligned_sequences = []
    for i, entry in enumerate(variant_sequences, start=1):
        variant_record = entry['record']
        collection_date = entry['date']
        print(f"Aligning sequence {i}/{len(variant_sequences)} from {collection_date} - ID: {variant_record.name}")
        
        variant_seq = str(variant_record.seq).upper()
        aligned_wuhan, aligned_variant = align_sequences(wuhan_seq, variant_seq)
        
        aligned_sequences.append({
            'date': collection_date,
            'aligned_wuhan': aligned_wuhan,
            'aligned_variant': aligned_variant,
        })
        
        print(f"Completed alignment for sequence {i}/{len(variant_sequences)} - ID: {variant_record.id}")

    print("All sequences aligned successfully.")
    # The aligned sequences can now be saved or passed to another function for further analysis

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python ParseData.py <VariantName>")
    else:
        variant_name = sys.argv[1]
        main(variant_name)
