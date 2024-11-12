import os
from datetime import datetime
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

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
    print("Initializing aligner...")
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
    
    # Set gap penalties
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Set end gap penalties to zero
    aligner.end_open_gap_score = 0
    aligner.end_extend_gap_score = 0
    
    print("Performing alignment...")
    alignments = aligner.align(wuhan_seq, variant_seq)
    
    # Get the best alignment
    alignment = alignments[0]
    
    # Extract the aligned sequences with gaps
    aligned_wuhan = alignment.target
    aligned_variant = alignment.query
    
    print("Alignment complete.")
    print(f"Aligned variant sequence is of length: {len(aligned_variant)}")
    print(f"Aligned Wuhan sequence is of length: {len(aligned_wuhan)}")
    return str(aligned_wuhan), str(aligned_variant)

def extract_custom_id(header):
    # Extracts the portion after "SARS-CoV-2/" up to the first comma
    if "SARS-CoV-2/" in header:
        start = header.find("SARS-CoV-2/") + len("SARS-CoV-2/")
        end = header.find(",", start)
        return header[start:end].strip()
    return header  # Fallback to entire header if pattern not found

def save_aligned_sequences(variant, aligned_sequences, aligned_wuhan):
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'aligned-sequences')
    os.makedirs(output_dir, exist_ok=True)

    # Save all aligned variant sequences to a single file
    variant_file_path = os.path.join(output_dir, f"{variant}_aligned_sequences.fasta")
    with open(variant_file_path, 'w') as f:
        for entry in aligned_sequences:
            f.write(f">Aligned Variant ID: {entry['custom_id']} | Date: {entry['date']}\n")
            f.write(entry['aligned_variant'] + "\n")
            f.write("#" * 80 + "\n")  # Delimiter between sequences

    print(f"Aligned variant sequences saved to {variant_file_path}")

    # Save the aligned Wuhan sequence separately
    wuhan_file_path = os.path.join(output_dir, 'Wuhan_aligned_sequence.fasta')
    with open(wuhan_file_path, 'w') as f:
        f.write(">Aligned Wuhan Reference Sequence\n")
        f.write(aligned_wuhan + "\n")

    print(f"Aligned Wuhan sequence saved to {wuhan_file_path}")

def main(variant):
    print(f"Starting analysis for variant: {variant}")
    
    wuhan_record = read_wuhan_sequence()
    wuhan_seq = str(wuhan_record.seq).upper()
    variant_sequences = read_variant_sequences(variant)
    
    print("Beginning alignment of each variant sequence with respect to the Wuhan reference...")
    print("This may take a while depending on the number of sequences.")
    print("-" * 100)
    
    aligned_sequences = []
    aligned_wuhan = None  # To store the aligned Wuhan reference sequence once
    for i, entry in enumerate(variant_sequences, start=1):
        variant_record = entry['record']
        collection_date = entry['date']
        custom_id = extract_custom_id(variant_record.description)
        print(f"Aligning sequence {i}/{len(variant_sequences)} from {collection_date} - ID: {custom_id}")
        
        variant_seq = str(variant_record.seq).upper()
        aligned_wuhan, aligned_variant = align_sequences(wuhan_seq, variant_seq)
        
        aligned_sequences.append({
            'date': collection_date,
            'aligned_variant': aligned_variant,
            'custom_id': custom_id,
        })
        
        print(f"Completed alignment for sequence {i}/{len(variant_sequences)} - ID: {custom_id}")

    # Save aligned sequences
    save_aligned_sequences(variant, aligned_sequences, aligned_wuhan)

    print("All sequences aligned and saved successfully.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python ParseData.py <VariantName>")
    else:
        variant_name = sys.argv[1]
        main(variant_name)
