from Bio import SeqIO
from io import StringIO

def extract_sequences(raw_data):
    sequences = []
    for entry in raw_data:
        with StringIO(entry) as handle:
            record = SeqIO.read(handle, "genbank")
            sequence_data = str(record.seq)
            collection_date = record.annotations.get("date", "Unknown")
            sequences.append((sequence_data, collection_date))
    return sequences