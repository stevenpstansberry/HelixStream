# AnalyzeData.py
from ParseData import main

def analyze_sequences():
    variant_name = 'Omicron'  
    sequences = main(variant_name)
    
    for idx, record in enumerate(sequences, start=1):
        print(f"Sequence {idx}: {record.id}")
        print(f"Length: {len(record.seq)}")
        print(f"Sequence:\n{record.seq}\n")

if __name__ == "__main__":
    analyze_sequences()