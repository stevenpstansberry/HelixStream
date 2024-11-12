# AlignData.py
from ParseData import main
from Bio import Align

#TODO use unix scripts to align sequences in ec2

def align_sequences():
    variant_name = 'Omicron'  
    sequences = main(variant_name)
    
    for idx, record in enumerate(sequences, start=1): # start at 1 to skip over the Wuhan sequence
        print(f"Sequence {idx}: {record.id}")
        #print(f"Length: {len(record.seq)}")
        #print(f"Sequence:\n{record.seq}\n")

    aligner = Align.PairwiseAligner()
    wuhanSequence = sequences[0].seq
    seq1 = sequences[1].seq




if __name__ == "__main__":
    align_sequences()