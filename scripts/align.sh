#!/bin/bash

# Example Usage:
# ./align.sh ../data/omicron-aggregated-sequences.fasta
# ./align.sh ../data/omicron-aggregated-sequences.fasta ../data/aligned-sequences/custom-output.txt

# Check if the input file name is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file> [output_file]"
  exit 1
fi

# Variables
INPUT_FILE="$1"
OUTPUT_FILE="${2:-../data/aligned-sequences/$(basename "$INPUT_FILE" .fasta)-aligned.txt}"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Input file $INPUT_FILE does not exist."
  exit 1
fi

# Create the output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Run Clustal Omega
clustalo -i "$INPUT_FILE" -o "$OUTPUT_FILE" --force

# Check if the alignment was successful
if [ $? -eq 0 ]; then
  echo "Alignment successful. Output saved to $OUTPUT_FILE"
else
  echo "Alignment failed."
  exit 1
fi