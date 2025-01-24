#!/bin/bash

# Example Usage:
# ./align.sh ../data/aggregated-sequences/omicron-aggregated-sequences.fasta
# ./align.sh ../data/aggregated-sequences/omicron-aggregated-sequences.fasta ../data/aligned-sequences/custom-output.txt

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
#!/bin/bash

# Example Usage:
# ./align.sh ../data/aggregated-sequences/omicron-aggregated-sequences.fasta
# ./align.sh ../data/aggregated-sequences/omicron-aggregated-sequences.fasta ../data/aligned-sequences/custom-output.txt

# Colors for messages
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Log file
LOG_FILE="../logs/alignment_$(date +%Y%m%d%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Trap interrupt signals for cleanup
trap 'echo -e "${RED}Alignment interrupted. Cleaning up...${NC}"; rm -f "$OUTPUT_FILE"; exit 1' INT

# Help message
if [ "$1" == "--help" ]; then
  echo -e "${BLUE}Usage: $0 <input_file> [output_file] [additional_clustalo_params]${NC}"
  echo "Align sequences in a FASTA file using Clustal Omega."
  echo
  echo "Arguments:"
  echo "  input_file               Path to the input FASTA file."
  echo "  output_file              Path to save the aligned output (optional)."
  echo "  additional_clustalo_params Parameters to pass to Clustal Omega (optional)."
  exit 0
fi

# Check if the input file name is provided
if [ -z "$1" ]; then
  echo -e "${RED}Usage: $0 <input_file> [output_file]${NC}"
  exit 1
fi

# Variables
INPUT_FILE="$1"
OUTPUT_FILE="${2:-../data/aligned-sequences/$(basename "$INPUT_FILE" .fasta)-aligned.txt}"
ADDITIONAL_PARAMS="$3"

# Check if Clustal Omega is installed
if ! command -v clustalo &> /dev/null; then
  echo -e "${RED}Error: Clustal Omega is not installed or not in the PATH.${NC}"
  exit 1
fi

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo -e "${RED}Input file $INPUT_FILE does not exist.${NC}"
  exit 1
fi

# Check if the input file has the correct extension
if [[ "$INPUT_FILE" != *.fasta ]]; then
  echo -e "${RED}Error: Input file must have a .fasta extension.${NC}"
  exit 1
fi

# Create the output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Start time for runtime calculation
START_TIME=$(date +%s)

# Inform the user that the alignment process is starting
echo -e "${BLUE}Starting alignment process for $INPUT_FILE...${NC}"
echo -e "This may take a while. Please be patient."

# Run Clustal Omega
clustalo -i "$INPUT_FILE" -o "$OUTPUT_FILE" --force $ADDITIONAL_PARAMS &

# Get the PID of the Clustal Omega process
CLUSTALO_PID=$!

# Periodically update the user about the alignment process
while kill -0 $CLUSTALO_PID 2> /dev/null; do
  echo -e "$(date): ${BLUE}Alignment is still in progress...${NC}"
  sleep 10  # Update every 10 seconds
done

# Check if the alignment was successful
if wait $CLUSTALO_PID; then
  END_TIME=$(date +%s)
  ELAPSED=$((END_TIME - START_TIME))
  echo -e "${GREEN}Alignment successful. Output saved to $OUTPUT_FILE.${NC}"
  echo -e "${GREEN}Total runtime: $ELAPSED seconds.${NC}"
else
  echo -e "${RED}Alignment failed.${NC}"
  exit 1
fi

# Create the output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Inform the user that the alignment process is starting
echo "Starting alignment process for $INPUT_FILE..."
echo "This may take a while. Please be patient."

# Run Clustal Omega
clustalo -i "$INPUT_FILE" -o "$OUTPUT_FILE" --force &

# Get the PID of the Clustal Omega process
CLUSTALO_PID=$!

# Periodically update the user about the alignment process
while kill -0 $CLUSTALO_PID 2> /dev/null; do
  echo "Alignment is still in progress..."
  sleep 60  # Update every 60 seconds
done

# Check if the alignment was successful
if wait $CLUSTALO_PID; then
  echo "Alignment successful. Output saved to $OUTPUT_FILE"
else
  echo "Alignment failed."
  exit 1
fi