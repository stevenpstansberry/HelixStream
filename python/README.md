# README.md

## Overview

This code sample is part a bioinformatics pipeline designed to fetch, process, align, and analyze SARS-CoV-2 variant sequences. It includes scripts for fetching sequences from the NCBI database, parsing and aggregating them, aligning sequences, and calculating evolutionary metrics. This is designed to be utulized by a backend API, but these scripts are able to be invoked manually. I've included an aggregated nucleotide sequence for the omicron variant since it can be somewhat tricky to use the Clustal Omega tool

## Prerequisites

- Python 3.x
- Biopython
- pandas
- numpy
- Clustal Omega (for sequence alignment)

## Setup

1. Install the required Python packages:
   ```sh
   cd python
   pip install -r requirements.txt
   ```
2. Ensure Clustal Omega is installed and available in your system's PATH.
   You can get this from [Clustal Omega](http://www.clustal.org/omega/)

   If you don't want to download Clustal Omega (it's a little tricky to use),
   you can align the sequences [here](https://www.ebi.ac.uk/jdispatcher/msa/clustalo) or use the provided aligned omicron sequences in the data directory

   Aligned sequences are required for analysis because they allow us to compare the sequences accurately by lining up similar regions. This helps in identifying differences and similarities between the sequences.

## Scripts

### 1. FetchData.py

This script fetches SARS-CoV-2 variant sequences from the NCBI Nucleotide database using the Entrez functionality from the BioPython library.

#### Usage

```sh
python FetchData.py [variantName]
```

- If `variantName` is provided, it fetches sequences for the specified variant.
- Otherwise, it fetches sequences for all predefined variants and the Wuhan reference sequen

#### Example

```py
python FetchData.py Omicron
```

### 2. ParseData.py

This script reads and processes SARS-CoV-2 genome sequences, including the Wuhan reference sequence and variant sequences.

#### Usage

```sh
python FetchData.py [variantName]
```

- If `variantName` is provided, it processes the specified variant.
- Otherwise, it processes all predefined variants.

#### Example

```sh
python ParseData.py Omicron
```

### 3. align.sh

This script aligns SARS-CoV-2 sequences using Clustal Omega.

#### Usage

```sh
./align.sh <input_file> [output_file]
```

- `<input_file>`: Path to the input FASTA file.
- `[output_file]` (optional): Path to the output file. If not provided, a default output file is created.

#### Example

```sh
./align.sh ../data/aggregated-sequences/omicron-aggregated-sequences.fasta
```

#### AnalyzeData.py

This script calculates evolutionary metrics for viral genome sequences, including percentage identity, substitutions, insertions, deletions, total mutations, mutation density per kilobase, and applies the Jukes-Cantor model to estimate evolutionary distances.
Note that this requires an aligned sequence

#### Usage

```sh
python AnalyzeData.py [variantName]
```

- If `variantName` is provided, it analyzes the specified variant.
- Otherwise, it analyzes all variants.

#### Example

```sh
python AnalyzeData.py Omicron
```
