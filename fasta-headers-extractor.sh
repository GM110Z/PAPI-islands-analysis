#!/bin/bash

# Input and output files
input_fasta=$1  # Replace with your input file
output_file="extracted_info.tsv"          # Replace with your desired output file

# Process the FASTA headers
grep '^>' "$input_fasta" | \
awk -F'[>: -]' '{print $2, $3, $4}' > "$output_file"

echo "Extraction complete. Results saved to $output_file"
