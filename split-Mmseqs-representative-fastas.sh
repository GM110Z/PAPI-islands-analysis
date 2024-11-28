#!/bin/bash

# Path to the multi-FASTA file
multi_fasta_file=$1/representatives.fasta_rep_seq.fasta  # Replace with your file path
output_folder=$1

# Split the multi-FASTA file
gawk -v output_folder="$output_folder" '
BEGIN {seq = ""; header = ""}
{
    if ($0 ~ /^>/) {
        if (seq != "") {
            # Extract the desired part of the header for filename
            split(header, parts, " ")
            file_name = output_folder "/" parts[1] ".fasta"
            print header "\n" seq > file_name
        }
        # Start a new sequence
        header = substr($0, 2)  # Remove ">"
        seq = ""
    } else {
        seq = seq $0
    }
}
END {
    # Write the last sequence
    if (seq != "") {
        split(header, parts, " ")
        file_name = output_folder "/" parts[1] ".fasta"
        print header "\n" seq > file_name
    }
}
' "$multi_fasta_file"

echo "Sequences split into separate files in folder: $output_folder"
