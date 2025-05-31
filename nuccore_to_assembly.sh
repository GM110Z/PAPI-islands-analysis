#!/bin/bash

# Convert input file to Unix format
awk '{ sub(/\r$/, ""); print }' "$1" > "$2"

input_file="$2"
output_file="$3"

if [ -z "$input_file" ] || [ -z "$output_file" ]; then
  echo "Usage: $0 <nc_list_file> <output_file>"
  exit 1
fi

# Empty the output file before starting
> "$output_file"

# Read all columns from the input file
while IFS=$'\t' read -r nc_id number; do
  # Get the first AssemblyAccession only
  assembly=$(esearch -db nucleotide -query "$nc_id" < /dev/null | \
    elink -target assembly | \
    esummary | \
    xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)
  
  printf "%s\t%s\t%s\n" "$nc_id" "$number" "$assembly" >> "$output_file"
done < "$input_file"


