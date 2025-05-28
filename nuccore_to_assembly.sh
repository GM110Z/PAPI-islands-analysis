#!/bin/bash

# Usage: ./nc_to_gcf.sh nc_list.txt output_file.tsv

input_file="$1"
output_file="$2"

if [ -z "$input_file" ] || [ -z "$output_file" ]; then
  echo "Usage: $0 <nc_list_file> <output_file>"
  exit 1
fi

# Empty the output file before starting
> "$output_file"

while IFS= read -r id; do
  printf "%s\t" "$id" >> "$output_file"
  esearch -db nucleotide -query "$id" < /dev/null | \
  elink -target assembly | \
  esummary | \
  xtract -pattern DocumentSummary -element AssemblyAccession >> "$output_file"
done < "$input_file"




