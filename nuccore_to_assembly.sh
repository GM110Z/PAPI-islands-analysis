#!/bin/bash

# Usage: ./nuccore2assembly.sh input_nc_list.txt cleaned_input.txt temp_output.tsv final_input.tsv final_output.txt

# Step 0: Convert input file to Unix format (strip CR)
awk '{ sub(/\r$/, ""); print }' "$1" > "$2"

input_file="$2"
output_file="$3"

if [ -z "$input_file" ] || [ -z "$output_file" ]; then
  echo "Usage: $0 <nc_list_file> <output_file>"
  exit 1
fi

# Set email environment for NCBI EDirect
export EMAIL=your.email@university.ac.uk  # <- set your real academic email

# Empty the output file before starting
> "$output_file"

# Function to get assembly accession with retries
get_assembly() {
  local nc_id="$1"
  local attempt=0
  local max_attempts=5
  local delay=5
  local result=""

  while [ "$attempt" -lt "$max_attempts" ]; do
    result=$(esearch -db nucleotide -query "$nc_id" < /dev/null | \
      elink -target assembly | \
      esummary | \
      xtract -pattern DocumentSummary -element AssemblyAccession | head -n 1)

    if [ -n "$result" ]; then
      echo "$result"
      return 0
    fi

    echo "Attempt $((attempt+1)) failed for $nc_id. Retrying in $delay sec..." >&2
    sleep "$delay"
    attempt=$((attempt+1))
    delay=$((delay+5))
  done

  echo "NA"
  return 1
}

# Read and process each line of the input file
while IFS=$'\t' read -r nc_id number; do
  # Clean carriage return
  nc_id=$(echo "$nc_id" | tr -d '\r')

  if [ -z "$nc_id" ]; then
    echo "Warning: empty NCID encountered. Skipping..." >&2
    continue
  fi

  # Get assembly with retries
  assembly=$(get_assembly "$nc_id")

  # Write result to output
  printf "%s\t%s\t%s\n" "$nc_id" "$number" "$assembly" >> "$output_file"

  # Optional logging for failures
  if [ "$assembly" = "NA" ]; then
    echo "FAILED: $nc_id" >> failed_queries.log
  fi

done < "$input_file"

# Final merge
awk 'BEGIN{FS=OFS="\t"} FNR==NR{map[$1]=$0; next} {print map[$1], $2, $3, $4}' "$3" "$4" > "$5"

