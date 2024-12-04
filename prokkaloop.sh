
#!/bin/bash
##edit header
# Directory containing your FASTA files
GENOME_DIR=$1

# Loop through each FASTA file in the directory
for fasta_file in "$GENOME_DIR"/*.fasta; do
  echo "Cleaning headers in $fasta_file..."

  # Clean headers by extracting only the accession number (up to the first space or colon)
  awk '/^>/ {split($1, arr, ":"); print arr[1]; next} {print}' "$fasta_file" > "${fasta_file}.tmp"

  # Overwrite the original file with the cleaned headers
  mv "${fasta_file}.tmp" "$fasta_file"
done

echo "All headers cleaned and files updated."


##now run prokka
echo "Annotating genomes now"                

# Define the directory containing genome files and the output directory
GENOME_DIR=$1       # Path to the directory containing genome files
OUTPUT_DIR=$1  # Path to save Prokka output

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each genome file in the directory
for genome_file in "$GENOME_DIR"/*.fasta; do
  # Extract the base name of the genome file (without path and extension)
  base_name=$(basename "$genome_file" .fasta)

  # Run Prokka with the genome file
  prokka --outdir "$OUTPUT_DIR/$base_name" --prefix "$base_name" "$genome_file"

  echo "Prokka annotation completed for $genome_file"
done

echo "All genome files processed."

