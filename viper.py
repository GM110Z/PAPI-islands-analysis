import os
import subprocess
import pandas as pd
import time
import sys

# Step 1: PSI-BLAST Processing Functions
def run_psiblast(query_file, db, output_file, evalue=0.01, num_threads=20):
    """Run PSI-BLAST for the given query against the database."""
    psiblast_cmd = [
        "psiblast",
        "-query", query_file,
        "-db", db,
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore stitle",
        "-num_threads", str(num_threads),
        "-out", output_file
    ]
    subprocess.run(psiblast_cmd, check=True)

def filter_psiblast_results(psiblast_output, coverage_threshold=0.4):
    """Filter PSI-BLAST results based on e-value and coverage."""
    filtered_hits = []
    
    # Load PSI-BLAST results into a DataFrame
    columns = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "stitle"]
    results = pd.read_csv(psiblast_output, sep="\t", names=columns)
    
    for _, row in results.iterrows():
        query_coverage = row['length'] / row['qlen']
        target_coverage = row['length'] / row['slen']

        if row['evalue'] <= 0.01 and query_coverage >= coverage_threshold and target_coverage >= coverage_threshold:
            filtered_hits.append(row)

    return pd.DataFrame(filtered_hits)

def resolve_similar_hits(filtered_hits):
    """Resolve similar hits by keeping the best scoring one based on bit score."""
    resolved_hits = []

    grouped = filtered_hits.groupby("qseqid")
    for _, group in grouped:
        group_sorted = group.sort_values(by="bitscore", ascending=False)
        best_hit = group_sorted.iloc[0]
        resolved_hits.append(best_hit)

    return pd.DataFrame(resolved_hits)

# Step 2: Main Process
def process_combined_proteome(proteome_file, output_dir, num_threads=20):
    """Process a combined proteome file against the VirulenceDatabase."""
    os.makedirs(output_dir, exist_ok=True)

    print(f"Processing {proteome_file}...")

    base_name = os.path.basename(proteome_file).split(".")[0]
    psiblast_output = os.path.join(output_dir, f"{base_name}_psiblast.txt")
    db = sys.argv[1]  # Use the correct database name

    start_time = time.time()
    run_psiblast(proteome_file, db, psiblast_output, num_threads=num_threads)
    print(f"PSI-BLAST completed in {time.time() - start_time:.2f} seconds.")

    print("Filtering results...")
    filtered_hits = filter_psiblast_results(psiblast_output)
    print(f"Filtering completed. {len(filtered_hits)} hits passed the thresholds.")

    print("Resolving similar hits...")
    resolved_hits = resolve_similar_hits(filtered_hits)
    print(f"Resolved similar hits. {len(resolved_hits)} unique hits retained.")

    filtered_results_file = os.path.join(output_dir, f"{base_name}_filtered_results.txt")
    resolved_hits.to_csv(filtered_results_file, sep="\t", index=False)

# Example usage
if __name__ == "__main__":
    proteome_file = sys.argv[2]  # Path to the combined proteome file
    output_dir = sys.argv[3]  # Directory to store PSI-BLAST results

    process_combined_proteome(proteome_file, output_dir, num_threads=8)
