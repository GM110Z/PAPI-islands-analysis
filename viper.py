#psiblast vs databases and combo with mmseqs clustering to avoid using replicate data in computation
import os
import subprocess
import pandas as pd
from Bio import SeqIO
import time

# Step 1: MMseqs Clustering - Assign Numeric Cluster IDs
def load_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t", header=None, names=["ClusterRep", "Protein_ID"])

def assign_numerical_ids(df):
    unique_clusters = df['ClusterRep'].unique()
    cluster_to_id = {cluster: idx + 1 for idx, cluster in enumerate(unique_clusters)}
    df['ClusterRep'] = df['ClusterRep'].map(cluster_to_id)
    return df, cluster_to_id

def save_updated_table(df, file_path):
    df.to_csv(file_path, sep="\t", index=False, header=False)

def process_mmseqs_clusters(mmseqs_file, output_file):
    """Process the MMseqs clustering file to generate a table with numeric cluster IDs."""
    clustering_df = load_clustering_table(mmseqs_file)
    updated_df, _ = assign_numerical_ids(clustering_df)
    save_updated_table(updated_df, output_file)
    print(f"MMseqs clustering table with numeric cluster IDs saved to {output_file}")

# Step 2: PSI-BLAST Processing Functions
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

def add_mmseqs_clusters(filtered_results_file, mmseqs_file, output_file):
    """Add MMseqs cluster information to the filtered PSI-BLAST results."""
    print("Adding MMseqs cluster information...")

    # Load filtered PSI-BLAST results and MMseqs cluster file
    filtered_results = pd.read_csv(filtered_results_file, sep="\t")
    mmseqs_clusters = pd.read_csv(mmseqs_file, sep="\t", names=["MMSEQ_CLUSTER","Protein_ID"])

    # Merge on qseqid and Protein_ID
    merged_results = filtered_results.merge(mmseqs_clusters, left_on="qseqid", right_on="Protein_ID", how="left")

    # Save the final merged results
    merged_results.to_csv(output_file, sep="\t", index=False)
    print(f"MMseqs cluster information added. Results saved to {output_file}")

# Step 3: Main Process
def process_combined_proteome(proteome_file, output_dir, mmseqs_file, num_threads=20):
    """Process a combined proteome file against the VirulenceDatabase."""
    os.makedirs(output_dir, exist_ok=True)

    # First, process the MMseqs clustering table to assign numeric IDs
    mmseqs_clusternumb_file = os.path.join(output_dir, "MMseqs-clusternumb.tsv")
    process_mmseqs_clusters(mmseqs_file, mmseqs_clusternumb_file)

    print(f"Processing {proteome_file}...")

    base_name = os.path.basename(proteome_file).split(".")[0]
    psiblast_output = os.path.join(output_dir, f"{base_name}_psiblast.txt")
    db = "VirulenceDatabase"  # Use the correct database name

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

    final_output_with_clusters = os.path.join(output_dir, f"{base_name}_final_results_with_clusters.txt")
    add_mmseqs_clusters(filtered_results_file, mmseqs_clusternumb_file, final_output_with_clusters)

# Example usage
if __name__ == "__main__":
    proteome_file = sys.argv[1] # Path to the combined proteome file
    output_dir =  sys.argv[2] # Directory to store PSI-BLAST results
    mmseqs_file =  sys.argv[3] # Path to the MMseqs cluster file

    process_combined_proteome(proteome_file, output_dir, mmseqs_file, num_threads=8)
