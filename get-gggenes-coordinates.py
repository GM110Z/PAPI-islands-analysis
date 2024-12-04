import csv
import sys
import os
import pandas as pd

# Function to load the clustering table (MMseqs data)
def load_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t")  # Assuming tab-separated format

# Function to assign numerical cluster IDs
def assign_numerical_ids(clustering_df):
    clustering_df['Mmseqs_cluster'] = clustering_df.groupby('cluster_id').ngroup()
    return clustering_df, clustering_df[['locus_tag', 'Mmseqs_cluster']]

# Function to save updated cluster information
def save_updated_table(updated_df, output_file):
    updated_df.to_csv(output_file, sep="\t", index=False)

# Function to load cluster data (from the updated clustering table)
def load_cluster_data(file_path):
    cluster_data = pd.read_csv(file_path, sep="\t")
    return dict(zip(cluster_data['locus_tag'], cluster_data['Mmseqs_cluster']))

# Function to load genomic data (from the prefiltered file)
def load_genomic_coordinates(input_file):
    genomic_data = pd.read_csv(input_file, sep="\t")  # Assuming tab-separated format
    return genomic_data

# Function to merge MMseqs cluster numbers with genomic data
def merge_with_cluster_data(genomic_data, cluster_data):
    genomic_data['Mmseqs_cluster'] = genomic_data['locus_tag'].map(cluster_data)
    return genomic_data

# Function to write the results to a file
def write_results(output_file, genomic_data):
    if genomic_data.empty:
        print("No genomic data to write.")
        return
    genomic_data.to_csv(output_file, sep="\t", index=False)

# Main execution flow
if __name__ == "__main__":
    # Paths to input and output files
    original_file_path = 'MMseqs-cluster.tsv'  # Path to your original MMseqs clustering file
    updated_file_path = 'MMseqs-clusternumb.tsv'  # Path to save updated table with cluster numbers

    # Step 1: Load and map cluster data
    clustering_df = load_clustering_table(original_file_path)
    updated_df, _ = assign_numerical_ids(clustering_df)
    save_updated_table(updated_df, updated_file_path)
    cluster_data = load_cluster_data(updated_file_path)

    # Part 2: Load genomic data (the prefiltered file)
    input_file = 'TableOfFeatures.txt'  # Path to your prefiltered file
    genomic_data = load_genomic_coordinates(input_file)
    if genomic_data.empty:
        print("Genomic data could not be loaded.")
        sys.exit(1)

    # Step 3: Merge MMseqs cluster numbers with genomic data
    genomic_data_with_clusters = merge_with_cluster_data(genomic_data, cluster_data)

    # Step 4: Write the result with MMseqs cluster numbers
    output_file = 'genomic_data_with_clusters.tsv'  # Output file with MMseqs cluster numbers
    write_results(output_file, genomic_data_with_clusters)

    print("Genomic data with clusters saved.")
    
    # Cleanup temporary files
    os.system("rm OutFile.gb") 
    os.system("rm prefiltered.txt")
