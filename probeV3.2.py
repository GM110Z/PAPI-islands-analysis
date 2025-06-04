import os
import sys
import glob
import subprocess as sp
import csv
import time
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Run HMMER search
def run_hmmer(model, fasta_file, output_file, threshold=10):
    print(f"Running HMMER on {fasta_file} with model {model}...")
    hmmsearch_cmd = f"hmmsearch -T {threshold} --incT {threshold} -o log --domtblout {output_file} {model} {fasta_file}"
    sp.run(hmmsearch_cmd, shell=True)

# Parse HMMER output to CSV
def parse_hmmer_output(hmm_output_file, csv_output_file):
    if not os.path.exists(hmm_output_file):
        print(f"Warning: {hmm_output_file} does not exist.")
        return
    with open(hmm_output_file, newline='') as input_file:
        for qresult in SearchIO.parse(input_file, 'hmmscan3-domtab'):
            hits = qresult.hits
            if hits:
                with open(csv_output_file, mode='w', newline='') as parsed_output:
                    writer = csv.writer(parsed_output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for hit in hits:
                        writer.writerow([hit.id, hit.description, hit.evalue, hit.bitscore, hit.seq_len])
                print(f"Hits found and written to {csv_output_file}")
            else:
                print(f"No hits found in {hmm_output_file}")

# Extract all protein IDs from CSV
def extract_protein_ids(csv_file):
    if not os.path.exists(csv_file):
        print(f"Warning: {csv_file} does not exist.")
        return []
    protein_ids = []
    with open(csv_file, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            protein_ids.append(row[0])
    return protein_ids

# Get CDS location from GenBank
def get_protein_location(genbank_file, protein_id):
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and protein_id in feature.qualifiers.get("protein_id", []):
                return int(feature.location.start), int(feature.location.end), record.id
    return None, None, None

# Extract GenBank region and translated CDS proteins
def extract_gbk_and_proteins(genbank_file, chrom_id, start, end, gbk_outfile, faa_outfile):
    for record in SeqIO.parse(genbank_file, "genbank"):
        if record.id != chrom_id:
            continue

        # Extract genomic region
        region = record[start:end]
        region.id = f"{record.id}_{start}_{end}"
        region.description = f"Extracted from {start} to {end}"

        # Write GenBank region
        with open(gbk_outfile, "a") as gbk_out:
            SeqIO.write(region, gbk_out, "genbank")
        print(f"Wrote GenBank region to {gbk_outfile}")

        # Extract translated proteins
        proteins = []
        for feature in record.features:
            if feature.type == "CDS" and feature.location.start >= start and feature.location.end <= end:
                qualifiers = feature.qualifiers
                if "translation" in qualifiers:
                    prot_seq = qualifiers["translation"][0]
                    prot_id = qualifiers.get("protein_id", [f"{record.id}_{feature.location.start}_{feature.location.end}"])[0]
                    desc = qualifiers.get("product", ["unknown"])[0]
                    proteins.append(SeqRecord(
                        Seq(prot_seq),
                        id=prot_id,
                        description=desc
                    ))

        if proteins:
            with open(faa_outfile, "a") as faa_out:
                SeqIO.write(proteins, faa_out, "fasta")
            print(f"Wrote {len(proteins)} proteins to {faa_outfile}")
        else:
            print("No proteins found in extracted region.")

# Main logic
def process_files(model1, model2, protein_directory, genbank_directory):
    fasta_files = glob.glob(f"{protein_directory}/*.fasta")
    genbank_files = glob.glob(f"{genbank_directory}/*.gb")
    genbank_dict = {os.path.splitext(os.path.basename(f))[0]: f for f in genbank_files}
    print(f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files.")

    for fasta_file in fasta_files:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        if base_name in genbank_dict:
            genbank_file = genbank_dict[base_name]
            print(f"Processing: {fasta_file} and {genbank_file}")

            run_hmmer(model1, fasta_file, "HMM_output-search1.txt")
            run_hmmer(model2, fasta_file, "HMM_output-search2.txt")

            parse_hmmer_output("HMM_output-search1.txt", "protein1.csv")
            parse_hmmer_output("HMM_output-search2.txt", "protein2.csv")

            protein_ids1 = extract_protein_ids("protein1.csv")
            protein_ids2 = extract_protein_ids("protein2.csv")

            if not protein_ids1 or not protein_ids2:
                print("One or both protein ID lists could not be extracted.")
            else:
                best_span = 0
                best_coords = None

                for pid1 in protein_ids1:
                    s1, e1, chrom1 = get_protein_location(genbank_file, pid1)
                    if not all([s1, e1, chrom1]):
                        continue

                    for pid2 in protein_ids2:
                        s2, e2, chrom2 = get_protein_location(genbank_file, pid2)
                        if not all([s2, e2, chrom2]):
                            continue

                        if chrom1 == chrom2:
                            span = max(e1, e2) - min(s1, s2)
                            dist = abs(s1 - e2) if s1 < s2 else abs(s2 - e1)

                            if dist <= 40000 and span > best_span:
                                best_span = span
                                best_coords = (chrom1, min(s1, s2), max(e1, e2))

                if best_coords:
                    chrom, start, end = best_coords
                    with open("Locus-coordinates.txt", "a") as f:
                        f.write(f"{chrom}\t{start}\t{end}\n")
                    print(f"Written coordinates: {chrom}, {start}-{end}")

                    extract_gbk_and_proteins(
                        genbank_file=genbank_file,
                        chrom_id=chrom,
                        start=start,
                        end=end,
                        gbk_outfile="Extracted-loci.gbk",
                        faa_outfile="Extracted-proteins.faa"
                    )
                else:
                    print("No valid model1-model2 pair found within distance threshold.")

            for file in ["HMM_output-search1.txt", "HMM_output-search2.txt", "protein1.csv", "protein2.csv"]:
                if os.path.exists(file):
                    os.remove(file)
            time.sleep(2.5)
        else:
            print(f"No matching GenBank file found for {fasta_file}.")

# Entry point
if __name__ == "__main__":
    model1 = sys.argv[1]
    model2 = sys.argv[2]
    protein_directory = sys.argv[3]
    genbank_directory = sys.argv[3]
    process_files(model1, model2, protein_directory, genbank_directory)
