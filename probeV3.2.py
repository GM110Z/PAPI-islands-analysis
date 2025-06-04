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

def run_hmmer(model, fasta_file, output_file, threshold=10):
    print(f"Running HMMER on {fasta_file} with model {model}...")
    hmmsearch_cmd = f"hmmsearch -T {threshold} --incT {threshold} -o log --domtblout {output_file} {model} {fasta_file}"
    sp.run(hmmsearch_cmd, shell=True)

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

def extract_protein_ids(csv_file):
    if not os.path.exists(csv_file):
        print(f"Warning: {csv_file} does not exist.")
        return []
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        return [row[0] for row in csv.reader(f)]

def get_protein_location(genbank_file, protein_id):
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and protein_id in feature.qualifiers.get("protein_id", []):
                return int(feature.location.start), int(feature.location.end), record.id
    return None, None, None

def extract_gbk_and_proteins(genbank_file, chrom_id, start, end, gbk_outfile, faa_outfile):
    for record in SeqIO.parse(genbank_file, "genbank"):
        if record.id != chrom_id:
            continue

        # GenBank region
        region = record[start:end]
        region.id = f"{record.id}_{start}_{end}"
        region.description = f"Extracted from {start} to {end}"
        with open(gbk_outfile, "w") as gbk_out:
            SeqIO.write(region, gbk_out, "genbank")
        print(f"Wrote GenBank region to {gbk_outfile}")

        # CDS proteins
        proteins = []
        for feature in record.features:
            if feature.type == "CDS" and feature.location.start >= start and feature.location.end <= end:
                qualifiers = feature.qualifiers
                if "translation" in qualifiers:
                    prot_seq = qualifiers["translation"][0]
                    prot_id = qualifiers.get("protein_id", [f"{record.id}_{feature.location.start}_{feature.location.end}"])[0]
                    desc = qualifiers.get("product", ["unknown"])[0]
                    proteins.append(SeqRecord(Seq(prot_seq), id=prot_id, description=desc))

        if proteins:
            with open(faa_outfile, "w") as faa_out:
                SeqIO.write(proteins, faa_out, "fasta")
            print(f"Wrote {len(proteins)} proteins to {faa_outfile}")
        else:
            print("No proteins found in region.")

def process_files(model1, model2, protein_dir, genbank_dir):
    fasta_files = glob.glob(f"{protein_dir}/*.fasta")
    genbank_files = glob.glob(f"{genbank_dir}/*.gb")
    genbank_dict = {os.path.splitext(os.path.basename(f))[0]: f for f in genbank_files}
    print(f"Found {len(fasta_files)} FASTA and {len(genbank_files)} GenBank files.")

    for fasta_file in fasta_files:
        base = os.path.splitext(os.path.basename(fasta_file))[0]
        genbank_file = genbank_dict.get(base)
        if not genbank_file:
            print(f"No GenBank file found for {base}")
            continue

        print(f"Processing: {fasta_file} and {genbank_file}")

        run_hmmer(model1, fasta_file, f"{base}_hmm1.txt")
        run_hmmer(model2, fasta_file, f"{base}_hmm2.txt")
        parse_hmmer_output(f"{base}_hmm1.txt", f"{base}_prot1.csv")
        parse_hmmer_output(f"{base}_hmm2.txt", f"{base}_prot2.csv")

        pids1 = extract_protein_ids(f"{base}_prot1.csv")
        pids2 = extract_protein_ids(f"{base}_prot2.csv")

        best_span = 0
        best_coords = None

        for pid1 in pids1:
            s1, e1, chrom1 = get_protein_location(genbank_file, pid1)
            if not all([s1, e1, chrom1]):
                continue

            for pid2 in pids2:
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
                f.write(f"{base}\t{chrom}\t{start}\t{end}\n")
            print(f"Wrote locus for {base}: {chrom}:{start}-{end}")

            gbk_out = f"{base}_locus.gbk"
            faa_out = f"{base}_proteins.faa"
            extract_gbk_and_proteins(genbank_file, chrom, start, end, gbk_out, faa_out)
        else:
            print(f"No valid locus found for {base}.")

        for tmp in [f"{base}_hmm1.txt", f"{base}_hmm2.txt", f"{base}_prot1.csv", f"{base}_prot2.csv"]:
            if os.path.exists(tmp):
                os.remove(tmp)

        time.sleep(2.5)

if __name__ == "__main__":
    model1 = sys.argv[1]
    model2 = sys.argv[2]
    protein_dir = sys.argv[3]
    genbank_dir = sys.argv[3]
    process_files(model1, model2, protein_dir, genbank_dir)

