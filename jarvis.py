#exclude proteins found in PADLOC from MMSeqs2 fasta output-script is called Jarvis

import sys
from Bio import SeqIO
import csv


#add list of ids from padloc to a list
list_of_accession = []
with open ('padloc-prots.txt', 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = '\t')
    for row in efetchin:
        list_of_accession.append(str(row[0]))

output_file = 'non-defence-plasmid-proteins.fasta'
with open(output_file, 'a') as output_handle:
    for record in SeqIO.parse("plasmid_rep_seq.fasta", 'fasta'):
        if not any(acc in record.name for acc in list_of_accession):
            SeqIO.write(record, output_handle, 'fasta')
