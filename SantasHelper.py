import sys
from Bio import SeqIO
import pandas as pd


# Load the operon results table with clusters
df2 = pd.read_csv(sys.argv[1], sep=",")  # Operon table with Protein and operon info

# Load PFAM domain table
df1 = pd.read_csv(sys.argv[2], sep="\t")  # PFAM table with domain info

# Load the PADLOC results
df_new = pd.read_csv(sys.argv[3], sep="\t")  # PADLOC output with system and target name

# Initialize merged DataFrame
merged_final_df = pd.DataFrame()

# Merge operon table with the PFAM domains table on Protein_ID and target name
merged_df = pd.merge(df2, df1, left_on='Protein_ID', right_on='target name', how='left')
final_df = merged_df[['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand',
                      'operon_number', 'product', 'accession',"ClusterRep",
                      'query name', 'E-value', 'description of target']]
merged_final_df = final_df

# Merge with PADLOC data
merged_final_df = pd.merge(merged_final_df, df_new[['system', 'target.name']],
                           left_on='Protein_ID', right_on='target.name',
                           how='left')
# Read the file
antidefense_df = pd.read_csv(sys.argv[4], sep="\t")  # Antidefense data with sys_beg

# Concatenate 'type' and 'subtype' columns
antidefense_df['Antidefense'] = antidefense_df['type'] + '_' + antidefense_df['subtype']

# Merge with the final dataframe
merged_final_df = pd.merge(merged_final_df, antidefense_df[['Antidefense', 'sys_beg']],
                            left_on='Protein_ID', right_on='sys_beg', how='left')

# Drop the redundant 'sys_beg' column
merged_final_df = merged_final_df.drop(columns=['sys_beg'])

# Merge with AMRFinder data if the file is provided

amrfinder_df = pd.read_csv("amrshort.txt", sep='\t')  # AMRFinder data with additional columns
    
amrfinder_df['AMRFinder'] = (
    amrfinder_df['Sequence name'] + '_' + 
    amrfinder_df['Scope'] + '_' + 
    amrfinder_df['Element type'] + '_' + 
    amrfinder_df['Element subtype'] + '_' + 
    amrfinder_df['Class'] + '_' + 
    amrfinder_df['Subclass']
)

amrfinder_df = amrfinder_df[['Protein identifier', 'AMRFinder']]

merged_final_df = pd.merge(merged_final_df, amrfinder_df,
                           left_on='Protein_ID', right_on='Protein identifier',how='left')

# Drop the redundant 'Protein identifier' column
merged_final_df = merged_final_df.drop(columns=['Protein identifier'])

# Save the merged data to the specified output file
merged_final_df.to_csv("Final-R-file.tsv", index=False, sep='\t')

