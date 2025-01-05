# compare padloc/amrfinder/defense finder chromosome hits with the original jarvis output to create a shorter 
#list of inputs for phorific when working with a lot of data

import pandas as pd

# Load Table 1 and Table 2 from their respective files
table1 = pd.read_csv("pagi3", sep="\t", header=None, names=["Accession","Start","Stop"])  # Adjust the separator if needed
table2 = pd.read_csv("pagi-padloc.txt", header=None,sep="\t", names=["Accession"])

# Filter rows in Table 1 where the Accession matches any in Table 2
matching_rows = table1[table1['Accession'].isin(table2['Accession'])]

# Save the matching rows to a new file
matching_rows.to_csv("matching_rows.csv", index=False)

print("Matching rows have been saved to 'matching_rows.csv'.")
