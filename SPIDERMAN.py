import requests
import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate
import csv
import sys


# Function to get KEGG entry details
def get_kegg_entry(kegg_id):
    url = f'http://rest.kegg.jp/get/{kegg_id}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return None

# Function to parse pathways from KEGG entry
def parse_pathways(entry):
    pathways = []
    for line in entry.split('\n'):
        if line.startswith('PATHWAY'):
            pathways.append(' '.join(line.split()[1:]))
    return pathways

# List of KEGG IDs to process/import from a file from blast Koala
kegg_ids = [ ]

with open (sys.argv[1], 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = '\t')
    for row in efetchin:
        kegg_ids.append(str(row[1]))

# Dictionary to store pathways for each KEGG ID
kegg_pathways = {}

# Fetch details and parse pathways for each KEGG ID
for kegg_id in kegg_ids:
    entry = get_kegg_entry(kegg_id)
    if entry:
        pathways = parse_pathways(entry)
        kegg_pathways[kegg_id] = pathways
    else:
        kegg_pathways[kegg_id] = []

# Prepare data for tabular format and plotting
table_data = []
pathway_names = []
pathway_counts = []

for kegg_id, pathways in kegg_pathways.items():
    if pathways:
        for pathway in pathways:
            table_data.append([kegg_id, pathway])
            pathway_names.append(pathway)
            pathway_counts.append(1)

# Create DataFrame for the table data
df = pd.DataFrame(table_data, columns=['KEGG ID', 'Pathway'])

# Save the table to a CSV file
csv_file = 'kegg_pathways.csv'
df.to_csv(csv_file, index=False)

# Count the occurrences of each pathway
pathway_counts_dict = {pathway: pathway_names.count(pathway) for pathway in set(pathway_names)}
pathway_names_unique = list(pathway_counts_dict.keys())
pathway_counts = list(pathway_counts_dict.values())
# Print the pathways in tabular format
print(tabulate(table_data, headers=['KEGG ID', 'Pathway'], tablefmt='grid'))

# Count the occurrences of each pathway
pathway_counts_dict = {pathway: pathway_names.count(pathway) for pathway in set(pathway_names)}
pathway_names_unique = list(pathway_counts_dict.keys())
pathway_counts = list(pathway_counts_dict.values())
# Print the pathways in tabular format
print(tabulate(table_data, headers=['KEGG ID', 'Pathway'], tablefmt='grid'))

# Plot all pathways on a pie chart
plt.figure(figsize=(10, 10))
plt.pie(pathway_counts, labels=pathway_names_unique, autopct='%1.1f%%', startangle=140, colors=plt.cm.Paired.colors)
plt.title('Pathways for All KEGG IDs')

# Save the plot to a file
plot_file = 'kegg_pathways.png'
plt.savefig(plot_file)

# Show the plot
plt.show()
