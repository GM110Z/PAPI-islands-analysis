#group-by-sizes
import pandas as pd

df = pd.read_csv('chromosome.tsv', delimiter='\t', header=None, names=['Sequence_ID', 'Size'])

# Sort the DataFrame by the 'Size' column
df = df.sort_values(by='Size')

# Define the size ranges for grouping (as specified)
size_ranges = [(20000, 40000), (41000, 60000), (61000, 80000), (81000, 150000)]

# Group the sequences based on size ranges
groups = {}
for low, high in size_ranges:
    # Filter rows based on size range
    group = df[(df['Size'] >= low) & (df['Size'] <= high)]
    groups[f'{low}-{high}'] = group

# Write each group to a separate file
for range_str, group in groups.items():
    filename = f'group_{range_str}.txt'  # Create a filename based on the range
    group.to_csv(filename, sep='\t', index=False, header=False)  # Write to file

    print(f"Written {filename}")
