#generate heatmap of traits (either AMRfinder or PADLOC results)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file
data = pd.read_csv('plasmid-defence.csv', sep=",")

# Create a pivot table for presence-absence
pivot_table = data.pivot_table(index='seqid', columns='system', aggfunc='size', fill_value=0)
pivot_table.to_csv('pivot_table_output.csv')


# Generate the heatmap
plt.figure(figsize=(10, 8))

# Define custom colormap: 0 values will be white, other values will use 'viridis'
cmap = sns.cm.rocket_r  # Choose any colormap you prefer for non-zero values
cmap.set_bad('white', 1.)  # Set color for 0 values (white with full transparency)

# Plot the heatmap with the custom colormap/make annot-True if you want the numbers in the heatmap
heatmap = sns.heatmap(pivot_table, cmap=cmap, cbar=True, annot=False, fmt="d", cbar_kws={'label': 'Presence-Absence'})

# Customize the color bar
colorbar = heatmap.collections[0].colorbar
colorbar.set_ticks([0, 1, 2, 3, 4])  # Adjust ticks according to your data
colorbar.set_ticklabels(['0', '1', '2', '3', '4'])  # Adjust tick labels

# Add labels and title
plt.xlabel('Class')
plt.ylabel('Contig id')
plt.title('Presence-Absence Heatmap')

# Save the heatmap to a file
plt.savefig('heatmap_with_legend.png',bbox_inches='tight', dpi=300)

# Display the heatmap
plt.show()
