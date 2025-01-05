import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, dendrogram

# Step 1: Gather all Hotspot*.txt files in the directory
files = [f for f in os.listdir() if f.startswith("Hotspot") and f.endswith(".txt")]

# Step 2: Initialize data structure to store presence/absence
hotspots = []
nuccore_set = set()

# Step 3: Process each file to collect data
for file in files:
    hotspot_id = os.path.splitext(file)[0]  # Extract hotspot name without extension
    hotspots.append(hotspot_id)
    
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t", header=None, names=["nuccore", "start", "end"])
    
    # Add nuccore IDs to the set
    nuccore_set.update(df["nuccore"].unique())

# Step 4: Create a DataFrame for the presence/absence matrix
nuccore_list = sorted(nuccore_set)
presence_absence = pd.DataFrame(index=nuccore_list, columns=hotspots, data=0)

# Step 5: Populate the presence/absence matrix
for file in files:
    hotspot_id = os.path.splitext(file)[0]  # Extract hotspot name without extension
    
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t", header=None, names=["nuccore", "start", "end"])
    
    # Mark presence in the matrix
    for nuccore in df["nuccore"].unique():
        presence_absence.at[nuccore, hotspot_id] = 1

# Step 6: Define a custom color palette
cmap = mcolors.ListedColormap(["#FFFFFF", "#800000"])  # White for absence, Dark red for presence

# Step 7: Perform hierarchical clustering (linkage matrix)
linkage_matrix = linkage(presence_absence, method='ward', metric='euclidean')

# Step 8: Create a side-by-side layout for the heatmap and dendrogram
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [5, 1]})

# Plot the heatmap
sns.heatmap(presence_absence, cmap=cmap, cbar_kws={"label": "Presence (1) / Absence (0)"},
            annot=False, fmt="d", linewidths=0, linecolor='none', 
            yticklabels=True, xticklabels=hotspots, ax=ax1)

# Plot the dendrogram using the same linkage matrix
dendrogram_info = dendrogram(
    linkage_matrix,
    labels=presence_absence.index,
    orientation='left',  # Vertical dendrogram (matching the heatmap)
    color_threshold=None,  # No color threshold
    ax=ax2
)

# Remove grid, labels, and borders from both axes
ax1.grid(False)  # Remove the grid from the heatmap
ax1.set_yticklabels([])  # Remove y-axis labels from heatmap
ax1.set_yticks([])  # Remove y-ticks
ax2.grid(False)  # Remove the grid from the dendrogram
ax2.set_xticklabels([])  # Remove x-axis labels from dendrogram
ax2.set_yticklabels([])  # Remove y-axis labels from dendrogram
ax2.set_xticks([])  # Remove x-ticks from dendrogram
ax2.set_yticks([])  # Remove y-ticks from dendrogram

# Customize plot title and axes
ax1.set_title("Hotspot Presence/Absence Heatmap", fontsize=14)
ax2.set_title("Nuccore Dendrogram", fontsize=14)

# Save the side-by-side plot as PDF
plt.tight_layout(pad=0)
plt.savefig("heatmap_and_dendrogram_side_by_side.pdf", dpi=300)
plt.show()

# Reorder presence/absence matrix to match dendrogram order
reordered_labels = [dendrogram_info['ivl'][i] for i in range(len(dendrogram_info['ivl']))]
presence_absence_reordered = presence_absence.loc[reordered_labels]

# Save the reordered presence/absence matrix as a CSV file
presence_absence_reordered.to_csv("presence_absence_matrix_reordered.csv")

# Write tree file metadata with labels
with open("tree_file_metadata_with_labels_reordered.txt", "w") as tree_file:
    tree_file.write("Linkage Matrix (Tree Data with Labels)\n")
    tree_file.write("-------------------------------------------------\n")
    
    # Iterate through the reordered dendrogram information and include labels
    for label in reordered_labels:
        tree_file.write(f"{label}\n")
    
    tree_file.write("Linkage matrix information:\n")
    tree_file.write("-------------------------------------------------\n")
    
    # Write out the linkage matrix, now in the order from the dendrogram
    for row in linkage_matrix:
        tree_file.write(f"{row[0]} {row[1]} {row[2]} {row[3]}\n")
