#heatmap presence absence for hotspot10-11--no dendrogram and hierarchical clustering not possible for 2 datasets only.

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Step 1: Gather all Hotspot*.txt files in the directory
files = [f for f in os.listdir() if f.startswith("Hotspot") and f.endswith(".txt")]

# Step 2: Initialize data structure to store presence/absence
hotspots = []
nuccore_set = set()

# Step 3: Process each file to collect data
for file in files:
    hotspot_id = os.path.splitext(file)[0]  # Extract hotspot name without extension
    hotspots.append(hotspot_id)
    
    df = pd.read_csv(file, sep="\t", header=None, names=["nuccore", "start", "end"])
    nuccore_set.update(df["nuccore"].unique())

# Step 4: Create presence/absence matrix
nuccore_list = sorted(nuccore_set)
presence_absence = pd.DataFrame(index=nuccore_list, columns=hotspots, data=0)

# Step 5: Populate the matrix
for file in files:
    hotspot_id = os.path.splitext(file)[0]
    df = pd.read_csv(file, sep="\t", header=None, names=["nuccore", "start", "end"])
    for nuccore in df["nuccore"].unique():
        presence_absence.at[nuccore, hotspot_id] = 1

# Step 6: Define custom colour map
cmap = mcolors.ListedColormap(["#FFFFFF", "#800000"])  # White = 0, Dark red = 1

# Step 7: Plot heatmap only (no dendrogram, no row labels, keep column labels)
plt.figure(figsize=(12, 6))
sns.heatmap(
    presence_absence,
    cmap=cmap,
    cbar_kws={"label": "Presence (1) / Absence (0)"},
    annot=False,
    fmt="d",
    linewidths=0,
    linecolor='none',
    yticklabels=False,            # REMOVE row labels
    xticklabels=hotspots          # KEEP column labels
)

plt.title("Hotspot Presence/Absence Heatmap", fontsize=14)
plt.xlabel("Hotspots")
plt.ylabel("")  # Remove y-axis label
plt.tight_layout()
plt.savefig("heatmap_only.pdf", dpi=300)
plt.show()
