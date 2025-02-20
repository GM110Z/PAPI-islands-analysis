import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm,PowerNorm


# Define the directory where files are stored
data_dir = "/Users/giusym/Desktop/Bioinformatics/Scripts to Plot Figures/"

# Get a list of all files
files = [f for f in os.listdir(data_dir) if f.endswith(".tsv")]

# Dictionary to store occurrence data for each island
data_dict = {}
all_defense_systems = set()

# Process each file
for file in files:
    file_path = os.path.join(data_dir, file)
    
    # Read the file (assume tab-separated or CSV format)
    df = pd.read_csv(file_path, sep="\t", header=0)  # Adjust delimiter if needed

    # Standardize column names
    df.columns = ["Defense system", "Occurrence"]

    # Extract island name from the file name
    island_name = os.path.splitext(file)[0]

    # Store occurrences in a dictionary
    data_dict[island_name] = dict(zip(df["Defense system"], df["Occurrence"]))

    # Collect unique defense systems
    all_defense_systems.update(df["Defense system"])

# Convert collected data into a DataFrame
all_defense_systems = sorted(all_defense_systems)
presence_absence_matrix = pd.DataFrame(0, index=all_defense_systems, columns=data_dict.keys())

# Fill in the matrix with occurrence values
for island, defenses in data_dict.items():
    for system, occurrence in defenses.items():
        presence_absence_matrix.at[system, island] = occurrence  # Store actual counts

# Create a masked version where 0s are fully white
mask = presence_absence_matrix == 0

# Increase figure width for better readability
plt.figure(figsize=(40, 30))  # Increased width to 18 inches

# Plot the heatmap with a log-normalized color scale
sns.heatmap(
    presence_absence_matrix, 
    cmap="Reds", 
    norm=PowerNorm(gamma=0.5, vmin=0.5, vmax=200),  # Log scale for better contrast
    linewidths=1.5,
    linecolor="black",
    cbar=True, 
    square=True, 
    mask=mask,
    cbar_kws={"shrink": 0.2}# White out true absences
)
plt.ylabel("Defense Systems")
plt.xticks(rotation=90,fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("heatmap.pdf",  dpi=300)
plt.show()
