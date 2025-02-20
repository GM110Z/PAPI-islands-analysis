# heatmap for islands with many defences
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,PowerNorm


def filter_and_plot_heatmap(input_csv, output_csv, heatmap_png):
    # Load the presence-absence matrix
    df = pd.read_csv(input_csv, index_col=0)
    
    # Count the number of different systems in each island (row-wise sum)
    island_counts = df.sum(axis=1)
    
    # Filter islands with more than 5 different systems
    filtered_df = df[island_counts > 5]
    
    # Save the new presence-absence matrix
    filtered_df.to_csv(output_csv)
    # Create a masked version where 0s are fully white
    mask = presence_absence_matrix == 0
    # Plot the heatmap
    plt.figure(figsize=(40, 20))
    sns.heatmap(
    presence_absence_matrix, 
    cmap="Reds", 
    norm=PowerNorm(gamma=0.5, vmin=0.5, vmax=200),  # Log scale for better contrast
    linewidths=1.5,
    linecolor="black",
    cbar=True, 
    mask=mask,
    cbar_kws={"shrink": 0.2}# White out true absences
)
    plt.ylabel("Defence Systems")
    plt.xticks(rotation=90,fontsize=16)
    plt.yticks(fontsize=16)
    
    # Save the heatmap
    plt.savefig(heatmap_png, dpi=300, bbox_inches='tight')
    plt.show()

# Example usage
filter_and_plot_heatmap("combined_defense_systems.csv", "filtered_matrix.csv", "heatmap2.pdf")
