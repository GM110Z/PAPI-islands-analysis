import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import sys

# Directory and file pattern
input_pattern = sys.argv[1]  #path for files to check with wildcard at the end with files extensions 
output_file = "overlap_results.csv"

def find_overlaps(df):
    overlaps = []
    total_pairs = len(df) * (len(df) - 1) // 2  # Total comparisons
    processed_pairs = 0

    print(f"Total comparisons to check: {total_pairs}")
    
    for i, row1 in df.iterrows():
        for j, row2 in df.iterrows():
            if i >= j:
                continue  # Avoid redundant comparisons and self-comparison
            
            processed_pairs += 1
            if processed_pairs % 100000 == 0:  # Print every 100,000 comparisons
                print(f"Processed {processed_pairs}/{total_pairs} comparisons...")

            if row1["Accession"] == row2["Accession"]:  # Check same accession
                # Check overlap in start and stop coordinates
                if max(row1["Start"], row2["Start"]) <= min(row1["Stop"], row2["Stop"]):
                    overlap_size = min(row1["Stop"], row2["Stop"]) - max(row1["Start"], row2["Start"]) + 1
                    overlaps.append({
                        "File1": row1["File"],
                        "Pathogenicity_Island1": row1["Pathogenicity_Island"],
                        "Accession1": row1["Accession"],
                        "Start1": row1["Start"],
                        "Stop1": row1["Stop"],
                        "File2": row2["File"],
                        "Pathogenicity_Island2": row2["Pathogenicity_Island"],
                        "Accession2": row2["Accession"],
                        "Start2": row2["Start"],
                        "Stop2": row2["Stop"],
                        "Overlap_Size": overlap_size
                    })
    print(f"Finished processing {processed_pairs}/{total_pairs} comparisons.")
    return overlaps

# Read all files matching the pattern
all_data = []
file_list = glob.glob(input_pattern)
print(f"Files found: {file_list}")

for file_path in file_list:
    file_name = os.path.basename(file_path)
    pathogenicity_island = file_name.split(".")[0]
    print(f"Processing file: {file_name}")
    df = pd.read_csv(file_path, sep="\t")  # Adjust delimiter if needed
    df["File"] = file_name
    df["Pathogenicity_Island"] = pathogenicity_island
    all_data.append(df)

# Combine data
combined_df = pd.concat(all_data, ignore_index=True)
print(f"Combined DataFrame:\n{combined_df.head()}")

# Ensure numeric types for Start and Stop columns
combined_df["Start"] = pd.to_numeric(combined_df["Start"], errors="coerce")
combined_df["Stop"] = pd.to_numeric(combined_df["Stop"], errors="coerce")

# Find overlaps
overlap_results = find_overlaps(combined_df)
print(f"Number of overlaps found: {len(overlap_results)}")

# Save results
if overlap_results:
    overlap_df = pd.DataFrame(overlap_results)
    overlap_df.to_csv(output_file, index=False)
    print(f"Overlap results saved to {output_file}")
else:
    print("No overlaps found.")

    
    
    
    
    
#plot these
# Load the data from the CSV file
df = pd.read_csv("overlap_results.csv", sep =",")
# Create a scatter plot
plt.figure(figsize=(10, 6))

# Loop through unique Pathogenicity_Island1 values to create separate scatter points
for key, grp in df.groupby(['Pathogenicity_Island1']):
    plt.scatter(grp['Pathogenicity_Island2'], grp['Overlap_Size'], label=key)

plt.title('Overlap Size between Pathogenicity Islands')
plt.xlabel('Pathogenicity Island 2')
plt.ylabel('Overlap Size')
plt.xticks(rotation=45)
plt.legend(title='Pathogenicity Island 1')
plt.grid()
plt.tight_layout()
plt.show()
