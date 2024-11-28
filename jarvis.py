
#Filter blast ouput


import pandas as pd

# Load your BLAST data into a DataFrame
df = pd.read_csv("newpapi.txt", sep="\t", header=None, 
                 names=["query_id", "subject_id", "identity", "alignment_length", 
                        "mismatches", "gap_opens", "q_start", "q_end", 
                        "s_start", "s_end", "evalue", "bit_score"])

# Clean any unwanted characters (e.g., non-breaking spaces) from 's_start' and 's_end' columns
df["s_start"] = df["s_start"].replace({r"\xa0": ""}, regex=True).str.strip()
df["s_end"] = df["s_end"].replace({r"\xa0": ""}, regex=True).str.strip()

# Convert 's_start' and 's_end' to integers, handle errors if conversion fails
df["s_start"] = pd.to_numeric(df["s_start"], errors="coerce")
df["s_end"] = pd.to_numeric(df["s_end"], errors="coerce")

# Step 1: Sort by subject and subject start position
df = df.sort_values(by=["subject_id", "s_start"])

# Step 2: Initialize variables to store regions
regions = []
current_region = None

# Step 3: Loop through the rows to group the hits into regions
for i, row in df.iterrows():
    subject_id = row["subject_id"]
    s_start = row["s_start"]
    s_end = row["s_end"]

    # Skip rows with invalid numeric values (NaN after coercion)
    if pd.isna(s_start) or pd.isna(s_end):
        continue

    # If it's the first hit, start a new region
    if current_region is None:
        current_region = {
            "subject_id": subject_id,
            "start": s_start,
            "end": s_end
        }
    else:
        # Check if the current hit is close enough to the last region
        if subject_id == current_region["subject_id"] and (s_start - current_region["end"] <= 110000):
            # Extend the current region if the current hit is close enough
            current_region["end"] = max(current_region["end"], s_end)
        else:
            # If the hits are far apart, save the current region and start a new one
            region_size = current_region["end"] - current_region["start"]
            if 20000 <= region_size <= 180000:  # Filter based on size
                regions.append(current_region)
            # Start a new region
            current_region = {
                "subject_id": subject_id,
                "start": s_start,
                "end": s_end
            }

# Final region check after loop finishes
if current_region is not None:
    region_size = current_region["end"] - current_region["start"]
    if 20000 <= region_size <= 180000:  # Filter based on size
        regions.append(current_region)

# Convert regions to a DataFrame for easy output as a table
regions_df = pd.DataFrame(regions)

# Convert start and end positions to integers (removing .0)
regions_df["start"] = regions_df["start"].astype(int)
regions_df["end"] = regions_df["end"].astype(int)

# Write the table to a CSV file
regions_df.to_csv("regions_output_table.csv", index=False)

# Optional: Print the DataFrame to see the output
print(regions_df)
