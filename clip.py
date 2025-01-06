#compare probe output with the blast table to ensure hotspots are in the PAPI's interval

import pandas as pd

# Load tables from files (adjust the paths as necessary)
table1_path = "newpapi1"  # Path to Table 1
table2_path = "Hotspot2-papi1.txt"  # Path to Table 2

# Read the tables into DataFrames
df1 = pd.read_csv(table1_path, sep="\t", header=None, names=["ID", "Start", "End", "Size"])
df2 = pd.read_csv(table2_path, sep="\t", header=None, names=["ID", "Start", "End"])

# Ensure numeric columns for start and end positions
df1["Start"] = pd.to_numeric(df1["Start"], errors="coerce")
df1["End"] = pd.to_numeric(df1["End"], errors="coerce")
df2["Start"] = pd.to_numeric(df2["Start"], errors="coerce")
df2["End"] = pd.to_numeric(df2["End"], errors="coerce")

# Ensure IDs are strings and clean up
df1["ID"] = df1["ID"].astype(str).str.strip()
df2["ID"] = df2["ID"].astype(str).str.strip()

# Function to compare ranges
def compare_ranges(df1, df2):
    results = []
    
    # Iterate through each row in Table 2
    for _, row2 in df2.iterrows():
        id2, start2, end2 = row2["ID"], row2["Start"], row2["End"]

        # Find matching rows in Table 1 by ID
        matching_rows = df1[df1["ID"] == id2]

        # Check if range from Table 2 is contained within any range from Table 1
        for _, row1 in matching_rows.iterrows():
            start1, end1 = row1["Start"], row1["End"]

            # Check if the range from Table 2 is fully contained in Table 1
            if start1 <= start2 and end1 >= end2:
                # If contained, add the corresponding row from Table 2 to results (without Size)
                results.append({
                    "ID": id2,
                    "Start": start2,
                    "End": end2
                })
                
    # Convert the results into a DataFrame and return
    return pd.DataFrame(results)

# Compare ranges between Table 1 and Table 2
result_df = compare_ranges(df1, df2)

# Output the result to a file and print the DataFrame
output_path = "contained_ranges_results.txt"
result_df.to_csv(output_path, sep="\t", index=False)

# Print the results
print("\nFinal Results:\n", result_df)
print(f"\nResults saved to: {output_path}")
