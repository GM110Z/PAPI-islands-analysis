#compare probe output with the blast table to ensure hotspots are in the PAPI's interval

import pandas as pd

# Load tables from files
# Replace 'table1.txt' and 'table2.txt' with your actual file paths
table1_path = "summary_sizes-small.txt"
table2_path = "Locus-coordinates.txt"


df1 = pd.read_csv(table1_path, sep="\t", header=None, names=["ID", "Start", "End"])
df2 = pd.read_csv(table2_path, sep="\t", header=None, names=["ID", "Start", "End"])

# Ensure numeric columns and clean IDs
df1["Start"] = pd.to_numeric(df1["Start"], errors="coerce")
df1["End"] = pd.to_numeric(df1["End"], errors="coerce")
df2["Start"] = pd.to_numeric(df2["Start"], errors="coerce")
df2["End"] = pd.to_numeric(df2["End"], errors="coerce")

df1["ID"] = df1["ID"].str.strip()
df2["ID"] = df2["ID"].str.strip()

# Function to compare ranges
def check_ranges(df1, df2):
    results = []
    for _, row2 in df2.iterrows():
        id2, start2, end2 = row2["ID"], row2["Start"], row2["End"]

        # Find matching rows in Table 1
        match = df1[
            (df1["ID"] == id2) & 
            (df1["Start"] <= start2) & 
            (df1["End"] >= end2)
        ]

        # Append result: True if match found, False otherwise
        results.append({
            "ID": id2,
            "Start": start2,
            "End": end2,
            "Included": not match.empty
        })

    return pd.DataFrame(results)

# Compare ranges
result_df = check_ranges(df1, df2)

# Save results to a file and print a summary
result_df.to_csv("comparison_results.txt", sep="\t", index=False)
print("\nFinal Results:\n", result_df
