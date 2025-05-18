import pandas as pd
from scipy.stats import binom_test, beta

# Read the data from a space- or tab-delimited file
df = pd.read_csv("hotspot_data.tsv", delim_whitespace=True)   #import the data stating proportion of cblaster genomes with antiphage or antidefence (column 1) vs total genomes found by cblaster (column 2)
# Background expectation
p_null = 0.05

# Calculate observed proportion and fold enrichment
df["Proportion"] = df["Genomes_with_system"] / df["Total_Genomes"]
df["Fold_Enrichment"] = df["Proportion"] / p_null

# Perform one-sided binomial test for enrichment (p > 0.05)
df["p_value"] = df.apply(
    lambda row: binom_test(row["Genomes_with_system"], row["Total_Genomes"], p=p_null, alternative='greater'),
    axis=1
)

# Calculate Bayesian posterior probability that true proportion > 0.05
df["Bayesian_P_p_gt_0.05"] = df.apply(
    lambda row: 1 - beta.cdf(p_null, row["Genomes_with_system"] + 1, row["Total_Genomes"] - row["Genomes_with_system"] + 1),
    axis=1
)

# Convert Bayesian probabilities to percentage
df["Bayesian_P_p_gt_0.05"] *= 100

# Print result
print(df)
