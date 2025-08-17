import pandas as pd
from scipy.stats import binom_test, beta

# Read the data from a space- or tab-delimited file
df = pd.read_csv("hotspot_data.tsv", delim_whitespace=True)
# (columns expected: Genomes_with_system, Total_Genomes)

# Background expectation
p_null = 0.05

# Calculate observed proportion and fold enrichment
df["Proportion"] = df["Genomes_with_system"] / df["Total_Genomes"]
df["Fold_Enrichment"] = df["Proportion"] / p_null

# One-sided binomial test for enrichment (p > 0.05)
df["p_value"] = df.apply(
    lambda row: binom_test(
        row["Genomes_with_system"],
        row["Total_Genomes"],
        p=p_null,
        alternative="greater"
    ),
    axis=1
)

# --- Bonferroni correction ---
m = df.shape[0]  # number of tests (hotspots)
alpha_bonf = 0.05 / m
df["p_value_Bonferroni"] = (df["p_value"] * m).clip(upper=1.0)
df["Significant_Bonferroni_(adj_alpha=0.05/m)"] = df["p_value_Bonferroni"] < 0.05
# (equivalently: raw p < alpha_bonf)
df["alpha_Bonferroni"] = alpha_bonf  # same value repeated, for reference

# Bayesian posterior probability that true proportion > 0.05
df["Bayesian_P_p_gt_0.05"] = df.apply(
    lambda row: 1 - beta.cdf(
        p_null,
        row["Genomes_with_system"] + 1,
        row["Total_Genomes"] - row["Genomes_with_system"] + 1
    ),
    axis=1
) * 100  # as %

# Print / save
print(df)
# df.to_csv("hotspot_results_with_bonferroni.tsv", sep="\t", index=False)
