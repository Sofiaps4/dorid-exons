import pandas as pd
import matplotlib.pyplot as plt

# Load the TSV table
df = pd.read_csv("seq_lengths.tsv", sep="\t", index_col=0)

# Extract gene lengths
gene_lengths = df.loc["MeanLength"]

# Drop the 'MeanLength' row to keep only species
species_df = df.drop("MeanLength")

# Calculate percentage of gene length for each species
percentage_df = (species_df.div(gene_lengths, axis=1)) * 100

# Flatten the values into a single series
all_percentages = percentage_df.values.flatten()

# Drop NaN or zero values
all_percentages = all_percentages[~pd.isna(all_percentages)]
all_percentages = all_percentages[all_percentages > 0]

# Count total valid entries
total_genes = len(all_percentages)

# Bin the values into ranges: 0–10%, 10–20%, ..., 90–100%
bins = range(0, 110, 10)
binned = pd.cut(all_percentages, bins=bins, right=False)
binned_counts = binned.value_counts().sort_index()

# Convert counts to percentages
binned_percentages = (binned_counts / total_genes) * 100

# Custom bar labels
labels = [f"{i}–{i+10}%" for i in range(0, 100, 10)]

# Plot the bar chart and save it
plt.figure(figsize=(6, 6))
plt.bar(labels, binned_percentages, color='#1f4e79', edgecolor='black', width=0.5)
plt.xlabel("Percentage of Gene Length Recovered", fontsize=12)
plt.ylabel("Percentage of Genes (%)", fontsize=12)
plt.title("Distribution of Gene Length Recovery Across Species", fontsize=14)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("barplot_percent_lengths.png", dpi=300)
plt.close()
