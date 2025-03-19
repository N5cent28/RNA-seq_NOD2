#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the gene counts from CSV
counts_df = pd.read_csv('gene_counts_comparison.csv')

# Get total mapped reads for each sample from summary files
def get_total_mapped_reads(summary_file):
    try:
        with open(summary_file, 'r') as f:
            for line in f:
                if 'Assigned' in line:
                    # Format: Assigned  nnnnnnnn  xx.xx%
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        return int(parts[1])
        return None
    except FileNotFoundError:
        print(f"Warning: Summary file {summary_file} not found")
        return None

# Try to get total mapped reads from summary files
srr19627924_mapped = get_total_mapped_reads('SRR19627924_gene_counts.txt.summary')
srr19627923_mapped = get_total_mapped_reads('SRR19627923_gene_counts.txt.summary')
srr19627925_mapped = get_total_mapped_reads('RNAseq1/gene_counts.txt.summary')

# If summary files aren't available, use the values from NOD2_expression.ipynb
if not srr19627924_mapped:
    srr19627924_mapped = 29743387  # From notebook
if not srr19627923_mapped:
    srr19627923_mapped = 29743387  # Using same value as first sample if not found
if not srr19627925_mapped:
    srr19627925_mapped = 29743387  # Using same value as first sample if not found

print(f"Total mapped reads:")
print(f"SRR19627924: {srr19627924_mapped}")
print(f"SRR19627923: {srr19627923_mapped}")
print(f"SRR19627925: {srr19627925_mapped}")

# Function to calculate FPKM
def calculate_fpkm(counts, length, total_mapped):
    return (counts * 1e9) / (length * total_mapped)

# Add FPKM values for each sample
counts_df['SRR19627924_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['SRR19627924_counts'], row['Length'], srr19627924_mapped),
    axis=1
)

counts_df['SRR19627923_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['SRR19627923_counts'], row['Length'], srr19627923_mapped),
    axis=1
)

counts_df['SRR19627925_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['SRR19627925_counts'], row['Length'], srr19627925_mapped),
    axis=1
)

# Calculate mean FPKM and standard deviation across replicates
counts_df['Mean_FPKM'] = counts_df[['SRR19627924_FPKM', 'SRR19627923_FPKM', 'SRR19627925_FPKM']].mean(axis=1)
counts_df['StdDev_FPKM'] = counts_df[['SRR19627924_FPKM', 'SRR19627923_FPKM', 'SRR19627925_FPKM']].std(axis=1)

# Sort by mean FPKM
counts_df = counts_df.sort_values('Mean_FPKM', ascending=False)

# Add expression level category
def get_expression_level(fpkm):
    if fpkm > 100:
        return "Very High"
    elif fpkm > 10:
        return "High"
    elif fpkm > 1:
        return "Moderate"
    else:
        return "Low"

counts_df['Expression_Level'] = counts_df['Mean_FPKM'].apply(get_expression_level)

# Round FPKM values for better readability
fpkm_columns = ['SRR19627924_FPKM', 'SRR19627923_FPKM', 'SRR19627925_FPKM', 'Mean_FPKM', 'StdDev_FPKM']
counts_df[fpkm_columns] = counts_df[fpkm_columns].round(2)

# Create a clean version of the table for display
display_df = counts_df[['Gene', 'Length', 
                       'SRR19627924_counts', 'SRR19627923_counts', 'SRR19627925_counts',
                       'SRR19627924_FPKM', 'SRR19627923_FPKM', 'SRR19627925_FPKM', 
                       'Mean_FPKM', 'StdDev_FPKM', 'Expression_Level']]

# Rename columns for cleaner display
display_df.columns = ['Gene', 'Length (bp)', 
                     'Counts Rep1', 'Counts Rep2', 'Counts Rep3',
                     'FPKM Rep1', 'FPKM Rep2', 'FPKM Rep3', 
                     'Mean FPKM', 'StdDev', 'Expression Level']

# Save to CSV
display_df.to_csv('gene_expression_comparison.csv', index=False)

# Print the table
print("\nGene Expression Comparison Across Replicates:")
print("=" * 100)
print(display_df.to_string(index=False))
print("\nResults saved to gene_expression_comparison.csv")

# Create bar plot of mean FPKM values with error bars
plt.figure(figsize=(14, 8))
plt.bar(counts_df['Gene'], counts_df['Mean_FPKM'], yerr=counts_df['StdDev_FPKM'], 
        capsize=5, color='skyblue', edgecolor='black')
plt.xticks(rotation=45, ha='right')
plt.yscale('log')  # Use log scale for better visibility
plt.ylabel('Mean FPKM (log scale)')
plt.xlabel('Gene')
plt.title('Mean Gene Expression Levels Across Three Replicates')
plt.tight_layout()
plt.savefig('gene_expression_comparison.png', dpi=300)
plt.close()

print("Plot saved to gene_expression_comparison.png")

# Create heatmap to show expression patterns across replicates
plt.figure(figsize=(10, 12))
heatmap_data = counts_df[['Gene', 'SRR19627924_FPKM', 'SRR19627923_FPKM', 'SRR19627925_FPKM']]
heatmap_data = heatmap_data.set_index('Gene')
heatmap_data.columns = ['Replicate 1', 'Replicate 2', 'Replicate 3']

# Apply log transformation for better visualization (add small value to avoid log(0))
log_data = np.log10(heatmap_data + 0.01)

sns.heatmap(log_data, cmap='viridis', annot=heatmap_data.round(1), fmt='.1f', 
            linewidths=0.5, cbar_kws={'label': 'log10(FPKM)'})
plt.title('Gene Expression (FPKM) Across Replicates')
plt.tight_layout()
plt.savefig('gene_expression_heatmap.png', dpi=300)
plt.close()

print("Heatmap saved to gene_expression_heatmap.png") 