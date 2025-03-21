#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the gene counts from CSV
try:
    counts_df = pd.read_csv('nlr_gene_counts_comparison.csv')
except FileNotFoundError:
    print("Error: nlr_gene_counts_comparison.csv not found.")
    print("Please run get_nlr_counts.py first.")
    exit(1)

# Get total mapped reads for each sample
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
run1_mapped = get_total_mapped_reads('NormoxiaRun1/gene_counts.txt.summary')
run2_mapped = get_total_mapped_reads('NormoxiaRun2/SRR19627924_gene_counts.txt.summary')
run3_mapped = get_total_mapped_reads('NormoxiaRun3/SRR19627923_gene_counts.txt.summary')

# If summary files aren't available, use default values
default_mapped = 29743387  # Default based on previous analysis
if not run1_mapped:
    run1_mapped = default_mapped
if not run2_mapped:
    run2_mapped = default_mapped
if not run3_mapped:
    run3_mapped = default_mapped

print(f"Total mapped reads:")
print(f"NormoxiaRun1: {run1_mapped}")
print(f"NormoxiaRun2: {run2_mapped}")
print(f"NormoxiaRun3: {run3_mapped}")

# Function to calculate FPKM
def calculate_fpkm(counts, length, total_mapped):
    # Handle cases where length is 0 or very small
    if length <= 10:  # Set a minimum length to avoid division by very small numbers
        return 0
    return (counts * 1e9) / (length * total_mapped)

# Add FPKM values for each sample
counts_df['NormoxiaRun1_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['NormoxiaRun1_counts'], row['Length'], run1_mapped),
    axis=1
)

counts_df['NormoxiaRun2_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['NormoxiaRun2_counts'], row['Length'], run2_mapped),
    axis=1
)

counts_df['NormoxiaRun3_FPKM'] = counts_df.apply(
    lambda row: calculate_fpkm(row['NormoxiaRun3_counts'], row['Length'], run3_mapped),
    axis=1
)

# Calculate mean FPKM and standard deviation across replicates
counts_df['Mean_FPKM'] = counts_df[['NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM']].mean(axis=1)
counts_df['StdDev_FPKM'] = counts_df[['NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM']].std(axis=1)

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
    elif fpkm > 0.1:
        return "Low"
    else:
        return "Very Low"

counts_df['Expression_Level'] = counts_df['Mean_FPKM'].apply(get_expression_level)

# Round FPKM values for better readability
fpkm_columns = ['NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM', 'Mean_FPKM', 'StdDev_FPKM']
counts_df[fpkm_columns] = counts_df[fpkm_columns].round(2)

# Create a clean version of the table for display
display_df = counts_df[['Gene', 'Length', 
                       'NormoxiaRun1_counts', 'NormoxiaRun2_counts', 'NormoxiaRun3_counts',
                       'NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM', 
                       'Mean_FPKM', 'StdDev_FPKM', 'Expression_Level']]

# Rename columns for cleaner display
display_df.columns = ['Gene', 'Length (bp)', 
                     'Counts Run1', 'Counts Run2', 'Counts Run3',
                     'FPKM Run1', 'FPKM Run2', 'FPKM Run3', 
                     'Mean FPKM', 'StdDev', 'Expression Level']

# Save to CSV
display_df.to_csv('nlr_gene_expression_comparison.csv', index=False)

# Print the table
print("\nNLR Gene Expression Comparison Across Replicates:")
print("=" * 100)
print(display_df.to_string(index=False))
print("\nResults saved to nlr_gene_expression_comparison.csv")

# Create bar plot of mean FPKM values with error bars for NLR genes only
plt.figure(figsize=(14, 8))

# Filter out housekeeping genes
nlr_only_df = counts_df[~counts_df['Gene'].isin(['actb1', 'actb2', 'b2m', 'tbp'])]

# Create the bar plot
plt.bar(nlr_only_df['Gene'], nlr_only_df['Mean_FPKM'], yerr=nlr_only_df['StdDev_FPKM'], 
        capsize=5, color='skyblue', edgecolor='black')
plt.xticks(rotation=45, ha='right')
plt.yscale('log')  # Use log scale for better visibility
plt.ylabel('Mean FPKM (log scale)')
plt.xlabel('Gene')
plt.title('NLR Gene Expression Levels Across Three Normoxia Replicates')
plt.tight_layout()
plt.savefig('nlr_gene_expression_comparison.png', dpi=300)
plt.close()

print("Plot saved to nlr_gene_expression_comparison.png")

# Create heatmap to show expression patterns across replicates
plt.figure(figsize=(12, 10))

# Prepare data for heatmap
heatmap_data = counts_df[['Gene', 'NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM']]
heatmap_data = heatmap_data.set_index('Gene')
heatmap_data.columns = ['Normoxia Run 1', 'Normoxia Run 2', 'Normoxia Run 3']

# Apply log transformation for better visualization (add small value to avoid log(0))
log_data = np.log10(heatmap_data + 0.01)

# Create clustered heatmap
sns.clustermap(log_data, cmap='viridis', annot=heatmap_data.round(2), fmt='.2f', 
               linewidths=0.5, figsize=(12, 10), 
               cbar_kws={'label': 'log10(FPKM + 0.01)'})
plt.savefig('nlr_gene_expression_heatmap.png', dpi=300)
plt.close()

print("Heatmap saved to nlr_gene_expression_heatmap.png")

# Create a comparison bar plot to show expression differences between NLR genes
plt.figure(figsize=(14, 10))

# Reshape data for grouped bar plot
nlr_only_df = counts_df[~counts_df['Gene'].isin(['actb1', 'actb2', 'b2m', 'tbp'])]
melted_df = pd.melt(nlr_only_df, 
                  id_vars=['Gene'], 
                  value_vars=['NormoxiaRun1_FPKM', 'NormoxiaRun2_FPKM', 'NormoxiaRun3_FPKM'],
                  var_name='Sample', value_name='FPKM')

# Rename sample values for better display
melted_df['Sample'] = melted_df['Sample'].replace({
    'NormoxiaRun1_FPKM': 'Normoxia Run 1',
    'NormoxiaRun2_FPKM': 'Normoxia Run 2',
    'NormoxiaRun3_FPKM': 'Normoxia Run 3'
})

# Create grouped bar plot
sns.barplot(x='Gene', y='FPKM', hue='Sample', data=melted_df)
plt.xticks(rotation=45, ha='right')
plt.title('NLR Gene Expression Across Normoxia Replicates')
plt.tight_layout()
plt.savefig('nlr_gene_expression_comparison_by_sample.png', dpi=300)
plt.close()

print("Comparison plot saved to nlr_gene_expression_comparison_by_sample.png") 