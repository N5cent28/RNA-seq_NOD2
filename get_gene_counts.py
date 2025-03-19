#!/usr/bin/env python3

import pandas as pd
import re

# Genes of interest based on NOD2_expression.ipynb
genes_of_interest = [
    'nod2', 'actb1', 'actb2', 'b2m', 'tbp', 'gapdh', 'rpl13a', 'il1b', 'tnfa', 'il6',
    'map1lc3a', 'map1lc3b', 'sqstm1', 'hif1aa', 'hif1ab'
]

def extract_gene_counts(filename):
    print(f"Extracting counts from {filename}...")
    
    # Dictionary to store gene counts and lengths
    gene_data = {gene: {'counts': 0, 'length': 0} for gene in genes_of_interest}
    
    with open(filename, 'r') as f:
        # Skip the header line
        next(f)
        
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 7:  # Ensure line has enough fields
                continue
            
            gene_id = parts[0].lower()
            
            # Check if this line is for a gene we're interested in
            for gene in genes_of_interest:
                if gene == gene_id:
                    try:
                        counts = int(parts[6])
                        # Calculate length from the coordinates in the file
                        # Format is like: NC_007112.7;NC_007112.7    6642;6892    6760;6955    -;-    635
                        lengths = parts[5]
                        
                        gene_data[gene]['counts'] = counts
                        gene_data[gene]['length'] = int(lengths)
                        break
                    except (ValueError, IndexError):
                        print(f"Error parsing line: {line}")
    
    return gene_data


# Check if three samples are present
try:
    srr19627925_data = extract_gene_counts('RNAseq1/gene_counts.txt')
    srr19627924_data = extract_gene_counts('SRR19627924_gene_counts.txt')
    srr19627923_data = extract_gene_counts('SRR19627923_gene_counts.txt')
    has_three_samples = True
except FileNotFoundError:
    print("Third replicate file not found in RNAseq1/gene_counts.txt")
    try:
        srr19627925_data = extract_gene_counts('SRR19627925_gene_counts.txt')
        has_three_samples = True
    except FileNotFoundError:
        print("Third replicate file not found as SRR19627925_gene_counts.txt")
        has_three_samples = False

# Create a comparison table
results = []
for gene in genes_of_interest:
    gene_entry = {
        'Gene': gene,
        'Length': srr19627924_data[gene]['length'],
        'SRR19627924_counts': srr19627924_data[gene]['counts'],
        'SRR19627923_counts': srr19627923_data[gene]['counts']
    }
    
    if has_three_samples:
        gene_entry['SRR19627925_counts'] = srr19627925_data[gene]['counts']
    
    results.append(gene_entry)

# Create dataframe and save to CSV
df = pd.DataFrame(results)
df.to_csv('gene_counts_comparison.csv', index=False)

# Print the results
print("\nGene Counts Comparison:")
print("=" * 80)
print(df.to_string(index=False))
print("\nResults saved to gene_counts_comparison.csv") 