#!/usr/bin/env python3

import pandas as pd
import re

# NLR genes of interest based on user request
nlr_genes_of_interest = [
    'nod1',  # Uniprot: X1WGQ4, replaced with A0A8M1RMD1
    'nod2',  # Uniprot: F8W3K2
    'nlrp1',  # Uniprot: A0A386CAB9
    'si:ch73-233m11.2',  # Uniprot: F1QZJ8
    'si:ch211-214c20.1',  # Uniprot: E9QC36, removed from Uniprot
    'b3dk63',  # ENSDARG00000095773 (Uniprot: B3DK63)
    'nlrc3l1',  # Uniprot: E7FBD8
    'nlrc6',  # Uniprot: A5PF24
    'si:ch211-195h23.3',  # Uniprot: B0V1H3
    
    # Adding housekeeping genes for comparison
    'actb1', 'actb2', 'b2m', 'tbp'
]

def extract_gene_counts(filename):
    print(f"Extracting counts from {filename}...")
    
    # Dictionary to store gene counts and lengths
    gene_data = {gene: {'counts': 0, 'length': 0} for gene in nlr_genes_of_interest}
    
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
            for gene in nlr_genes_of_interest:
                if gene == gene_id:
                    try:
                        counts = int(parts[6])
                        # Length is the sum of all exon lengths
                        lengths = parts[5]
                        
                        gene_data[gene]['counts'] = counts
                        gene_data[gene]['length'] = int(lengths)
                        break
                    except (ValueError, IndexError):
                        print(f"Error parsing line: {line}")
    
    return gene_data


# Extract gene counts from all three samples
try:
    run1_data = extract_gene_counts('NormoxiaRun1/gene_counts.txt')
    run2_data = extract_gene_counts('NormoxiaRun2/SRR19627924_gene_counts.txt')
    run3_data = extract_gene_counts('NormoxiaRun3/SRR19627923_gene_counts.txt')
    
    # Create a comparison table
    results = []
    for gene in nlr_genes_of_interest:
        gene_entry = {
            'Gene': gene,
            'Length': run1_data[gene]['length'],
            'NormoxiaRun1_counts': run1_data[gene]['counts'],
            'NormoxiaRun2_counts': run2_data[gene]['counts'],
            'NormoxiaRun3_counts': run3_data[gene]['counts']
        }
        
        results.append(gene_entry)

    # Create dataframe and save to CSV
    df = pd.DataFrame(results)
    df.to_csv('nlr_gene_counts_comparison.csv', index=False)

    # Print the results
    print("\nNLR Gene Counts Comparison:")
    print("=" * 80)
    print(df.to_string(index=False))
    print("\nResults saved to nlr_gene_counts_comparison.csv")
    
except FileNotFoundError as e:
    print(f"Error: {e}")
    print("Please make sure all three gene count files are available") 