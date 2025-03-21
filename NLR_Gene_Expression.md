# NLR Gene Expression Analysis in Zebrafish

<head>
<style>
.code-container {
  position: relative;
  margin-bottom: 1em;
}
.code-block {
  background-color: #f5f5f5;
  border-radius: 4px;
  padding: 1em;
  overflow-x: auto;
  margin: 0;
}
.copy-button {
  position: absolute;
  top: 5px;
  right: 5px;
  padding: 4px 8px;
  background-color: #f1f1f1;
  border: 1px solid #ccc;
  border-radius: 4px;
  font-size: 12px;
  cursor: pointer;
  opacity: 0.7;
  transition: opacity 0.3s;
}
.copy-button:hover {
  opacity: 1;
}
.copy-success {
  color: green;
}
</style>
<script>
document.addEventListener('DOMContentLoaded', function() {
  const codeBlocks = document.querySelectorAll('pre code');
  
  codeBlocks.forEach((codeBlock, index) => {
    const container = document.createElement('div');
    container.className = 'code-container';
    const parent = codeBlock.parentNode;
    parent.parentNode.insertBefore(container, parent);
    container.appendChild(parent);
    
    const copyButton = document.createElement('button');
    copyButton.className = 'copy-button';
    copyButton.textContent = 'Copy';
    copyButton.addEventListener('click', function() {
      const textToCopy = codeBlock.textContent;
      navigator.clipboard.writeText(textToCopy)
        .then(() => {
          copyButton.textContent = 'Copied!';
          copyButton.classList.add('copy-success');
          setTimeout(() => {
            copyButton.textContent = 'Copy';
            copyButton.classList.remove('copy-success');
          }, 2000);
        })
        .catch(err => {
          console.error('Failed to copy text: ', err);
        });
    });
    
    container.appendChild(copyButton);
  });
});
</script>
</head>

## Overview

This document describes the analysis of the expression levels of various NLR (Nucleotide-binding domain and Leucine-rich Repeat containing) genes in zebrafish RNA-seq data from three normoxia conditions. The analysis focuses on target and counter-selection NLR genes of interest.

## Genes Analyzed

### Target & Counter-selection Genes
- NOD1 (Uniprot: X1WGQ4, replaced with A0A8M1RMD1)
- NOD2a (Uniprot: F8W3K2)
- NLRP1 (Uniprot: A0A386CAB9)
- si:ch73-233m11.2 (Uniprot: F1QZJ8)
- si:ch211-214c20.1 (Uniprot: E9QC36, removed from Uniprot)
- ENSDARG00000095773 (Uniprot: B3DK63)

### Counter-selection Only Genes
- NLRC3L1 (Uniprot: E7FBD8)
- NLRC6 (Uniprot: A5PF24)
- si:ch211-195h23.3 (Uniprot: B0V1H3)

### Housekeeping Genes (for comparison)
- actb1, actb2 (β-actin)
- b2m (β-2-microglobulin)
- tbp (TATA-box binding protein)

## Methods

### Data Sources
Expression data was extracted from RNA-seq datasets of zebrafish under normoxic conditions:
- Normoxia Run 1 (SRR19627925)
- Normoxia Run 2 (SRR19627924)
- Normoxia Run 3 (SRR19627923)

### Analysis Workflow

1. **Gene Count Extraction**:
   - Raw gene counts were extracted from featureCounts output files for each replicate.
   - The counts were saved in a CSV file for further analysis.

2. **FPKM Calculation**:
   - FPKM (Fragments Per Kilobase Million) values were calculated using the formula:
     ```
     FPKM = (Counts * 10^9) / (Gene Length * Total Mapped Reads)
     ```
   - FPKM values were calculated for each gene in each replicate.
   - Mean FPKM and standard deviation were calculated across the three replicates.

3. **Expression Level Classification**:
   - Genes were classified into expression categories based on their mean FPKM values:
     - Very High: FPKM > 100
     - High: 10 < FPKM ≤ 100
     - Moderate: 1 < FPKM ≤ 10
     - Low: 0.1 < FPKM ≤ 1
     - Very Low: FPKM ≤ 0.1

4. **Visualization**:
   - Bar plots showing mean FPKM values with error bars
   - Heatmap showing expression patterns across replicates
   - Grouped bar plot comparing expression levels between NLR genes across replicates

## Running the Analysis

The analysis consists of two main scripts:

1. **get_nlr_counts.py**: Extracts gene counts from the featureCounts output files
   ```bash
   python get_nlr_counts.py
   ```

2. **calculate_nlr_fpkm.py**: Calculates FPKM values and generates visualizations
   ```bash
   python calculate_nlr_fpkm.py
   ```

## Results

The analysis produces the following output files:
- **nlr_gene_counts_comparison.csv**: Raw gene counts across replicates
- **nlr_gene_expression_comparison.csv**: FPKM values and expression statistics
- **nlr_gene_expression_comparison.png**: Bar plot of mean FPKM values
- **nlr_gene_expression_heatmap.png**: Clustered heatmap of expression patterns
- **nlr_gene_expression_comparison_by_sample.png**: Grouped bar plot comparing expression across replicates

## Interpretation

This analysis provides insights into the expression levels of NLR genes in zebrafish under normoxic conditions. The results can be used to:

1. Identify which NLR genes are expressed at detectable levels
2. Compare expression levels between different NLR genes
3. Assess the consistency of expression across replicates
4. Select genes for further functional studies based on their expression patterns

The expression levels of NLR genes can inform target selection for downstream applications such as CRISPR knockout, antibody development, or functional characterization studies.

## Limitations and Considerations

- Gene identification in the gene counts files is based on gene symbols, which may not capture all transcript variants
- Some genes might be annotated differently in the genome than their UniProt identifiers
- FPKM values provide relative expression levels but should be validated with additional methods (qPCR, protein detection)
- Expression under normoxic conditions may differ from expression under stimulated conditions (e.g., immune challenge) 