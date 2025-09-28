# Analysis Workflow Documentation

## Overview

This document describes the computational analysis workflow for the COL4A5 splice modulation study in X-linked Alport syndrome kidney organoids.

## Analysis Pipeline

### 1. Bulk RNA-seq Analysis

**Input:** Raw FASTQ files → BAM files (alignment)
**Script:** `scripts/bulk-rnaseq and proteomics integrative analysis.R`
**Notebook:** `bulk_rnaseq/DEG_analysis.Rmd`

#### Steps:
1. **Quality Control**
   - FastQC for raw data quality assessment
   - Trimming and filtering if needed

2. **Alignment and Quantification**
   - STAR alignment to human reference genome
   - FeatureCounts for gene expression quantification

3. **Differential Expression Analysis**
   - DESeq2 for statistical analysis
   - Comparisons:
     - Time course: d22 vs d32 vs d42
     - Genotypes: MT1/MT2 vs WT
     - Treatment: PMO vs SC

4. **Functional Enrichment**
   - GO (Biological Process, Cellular Component, Molecular Function)
   - KEGG pathway analysis
   - Reactome pathway analysis

### 2. Single-cell RNA-seq Analysis

**Input:** 10X Genomics data
**Scripts:** 
- `scripts/singleCell_early_versus_late.R`
- `scripts/singleCell_mutant_versus_control.R`

#### Steps:
1. **Preprocessing**
   - Cell filtering (nFeature_RNA, nCount_RNA, percent.mt)
   - Normalization and scaling
   - Dimensionality reduction (PCA, UMAP)

2. **Clustering and Annotation**
   - FindNeighbors and FindClusters
   - Cell type annotation using reference datasets
   - Marker gene identification

3. **Differential Expression**
   - FindMarkers for cluster-specific genes
   - FindAllMarkers for comprehensive analysis

### 3. Proteomics Analysis

**Input:** Mass spectrometry data (peptide.txt files)
**Script:** `scripts/msqrob2.R`

#### Steps:
1. **Data Preprocessing**
   - Peptide filtering and quality control
   - Missing value imputation
   - Log2 transformation

2. **Protein Quantification**
   - Peptide-to-protein summarization
   - Normalization (median centering)

3. **Statistical Analysis**
   - msqrob2 for differential protein expression
   - Multiple testing correction (FDR)

4. **Functional Annotation**
   - Matrix protein classification (NABA database)
   - Membrane protein identification
   - Pathway enrichment analysis

### 4. Alternative Splicing Analysis

**Input:** BAM files from bulk RNA-seq
**Script:** `scripts/run_FRASER.R`

#### Steps:
1. **Junction Counting**
   - Extract splice junctions from BAM files
   - Count junction reads per sample

2. **FRASER Analysis**
   - Fit negative binomial model
   - Identify differential splicing events
   - Multiple testing correction

3. **Results Interpretation**
   - Splicing event annotation
   - Functional impact assessment

## File Organization

### Input Data Structure
```
data/
├── bulk_rnaseq/
│   ├── bam_files/          # Aligned BAM files
│   ├── fastq_files/        # Raw sequencing data
│   └── metadata/           # Sample information
├── singlecell/
│   ├── dataset1/           # Early vs Late comparison
│   └── dataset2/           # Mutant vs Control comparison
└── proteomics/
    ├── peptide_files/      # Mass spec peptide data
    └── protein_databases/  # Reference databases
```

### Output Structure
```
results/
├── bulk_rnaseq/
│   ├── deg_results/        # Differential expression results
│   ├── enrichment/         # GO/KEGG enrichment results
│   └── visualizations/     # Plots and figures
├── singlecell/
│   ├── clustering/         # Cell type clustering results
│   ├── markers/           # Cell type markers
│   └── differential/      # scRNA-seq DE results
├── proteomics/
│   ├── protein_quant/     # Protein quantification
│   ├── differential/     # Differential protein expression
│   └── functional/       # Functional annotation
└── splicing/
    ├── junction_counts/   # Splicing junction data
    └── fraser_results/   # FRASER analysis results
```

## Quality Control Metrics

### Bulk RNA-seq
- Total reads per sample
- Mapping rate
- Duplication rate
- Gene detection rate

### Single-cell RNA-seq
- Cells per sample
- Genes per cell
- Mitochondrial gene percentage
- Doublet rate

### Proteomics
- Peptide identification rate
- Protein coverage
- Missing value percentage
- CV (coefficient of variation)

## Statistical Parameters

### Differential Expression
- **Bulk RNA-seq:** DESeq2, FDR < 0.05, |log2FC| > 1
- **Single-cell:** FindMarkers, p_val_adj < 0.05
- **Proteomics:** msqrob2, FDR < 0.05, |log2FC| > 1

### Enrichment Analysis
- **GO:** FDR < 0.05, minimum 5 genes
- **KEGG:** FDR < 0.05, minimum 3 genes
- **Reactome:** FDR < 0.05, minimum 3 genes

## Reproducibility

### Environment
- R version: 4.2.0+
- Bioconductor: 3.15+
- All package versions specified in `requirements.txt`

### Random Seeds
- Set seed for reproducible results
- Document seed values in scripts

### Version Control
- Git for code versioning
- Zenodo for data archiving
- DOI for dataset citation

## Troubleshooting

### Common Issues
1. **Memory requirements:** Large datasets may require high-memory systems
2. **Package conflicts:** Use renv or conda for environment management
3. **File paths:** Ensure relative paths are correct
4. **Missing dependencies:** Check all required packages are installed

### Performance Optimization
- Use parallel processing where possible
- Consider using HPC clusters for large datasets
- Optimize memory usage for single-cell analysis
