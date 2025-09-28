# Therapeutic splice modulation of COL4A5 reinstates collagen IV assembly in an organoid model of X-linked Alport syndrome

**Authors:** Hassan Saei, Bruno Estebe, Nicolas Gaudin, Mahsa Esmailpour, Julie Haure, Olivier Gribouval, Christelle Arrondel, Vincent Moriniere, Pinyuan Tian, Rachel Lennon, Corinne Antignac, Geraldine Mollet*, Guillaume Dorval*

*These authors equally supervised the work.

This repository contains the computational analysis code and datasets for a study investigating therapeutic splice modulation of COL4A5 in X-linked Alport syndrome using kidney organoid models. The analysis includes bulk RNA-seq, single-cell RNA-seq, and proteomics data.

## Data Availability

- **Complete analysis datasets** are available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16585561.svg)](https://doi.org/10.5281/zenodo.16585561)
- **Raw and processed files** are deposited in NCBI GEO:
  - [GSE281080](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281080) - Bulk RNA-seq data
  - [GSE181081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181081) - Single-cell RNA-seq data

- **Preprint:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.10.658776v1)

## Repository Structure

```
├── scripts/                          # Analysis scripts
│   ├── bulk-rnaseq and proteomics integrative analysis.R
│   ├── msqrob2.R                    # Proteomics analysis (msqrob2)
│   ├── run_FRASER.R                 # Alternative splicing analysis
│   ├── singleCell_early_versus_late.R
│   ├── singleCell_mutant_versus_control.R
│   └── sample_table.tsv             # Sample metadata
├── bulk_rnaseq/                      # Bulk RNA-seq results
│   ├── DEG_analysis.Rmd             # Main analysis notebook
│   ├── DEG_analysis.pdf             # Analysis report
│   └── Listes_Res_MultiTests/       # Differential expression results
└── README.md
```

## Quick Start
### Prerequisites

```r
# Install required R packages
install.packages(c("shiny", "dplyr", "ggplot2"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("msqrob2", "QFeatures", "limma", "clusterProfiler"))
```

### Interactive Analysis
#### 1. Bulk RNA-seq Analysis

```r
# Run locally
library(bulkAnalyseR)
library(shiny)
runApp("shinyApp_GSE281080/")
```
[Access online](https://hassansaei.shinyapps.io/shiny_gse281080/)


#### 2. Single-cell Analysis
```r
# Dataset 1: Early vs Late comparison
runApp("shinyApp_dataset1/")

# Dataset 2: Mutant vs Control comparison  
runApp("shinyApp_dataset2/")
```
[Access online dataset 1](https://hassansaei.shinyapps.io/shinyapp_dataset1/)

[Access online dataset 2](https://hassansaei.shinyapps.io/shinyapp_dataset2/)

## Analysis Workflow

### 1. Bulk RNA-seq Analysis
- **Script:** `scripts/bulk-rnaseq and proteomics integrative analysis.R`
- **Notebook:** `bulk_rnaseq/DEG_analysis.Rmd`
- **Results:** Differential expression analysis comparing:
  - Time points: Day 22 vs Day 32 vs Day 42
  - Conditions: MT1/MT2 vs WT controls

### 2. Single-cell RNA-seq Analysis
- **Early vs Late:** `scripts/singleCell_early_versus_late.R`
- **Mutant vs Control:** `scripts/singleCell_mutant_versus_control.R`
- **Interactive apps:** Available as Shiny applications

### 3. Proteomics Analysis
- **Script:** `scripts/msqrob2.R`
- **Method:** msqrob2 for differential protein expression
- **Focus:** Extracellular matrix and membrane proteins

### 4. Alternative Splicing Analysis
- **Script:** `scripts/run_FRASER.R`
- **Results:** Available in Zenodo (Fraser_out folder)

## Experimental Design

### Sample Groups
- **Time course:** Day 22 (early), Day 32 (mid), Day 42 (late)
- **Genotypes:** 
  - WT: Wild-type controls
  - MT1: Mutation 1 (PMO vs SC treatment)
  - MT2: Mutation 2 (PMO vs SC treatment)

### Data Types
- **Bulk RNA-seq:** Illumina NovaSeq 6000
- **Single-cell RNA-seq:** 10X Genomics
- **Proteomics:** Mass spectrometry (msqrob2 analysis)

## Key Results
- Differential expression analysis across time points and genotypes
- Functional enrichment analysis (GO, KEGG, Reactome)
- Alternative splicing analysis using FRASER
- Interactive visualization tools for data exploration

## Citation
If you use any code or workflows from this repository, please cite:

```bibtex
@article{saei2025therapeutic,
  title={Therapeutic splice modulation of COL4A5 reinstates collagen IV assembly in an organoid model of X-linked Alport syndrome},
  author={Saei, Hassan and Estebe, Bruno and Gaudin, Nicolas and others},
  journal={[bioRxiv]},
  year={2025},
  doi={[https://doi.org/10.1101/2025.06.10.658776]}
}
```

## Contact
- **Email:** hassan.saeiahan@gmail.com
- **Repository:** [GitHub](https://github.com/[username]/Saei_JCI_Insight_2025)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

