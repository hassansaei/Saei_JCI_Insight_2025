## Script to analyze bulk RNA-seq dataset GSE281080
## Comparing organoids at different time points of culture (d22, d32 and d42)
# loading libraries

#BiocManager::install("clusterProfiler")
rm(list=ls())
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(radiant)
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library("readxl")
library(data.table)
library(pheatmap)

# convert ENSG to gene symbol

str(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))
ens2sym_entrez <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86), columns = c("SYMBOL"))

# Read bulk RNA-seq results (d32_vs_d22_f1.2 and d42_vs_d22_f1.2)

DEG_32 <- read.table("../0_NGS2024_7984/Listes_Res_MultiTests/(d32)_vs_(d22)_f1.2_(7369).txt", sep = "\t", header = T)
DEG_42 <- read.table("../0_NGS2024_7984/Listes_Res_MultiTests/(d42)_vs_(d22)_f1.2_(11925).txt", sep = "\t", header = T)
NABA <- read.table("NABA_MATRISOME.v2023.2.txt", sep = "\t", header = T)

# function to clean, annotate and extract matrix genes
process_DEG <- function(df, naba_symbols) {
  df %>%
    filter(!is.na(pval_edgeR), !is.na(pval_Voom), !is.na(pval_DEseq2)) %>%
    mutate(Category = ifelse(Symbol %in% naba_symbols, "Matrix", "No")) %>%
    filter(Category == "Matrix") %>%
    pull(Symbol)
}

# process both datasets
genes_d32 <- process_DEG(DEG_32, NABA$Symbol)
genes_d42 <- process_DEG(DEG_42, NABA$Symbol)

# generating venndiagram comparing differential matrix encoding genes at day 42 and day 32 
library("ggVennDiagram")
list <- list(genes_d32, genes_d42)
p <- ggVennDiagram(list, label_alpha = 0, label_size = 7,
              category.names = c("Day 32", "Day 42")) + ggplot2::scale_fill_gradient(low="blue",high = "yellow") + ggplot2::coord_flip()

# make a folder to save figures
dir.create("Figures")

# save veendiagram as png
png(filename = "Figures/Venn_matrix_genes.png", width = 150, height = 100, res = 300, units = "mm")
plot(p)
dev.off()

# extrcat common genes between d42_vs_22 and d32_vs_d22 comparisons
common_genes <- intersect(genes_d32, genes_d42)
write.table(common_genes, file = "Figures/Common_DE_matrix_encoding_genes.txt", row.names = F, col.names = "Common_genes")


# MDS plot
sampleTable <- read.table("../0_NGS2024_7984/DataNorm/DataNormDESeq2.txt", header = T)
sampleTable <- sampleTable[,c(1:9)]
group <- factor(c(rep("d22", 3), rep("d32", 3), rep("d42", 3)))
colors <- c("darkred", "darkblue", "darkgreen")[group]
limma::plotMDS(sampleTable, col = colors, cex = 0.5)
mds <- plotMDS(sampleTable, plot = FALSE)
plot(mds$x, mds$y, col = colors, pch = 19, cex = 1.5, xlab = "Dimension 1 (92%)", ylab = "Dimension 2 (6%)", main = "")
legend("topleft", legend = colnames(sampleTable), col = colors, pch = 19, cex = 0.8, bty = "n")

png(filename = "Figures/PlotMDS_bulk_rnaseq.png", width = 140, height = 140, res = 300, units = "mm")
plot(mds$x, mds$y, col = colors, pch = 19, cex = 1.5, xlab = "Dimension 1 (92%)", ylab = "Dimension 2 (6%)", main = "")
legend("topleft", legend = colnames(sampleTable), col = colors, pch = 19, cex = 0.8, bty = "n")
dev.off()


# Volcano plot comparing d42vsd22 and d32vsd22
# SETTINGS
deg_file <- "../0_NGS2024_7984/Analyse_DEseq2/(d42)_vs_(d22)/res_DEseq2.0_(d42)_vs_(d22).xlsx"
ratio_col <- "RatiosMoys_(d42)_vs_(d22)"
threshold_pvalue <- 0.05
threshold_lfc <- 1.5

highlight_genes <- c("COL4A1","COL4A2","COL4A3","COL4A4","COL4A5","COL4A6",
                     "LAMA1","LAMB1","LAMB2","LAMC1","AGRN","LAMA5")

# load data
DEG <- read_excel(deg_file) %>%
  mutate(
    logRatio = log2(.data[[ratio_col]]),
    logPval  = -log10(pval_BH)
  )

# intersection of matrix genes (assumes genes_d32, genes_d42 exist)
common <- intersect(genes_d32, genes_d42)

# categorize
DEG <- DEG %>%
  mutate(
    Category = ifelse(Symbol %in% common, "Matrix", "No"),
    significant = ifelse(pval_BH < threshold_pvalue & abs(logRatio) > threshold_lfc, "yes", "no"),
    category = case_when(
      pval_BH < threshold_pvalue & logRatio > threshold_lfc  ~ "Up",
      pval_BH < threshold_pvalue & logRatio < -threshold_lfc ~ "Down",
      TRUE                                                  ~ "Non-sig"
    ),
    highlight = case_when(
      category == "Up"    & Symbol %in% common ~ "Up_highlight",
      category == "Down"  & Symbol %in% common ~ "Down_highlight",
      category == "Non-sig" & Symbol %in% common ~ "Non-sig_highlight",
      category == "Up"    ~ "Up",
      category == "Down"  ~ "Down",
      TRUE                ~ "Non-sig"
    ),
    matrix = case_when(
      highlight == "Up_highlight"   & Symbol %in% highlight_genes ~ "Up_matrix",
      highlight == "Down_highlight" & Symbol %in% highlight_genes ~ "Down_matrix",
      highlight == "Non-sig_highlight" ~ "Non-sig_highlight",
      highlight == "Up"             ~ "Up",
      highlight == "Down"           ~ "Down",
      TRUE                          ~ "Non-sig"
    )
  )

# Optional: Remove invalid padj rows
DEG <- DEG %>% filter(padj >= 3.087506e-300)

# COLORS
color_map <- c(
  "Up_highlight"       = "darkred",
  "Down_highlight"     = "blue",
  "Non-sig_highlight"  = "gray80",
  "Up"                 = "lightcoral",
  "Down"               = "lightblue",
  "Non-sig"            = "gray80"
)


# Draw the volcano plot
volcano <- ggplot(DEG, aes(x= logRatio, y = -log10(pval_BH) , color=category)) +
  geom_point(alpha=0.6, size=1.5) +
  theme_minimal() +
  theme(legend.position="top", axis.line = element_line(color = "black")) +
  scale_color_manual(values=color_map) +
  labs(title="Volcano Plot",
       x="Log2(ratio)",
       y="-Log10 P-value",
       color="Volcano CTX vs CTRL") +
  geom_hline(yintercept=-log10(threshold_pvalue), linetype="dashed", color = "black") +
  geom_vline(xintercept=c(-threshold_lfc, threshold_lfc), linetype="dashed", color = "black")


# Print the plot
volcano <- volcano +
  geom_text_repel(data = subset(DEG, highlight %in% c("Up_highlight", "Down_highlight")),
                  aes(label = Symbol),
                  size = 3,
                  max.overlaps = 20)

print(volcano)


# Enrichment
# SETTINGS 
enrich_file <- "../0_NGS2024_7984/AnalyseGO/GO_All_Venn_(d32)_vs_(d22)_f1.5_(3360).txt"
qval_threshold <- 0.05
top_n_terms <- 3  # number of top terms per ontology

# FUNCTION
process_enrichment <- function(file, q_threshold = 0.05, top_n = 3) {
  enrich <- fread(file, header = TRUE, sep = "\t") %>%
    filter(qvalue <= q_threshold)
  
  # Split by ontology and select top terms
  df_list <- lapply(c("BP", "CC", "MF"), function(onto) {
    enrich %>%
      filter(ONTOLOGY == onto) %>%
      slice_head(n = top_n)
  })
  
  # Combine and add transformed q-values
  df_combined <- bind_rows(df_list) %>%
    mutate(log10_Adjusted = -log10(qvalue))
  
  return(df_combined)
}


df_combined <- process_enrichment(enrich_file, q_threshold = qval_threshold, top_n = top_n_terms)

# Barlpot for enriched terms
p <- ggplot(df_combined, aes(x = log10_Adjusted, y = reorder(Description, log10_Adjusted), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("BP" = "blue3", "CC" = "coral", "MF" = "azure3")) +
  labs(
    title = "",
    x = "-log10(q-value)",
    y = "Enrichments",
    fill = "Databases"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
  )


png(filename = "Figures/Enrichment_d32vsd22.png", width = 200, height = 150, res = 300, units = "mm")
plot(p)
dev.off()

# add gene names to the normalized matrix
sampleTable <- read.table("../0_NGS2024_7984/DataNorm/DataNormDESeq2.txt", header = T)
sampleTable <- sampleTable[,c(1:9)]
sampleTable$GENEID <- rownames(sampleTable)
sampleTable <- merge(sampleTable, ens2sym_entrez, by="GENEID")
sampleTable <- sampleTable %>% 
  distinct(SYMBOL, .keep_all = TRUE)
rownames(sampleTable) <- sampleTable$SYMBOL
sampleTable <- sampleTable[,-c(1)]
sampleTable <- sampleTable[,-10]
sampleTable$Gene.names <- rownames(sampleTable)
sampleTable <- sampleTable[,c("Gene.names", "d22_R1", "d22_R2","d22_R3", "d32_R1", "d32_R2", "d32_R3", "d42_R1", "d42_R2","d42_R3")]
#write.table(sampleTable, file = "../0_NGS2024_7984/DataNorm/NormDataDESeq2.annot.txt", row.names = F)

# extrcat cell type specific DEGs 
df1 <- DEG_d42
df2 <- read.table("MarkerGenes_NGS2024_7386.txt", header = T)
colnames(df2) <- c("p_val", "avg_log2FC","pct.1",      "pct.2",    "p_val_adj",  "cluster",    "Symbol")
df2 <- df2[df2$p_val_adj <= 0.05,]
merged_df <- merge(df1, df2, by = "Symbol")
merged_df <- merged_df[order(merged_df$Symbol, -merged_df$avg_log2FC), ]
merged_df <- merged_df[!duplicated(merged_df$Symbol), ]
merged_df <- merged_df[merged_df$cluster == "Podocytes",]
merged_df$pval_DEseq2 <- gsub(",", ".", merged_df$pval_DEseq2)
merged_df$pval_DEseq2 <- as.numeric(merged_df$pval_DEseq2)

deg_count_by_celltype <- table(merged_df$cluster)
deg_count_df <- as.data.frame(deg_count_by_celltype)

# Sort the data frame by DEG count (descending order)
deg_count_df <- deg_count_df[order(-deg_count_df$Freq), ]

# Plot the vertical bar chart, sorted by DEG count, colored based on count
ggplot(deg_count_df, aes(x = reorder(Var1, -Freq), y = Freq, fill = Freq)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color based on DEG count
  labs(title = "Number of DEGs by Cell Type", x = "Cell Type", y = "DEG Count") +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none") 


p <- ggplot(deg_count_df, aes(x = reorder(Var1, Freq), y = Freq, fill = Freq)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), hjust = -0.3, color = "black", size = 4) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color based on DEG count
  labs(title = "", x = "Cell Types", y = "Number of differentially regulated genes (d42vs22)") +
  theme_minimal() +
  coord_flip() +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black"),  # Add lines to x and y axes
    axis.text = element_text(color = "black", size = 10),  # Make axis values black
    legend.position = "none"  # Hide the legend as colors represent counts
  )


png(filename = "DEG.d42vsd22.celltype.png", width = 280, height = 150, res = 300, units = "mm")
plot(p)
dev.off()


# Setting
rna_files <- list(
  d32 = list(file = "../0_NGS2024_7984/Listes_Res_MultiTests/(d32)_vs_(d22)_f1.2_(7369).txt",
             ratio_col = "RatiosMoys_.d32._vs_.d22."),
  d42 = list(file = "../0_NGS2024_7984/Listes_Res_MultiTests/(d42)_vs_(d32)_f1.2_(8892).txt",
             ratio_col = "RatiosMoys_.d42._vs_.d32.")
)

pro_files <- list(
  d34 = "../Proteomics.perseus.d34_vs_d22_DEPs.txt",
  d38 = "../Proteomics.perseus.d38_vs_d22_DEPs.txt"
)

# functions for data preprocessing > upset plot
process_rna <- function(file, ratio_col) {
  read.table(file, header = TRUE, sep = "\t") %>%
    select(1:3, 5) %>%              # Keep relevant cols
    select(-1) %>%
    mutate(
      log2FC = log2(.data[[ratio_col]]),
      pval_DEseq2 = as.numeric(gsub(",", ".", pval_DEseq2)),
      Pvalue = -log(pval_DEseq2)
    ) %>%
    filter(!is.na(Pvalue), is.finite(Pvalue)) %>%
    select(Symbol = Symbol, Pvalue, LogFC = log2FC)
}

process_proteomics <- function(file) {
  fread(file) %>%
    filter(Significant == "+") %>%
    select(Pvalue = 2, LogFC = 3, SYMBOL = 7) %>%
    select(SYMBOL, Pvalue, LogFC)
}

# RUN
rna_list <- lapply(rna_files, function(x) process_rna(x$file, x$ratio_col))
pro_list <- lapply(pro_files, process_proteomics)

# Example: Compare d32 RNA vs d34 Proteomics
DEG_genes <- rna_list$d42$Symbol
pro_genes <- pro_list$d38$SYMBOL

# Remove duplicates
DEG_genes <- unique(DEG_genes)
pro_genes <- unique(pro_genes)

# Create gene list for UpSet
gene_list <- list(
  RNAseq = DEG_genes,
  Proteomics = pro_genes
)

# Find common genes
common_genes <- intersect(gene_list$RNAseq, gene_list$Proteomics)
length(common_genes)


#setting colors
#this can also be done with hexadecimal
main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")

mb_ratio1 <- c(0.55,0.45)
set_vars <- c("RNAseq", "Proteomics")
text_scale_options1 <- c(1, 1, 1, 1, 0.75, 1)
text_scale_options2 <- c(1.3, 1.3, 1, 1, 2, 0.75)
text_scale_options3 <- c(1.5, 1.25, 1.25, 1, 2, 1.5)

p <- upset(fromList(gene_list), 
      sets = set_vars,
      mb.ratio = mb_ratio1, 
      mainbar.y.label = "Intersection Size", 
      order.by = "freq",
      show.numbers = TRUE,
      point.size = 4, 
      line.size = 1,
      text.scale=text_scale_options3,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col)


png(filename = "Figures/Upset.d38_vs_d32.perseus.DEP.png", width = 150, height = 160 , units = "mm", res = 300)
p
dev.off()


# Enrichment analysis of the shared genes and proteins 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

common_rna_pro <- intersect(gene_list$RNAseq, gene_list$Proteomics)
original_gene_list <- common_rna_pro

gene_sets_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_sets_GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets_GO_CC <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
gene_sets_GO_MF <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
gene_sets_Reactome <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")

gene_sets_hallmark <- gene_sets_hallmark %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_BP <- gene_sets_GO_BP %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_CC <- gene_sets_GO_CC %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_MF <- gene_sets_GO_MF %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_Reactome <- gene_sets_Reactome %>%
  dplyr::select(gs_name, gene_symbol)


# FUNCTION
process_enrichment <- function(gene_list, term2gene, top_n = 6) {
  # Run enrichment
  enrich_res <- enricher(
    gene = gene_list,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # Convert to dataframe and curate GeneRatio
  df <- as.data.frame(enrich_res)
  
  if (nrow(df) > 0) {
    split_values <- strsplit(df$GeneRatio, "/")
    numerator <- sapply(split_values, function(x) as.numeric(x[1]))
    denominator <- sapply(split_values, function(x) as.numeric(x[2]))
    df$GR <- numerator / denominator
    
    # Order and keep top N
    df <- df %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n) %>%
      arrange(Count)
    
    # Keep factor levels for plotting
    df$ID <- factor(df$ID, levels = df$ID)
  }
  
  return(df)
}

# RUN
df_CC <- process_enrichment(original_gene_list, gene_sets_GO_CC)
df_BP <- process_enrichment(original_gene_list, gene_sets_GO_BP)
df_MF <- process_enrichment(original_gene_list, gene_sets_GO_MF)
df_Reactome <- process_enrichment(original_gene_list, gene_sets_Reactome)

df_CC$source <- "GO:CC"
df_BP$source <- "GO:BP"
df_MF$source <- "GO:MF"
df_Reactome$source <- "Reactome"

#df_MF$ID <- gsub("_CONSTITUENT_CONFERRING_TENSILE_STRENGTH", "", df_MF$ID)
combined_df <- bind_rows(df_CC, df_BP, df_MF)
combined_df$ID <- gsub("GOCC_|GOBP_|GOMF_|REACTOME_|KEGG_", "", combined_df$ID)
combined_df$ID <- gsub("_", " ", combined_df$ID)
combined_df$ID  <- tolower(combined_df$ID ) # Lowercase all
combined_df$ID <- paste0(toupper(substr(combined_df$ID , 1, 1)), substr(combined_df$ID , 2, nchar(combined_df$ID))) # only first letter be uppercase
combined_df <- combined_df %>% arrange(p.adjust)
combined_df <- combined_df %>% arrange(Count)
combined_df$ID <- factor(combined_df$ID, levels = combined_df$ID)
combined_df$ID <- factor(combined_df$ID, levels = combined_df$ID[order(combined_df$GR)])

# Save the table
write.table(combined_df, file = "Figures/Enrichment.830.intersect.txt", row.names = F)


png(filename = "Figures/Integrated_Enrichment.830.intersect.png", width = 200, height = 200, units = "mm", res = 300)
ggplot(combined_df, aes(x =GR, y = ID, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("GO:CC" = "navy", "GO:BP" = "lightblue", "GO:MF" = "darkred", "Reactome" = "gray")) +  # Colors for different data frames
  scale_alpha_continuous(range = c(0.3, 1)) +  # Intensity of bar based on p-value
  labs(x = "GeneRatio", y = "", fill = "Databases", alpha = "p-value Intensity (-log10)") +
  theme_bw()  + theme(legend.position = "top", 
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10, colour = 'black'),
                      axis.text.y = element_text(face = "bold", size = 12, colour = "black"))
dev.off()


# Heatmap for enrichment
# RNAseq and proteomics matrices
rna_data <- read.table("../0_NGS2024_7984/DataNorm/NormDataDESeq2.annot.txt", header = T)
rna_data <- rna_data[,c("Gene.names", "d22_R1", "d22_R2", "d22_R3", "d32_R1", "d32_R2","d32_R3",
                        "d42_R1", "d42_R2", "d42_R3")]
colnames(rna_data) <- c("Gene.names", "Early_R1", "Early_R2", "Early_R3", "Mid_R1", "Mid_R2", "Mid_R3",
                        "Late_R1", "Late_R2", "Late_R3")
rna_data <- rna_data %>% 
  distinct(Gene.names, .keep_all = TRUE)
rownames(rna_data) <- rna_data[,1]
rna_data <- rna_data[rowSums(rna_data > 10) > 0, ]

proteomics_data <- read.table("../0_NGS2024_7984/DataNorm/Proteomics.timeseries.merged.txt", sep = "\t", header = T)
colnames(proteomics_data)
proteomics_data <- proteomics_data[,c("Gene.names", "d22_R1", "d22_R2", "d22_R3",
                                      "d32_R1", "d32_R2", 
                                      "d32_R3", "d38_R1", "d38_R2", 
                                      "d38_R3")]
colnames(proteomics_data) <- c("Gene.names", "Early_R1", "Early_R2", "Early_R3", "Mid_R1", "Mid_R2", "Mid_R3",
                               "Late_R1", "Late_R2", "Late_R3")
proteomics_data <- proteomics_data %>% 
  distinct(Gene.names, .keep_all = TRUE)
rownames(proteomics_data) <- proteomics_data[,1]
rna_data <- rna_data[,-1]
proteomics_data <- proteomics_data[,-1]
  

#common_genes <- intersect(rownames(rna_data), rownames(proteomics_data))
#rna_data <- rna_data[rownames(rna_data) %in% common_genes, ]
#proteomics_data <- proteomics_data[rownames(proteomics_data) %in% common_genes, ]

# Split the gene names in df1
df1_split <- strsplit(as.character(combined_df$geneID), "/")
df1_split <- strsplit(as.character(df_combined$geneID), "/")

# Function to check if a gene from df2 is in any row of df1
gene_in_df1 <- function(gene, gene_list) {
  any(sapply(gene_list, function(x) gene %in% x))
}

# Prepare RNA-seq and proteomics to make a single scaled dataframe
# Identify common genes before binding dfs together
common_genes <- intersect(rownames(rna_data), rownames(proteomics_data))

# Subset and reorder the bulk_matrix and proteomics_matrix based on common_genes
bulk_matrix_aligned <- rna_data[match(common_genes, rownames(rna_data)), ]
bulk_matrix_aligned <- bulk_matrix_aligned[,-c(4:6)]
proteomics_matrix_aligned <- proteomics_data[match(common_genes, rownames(proteomics_data)), ]
proteomics_matrix_aligned <- proteomics_matrix_aligned[,-c(7:9)]
# Ensure both matrices now have the same gene order
all(rownames(bulk_matrix_aligned) == rownames(proteomics_matrix_aligned))

# Scale both dfs separately
bulk_matrix_aligned <- as.data.frame(t(apply(bulk_matrix_aligned,1, scale)))
colnames(bulk_matrix_aligned) <- c("Early_R1", "Early_R2", "Early_R3", "Late_R1", "Late_R2", "Late_R3")
proteomics_matrix_aligned <- as.data.frame(t(apply(proteomics_matrix_aligned,1, scale)))
colnames(proteomics_matrix_aligned) <- c("Early_P1", "Early_P2", "Early_P3", "Late_P1", "Late_P2", "Late_P3")
combined_matrix <- cbind(bulk_matrix_aligned, proteomics_matrix_aligned)
combined_matrix <- combined_matrix[-15,]

# Filter df2 based on the presence of genes in df1
## Extracellular matrix structural constituent
df2_filtered <- combined_matrix[sapply(rownames(combined_matrix), gene_in_df1, gene_list = df1_split[[3]]), ]

# add annotation to the heatmap (timepoints and data sources)
annotation_col <- data.frame(Source = c(rep("RNA-seq", ncol(bulk_matrix_aligned)), 
                                        rep("Proteomics", ncol(proteomics_matrix_aligned))), 
                             Timepoints = c(rep("Early", 3), rep("Late", 3), rep("Early", 3), rep("Late", 3)))

rownames(annotation_col) <- colnames(df2_filtered)

ann_colors <- list(
  Source = c("RNA-seq" = "#1b9e77", "Proteomics" = "#d95f02"),
  Timepoints = c("Early" = "#7570b3", "Late" = "#e7298a")
)

p <- pheatmap(df2_filtered, cluster_rows = TRUE, cluster_cols = F, 
         annotation_col = annotation_col, 
         show_rownames = T, 
         color = brewer.pal(n = 9, name = "Reds"), annotation_colors = ann_colors,show_colnames = F,
         main = "Extracellular matrix structural constituent")

# flip the heatmap orientation
p <- pheatmap(t(df2_filtered),                # transpose the matrix
              cluster_rows = FALSE,           # flip clustering accordingly
              cluster_cols = TRUE,            
              annotation_row = annotation_col,
              annotation_colors = ann_colors,
              show_rownames = F,
              color = brewer.pal(n = 9, name = "Reds"),
              main = "Extracellular matrix structural constituent")



png(filename = "Figures/Heatmap_Extracellular matrix structural constituent.png", width = 150, height = 180, res = 300, units = "mm")
p
dev.off()


# Heatmap for GBM genes in bulk RNA-seq
sampleTable <- read.table("0_NGS2024_7984/DataNorm/DataNormDESeq2.txt", header = T)
sampleTable <- sampleTable[,c(1:9)]
sampleTable$GENEID <- rownames(sampleTable)
sampleTable <- merge(sampleTable, ens2sym_entrez, by="GENEID")
sampleTable <- sampleTable %>% 
  distinct(SYMBOL, .keep_all = TRUE)
rownames(sampleTable) <- sampleTable$SYMBOL
sampleTable <- sampleTable[,-c(1)]

keep <- c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "LAMA1", "LAMB2", "LAMA5", "LOXL2")
sampleTable <- subset(sampleTable, SYMBOL %in% keep)
sampleTable <- sampleTable[,-10]
df2_filtered <- sampleTable


mat_z <- t(apply(df2_filtered,1, scale))
colnames(mat_z) <- c("d22_R1", "d22_R2", "d22_R3", "d32_R1", "d32_R2", "d32_R3", "d42_R1", "d42_R2", "d42_R3")
my_sample_col <- data.frame(Timepoints = rep(c("Early", "Mid", "Late"), c(3,3,3)))
row.names(my_sample_col) <- colnames(mat_z)

ann_colors <- list(
  Timepoints = c("Early" = "#7570b3", "Mid" = "#1b9e77","Late" = "#e7298a")
)

p <- pheatmap(mat_z,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = my_sample_col,
         annotation_colors = ann_colors, cluster_rows = TRUE, cluster_cols = F, show_colnames = F,
         fontsize = 14, color = brewer.pal(n = 9, name = "Reds"))

png(filename = "d22vsd38/Figures/Heatmap_GBM_collagens.png", width = 200, height = 80, res = 300, units = "mm")
p
dev.off()


library(ggrepel)
pro <- pro %>%
  mutate(
    highlight = case_when(
      Pvalue > 1 & LogFC > 0.1 ~ "Up",
      Pvalue > 1 & LogFC < -0.1 ~ "Down",
      TRUE ~ "Not_sig"
    )
  )

color_map <- c(
  "Up" = "darkred",
  "Not_sig" = "azure4",
  "Down" = "blue"
)

volcano <- ggplot(pro,
                  aes(x = LogFC, y = Pvalue, colour= highlight)) +
  geom_point(cex = 1.0, stroke = 0.6, color = "black") +
  geom_point(aes(fill = highlight), cex = 1.0, stroke = 0.5) +
  theme_minimal() +
  labs(x="LogFC",
       y="-Log10(p-value)") +
  scale_color_manual(values=color_map) +
  theme(legend.position="top", axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(face = "bold", size = 10)) +
  geom_hline(yintercept=1, linetype="dashed", color = "gray") +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color = "gray")

volcano

png(filename = "Volcano.d38_vs_d22.perseus.DEP.png", width = 100, height = 100, units = "mm", res = 300)
plot(volcano)
dev.off()
