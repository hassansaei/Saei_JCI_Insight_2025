# =============================================================================
# Proteomics Data Analysis using msqrob2
#
# This script performs comprehensive proteomics data analysis including:
# - Data preprocessing and quality control
# - Peptide-level and protein-level quantification
# - Statistical analysis for differential expression
# - Functional enrichment analysis
# - Visualization (MDS plots, volcano plots, enrichment plots)
# 
# Author: Hassan Saei
# 
# Dependencies: msqrob2, QFeatures, limma, clusterProfiler, ggplot2
# =============================================================================

# install required packages (uncomment if needed)
#BiocManager::install("msqrob2")
#BiocManager::install("msdata")
#BiocManager::install("MSnbase")
#devtools::install_github("bartongroup/proteusLabelFree")
#devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes = FALSE)
#BiocManager::install("statomics/msqrob2gui")
#remotes::install_github("statOmics/msqrob2gui@master", force = T)

# load required libraries
library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(gridExtra)
library(proteusLabelFree)
library(dplyr)

# set working directory and read peptide data
# note: change path according to your data location
setwd("../1st Project/MQ/ECM/")
peptidesFile <- "ECM_peptides.txt"
peptidesFile <- "F1_peptides.txt"  # using soluble fraction dataset for analysis

# read peptide data from file
pep <- read.table(peptidesFile, sep = "\t", header = T)

# data preprocessing: clean protein identifiers and filter data
pep <- pep %>%
  mutate(Proteins = ifelse(Proteins == "", NA, Proteins)) %>%
  #mutate(Intensity.1C3_d37_ECM_R1 = ifelse(Intensity.1C3_d37_ECM_R1 == "", NA, Intensity.1C3_d37_ECM_R1)) %>%
  filter(Proteins != "NA")
  #filter(Intensity.1C3_d37_ECM_R1 != "NA")

# ECM
#pep <- pep[,-c(87:91, 95,96,97)]
# F1
colnames(pep)
pep <- pep[,-c(79:81, 85:87)]  # remove unwanted columns
pep <- pep[,-c(63:65, 69:71)]  # remove additional unwanted columns
pep <- pep[,-c(98:103)]        # remove more unwanted columns
pep <- pep[,-c(95:100)]        # remove final set of unwanted columns

# simplify protein identifiers by taking only the first protein in each group
pep <- pep %>%
  mutate(Proteins = sapply(str_split(Proteins, ";"), function(x) paste(x[1], collapse = ";")))

# check column names and count missing values
colnames(pep)

# function to count blank values in each column
count_blanks <- function(x) {
  sum(x == "")
}
blank_counts <- sapply(pep, count_blanks)
print(blank_counts)

# identify intensity columns for QFeatures object creation
ecols <- grep("Intensity\\.", names(pep))

# create QFeatures object for proteomics analysis
pe <- readQFeatures(assayData = pep, fnames = 1,
                    ecol = ecols, name = "peptideRaw", sep = "\t")

# handle missing values and zero intensities
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)
pe <- zeroIsNA(pe, "peptideRaw")  # convert zeros to NA for proper missing value handling

# visualize missing value patterns
MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

# apply log2 transformation to intensity data
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")

# plot density distributions after log transformation
limma::plotDensities(assay(pe[["peptideLog"]]), legend = FALSE)
legend("topright", legend = colnames(assay(pe[["peptideLog"]])), 
       col = 1:ncol(assay(pe[["peptideLog"]])), lty = 1, cex = 0.6)


# handle overlapping protein groups by selecting smallest unique groups
Protein_filter <- rowData(pe[["peptideLog"]])$Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins)
pe <- pe[Protein_filter,,]

# filter out reverse hits and potential contaminants
pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")

# filter peptides with at least 2 non-zero values across samples
pe <- filterFeatures(pe, ~ nNonZero >= 2)
nrow(pe[["peptideLog"]])  # check number of peptides after filtering


# apply median normalization to reduce technical variation
pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)

# visualize peptide distributions after normalization
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distributions after normalization", 
        ylab = "intensity", las=2)

# save boxplot to file
par(mar = c(10, 4, 4, 2) + 0.1)
png(filename = "ECM/1_Analysis_1/BoxPlot_afterfilteration.png", 
    width = 100, height = 100, res = 300, units = "mm")
boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distributions after normalization",
        ylab = "intensity", cex.axis = 0.3, las = 2, cex.lab=0.4, cex.main = 0.5)
dev.off()

# plot density distributions after normalization
limma::plotDensities(assay(pe[["peptideNorm"]]), legend = FALSE)
legend("topright", legend = colnames(assay(pe[["peptideLog"]])), 
       col = 1:ncol(assay(pe[["peptideLog"]])), lty = 1, cex = 0.6)


# define experimental groups for analysis
# note: adjust group assignments based on your experimental design
group <- factor(c(rep("MT1.PMO", 3), rep("MT1.SC", 3), rep("MT2.PMO", 3), rep("MT2.SC", 3), rep("WT", 4)))
colors <- c("darkred", "darkblue", "darkgreen", "black")[group]

# perform multidimensional scaling (MDS) analysis
limma::plotMDS(assay(pe[["peptideNorm"]]), col = colors, cex = 0.5)

# create enhanced MDS plot
mds <- plotMDS(assay(pe[["peptideNorm"]]), plot = FALSE)
plot(mds$x, mds$y, col = colors, pch = 19, cex = 1.5, 
     xlab = "Dimension 1 (39%)", ylab = "Dimension 2 (33%)", main = "MDS Plot")
legend("topleft", legend = colnames(assay(pe[["peptideLog"]])), 
       col = colors, pch = 19, cex = 0.8, bty = "n")

# save MDS plot to file
png(filename = "PlotMDS_F1_A1.png", width = 140, height = 140, res = 300, units = "mm")
plot(mds$x, mds$y, col = colors, pch = 19, cex = 1.5, 
     xlab = "Dimension 1 (39%)", ylab = "Dimension 2 (33%)", main = "")
dev.off()

# aggregate peptides to protein level for statistical analysis
pe <- aggregateFeatures(pe, i = "peptideNorm", na.rm = TRUE, name = "protein", fcol= "Proteins")

# perform MDS analysis at protein level
limma::plotMDS(assay(pe[["protein"]]), col = colors, cex = 0.8)
mds <- plotMDS(assay(pe[["protein"]]), plot = FALSE)
plot(mds$x, mds$y, col = colors, pch = 19, cex = 1.5, 
     xlab = "Dimension 1 (56%)", ylab = "Dimension 2 (25%)", 
     main = "MDS Plot after aggregation")
legend("topleft", legend = colnames(assay(pe[["peptideLog"]])), 
       col = colors, pch = 19, cex = 0.8, bty = "n")

# set up experimental design for statistical modeling
# note: adjust condition assignments based on your experimental groups
condition <- data.frame(rownames(colData(pe)))
condition$condition <- c("G3", "G3", "G3", "G2", "G2", "G1", "G1", "G1", "G4", "G4", "G4")
colData(pe)$condition <- as.factor(condition$condition)

# fit msqrob2 model for differential expression analysis
pe <- msqrob(object = pe, i = "protein", formula = ~condition, overwrite=TRUE)

# check model coefficients
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

# visualize experimental design matrix
library(ExploreModelMatrix)
VisualizeDesign(colData(pe),~condition)$plotlist[[1]]

# perform differential expression analysis
# get coefficient names from the fitted model
coefs <- names(getCoef(rowData(pe[["protein"]])$msqrobModels[[1]]))[-1] 

# generate all pairwise combinations of coefficients
coef_comb_data <- expand.grid(coefs, coefs) |>
  filter(Var1 != Var2)

# remove duplicate combinations (same comparison in different order)
indx <- !duplicated(t(apply(coef_comb_data, 1, sort)))

# format contrast names
coef_comb <- coef_comb_data[indx, ] |>
  unite(contrast, Var1, Var2, sep = " - ") |>
  pull(var = contrast, name = NULL)

# create contrast list for hypothesis testing
contrast_names <- c(coefs, coef_comb)
contrast_list <- paste0(c(coefs, coef_comb), "=0")

# make contrasts matrix
contrasts <- makeContrast(
  contrast_list,
  parameterNames = coefs
)

# perform hypothesis testing for specific contrast (conditionG2 vs reference)
contrasts <- makeContrast("conditionG2=0", parameterNames = c("conditionG2"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = contrasts, overwrite = TRUE)

# extract differentially expressed proteins (DEPs) with adjusted p-value < 0.05
DEG_d34d22 <- rowData(pe[["protein"]])$"conditionG2"
DEG_d34d22 <- rowData(pe[["protein"]])$"conditionG2" %>%
  filter(adjPval < 0.05)

# note: the following lines appear to be for a different comparison
# DEG_d38d22 <- rowData(pe[["protein"]])$"conditiond37"
# DEG_d38d22 <- rowData(pe[["protein"]])$"conditiond37" %>%
#   filter(adjPval < 0.05)

# merge with gene names for annotation
match <- pep[,c("Proteins", "Gene.names")]
DEG_d34d22$Proteins <- rownames(DEG_d34d22)
# DEG_d38d22$Proteins <- rownames(DEG_d38d22)
DEG_d34d22 <- merge(DEG_d34d22, match, by="Proteins")
# DEG_d38d22 <- merge(DEG_d38d22, match, by="Proteins")

# clean gene names by taking only the first gene symbol
DEG_d34d22 <- DEG_d34d22 %>%
  mutate(Gene.names = sapply(str_split(Gene.names, ";"), `[`, 1))

# DEG_d38d22 <- DEG_d38d22 %>%
#   mutate(Gene.names = sapply(str_split(Gene.names, ";"), `[`, 1))

# remove duplicate genes, keeping the first occurrence
DEG_d34d22 <- DEG_d34d22 %>%
  distinct(Gene.names, .keep_all = TRUE)

# DEG_d38d22 <- DEG_d38d22 %>%
#   distinct(Gene.names, .keep_all = TRUE)


# rename columns for better readability
colnames(DEG_d34d22) <- c("ProteinIDs", "logFC", "se", "df", "t", "pval", "adjPval", "Symbol")
# colnames(DEG_d38d22) <- c("ProteinIDs", "logFC", "se", "df", "t", "pval", "adjPval", "Symbol")

# load annotation databases for functional classification
NABA <- read.table("../../../Alport paper/11_Bulk_RNASEQ/NABA_MATRISOME.v2023.2.txt", sep = "\t", header = T)
BM <- read.table("proteins_Membrane.csv", sep = "\t", header = T)
colnames(BM) <- c("ProteinIDs", "Description")

# annotate proteins with matrix and membrane classifications
DEG_d34d22 <- DEG_d34d22 %>%
  mutate(Matrix = ifelse(Symbol %in% NABA$Symbol, "Matrix", "Other")) %>%
  mutate(BM = ifelse(ProteinIDs %in% BM$ProteinIDs, "BM", "Other"))

# DEG_d38d22 <- DEG_d38d22 %>%
#   mutate(Matrix = ifelse(Symbol %in% NABA$Symbol, "Matrix", "Other")) %>%
#   mutate(BM = ifelse(ProteinIDs %in% BM$ProteinIDs, "BM", "Other"))

# save results to files
write.table(DEG_d34d22, file = "Proteomics_1C3vsCTRL_ECM_DEP.txt", sep = "\t", row.names = F)
# write.table(DEG_d38d22, file = "Proteomics_d38vsd22_F1_DEP.txt", sep = "\t", row.names = F)


# identify matrix proteins for comparison
matrix_d34 <- DEG_d34d22$Matrix == "Matrix"
# matrix_d38 <- DEG_d38d22$Matrix == "Matrix"

genes_d34 <- DEG_d34d22 %>%
  filter(matrix_d34) %>%
  pull(Symbol)

# genes_d38 <- DEG_d38d22 %>%
#   filter(matrix_d38) %>%
#   pull(Symbol)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(Genes1 = genes_d34, Genes2 = genes_d38),
  category.names = c("Day 32", "Day 42"),
  filename = NULL, # Set to NULL to display in RStudio plot window
  output = TRUE,
  fill = c("skyblue", "lightgreen"),
  alpha = 0.5, # Transparency
  cat.col = c("black", "black"),
  cat.cex = 1.5, # Category label size
  cex = 2, # Font size for numbers
  lwd = 3, # Line width
  col = "black" # Line color
)
grid.draw(venn.plot)

# perform functional enrichment analysis
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# get gene list for enrichment analysis
original_gene_list <- DEG_d34d22$Symbol

# load gene set databases for enrichment analysis
gene_sets_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_sets_GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets_GO_CC <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
gene_sets_GO_MF <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
gene_sets_KEGG <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
gene_sets_Reactome <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")

# prepare gene sets for enrichment analysis
gene_sets_hallmark <- gene_sets_hallmark %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_BP <- gene_sets_GO_BP %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_CC <- gene_sets_GO_CC %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_GO_MF <- gene_sets_GO_MF %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_KEGG <- gene_sets_KEGG %>%
  dplyr::select(gs_name, gene_symbol)

gene_sets_Reactome <- gene_sets_Reactome %>%
  dplyr::select(gs_name, gene_symbol)

# perform gene ontology enrichment analysis
CC <- as.data.frame(enricher(gene = original_gene_list,
                             TERM2GENE = gene_sets_GO_CC, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
BP <- as.data.frame(enricher(gene = original_gene_list, TERM2GENE = gene_sets_GO_BP,
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))
MF <- as.data.frame(enricher(gene = original_gene_list, TERM2GENE = gene_sets_GO_MF,
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05))

# combine enrichment results from different GO categories
df_combined <- rbind(CC, BP, MF)
df_combined$log10_Adjusted <- -log10(df_combined$qvalue)
df_combined <- df_combined %>%
  mutate(ONTOLOGY = case_when(
    grepl("^GOMF_", Description) ~ "MF",
    grepl("^GOCC_", Description) ~ "CC",
    grepl("^GOBP_", Description) ~ "BP",
    TRUE ~ NA_character_
  ),
  Description = gsub("GOMF_|GOCC_|GOBP_", "", Description))

# save enrichment results
write.table(df_combined, file = "GO_d38vsd22_F1.txt", sep = "\t")

# select top 4 enriched terms from each category for visualization
df_BP <- BP[1:4,]
df_CC <- CC[1:4,]
df_MF <- MF[1:4,]

# combine top enriched terms for plotting
df_combined <- rbind(df_BP, df_CC, df_MF)
df_combined$ID <- gsub("GOMF_|GOCC_|GOBP_", "", df_combined$ID)
df_combined$log10_Adjusted <- -log10(df_combined$qvalue)
df_combined <- df_combined %>%
  mutate(ONTOLOGY = case_when(
    grepl("^GOMF_", Description) ~ "MF",
    grepl("^GOCC_", Description) ~ "CC",
    grepl("^GOBP_", Description) ~ "BP",
    TRUE ~ NA_character_
  ),
  Description = gsub("GOMF_|GOCC_|GOBP_", "", Description))

# create enrichment bar plot
p <- ggplot(df_combined, aes(x = log10_Adjusted, y = reorder(ID, log10_Adjusted), fill = ONTOLOGY)) +
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
    plot.title = element_text(size = 16, face = "bold")
  )

# save enrichment plot
png(filename = "Enrichment_d38vsd22.png", width = 200, height = 150, res = 300, units = "mm")
plot(p)
dev.off()

# create volcano plot with functional annotations
DEG_d34d22 <- DEG_d34d22 %>%
  mutate(
    highlight = case_when(
      pval < 0.05 & grepl("BM", BM) ~ "BM",
      pval < 0.05 & grepl("Matrix", Matrix) ~ "Matrix",
      pval < 0.05 ~ "Sig",
      TRUE ~ "Not_sig"
    )
  )

DEG_d38d22 <- DEG_d38d22 %>%
  mutate(
    highlight = case_when(
      pval < 0.05 & grepl("BM", BM) ~ "BM",
      pval < 0.05 & grepl("Matrix", Matrix) ~ "Matrix",
      pval < 0.05 ~ "Sig",
      TRUE ~ "Not_sig"
    )
  )

color_map <- c(
  "Sig" = "lightblue",
  "Not_sig" = "black",
  "Matrix" = "blue",
  "BM" = "darkred"
)

# create volcano plot
library(ggrepel)
volcano <- ggplot(DEG_d34d22,
  aes(x = logFC, y = -log10(pval), colour= highlight)) +
  geom_point(cex = 2.0, stroke = 0.6, color = "black") +
  geom_point(aes(fill = highlight), cex = 2.0, stroke = 0.5) +
  theme_minimal() +
  labs(x="Log2(fold change)",
       y="-Log10(p-value)") +
  scale_color_manual(values=color_map) +
  theme(legend.position="top", axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "gray") +
  geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed", color = "gray")

# add gene labels for significant matrix and membrane proteins
volcano <- volcano +
  geom_text_repel(data = subset(DEG_d34d22, highlight %in% c("BM", "Matrix") & abs(logFC) > 1.2),
                  aes(label = Symbol),
                  size = 3,
                  max.overlaps = 20)
volcano

# save volcano plot
png(filename = "Volcano_F1_A1_d38vsd22.png", width = 180, height = 150, res = 300, units = "mm")
plot(volcano)
dev.off()
