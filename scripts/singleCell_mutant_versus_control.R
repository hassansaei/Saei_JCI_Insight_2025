################################################################################
# Script to analyze single-cell data obtained from comparing COL4A5 mutated organoids with isogenic controls
# Author: Hassan Saei
# Contact: hassan.saeiahan@gmail.com
# Integrated Seurat object with cell type annotations (.h5ad file is available in the Zenodo repository)
################################################################################

# load required libraries
library(dplyr)          # data manipulation and transformation
library(Seurat)         # single-cell analysis toolkit
library(reticulate)     # interface to Python
library(anndata)        # Python anndata interface
library(patchwork)      # plot composition and layout
library(reshape2)       # data reshaping functions
library(RColorBrewer)   # color palettes for plots
library(ggplot2)        # grammar of graphics plotting
library(SoupX)          # ambient RNA removal
library(openxlsx)       # Excel file reading/writing
library(scCustomize)    # custom Seurat functions
library(dittoSeq)       # single-cell visualization

# define sample identifiers for the organoid datasets
id <- c("2G9_B1", "2G9_B2", "1C3_B1", "2F10_B1")

# function to identify outliers using median absolute deviation (MAD)
# obj: Seurat object
# metric: column name in metadata to evaluate
# nmads: number of MADs from median to consider as outlier
mad_outlier <- function(obj, metric, nmads){
  M <- obj@meta.data[[metric]]                    # extract metric values
  median_M <- median(M, na.rm = TRUE)             # calculate median
  mad_M <- mad(M, na.rm = TRUE)                   # calculate median absolute deviation
  outlier <- (M < (median_M - nmads * mad_M)) | (M > (median_M + nmads * mad_M))  # identify outliers
  return(outlier)
}

# set Seurat object assay version to v3 for compatibility
options(Seurat.object.assay.version = "v3")

# function to preprocess individual samples
# id: sample identifier string
preprocess <- function(id){
  path <- paste0("0_rawData/", id, "/filtered_feature_bc_matrix/")  # construct path to 10X data
  obj <-  Read10X(data.dir = path)                                    # read 10X genomics data
  obj <- CreateSeuratObject(counts = obj, project = id, min.cells = 0, min.features = 200)  # create Seurat object
  obj$sample_id <- id                                                 # add sample identifier to metadata
  
  # calculate quality control metrics
  obj$log1p_total_counts <- log1p(obj@meta.data$nCount_RNA)           # log-transformed total UMI counts
  obj$log1p_n_genes_by_counts <- log1p(obj@meta.data$nFeature_RNA)   # log-transformed gene counts
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # mitochondrial gene percentage
  obj[["percent.rbp"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")  # ribosomal gene percentage
  
  # outlier detection (commented out - will be applied later)
  #bool_vector <- !mad_outlier(obj, 'log1p_total_counts', 5) & !mad_outlier(obj, 'log1p_n_genes_by_counts', 5) & !mad_outlier(obj, 'percent.mt', 3)
  #obj <- subset(obj, cells = which(bool_vector))
  
  return(obj)
}


# create list of Seurat objects by preprocessing each sample
data_list <- sapply(id, preprocess)

# print dimensions of each dataset before filtering
for (i in 1:length(data_list)){
  print(dim(data_list[[i]]))
}


# function to create and save violin plots for quality control metrics
createAndSaveVlnPlot <- function(data, filename) {
  png(filename, width = 952, height = 631)  # create PNG file with specified dimensions
  p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4)  # create violin plot
  print(p)                                                                                                  # display plot
  dev.off()                                                                                                 # close graphics device
}

# generate violin plots for all samples to assess quality control metrics
for (i in 1:length(data_list)) {
  filename <- paste0("VlnPlot_AF_", i, ".png")  # create filename for each sample
  createAndSaveVlnPlot(data_list[[i]], filename)  # generate and save plot
  cat("Saved plot", i, "as", filename, "\n")      # print confirmation message
}


# function to filter out low-quality cells using MAD-based outlier detection
preprocess_2 <- function(obj){
  
  # identify outliers based on multiple QC metrics
  bool_vector <- !mad_outlier(obj, 'log1p_total_counts', 5) & !mad_outlier(obj, 'log1p_n_genes_by_counts', 5) & !mad_outlier(obj, 'percent.mt', 3)
  obj <- subset(obj, cells = which(bool_vector))  # keep only non-outlier cells
  
  return(obj)
}

# apply cell filtering to all samples
data_list <- sapply(data_list, preprocess_2)

# print dimensions of each dataset after filtering
for (i in 1:length(data_list)){
  print(dim(data_list[[i]]))
}

# perform SCTransform normalization and dimensionality reduction for each sample individually
for (i in 1:length(data_list)) {
  data_list[[i]] <- SCTransform(data_list[[i]], vst.flavor = "v2", verbose = FALSE, vars.to.regress = "percent.mt")  # normalize data and regress out mitochondrial genes
  data_list[[i]] <- RunPCA(data_list[[i]], features = VariableFeatures(object = data_list[[i]]), npcs = 50)          # principal component analysis
  data_list[[i]] <- FindNeighbors(data_list[[i]], dims = 1:50, reduction = "pca", k.param = 20)                     # build k-nearest neighbor graph
  data_list[[i]] <- FindClusters(data_list[[i]], resolution = 0.8, cluster.name = "unintegrated_clusters")           # perform clustering
  data_list[[i]] <- RunUMAP(data_list[[i]], dims = 1:50, reduction = "pca", reduction.name = "umap_unintegrated_0.8") # UMAP embedding
  data_list[[i]] <- RunTSNE(data_list[[i]], reduction = "pca", dims = 1:50)                                          # t-SNE embedding
}


# visualize UMAP and gene expression for sample 2 (1C3_B1)
DimPlot(data_list[[2]], reduction = "umap_unintegrated_0.8", label = TRUE, pt.size = 1)  # corrected reduction name
VlnPlot(data_list[[2]], features = c("GATA3", "NPHS2"))  # plot expression of kidney development markers

# ScType automated cell type annotation
# load required libraries for cell type annotation
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")  # load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")      # load scoring function

# define database and tissue type for annotation
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue = "Kidney"  # tissue type for kidney-specific cell type annotation

# prepare gene sets for kidney cell types
gs_list = gene_sets_prepare(db_, tissue)

# perform automated cell type annotation using ScType for each sample
for (i in 1:length(data_list)){
  # calculate cell type scores using positive and negative gene sets
  es.max = sctype_score(scRNAseqData = data_list[[i]][["SCT"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # calculate cluster-level scores for each cell type
  cL_resutls = do.call("rbind", lapply(unique(data_list[[i]]@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data_list[[i]]@meta.data[data_list[[i]]@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data_list[[i]]@meta.data$seurat_clusters==cl)), 10)
  }))
  
  # select top scoring cell type for each cluster
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"  # assign unknown for low-confidence predictions
  
  # initialize cell type identity column
  data_list[[i]]@meta.data$ScType_ident = ""
  
  # assign cell type identities to cells based on cluster annotations
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    data_list[[i]]@meta.data$ScType_ident[data_list[[i]]@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
}

# function to create and save UMAP plots colored by cell type annotation
createAndSaveVlnPlot <- function(data, filename) {
  png(filename, width = 1024, height = 631)  # create PNG file with specified dimensions
  p <- DimPlot(data, reduction = "umap_unintegrated_0.8", pt.size = 1.5)  # create UMAP plot
  print(p)                                                                  # display plot
  dev.off()                                                                 # close graphics device
}

# generate UMAP plots colored by ScType annotations for all samples
for (i in 1:length(data_list)) {
  filename <- paste0("uMAP_nolabel_res0.8_", i, ".png")  # create filename for each sample
  Idents(data_list[[i]]) <- "ScType_ident"               # set cell type as identity
  createAndSaveVlnPlot(data_list[[i]], filename)          # generate and save plot
  cat("Saved plot", i, "as", filename, "\n")              # print confirmation message
}


# save the processed data list for later use
saveRDS(data_list, file = "data_list_sct.rds")

################################################################################

# load pre-integrated data (alternative to running integration)
data_list <- readRDS(file = "NGS2024_7386/2_Integrated/Data_integrated.sct.cca.rds")

# select integration features and prepare for integration
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3500)  # select highly variable features
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)  # prepare SCT integration

# find integration anchors and integrate data
data.anchors.cca <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = features, reduction = "cca")
integrated.cca <- IntegrateData(anchorset = data.anchors.cca, normalization.method = "SCT")

# load the final integrated dataset
integrated.cca <- readRDS("integrated.sct.cca.a2.v2.rds")

# switch to RNA assay for additional processing
DefaultAssay(integrated.cca) <- "RNA"
integrated.cca <- NormalizeData(object = integrated.cca)  # normalize RNA assay
integrated.cca <- ScaleData(integrated.cca)              # scale RNA assay

# switch to integrated assay for dimensionality reduction and clustering
DefaultAssay(integrated.cca) <- "integrated"
integrated.cca <- RunPCA(integrated.cca, npcs = 50, verbose = FALSE)  # principal component analysis
integrated.cca <- FindNeighbors(integrated.cca, dims = 1:50, reduction = "pca", k.param = 20)  # build k-nearest neighbor graph
integrated.cca <- FindClusters(integrated.cca, resolution = 0.5, cluster.name = "cca_integrated_res0.5")  # clustering at resolution 0.5
integrated.cca <- FindClusters(integrated.cca, resolution = 1.6, cluster.name = "cca_integrated_res1.6")  # clustering at resolution 1.6
integrated.cca <- RunUMAP(integrated.cca, dims = 1:50, reduction = "pca", reduction.name = "umap.cca")  # UMAP embedding
integrated.cca <- RunTSNE(integrated.cca, reduction = "pca", dims = 1:50)  # t-SNE embedding


## CELL TYPE LABEL TRANSFER SECTION
# transfer cell type labels from reference dataset to query dataset
ref <- readRDS("integrated.sct.cca.a1.v2.rds")  # load reference dataset with known cell types
query <- integrated.cca                          # use current integrated dataset as query

# find highly variable genes for transfer
features <- SelectIntegrationFeatures(
  object.list = list(ref, query),
  nfeatures = 3000
)

# ensure features are present in both datasets
features <- intersect(features, rownames(ref))
features <- intersect(features, rownames(query))
                      
common_genes <- intersect(rownames(ref), rownames(query))  # find common genes between datasets

# set default assay to SCT for both datasets
DefaultAssay(ref) <- 'SCT'
DefaultAssay(query) <- 'SCT'

# find transfer anchors between reference and query
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  normalization.method = "SCT",
  dims = 1:30,
  features = features
)

# define the column containing cell type labels in reference
label_col <- "CellType_V3"

# transfer cell type labels from reference to query
pred <- TransferData(
  anchorset = anchors,
  refdata = ref[[label_col, drop = TRUE]],
  dims = 1:30
)

# add predictions to query dataset
query <- AddMetaData(query, pred)
query$predicted_label <- query$predicted.id

# examine prediction quality
head(query@meta.data[, c("predicted_label", "prediction.score.max")])

# count cells per predicted label
table(query$predicted_label)

# filter low-confidence assignments using threshold
threshold <- 0.5
query$CellType_v3 <- ifelse(
  query$prediction.score.max >= threshold, query$predicted_label, "Unassigned"
)
table(query$CellType_v3)


# visualize cell type annotations on UMAP
Idents(query) <- "CellType_v3"
p <- DimPlot_scCustom(query, reduction = "umap.cca", label = T, pt.size = 0.7,  # corrected reduction name
        label.size = 4, colors_use = "stepped", split_seurat = TRUE) + NoLegend()

# save UMAP plot with cell type annotations
png(filename = "uMAP.annotation.WT_vs_MT.png", width = 250, height = 200, res = 300, units = "mm")
plot(p)
dev.off()

# update integrated dataset with transferred labels
integrated.cca <- query

# MANUAL CELL TYPE ANNOTATION BASED ON CLUSTERS
# define cluster-to-cell-type mapping for kidney organoid cell types
cluster_mapping <- list(
  'Tubular epithelial cells' = c(20,1,21,11,36,9,18,26,16),  # kidney tubule clusters
  'Podocytes' = c(6,2,14,8,31,24,10),                        # glomerular podocyte clusters
  'Other' = c(22,33,27,29,17,19,28,34,15, 35, 30),          # unclassified clusters
  'Stroma' = c(7,5,3,0,4,32,13,25,12,23)                    # stromal/mesenchymal clusters
)

# create empty vector to store segment annotations
integrated.cca@meta.data$segments <- NA

# assign cell type annotations based on cluster membership
for (cluster_name in names(cluster_mapping)) {
  cluster_numbers <- cluster_mapping[[cluster_name]]
  integrated.cca@meta.data$segments[
    integrated.cca@meta.data$cca_integrated_res1.6 %in% cluster_numbers  # corrected column name
  ] <- cluster_name
}

# SAMPLE GROUP ANNOTATION
# define sample-to-group mapping for experimental conditions
cluster_mapping <- list(
  'Control' = c("2G9_B1", "2G9_B2"),    # control samples
  'MT_severe' = c("2F10_B1"),            # severe mutant sample
  'MT_moderate' = c("1C3_B1")            # moderate mutant sample
)

# create empty vector to store group annotations
integrated.cca@meta.data$group2 <- NA

# assign group annotations based on sample identity
for (cluster_name in names(cluster_mapping)) {
  cluster_numbers <- cluster_mapping[[cluster_name]]
  integrated.cca@meta.data$group2[
    integrated.cca@meta.data$orig.ident %in% cluster_numbers
  ] <- cluster_name
}


# DIFFERENTIAL EXPRESSION ANALYSIS
# find marker genes for each cluster
integrated.cca.markers <- FindAllMarkers(integrated.cca, only.pos = TRUE)  # find positive markers only
integrated.cca.markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  # filter for significant markers


# GENE EXPRESSION VISUALIZATION
DefaultAssay(integrated.cca) <- "SCT"

# visualize podocyte-specific markers
VlnPlot(integrated.cca, features = c("TGFB1"), cols = c("lightgray", "darkblue"))  # TGF-beta1 expression
VlnPlot(integrated.cca, features = c("DDN"))                                      # dendrin expression
FeaturePlot(integrated.cca, features = c("TGFB1"), cols = c("lightgray", "darkblue"))  # spatial expression
# visualize podocyte precursor markers
VlnPlot(integrated.cca, features = c("CTGF", "OLFM3", "MAFB", "NPHS1"))  # connective tissue growth factor, olfactomedin, MAF bZIP, nephrin
FeaturePlot(integrated.cca, features = c("CTGF", "OLFM3", "MAFB", "NPHS1"), cols = c("lightgray", "darkblue"))

# visualize nephron progenitor markers
VlnPlot(integrated.cca, features = c("DAPL1", "LYPD1", "SIX1", "SIX2","CRABP2"))  # death-associated protein, LY6/PLAUR, SIX homeobox genes
FeaturePlot(integrated.cca, features = c("DAPL1", "LYPD1", "SIX1", "SIX2","CRABP2"), cols = c("lightgray", "darkblue"))

# visualize epithelial markers
VlnPlot(integrated.cca, features = c("PAX2", "PAX8", "KRT19", "EPCAM","LRP2"), ncol = 4)  # paired box genes, keratin, epithelial cell adhesion, LDL receptor
FeaturePlot(integrated.cca, features = c("PAX2", "PAX8", "KRT19", "EPCAM","LRP2"), cols = c("lightgray", "darkblue"))

# visualize loop of Henle markers
VlnPlot(integrated.cca, features = c("ESRRG", "POU3F3", "SLC12A1"), ncol = 3)  # estrogen-related receptor, POU class 3, sodium-potassium-chloride cotransporter
FeaturePlot(integrated.cca, features = c("ESRRG", "POU3F3", "SLC12A1"), cols = c("lightgray", "darkblue"))

# visualize distal progenitor markers
VlnPlot(integrated.cca, features = c("EPCAM", "EMX2", "SPP1", "MAL", "PAX2", "GATA3", "TFAP2A"), ncol = 4)  # epithelial, empty spiracles, osteopontin, myelin, paired box, GATA binding, transcription factor
FeaturePlot(integrated.cca, features = c("EPCAM", "EMX2", "SPP1", "MAL", "PAX2", "GATA3", "TFAP2A"), cols = c("lightgray", "darkblue"))

# visualize proximal tubule precursor markers
VlnPlot(integrated.cca, features = c("IGFBP7", "FXYD2", "CDH6", "HNF1B", "HNF4G", "HNF4A", "SLC3A1"), ncol = 4)  # insulin-like growth factor binding, FXYD domain, cadherin, hepatocyte nuclear factors, solute carrier
FeaturePlot(integrated.cca, features = c("IGFBP7", "FXYD2", "CDH6", "HNF1B", "HNF4G", "HNF4A", "SLC3A1"),
            cols = c("lightgray", "darkblue"), max.cutoff = 1.5, min.cutoff = 0.5)
# visualize claudin tight junction markers
VlnPlot(integrated.cca, features = c("CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN7", "CLDN8"), ncol = 4)  # removed duplicate CLDN4
FeaturePlot(integrated.cca, features = c("CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN7", "CLDN8"), cols = c("lightgray", "darkblue"))

# visualize collagen markers
VlnPlot(integrated.cca, features = c("COL4A5", "COL4A6"), ncol = 2)  # type IV collagen alpha chains
FeaturePlot(integrated.cca, features = c("COL1A1"),
            cols = c("lightgray", "darkblue"), reduction = "tsne", split.by = "group")  # type I collagen
# Muscle progenitor
VlnPlot(integrated.cca, features = c("MYOG", "MYOD1"))
FeaturePlot(integrated.cca, features = c("MYOG", "MYOD1", "SIX1", "SIX2"), cols = c("lightgray", "darkblue"))
# Neural progenitor
VlnPlot(integrated.cca, features = c("HES6", "STMN2"))
FeaturePlot(integrated.cca, features = c("HES6", "STMN2"), cols = c("lightgray", "darkblue"))
# Glial
VlnPlot(integrated.cca, features = c("FABP7", "TTYH1","SOX2"))
FeaturePlot(integrated.cca, features = c("FABP7", "TTYH1","SOX2"), cols = c("lightgray", "darkblue"))
# Cell cycle (Mitotic cell cycle process)
VlnPlot(integrated.cca, features = c("CENPF", "HMGB2", "UBE2C","HIST1H4C", "PCLAF", "TYMS"))
FeaturePlot(integrated.cca, features = c("CENPF", "HMGB2", "UBE2C","HIST1H4C", "PCLAF", "TYMS"), cols = c("lightgray", "darkblue"))
# Mesangial cells
VlnPlot(integrated.cca, features = c("LAMC3", "NKD1", "TNC", "PHACTR3", "PDGFRB"))
FeaturePlot(integrated.cca, features = c("LAMC3", "NKD1", "TNC", "PHACTR3", "PDGFRB"), cols = c("lightgray", "darkblue"))
# Myofibroblast
VlnPlot(integrated.cca, features = c("COL4A1", "COL1A2", "POSTN", "COL3A1", "COL6A3", "COL14A1", "OGN", "GREM2", "FRZB"))
FeaturePlot(integrated.cca, features = c("COL4A1", "COL1A2", "POSTN", "COL3A1", "COL6A3", "COL14A1", "OGN", "GREM2", "FRZB"), cols = c("lightgray", "darkblue"))
FeaturePlot(integrated.cca, features = c("ACTA2", "POSTN", "COL1A1", "COL1A2"), cols = c("lightgray", "darkblue"), min.cutoff = 2.5, max.cutoff = 4)
# Endothelium
VlnPlot(integrated.cca, features = c("CLDN5", "PECAM1", "KDR", "GNG11", "CALM1"))

# create bar plots showing cell type proportions
p <- dittoBarPlot(integrated.cca, "CellType_v3", group.by = "group2", ylab = "Percentage of cell types")  # corrected group.by
p <- dittoBarPlot(integrated.cca, "segments", group.by = "group2", ylab = "Percentage of cell types",  # corrected group.by
                  color.panel = c("red", "deeppink", "gray27", "lightgray"))


# create dot plot for gene expression comparison
Idents(integrated.cca) <- "group2"  # corrected identity
DefaultAssay(integrated.cca) <- "SCT"

genes <- c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "LAMA1", "LAMA5", "LAMB1", "LAMB2")
genes <- c("COL4A1", "COL4A2", "COL4A5", "COL4A6", "LAMA1", "LAMA2", "LAMB1", "HSPG2", "NID2", "ITGB3", "MAGI1", "SLC9A3R2",
           "ARHGAP24", "NEXN", "DIAPH3", "MYO10", "TNS1", "TNS3", "RND3", "SORBS2", "TIMP3", "ADAMTS1", "HTRA1", "FBLN1", "FBLN2", "VCAN")

p <- DotPlot(integrated.cca, features= genes,
        dot.scale=10, group.by ="group", cols="RdBu") + theme(axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1,face = "bold"),
                                                                axis.text.y = element_text(size = 14, vjust = 0.5, hjust=1, face = "italic"), axis.title = element_blank()) + coord_flip()
png(filename = "DotPlot.MT_vs_WT.png", width = 150,height = 150, units = "mm", res = 300)
plot(p)
dev.off()


# DIFFERENTIAL EXPRESSION ANALYSIS
# define sample groups for comparison
cluster_mapping <- list(
  'Control' = c("2G9_B1", "2G9_B2"),
  'Severe' = c("2F10_B1"),
  'Moderate' = c("1C3_B1")
)

# create sample annotation column
integrated.cca@meta.data$sample <- NA

# assign sample groups based on orig.ident
for (cluster_name in names(cluster_mapping)) {
  cluster_numbers <- cluster_mapping[[cluster_name]]
  integrated.cca@meta.data$sample[
    integrated.cca@meta.data$orig.ident %in% cluster_numbers
  ] <- cluster_name
}

# create combined condition annotation
integrated.cca$condition2 <- paste(integrated.cca@meta.data$sample, integrated.cca@meta.data$segments, sep = ".")

# perform differential expression analysis
DefaultAssay(integrated.cca) <- "RNA"
Idents(integrated.cca) <- "condition2"

# compare severe mutant podocytes vs control podocytes
de.severe <- FindMarkers(integrated.cca, ident.1 = "Severe.Podocytes",
                       ident.2 = "Control.Podocytes", verbose = FALSE)

# compare moderate mutant podocytes vs control podocytes
de.moderate <- FindMarkers(integrated.cca, ident.1 = "Moderate.Podocytes",
                       ident.2 = "Control.Podocytes", verbose = FALSE)


ribo.genes <- read.table(file = "Ribosomal_Genes.txt", col.names = "ribo_genes")
de.severe$ribo <- rownames(de.severe) %in% ribo.genes$ribo_genes
de.moderate$ribo <- rownames(de.moderate) %in% ribo.genes$ribo_genes
de.severe <- de.severe[de.severe$ribo == FALSE, -6]
de.moderate <- de.moderate[de.moderate$ribo == FALSE, -6]

de.severe <- de.severe %>%
  filter(p_val_adj < 0.000001 & abs(avg_log2FC) > 0.5)
de.moderate <- de.moderate %>%
  filter(p_val_adj < 0.000001 & abs(avg_log2FC) > 0.5)

de.severe$Symbol <- rownames(de.severe)
de.moderate$Symbol <- rownames(de.moderate)

# Find common and unique genes
common_genes <- intersect(de.severe$Symbol, de.moderate$Symbol)
unique_genes_de.severe <- setdiff(de.severe$Symbol, de.moderate$Symbol)
unique_genes_de.moderate <- setdiff(de.moderate$Symbol, de.severe$Symbol)
write.table(common_genes, file = "Common_DEG_severe_moderate_XLAS.txt", sep = "\t")

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    df1 = de.severe$Symbol,
    df2 = de.moderate$Symbol
  ),
  category.names = c("Severe", "Moderate"),
  filename = NULL,
  fill = c("blue", "red"),
  output = TRUE
)
grid.draw(venn.plot)

# write differential expression results to files
write.table(de.severe, file = "DEG_2F10vsCTRL_Podocytes.txt", sep = "\t", row.names = F)  # corrected variable name
write.table(de.moderate, file = "DEG_1C3vsCTRL_Podocytes.txt", sep = "\t", row.names = F)  # corrected variable name

# visualize RORB expression across conditions
DefaultAssay(integrated.cca) <- "RNA"
VlnPlot(integrated.cca, features = c("RORB"), idents = c("Control", "Severe", "Moderate"), group.by = "group2")  # corrected identities and group.by 

# GENE COUNT MATRIX GENERATION FOR BULK-LIKE ANALYSIS
# subset podocytes and tubular cells
podo <- subset(integrated.cca, subset = segments == "Podocytes")
tubules <- subset(integrated.cca, subset = segments == "Tubular epithelial cells")

# aggregate expression by group for podocytes
df_macro <- AggregateExpression(podo, group.by = c("group2"),  # corrected group.by
                                assays = "RNA", slot = "counts", return.seurat = F)

# aggregate expression by group for tubular cells
df_macro_tubules <- AggregateExpression(tubules, group.by = c("group2"),  # corrected group.by
                                assays = "RNA", slot = "counts", return.seurat = F)

# aggregate expression by group for all cells
df_macro_all <- AggregateExpression(integrated.cca, group.by = c("group2"),  # corrected group.by
                                assays = "RNA", slot = "counts", return.seurat = F)


df_macro  <- as.matrix(df_macro$RNA)
df_macro <- as.data.frame(df_macro)
names(df_macro) <- c("Control_1", "Control_2", "MT_moderate", "MT_severe")
write.table(df_macro, file = "CountTable_tubules_raw.txt", sep = "\t")

colData <- data.frame(samples = colnames(df_macro))
colData$condition <- ifelse(colData$samples %in% c("Control_1", "Control_2"), "Control", "Mutant")
colData <- as.data.frame(colData)
rownames(colData) <- colData$samples

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = df_macro, colData = colData, design = ~condition )
# Filter dds
keep.1 <- rowSums(counts(dds)) >= 10
dds <- dds[keep.1,]
dds  <- DESeq(dds)
dim(dds)

normCount <- counts(dds, normalized=T)
normCount <- as.data.frame(normCount)
normCount$Symbol <- rownames(normCount)
write.csv(normCount, "Normalized_Podo_DEseq2.csv")
normCount <- read.csv("A2_WTvsMT/Normalized_Podo_DEseq2.csv")


resultsNames(dds)
res <- results(dds, name = "condition_Mutant_vs_Control", alpha = 0.05)
res <- as.data.frame(res)
res$Symbol <- rownames(res)
res <- res %>%
  filter(padj< 0.05)
write.table(res, file = "DEG.podocyte_DESeq2.txt", sep = "\t", row.names = F)
write.table(res, file = "DEG.tubules_DESeq2.txt", sep = "\t", row.names = F)

keep <- c("CDKN1C",
          "NRP2",
          "PTPRR",
          "SEMA3A",
          "PDE1A",
          "SLC7A11",
          "PTPRK",
          "NID2",
          "RRAS",
          "PTTG1",
          "ADAMTS1",
          "TIMP3",
          "E2F3",
          "SLIT3",
          "S100A10",
          "MAP3K5",
          "CA12",
          "SVIL",
          "ARNT2",
          "JAG1",
          "IL1R1",
          "PLK2",
          "CAV1",
          "ABCC9",
          "ANK3",
          "BMP4",
          "INPP4B",
          "COL1A2",
          "COL4A2")

keep <- c("ANXA1", "COL13A1", "LAMA2", "LAMA1", "COL23A1", "HTRA1",
          "HSPG2", "COL1A2", "IGFBPL1", "COL4A2", "SMOC1",
          "HMCN1","S100A10", "S100A7", "LAMA2", "COL4A2", "LAMA1", "SMOC1",
          "COL23A1", "HMCN1", "HSPG2", "COL1A2", "COL4A2", "LAMA2",
          "LAMA1", "ITGA8", "HSPG2")

normCount <- as.data.frame(normCount)
normCount$Symbol <- rownames(normCount)
normCount <- subset(normCount, X %in% keep)
rownames(normCount) <- normCount$X
normCount <- normCount[,-1]
normCount <- as.matrix(normCount)
mat.z <- t(apply(normCount, 1, scale))
colnames(mat.z) <- c("Control_1", "Control_2", "MT_moderate", "MT_severe")
mat.z <- mat.z[,c("MT_moderate", "MT_severe","Control_1", "Control_2")]
library(pheatmap)
pheatmap(mat.z,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete", fontsize_row = 10,
         main = "", fontsize = 10, color = colorRampPalette(c("blue", "white", "red"))(100))

png(filename = "Heatmap_Podo.png", width = 200, height = 150, res = 300, units = "mm")
pheatmap(mat.z,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "", fontsize = 12, color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Add cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
integrated.cca <- CellCycleScoring(integrated.cca, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(integrated.cca, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
integrated.cca <- RunPCA(integrated.cca, features = c(s.genes, g2m.genes))
DimPlot(integrated.cca)
Idents(integrated.cca) <- "Phase"
DimPlot(integrated.cca, reduction = "tsne", label = F, pt.size = 0.6)

# volcano plot DESeq2
library(ggplot2)

# Define significance and log2 fold change thresholds, for example:
threshold_pvalue <- 0.05
threshold_lfc <- 0.5

# Create a significance column
res$significant <- ifelse(res$padj < threshold_pvalue & abs(res$log2FoldChange) > threshold_lfc, "yes", "no")
res$category <- with(res, ifelse(padj < threshold_pvalue & log2FoldChange > threshold_lfc, "up", 
                                     ifelse(padj < threshold_pvalue & log2FoldChange < -threshold_lfc, "down", "non-significant")))

res <- as.data.frame(res)
res$SYMBOL <- rownames(res)
# Draw the volcano plot
volcano <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=category)) +
  geom_point(alpha=0.6, size=1.5) +
  theme_minimal() +
  theme(legend.position="top") +
  scale_color_manual(values=c("up"="red", "down"="blue", "non-significant"="gray")) +
  labs(title="Volcano Plot",
       x="Log2 Fold Change",
       y="-Log10 P-value",
       color="Volcano CTX vs CTRL") +
  geom_hline(yintercept=-log10(threshold_pvalue), linetype="dashed", color = "blue") +
  geom_vline(xintercept=c(-threshold_lfc, threshold_lfc), linetype="dashed", color = "blue")

# Calculate the counts of upregulated and downregulated genes
upregulated_count <- sum(res$padj < threshold_pvalue & res$log2FoldChange > 0.5)
downregulated_count <- sum(res$padj < threshold_pvalue & res$log2FoldChange < -0.5)


volcano <- volcano + labs(title = paste0("Upregulated: ", upregulated_count, " & Downregulated: ", downregulated_count))
# Print the plot
print(volcano)

# Enrichment
enrich <- read.table("A2_WTvsMT/Enichment.txt", sep = "\t", header = T)
enrich$log10_Adjusted <- -log10(enrich$Adjusted.P.value)
p <- ggplot(enrich, aes(x = log10_Adjusted, y = reorder(Term, log10_Adjusted), fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Pathways" = "darkred", "GO-BP" = "lightblue", "GO-CC" = "navy")) +
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

png(filename = "Enrichment_Podo.2.png", width = 250, height = 150, res = 300, units = "mm")
plot(p)
dev.off()

# Heatmap for enrichment
# Step 1: Separate genes for each ontology term
# Split Genes column and create a long-format data frame with one gene per row
library(tidyverse)
library(pheatmap)
df_long <- enrich %>%
  separate_rows(Genes, sep = ";") %>%
  distinct()

# Step 2: Create a binary matrix for ontology terms and genes
# Spread the data to have a binary matrix of Terms vs. Genes
binary_matrix <- df_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Genes, values_from = value, values_fill = 0) %>%
  column_to_rownames("Term")

binary_matrix[is.na(binary_matrix)] <- 0
binary_matrix <- binary_matrix[,-c(1:7)]
gene_annotation <- data.frame(Ontology = rownames(binary_matrix))
rownames(gene_annotation) <- rownames(binary_matrix)

# Convert binary matrix to long format for ggplot
binary_long <- as.data.frame(binary_matrix) %>%
  rownames_to_column("Term") %>%
  pivot_longer(-Term, names_to = "Gene", values_to = "Presence")

# Step 3: Plot heatmap with ggplot2, adding borders and custom colors
p <- ggplot(binary_long, aes(x = Gene, y = Term, fill = as.factor(Presence))) +
  geom_tile(color = "black", size = 0.5) +  # Add border color and adjust thickness
  scale_fill_manual(values = c("0" = "white", "1" = "darkred")) +  # Custom color for binary values
  labs(title = "Gene Presence by Ontology Term", x = "Genes", y = "Ontology") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(colour = "black",angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 10),
    axis.text.y = element_text(colour = "black", face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()  # Remove grid lines to emphasize borders
  )

png(filename = "Erichment.heatmap.png", width = 450, height = 100, res = 300, units = "mm")
plot(p)
dev.off()


library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

res <- res %>%
  filter(padj< 0.05)
original_gene_list <- res$Symbol

gene_sets_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_sets_GO_BP <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets_GO_CC <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
gene_sets_GO_MF <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
gene_sets_KEGG <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
gene_sets_Reactome <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")

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

CC <- enricher(gene = original_gene_list, TERM2GENE = gene_sets_GO_CC, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
BP <- enricher(gene = original_gene_list, TERM2GENE = gene_sets_GO_BP, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
MF <- enricher(gene = original_gene_list, TERM2GENE = gene_sets_GO_MF, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
KEGG <- enricher(gene = original_gene_list, TERM2GENE = gene_sets_KEGG, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
React <- enricher(gene = original_gene_list, TERM2GENE = gene_sets_Reactome, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

#write.table(SMZ@result, file = "RNASeqProject_Mehrdad/DESeq2_GO_CC_TMPvsCTRL.txt", sep = "\t")
#write.table(SMZ@result, file = "RNASeqProject_Mehrdad/DESeq2_GO_BP_TMPvsCTRL.txt", sep = "\t")
#write.table(SMZ@result, file = "RNASeqProject_Mehrdad/DESeq2_GO_MF_SMZvsCTRL.txt", sep = "\t")
#write.table(SMZ@result, file = "RNASeqProject_Mehrdad/DESeq2_KEGG_SMZvsCTRL.txt", sep = "\t")
write.table(SMZ@result, file = "RNASeqProject_Mehrdad/DESeq2_reactome_TMPvsCTRL.txt", sep = "\t")

# curate GeneRatio information

SMZ <- as.data.frame(SMZ@result)
split_values <- strsplit(SMZ$GeneRatio, "/")
numerator <- sapply(split_values, function(x) as.numeric(x[1]))
denominator <- sapply(split_values, function(x) as.numeric(x[2]))
SMZ$GR <- numerator / denominator
SMZ <- SMZ %>% arrange(p.adjust)
SMZ <- SMZ[c(1:10),]
SMZ <- SMZ %>% arrange(Count)
SMZ$ID <- factor(SMZ$ID, levels = SMZ$ID)

ggplot(SMZ, aes(x=GR, y= ID, color=p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("GeneRatio") + 
  theme(axis.text.y = element_text(face = "bold", size = 12), axis.text.x = element_text(face = "bold", size = 12)) +
  ggtitle("GO_CC enrichment: SMZ vs CTRL") + guides(colour = guide_colourbar(order = 1))


# GSEA analysis using ClusterProfiler

DEseq <- read.table("RNASeqProject_Mehrdad/3_DEG_analysis_limma_edgeR_DESeq2/DESeq2_CTXvsCTRL_result_NotFiltered.txt", sep = "\t", header = T)

# we want the log2 fold change from DESeq2

original_gene_list <- res$log2FoldChange

# name the vector
names(original_gene_list) <- res$Symbol
# omit any NA values 
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
keytypes(org.Hs.eg.db)

gse <- gseGO(geneList= gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = FALSE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")


write.table(gse@result, file = "GSEA_merged_DMSOvsCTRL_GO_CC.txt", sep = "\t")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
ridgeplot(gse) + labs(x = "enrichment distribution")
# create GSEA plot for specific gene set
# note: 'result' variable is undefined - should be gse@result
gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)  # corrected variable reference

# KEGG Enrichment
#https://www.genome.jp/kegg/catalog/org_list.html
# convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = res[res$Symbol %in% dedup_ids$SYMBOL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


write.table(kk2@result, file = "GSEA_merged_DMSOvsCTRL_KEGG.txt", sep = "\t")

# DOT PLOT VISUALIZATION
# create combined group and segment annotation
integrated.cca$groupSegment <- paste(integrated.cca@meta.data$group2, integrated.cca@meta.data$segments, sep = ".")  # corrected group reference
Idents(integrated.cca) <- "group2"  # corrected identity
DefaultAssay(integrated.cca) <- "SCT"

# create dot plot for laminin and collagen genes
p <- DotPlot(integrated.cca, features=c("LAMB1", "LAMB2", "LAMA5","LAMA1", "COL4A6", "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1"),
             dot.scale=8, group.by ="group2", cols="RdBu") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1,face = "italic"),  # corrected group.by
                                                                       axis.text.y = element_text( vjust = 0.5, hjust=1, face = "bold"), axis.title = element_blank())

# create dot plot for combined group and segment
p <- DotPlot(integrated.cca, features=c("LAMB1", "LAMB2", "LAMA5","LAMA1", "COL4A6", "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1"),
             dot.scale=8, group.by ="groupSegment", cols="RdBu") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1,face = "italic"),
                                                                  axis.text.y = element_text( vjust = 0.5, hjust=1, face = "bold"), axis.title = element_blank())


# create dot plot for collagen and laminin genes by group
p <- DotPlot(integrated.cca, features=c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6", "LAMB1", "LAMA1", "LAMB2", "LAMA5"),
             dot.scale=10, group.by ="group2", cols="RdBu") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1,face = "italic", size = 14),  # corrected group.by
                                                                        axis.text.y = element_text( vjust = 0.5, hjust=1,face = "bold", size = 14), axis.title = element_blank())

# create dot plot for cell type markers
p <- DotPlot(integrated.cca, features=c("CENPF", "HMGB2", "UBE2C","HIST1H4C", "PCLAF",
                                        "TYMS","LAMC3", "NKD1", "TNC", "PHACTR3", "PDGFRB",
                                        "COL4A1", "COL1A2", "POSTN", "COL3A1", "COL6A3",
                                        "COL14A1", "OGN", "GREM2", "FRZB","CLDN5", "PECAM1", "KDR"),
             dot.scale=10, group.by ="CellType_v3", cols="RdBu") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1, size = 14),  # corrected group.by
                                                                         axis.text.y = element_text( vjust = 0.5, hjust=1,face = "italic", size = 14), axis.title = element_blank())
png(filename = "DotPlot_mutvsctrl.2.png", width = 300,height = 200, units = "mm", res = 300)
plot(p)
dev.off()

# VOLCANO PLOT FOR DIFFERENTIAL EXPRESSION ANALYSIS
library(data.table)
DEG <- fread("A2_WTvsMT/DEG.podocyte_DESeq2.txt")  # load differential expression results
DEG <- DEG[DEG$significant == "yes",]  # filter for significant genes

library(ggrepel)  # for text labels in plots
DEG <- DEG %>%
  mutate(
    highlight = case_when(
      padj < 0.05 & log2FoldChange > 0.2 ~ "Up",      # upregulated genes
      padj < 0.05 & log2FoldChange < -0.2 ~ "Down",   # downregulated genes
      TRUE ~ "Not_sig"                                # non-significant genes
    )
  )

color_map <- c(
  "Up" = "darkred",
  "Not_sig" = "azure4",
  "Down" = "blue"
)

DEG$log10Pval <- -log10(DEG$padj)

volcano <- ggplot(DEG,
                  aes(x = log2FoldChange, y = log10Pval, colour= highlight)) +
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

# MA PLOT FOR DIFFERENTIAL EXPRESSION VISUALIZATION
data <- DEG  # use the DEG data

# define significance and fold change thresholds
padj_threshold <- 0.05
log2fc_bold_threshold <- 1

# filter data for significant genes
significant_data <- subset(data, padj < padj_threshold)
significant_data$color <- ifelse(significant_data$log2FoldChange > 0, "darkred", "darkblue")  # color by direction

# identify genes with high fold change for labeling
label_data <- subset(significant_data, abs(log2FoldChange) > 2)

# create MA plot (mean expression vs log2 fold change)
ggplot(significant_data, aes(x = baseMean, y = log2FoldChange, label = SYMBOL)) +
  geom_point(aes(color = color), size = 3, alpha = 0.7) +
  scale_x_log10() +  # log scale for mean expression
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # zero line
  geom_text_repel(data = label_data, aes(label = SYMBOL), fontface = "bold", size = 3, 
                  box.padding = 0.5, max.overlaps = Inf) +  # label high fold change genes
  labs(x = "Mean Expression (baseMean)", y = "log2 Fold Change",
       title = "MA Plot of Significant Differentially Expressed Genes (Dark Red for Up, Dark Blue for Down)") +
  theme_minimal() +
  scale_color_identity()


# Make shineyCell ui for this dataset
library(ShinyCell)
DefaultAssay(integrated.cca) <- "RNA"

integrated.cca <- NormalizeData(object = integrated.cca,  verbose = FALSE)
integrated.cca <- FindVariableFeatures(object = integrated.cca, nfeatures = 3500, verbose = FALSE, selection.method = 'vst')
integrated.cca <- ScaleData(object = integrated.cca)


integrated.cca$CellType_v3 <- factor(integrated.cca$CellType_v3)
integrated.cca$predicted_label <- factor(integrated.cca$predicted_label)
integrated.cca$CellType_v2 <- factor(integrated.cca$CellType_v2)
integrated.cca$Annotation_v1 <- factor(integrated.cca$Annotation_v1)
integrated.cca$segments <- factor(integrated.cca$segments)

unique(Idents(integrated.cca))

scConf = createConfig(integrated.cca)
integrated.cca <- UpdateSeuratObject(integrated.cca)
makeShinyApp(integrated.cca, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell Dataset X-linked Alport syndrome Organoids")

# save the final integrated Seurat object with all annotations
saveRDS(integrated.cca, file = "integrated.sct.cca.a2.v2.rds")
sessionInfo()
