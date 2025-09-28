################################################################################
# Single-cell RNA-seq Analysis: Early vs Late Kidney Organoids
# Author: Hassan Saei
# Contact: hassan.saeiahan@gmail.com
# Purpose: Compare gene expression patterns between day 22 (early) and day 38 (late) organoids
# Integrated Seurat object with cell type annotations (.h5ad file is available in the Zenodo repository)
################################################################################

# loading libraries
library(dplyr)
library(Seurat)
library(reticulate)
library(anndata)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(SoupX)
library(openxlsx)
library(scCustomize)
library(dittoSeq)

options(Seurat.object.assay.version = "v3")

################################################################################
# load pre-processed data list for integration
data_list <- readRDS(file = "NGS2024_7386/2_Integrated/Data_integrated.sct.cca.rds")

# select integration features for CCA integration (3500 most variable features)
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3500)

# prepare data for SCT integration
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)

# find integration anchors using CCA method
data.anchors.cca <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = features, reduction = "cca")

# integrate datasets using SCT normalization
integrated.cca <- IntegrateData(anchorset = data.anchors.cca, normalization.method = "SCT")

# set integrated assay as default and perform dimensionality reduction
DefaultAssay(integrated.cca) <- "integrated"
integrated.cca <- RunPCA(integrated.cca, npcs = 50, verbose = FALSE)
integrated.cca <- FindNeighbors(integrated.cca, dims = 1:50, reduction = "pca", k.param = 20)
# find clusters at two different resolutions for comparison
integrated.cca <- FindClusters(integrated.cca, resolution = 0.5, cluster.name = "cca_integrated_res0.5")
integrated.cca <- FindClusters(integrated.cca, resolution = 1.2, cluster.name = "cca_integrated_res1.2")
# generate UMAP and t-SNE embeddings
integrated.cca <- RunUMAP(integrated.cca, dims = 1:50, reduction = "pca", reduction.name = "umap.cca.res")
integrated.cca <- RunTSNE(integrated.cca, reduction = "pca", dims = 1:50)

# Load the final integrated object (this overwrites the above integration steps)
integrated.cca <- readRDS("integrated.sct.cca.a1.v2.rds")

# cell-type annotation and renaming
# rename cell types for consistency and clarity
# original names contain marker information in parentheses
integrated.cca@meta.data$CellType_V3 <- integrated.cca@meta.data$CellType_V2
integrated.cca@meta.data$CellType_V3 <- as.character(integrated.cca@meta.data$CellType_V3)
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Myofibroblast_I (ACTA2 +)"] <- "Myofibroblast_I"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Nephron/Podocyte progenitors (CTGF+|DAPL1+)"] <- "Podocyte progenitors"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Proliferative stroma_II"] <- "Proliferative stroma"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Fibroblast_I (PDGFRA+)"] <- "Fibroblast_I"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Proliferative epithelium progenitor"] <- "Proliferative epithelium"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "LOH-like (SLC12A1+)"] <- "LOH-like"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Fibroblast_II (COL14A1+)"] <- "Fibroblast_II"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Glial (FABP7+)"] <- "Glial"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Proliferative stroma_I"] <- "Proliferative stroma"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "PT progenitors (PEC-like)"] <- "PEC-like"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "DIFF_fibroblasts (OGN+)"] <- "Fibroblast_III"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Distal tubules (GATA3+)"] <- "DCT"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "PT (HNF4A+|HNF4G+)"] <- "PT"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "DEV_PT_II (SLC3A1+)"] <- "DEV_PT_II"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Myofibroblast_I (ACTA2+)"] <- "Myofibroblast_I"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "DEV_podocyte progenitors"] <- "Podocyte progenitors"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Endothelial cells"] <- "Endothelial progenitors"
integrated.cca@meta.data$CellType_V3[integrated.cca@meta.data$CellType_V3 == "Glial (TTYH1+)"] <- "Glial"

Idents(integrated.cca) <- "CellType_V3"
p <- DimPlot_scCustom(integrated.cca, reduction = "umap.cca.res", label = T, pt.size = 0.7,
                      label.size = 4, colors_use = "stepped", split_seurat = TRUE)+ NoLegend()

# save annotated UMAP
png(filename = "uMAP.annotation.WT_vs_MT.2.png", width = 250, height = 200, res = 300, units = "mm")
plot(p)
dev.off()

# find marker genes in each cluster
integrated.cca.markers <- FindAllMarkers(integrated.cca, only.pos = TRUE)
integrated.cca.markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.table(integrated.cca.markers, file = "Markers_CellTypes.txt", sep = "\t", row.names = F)


# visualization and marker analysis
# set SCT assay as default for visualization
DefaultAssay(integrated.cca) <- "SCT"

# podocyte markers - key markers for mature podocytes
VlnPlot(integrated.cca, features = c("NPHS1", "NPHS2", "DDN", "WT1", "PODXL"), ncol = 2)
FeaturePlot(integrated.cca, features = c("NPHS1", "NPHS2", "PODXL", "DDN"), cols = c("lightgray", "darkblue"))
# collagen IV isoforms - basement membrane components
VlnPlot(integrated.cca, features = c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6"), ncol = 2)
FeaturePlot(integrated.cca, features = c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6"), cols = c("lightgray", "darkblue"))

# podocyte precursor markers - early podocyte development
VlnPlot(integrated.cca, features = c("CTGF", "OLFM3", "MAFB", "NPHS1"))
FeaturePlot(integrated.cca, features = c("CTGF", "OLFM3", "MAFB", "NPHS1"), cols = c("lightgray", "darkblue"))

# nephron progenitors - SIX1/SIX2 positive progenitor cells
VlnPlot(integrated.cca, features = c("DAPL1", "LYPD1", "SIX1", "SIX2","CRABP2"))
FeaturePlot(integrated.cca, features = c("DAPL1", "LYPD1", "SIX1", "SIX2","CRABP2"), cols = c("lightgray", "darkblue"))

# epithelial markers - general epithelial cell markers
VlnPlot(integrated.cca, features = c("PAX2", "PAX8", "KRT19", "EPCAM","LRP2"), ncol = 4)
FeaturePlot(integrated.cca, features = c("PAX2", "PAX8", "KRT19", "EPCAM","LRP2"), cols = c("lightgray", "darkblue"))

# loop of Henle (LOH) markers - distal nephron segment
VlnPlot(integrated.cca, features = c("ESRRG", "POU3F3", "SLC12A1"), ncol = 3)
FeaturePlot(integrated.cca, features = c("ESRRG", "POU3F3", "SLC12A1"), cols = c("lightgray", "darkblue"))

# distal progenitor markers - developing distal tubule
VlnPlot(integrated.cca, features = c("EPCAM", "EMX2", "SPP1", "MAL", "PAX2", "GATA3", "TFAP2A"), ncol = 4)
FeaturePlot(integrated.cca, features = c("EPCAM", "EMX2", "SPP1", "MAL", "PAX2", "GATA3", "TFAP2A"), cols = c("lightgray", "darkblue"))

# proximal tubule precursor markers - developing proximal tubule
VlnPlot(integrated.cca, features = c("IGFBP7", "FXYD2", "CDH6", "HNF1B", "HNF4G", "HNF4A", "SLC3A1"), ncol = 4)
FeaturePlot(integrated.cca, features = c("IGFBP7", "FXYD2", "CDH6", "HNF1B", "HNF4G", "HNF4A", "SLC3A1"),
            cols = c("lightgray", "darkblue"), max.cutoff = 1.5, min.cutoff = 0.5)

# mesangial cell markers - glomerular mesangial cells
VlnPlot(integrated.cca, features = c("HOPX"), ncol = 1)
FeaturePlot(integrated.cca, features = c("HOPX", "AGTR1"), cols = c("lightgray", "darkblue"))
FeaturePlot(integrated.cca, features = c("HOPX", "PDGFRB", "COL1A1", "ACTA2", "ITGA8"), cols = c("lightgray", "darkblue"))

# claudins - tight junction proteins
VlnPlot(integrated.cca, features = c("CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN7", "CLDN8"), ncol = 4)
FeaturePlot(integrated.cca, features = c("CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN7", "CLDN8"), cols = c("lightgray", "darkblue"))

# muscle progenitor markers - myogenic differentiation
VlnPlot(integrated.cca, features = c("MYOG", "MYOD1"))
FeaturePlot(integrated.cca, features = c("MYOG", "MYOD1", "SIX1", "SIX2"), cols =c("lightgray", "darkblue"))

# neural progenitor markers - neurogenic differentiation
VlnPlot(integrated.cca, features = c("HES6", "STMN2"))
FeaturePlot(integrated.cca, features = c("HES6", "STMN2"), cols = c("lightgray", "darkblue"))

# glial markers - glial cell differentiation
VlnPlot(integrated.cca, features = c("FABP7", "TTYH1","SOX2"))
FeaturePlot(integrated.cca, features = c("FABP7", "TTYH1","SOX2"), cols = c("lightgray", "darkblue"))

# cell cycle markers - mitotic cell cycle process
VlnPlot(integrated.cca, features = c("CENPF", "HMGB2", "UBE2C","HIST1H4C", "PCLAF", "TYMS"))
FeaturePlot(integrated.cca, features = c("CENPF", "HMGB2", "UBE2C","HIST1H4C", "PCLAF", "TYMS"), cols = c("lightgray", "darkblue"))

# mesangial cell markers - glomerular mesangial cells
VlnPlot(integrated.cca, features = c("LAMC3", "NKD1", "TNC", "PHACTR3", "PDGFRB", "HOPX"))
FeaturePlot(integrated.cca, features = c("LAMC3", "NKD1", "TNC", "PHACTR3", "PDGFRB"), cols = c("lightgray", "darkblue"))

# myofibroblast markers - activated fibroblasts with contractile properties
VlnPlot(integrated.cca, features = c("COL4A1", "COL1A2", "POSTN", "COL3A1", "COL6A3", "COL14A1", "OGN", "GREM2", "FRZB"))
FeaturePlot(integrated.cca, features = c("COL4A1", "COL1A2", "POSTN", "COL3A1", "COL6A3", "COL14A1", "OGN", "GREM2", "FRZB"), cols = c("lightgray", "darkblue"))
FeaturePlot(integrated.cca, features = c("ACTA2", "POSTN", "COL1A1", "COL1A2"), cols = c("lightgray", "darkblue"), min.cutoff = 2.5, max.cutoff = 4)

# endothelial markers - vascular endothelial cells
VlnPlot(integrated.cca, features = c("CLDN5", "PECAM1", "KDR", "GNG11", "CALM1"))


#polychrome_pal <- DiscretePalette_scCustomize(num_colors = 25, palette = "polychrome")
#p <- DimPlot_scCustom(integrated.cca, reduction = "umap.cca.res1.2", label = F,
#                      pt.size = 0.4, label.size = 2.5, colors_use = polychrome_pal, repel = T, label.box = T) + NoLegend()


# threshold-based cell highlighting
# highlight cells expressing COL4 genes above specific thresholds
# thresholds determined based on expression distribution analysis
thresholds <- c(COL4A1=1, COL4A2=1, COL4A3=0.5, COL4A4=0.5, COL4A5=1, COL4A6=1)

# build list of cells above threshold for each gene
cells_list <- lapply(names(thresholds), function(g){
  vals <- FetchData(integrated.cca, vars = g)[,1]
  colnames(integrated.cca)[vals > thresholds[g]]
})
names(cells_list) <- names(thresholds)

pdf("UMAP_highlight_COL4_all.pdf", width = 6, height = 5)  # size as you like
for(g in names(cells_list)){
  p <- Cell_Highlight_Plot(
    seurat_object = integrated.cca,
    cells_highlight = list(g = cells_list[[g]]),
    highlight_color = "navy",
    pt.size = 0.1
  ) + ggtitle(g)
  print(p)
}
dev.off()


# cluster annotation and grouping
# define main cluster groups based on biological function
# cluster numbers correspond to Seurat clusters at resolution 1.2
cluster_mapping <- list(
  'Tubules' = c(5,14,3,11,23,16,21,19,22,17,31,28,32),
  'Podocytes' = c(1,6,26,8),
  'Off-target' = c(24,13,12,30,27,10),
  'Stroma' = c(0,2,4,15,25,9,18,7,20, 29)
)

# create an empty vector to store the cluster names
integrated.cca@meta.data$Mainclusters <- NA

# assign cluster names based on the cluster numbers
for (cluster_name in names(cluster_mapping)) {
  cluster_numbers <- cluster_mapping[[cluster_name]]
  integrated.cca@meta.data$Mainclusters[
    integrated.cca@meta.data$cca_integrated_res1.2 %in% cluster_numbers
  ] <- cluster_name
}

# highlight some clusters
Idents(integrated.cca) <- "Mainclusters"
p <- Cluster_Highlight_Plot(seurat_object = integrated.cca, cluster_name = c('Tubules', 'Stroma', 'Off-target', 'Podocytes'),
                            highlight_color = c("darkred","forestgreen", "lightgray", "darkblue"))

png(filename = "uMAP_highlight_mainclusters.png", width = 150, height = 120, units = "mm", res = 300)
p
dev.off()

# time point analysis
# define time point groups based on organoid culture days
# early: day 25, Late: day 38
cluster_mapping <- list(
  'early' = c("05_d25"),
  'late' = c("05_d38")
)

# create an empty vector to store the cluster names
integrated.cca@meta.data$group <- NA

# assign cluster names based on the cluster numbers
for (cluster_name in names(cluster_mapping)) {
  cluster_numbers <- cluster_mapping[[cluster_name]]
  integrated.cca@meta.data$group[
    integrated.cca@meta.data$orig.ident %in% cluster_numbers
  ] <- cluster_name
}

# coparitive analysis and VISUALvisualization
# generate stacked bar plots showing cell type composition by time point
p <- dittoBarPlot(integrated.cca, "CellType_V3", group.by = "group", ylab = "Percentage of cell types")
p <- dittoBarPlot(integrated.cca, "Mainclusters", group.by = "group", ylab = "Percentage of cell types")

# create combined day-cluster identifier for comparative analysis
integrated.cca$DayCluster <- paste(integrated.cca@meta.data$Mainclusters, integrated.cca@meta.data$group, sep = " ")
Idents(integrated.cca) <- "DayCluster"
DefaultAssay(integrated.cca) <- "SCT"
p <- DotPlot(integrated.cca, features=c("LAMB1", "LAMB2", "LAMA5","LAMA1", "COL4A6", "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1"),
        dot.scale=10, group.by ="DayCluster", cols="RdBu") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1,face = "bold"),
                                                                axis.text.y = element_text( vjust = 0.5, hjust=1, face = "italic"), axis.title = element_blank()) + coord_flip()

Clustered_DotPlot(seurat_object = integrated.cca, features = c("LAMB1", "LAMB2", "LAMA5","LAMA1", "COL4A6", "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1"),
                  k = 8,
                  show_ident_legend = FALSE, show_ident_colors = F, cluster_ident = T, flip = F,
                  colors_use_exp = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")), x_lab_rotate = 90, row_label_fontface = "italic")




# final annotation dotplot
# generate comprehensive dotplot showing key markers for all cell types
Idents(integrated.cca) <- "CellType_V3"

# define comprehensive marker gene set for cell type annotation
# includes markers for: podocytes, nephron progenitors, epithelial cells, 
# stromal cells, muscle/neural progenitors, and endothelial cells
features <- c("NPHS1", "NPHS2", "DDN", "DAPL1","LYPD1", "CLDN1", "PAX8", "PAX2","SLC3A1",
           "LRP2" ,"SLC6A13","APOE","HNF4A", "HNF4G", "CDH17","EPCAM", "KRT19","CDH1","CXCL3", "MAL", "GATA3", "CALB1","DEPDC1B","UBE2C", "SCLT1","COL1A1","ACTA2","COL3A1", "SERPINE1", "TNC","FN1","LOX","ITGA8",
           "PDGFRL", "SFRP2","COL14A1", "PLAT", "PDGFRB", "FRZB",
           "MYOG","HES6", "PECAM1", "MSX1")

integrated.cca@meta.data$CellType_V3 <- factor(integrated.cca@meta.data$CellType_V3, 
                                               levels=c("Podocytes",
                                                        "Podocyte progenitors",
                                                        "PEC-like",
                                                        "PT",
                                                        "DEV_PT_I",
                                                        "DEV_PT_II",
                                                        "DEV_TUB",
                                                        "Proliferative epithelium",
                                                        "LOH-like",
                                                        "DCT",
                                                        "Neural progenitors",
                                                        "Glial",
                                                        "Muscle progenitors",
                                                        "Proliferative stroma",
                                                        "Myofibroblast_I",
                                                        "Myofibroblast_II",
                                                        "Fibroblast_I",
                                                        "Fibroblast_II",
                                                        "Fibroblast_III", 
                                                        "Endothelial progenitors"
                                                        ))

png(filename = "DotPlot_cell_markers.update.png", width = 200,height = 200, units = "mm", res = 300)
Clustered_DotPlot(seurat_object = integrated.cca, features = features, k = 8,
                  show_ident_legend = FALSE, show_ident_colors = F, cluster_ident = T, flip = F,
                  colors_use_exp = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")), x_lab_rotate = 90, row_label_fontface = "italic")
dev.off()

# Make ShinyCell ui for this datase
library(Seurat)
library(ShinyCell)

DefaultAssay(integrated.cca) <- "RNA"

integrated.cca <- NormalizeData(object = integrated.cca,  verbose = FALSE)
integrated.cca <- FindVariableFeatures(object = integrated.cca, nfeatures = 3500, verbose = FALSE, selection.method = 'vst')
integrated.cca <- ScaleData(object = integrated.cca)

integrated.cca$CellType_V1 <- factor(integrated.cca$CellType_V1)
integrated.cca$CellType_V2 <- factor(integrated.cca$CellType_V2)
integrated.cca$CellType_V3 <- factor(integrated.cca$CellType_V3)

Idents(integrated.cca) <- "CellType_V3"
unique(Idents(integrated.cca))
typeof(integrated.cca$CellType_V3)

scConf = createConfig(integrated.cca)

integrated.cca <- UpdateSeuratObject(integrated.cca)
makeShinyApp(integrated.cca, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell Dataset Kidney Organoids (Early versus Late)")

# Save final object
saveRDS(integrated.cca, file = "integrated.sct.cca.a1.v2.rds")
sessionInfo()