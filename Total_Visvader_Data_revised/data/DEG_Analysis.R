# Single-cell RNA-seq - DE Gene Analysis

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_subFinal2.RData")

# Marker ID Setup
DefaultAssay(seurat_subFinal) <- "RNA"

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_subFinal,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

write.csv(markers, "markersFinal.csv")

# Make clustered Heatmap
dim(markers)
head(markers)
table(table(markers$gene))
markers_single <- markers[markers$gene %in% names(table(markers$gene))[table(markers$gene) == 1],]
dim(markers_single)
table(table(markers_single$gene))
table(markers_single$cluster)
top50 <- markers_single %>% group_by(cluster) %>% top_n(50, avg_log2FC)
dim(top50)
pdf("Heatmap_seurat_subFinal.pdf")
DoHeatmap(
  object = seurat_sub, 
  features = top50$gene)
dev.off()

# Extract top 50 markers per cluster
grouped_markers <- grouped_df(markers, "cluster")
top50FC <- slice_max(.data = grouped_markers, order_by = avg_log2FC, n = 50)
write.csv(top50FC, "top50.csv")

# Find Markers defining subtypes
Idents(object = seurat_subFinal) <- "subtype"
markersHR <- FindMarkers(seurat_subFinal, ident.1 = "HR", ident.2 = c("HER2", "TNBC"))
markersHR <- rownames_to_column(markersHR, var = "gene")
write.csv(markersHR, "markersHR.csv")

markersHER2 <- FindMarkers(seurat_subFinal, ident.1 = "HER2", ident.2 = c("HR", "TNBC"))
markersHER2 <- rownames_to_column(markersHER2, var = "gene")
write.csv(markersHER2, "markersHER2.csv")

markersTNBC <- FindMarkers(seurat_subFinal, ident.1 = "TNBC", ident.2 = c("HR", "HER2"))
markersTNBC <- rownames_to_column(markersTNBC, var = "gene")
write.csv(markersTNBC, "markersTNBC.csv")

Idents(object = seurat_subFinal) <- "RNA_snn_res.1.9"

# Plasma cell plots
pdf("Plasmacell_plots.pdf")
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD38", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "SDC1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "TNFRSF17", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "MS4A1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD27", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CXCR4", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IGHD", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "PRDM1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IRF4", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "XBP1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD79A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD19", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# Memory B plots
pdf("MemoryB_plots.pdf")
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD80", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD86", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CR2", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD40", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "FAS", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "TNFRSF13B", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "PTPRJ", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "MS4A1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD27", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "PAX5", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "POU2AF1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "SPIB", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# Regulatory B plots
pdf("RegulatoryB_plots.pdf")
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IL10", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IL12A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "EBI3", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD5", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD1D", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD24", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD40", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IGHM", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CD38", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "CR2", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# Pre-chemo B plots
pdf("Pre_chemo_signature_plots.pdf")
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "ISG15", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "MX1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IFI6", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "LY6E", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "WARS", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "GZMB", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "MKI67", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "IL10", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "HMGB2", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_subFinal, reduction = "tsne", features = "HIST1H4C", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()
