# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_clustered2.RData")

# Major TME populations
pdf("Immune.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "PTPRC", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

pdf("Fibroblasts.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "ACTA2", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "PDPN", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "PDGFRA", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

pdf("Vasculature.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "PECAM1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "LYVE1", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

pdf("Epithelial.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "KRT8",  order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "KRT18", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "KRT19", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "KRT17", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "KRT5",  order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "EPCAM", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# Monocytes/Macrophages
pdf("Monocytes_Macrophages.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD14", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "LYZ", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "FCGR3A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "MS4A7", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "ITGAM", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CX3CR1", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# cDCs
pdf("cDCs.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "FCER1A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CST3", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "ITGAX", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# pDCs
pdf("pDCs.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IL3RA", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "GZMB", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "SERPINF1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "ITM2C", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "ITGAX", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# B cells
pdf("Bcells.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "MS4A1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD19", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD79A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD79B", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "SDC1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "JCHAIN", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "PRDM1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IRF4", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHD", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHM", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHG1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHG2", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHG3", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHG4", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHA1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHA2", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "IGHE", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# T cells
pdf("Tcells.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD3E", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD3D", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD3G", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD4", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD8A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD8B", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "FOXP3", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()

# NK cells
pdf("NKcells.pdf")
FeaturePlot(seurat_clustered, reduction = "tsne", features = "NCAM1", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "CD3E", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "FCGR3A", order = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(seurat_clustered, reduction = "tsne", features = "NCR1", order = TRUE, min.cutoff = 'q10', label = TRUE)
dev.off()