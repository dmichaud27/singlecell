# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_filtered.RData")

# Normalize the counts
seurat_normalized <- NormalizeData(filtered_seurat)

# Identify the most variable genes
seurat_normalized <- FindVariableFeatures(seurat_normalized, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_normalized <- ScaleData(seurat_normalized)

# Check quartile values
summary(seurat_normalized@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_normalized@meta.data$mitoFr <- cut(seurat_normalized@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.03707, 0.06271, 0.10069, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

# Run PCA
seurat_normalized <- RunPCA(seurat_normalized, features = VariableFeatures(seurat_normalized), nfeatures.print = 10, ndims.print = 1:10)
pdf("PCA_seurat_normalized.pdf")
PCAPlot(seurat_normalized, group.by = "sample")
dev.off()

# Load cell cycle markers
load("cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_normalized, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes,
                                 set.ident = TRUE)

# Rescaling on regressed mito Ratio
seurat_phase <- ScaleData(seurat_phase, vars.to.regress = "mitoRatio")

# Test regression of PCA
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), nfeatures.print = 10, ndims.print = 1:10)
pdf("PCA_seurat_phase_postregression.pdf")
DimPlot(seurat_phase, group.by = "sample")
DimPlot(seurat_phase, split.by = "Phase", group.by = "sample")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "sample")
DimPlot(seurat_phase, group.by = "subtype")
DimPlot(seurat_phase, split.by = "Phase", group.by = "subtype")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "subtype")
dev.off()

# Run tSNE
seurat_phase <- RunTSNE(seurat_phase, 
                        dims = 1:40,
                        reduction = "pca")
                        
pdf("tsne_seurat_phase_postregression.pdf")                        
DimPlot(seurat_phase, group.by = "sample")
DimPlot(seurat_phase, split.by = "Phase", group.by = "sample")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "sample")
DimPlot(seurat_phase, group.by = "subtype")
DimPlot(seurat_phase, split.by = "Phase", group.by = "subtype")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "subtype")
dev.off()

# Run UMAP
seurat_phase <- RunUMAP(seurat_phase, 
                             dims = 1:40,
                             reduction = "pca")

pdf("umap_seurat_phase_postregression.pdf")                        
DimPlot(seurat_phase, group.by = "sample")
DimPlot(seurat_phase, split.by = "Phase", group.by = "sample")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "sample")
DimPlot(seurat_phase, group.by = "subtype")
DimPlot(seurat_phase, split.by = "Phase", group.by = "subtype")
DimPlot(seurat_phase, split.by = "mitoFr", group.by = "subtype")
dev.off()

# Save seurat object
save(seurat_phase, file = "seurat_phase.RData")