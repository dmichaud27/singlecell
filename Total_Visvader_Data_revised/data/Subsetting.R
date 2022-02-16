# Single-cell RNA-seq - Subsetting & Reanalysis

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_clustered2.RData")

# Subset data to get rid of unwanted clusters
seurat_sub <- subset(seurat_clustered, idents = c(13, 29))
levels(seurat_sub)

# Count number of cells in existing clusters
n_cells_sub <- FetchData(seurat_sub, 
                     vars = c("ident", "orig.ident")) %>%
dplyr::count(ident, orig.ident) %>%
tidyr::spread(ident, n)

write.csv(n_cells_sub, file = "n_cells_sub.csv")

# Save subsetted data
save(seurat_sub, file = "seurat_sub.RData")

# Reanalyze data
seurat_sub <- NormalizeData(seurat_sub)

seurat_sub <- FindVariableFeatures(seurat_sub, 
                                         selection.method = "vst",
                                         nfeatures = 2000, 
                                         verbose = FALSE)

# Scale subset data
seurat_sub <- ScaleData(seurat_sub, vars.to.regress = "mitoRatio")
seurat_sub <- RunPCA(seurat_sub, features = VariableFeatures(seurat_sub), nfeatures.print = 10, ndims.print = 1:10)
pdf("PCA_seurat_sub_postregression.pdf")
DimPlot(seurat_sub, group.by = "sample", reduction = "pca")
DimPlot(seurat_sub, group.by = "subtype", reduction = "pca")
DimPlot(seurat_sub, split.by = "Phase", group.by = "sample", reduction = "pca")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "sample", reduction = "pca")
DimPlot(seurat_sub, split.by = "Phase", group.by = "subtype", reduction = "pca")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "subtype", reduction = "pca")

dev.off()

# tSNE subset data
seurat_sub <- RunTSNE(seurat_sub, dims = 1:40, reduction = "pca")
pdf("tsne_seurat_sub_postregression.pdf")                        
DimPlot(seurat_sub, group.by = "sample", reduction = "tsne")
DimPlot(seurat_sub, group.by = "subtype", reduction = "tsne")
DimPlot(seurat_sub, split.by = "Phase", group.by = "sample", reduction = "tsne")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "sample", reduction = "tsne")
DimPlot(seurat_sub, split.by = "Phase", group.by = "subtype", reduction = "tsne")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "subtype", reduction = "tsne")
dev.off()

# UMAP subset data
pdf("umap_seurat_sub_postregression.pdf")                        
DimPlot(seurat_sub, group.by = "sample", reduction = "umap")
DimPlot(seurat_sub, group.by = "subtype", reduction = "umap")
DimPlot(seurat_sub, split.by = "Phase", group.by = "sample", reduction = "umap")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "sample", reduction = "umap")
DimPlot(seurat_sub, split.by = "Phase", group.by = "subtype", reduction = "umap")
DimPlot(seurat_sub, split.by = "mitoFr", group.by = "subtype", reduction = "umap")
dev.off()

# Determine the K-nearest neighbor graph
seurat_sub <- FindNeighbors(seurat_sub, 
                              dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_sub <- FindClusters(object = seurat_sub,
                             resolution = seq(0.1, 2, 0.1))

sapply(grep("res",colnames(seurat_sub@meta.data),value = TRUE),
       function(x) length(unique(seurat_sub@meta.data[,x])))

# Assign identity of clusters
Idents(object = seurat_sub) <- "RNA_snn_res.0.5"

# Plot the tSNE
pdf("tsne_seurat_sub_clustered.pdf")
DimPlot(seurat_sub,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)
DimPlot(seurat_sub,
        reduction = "tsne",
        group.by = "subtype")
dev.off()

# Plot the UMAP
pdf("umap_seurat_sub_clustered.pdf")
DimPlot(seurat_sub,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
DimPlot(seurat_sub,
        reduction = "umap",
      group.by = "subtype")
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_sub) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_sub <- NormalizeData(seurat_sub, verbose = FALSE)

# Save subsetted data
save(seurat_sub, file = "seurat_sub2.RData")