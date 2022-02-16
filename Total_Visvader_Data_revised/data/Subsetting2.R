# Single-cell RNA-seq - Subsetting & Reanalysis

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_sub2.RData")

# Subset data to get rid of unwanted clusters
seurat_subFinal <- subset(seurat_sub, idents = c(24,17,30,31,28), invert = TRUE)
levels(seurat_subFinal)

# Count number of cells in existing clusters
n_cells_sub2 <- FetchData(seurat_subFinal, 
                     vars = c("ident", "orig.ident")) %>%
dplyr::count(ident, orig.ident) %>%
tidyr::spread(ident, n)

write.csv(n_cells_sub2, file = "n_cells_subFinal.csv")

# Save subsetted data
save(seurat_subFinal, file = "seurat_subFinal.RData")

# Reanalyze data
seurat_subFinal <- NormalizeData(seurat_subFinal)

seurat_subFinal <- FindVariableFeatures(seurat_subFinal, 
                                         selection.method = "vst",
                                         nfeatures = 2000, 
                                         verbose = FALSE)

# Scale subset data
seurat_subFinal <- ScaleData(seurat_subFinal, vars.to.regress = "mitoRatio")
seurat_subFinal <- RunPCA(seurat_subFinal, features = VariableFeatures(seurat_subFinal), nfeatures.print = 10, ndims.print = 1:10)
pdf("PCA_seurat_subFinal_postregression.pdf")
DimPlot(seurat_subFinal, group.by = "sample", reduction = "pca")
DimPlot(seurat_subFinal, group.by = "subtype", reduction = "pca")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "sample", reduction = "pca")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "sample", reduction = "pca")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "subtype", reduction = "pca")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "subtype", reduction = "pca")
dev.off()

# tSNE subset data
seurat_subFinal <- RunTSNE(seurat_subFinal, dims = 1:40, reduction = "pca")
pdf("tsne_seurat_subFinal_postregression.pdf")                        
DimPlot(seurat_subFinal, group.by = "sample", reduction = "tsne")
DimPlot(seurat_subFinal, group.by = "subtype", reduction = "tsne")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "sample", reduction = "tsne")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "sample", reduction = "tsne")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "subtype", reduction = "tsne")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "subtype", reduction = "tsne")
dev.off()

# UMAP subset data
pdf("umap_seurat_subFinal_postregression.pdf")                        
DimPlot(seurat_subFinal, group.by = "sample", reduction = "umap")
DimPlot(seurat_subFinal, group.by = "subtype", reduction = "umap")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "sample", reduction = "umap")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "sample", reduction = "umap")
DimPlot(seurat_subFinal, split.by = "Phase", group.by = "subtype", reduction = "umap")
DimPlot(seurat_subFinal, split.by = "mitoFr", group.by = "subtype", reduction = "umap")
dev.off()

# Determine the K-nearest neighbor graph
seurat_subFinal <- FindNeighbors(seurat_subFinal, 
                              dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_subFinal <- FindClusters(object = seurat_subFinal,
                             resolution = seq(0.1, 2, 0.1))

sapply(grep("res",colnames(seurat_subFinal@meta.data),value = TRUE),
       function(x) length(unique(seurat_subFinal@meta.data[,x])))

# Assign identity of clusters
Idents(object = seurat_subFinal) <- "RNA_snn_res.1.9"

# Plot the tSNE
pdf("tsne_seurat_subFinal_clustered.pdf")
DimPlot(seurat_subFinal,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)
DimPlot(seurat_subFinal,
        reduction = "tsne",
        group.by = "subtype")
DimPlot(seurat_subFinal,
        reduction = "tsne",
        label = TRUE,
        label.size = 6,
        split.by = "subtype")
dev.off()

# Plot the UMAP
pdf("umap_seurat_subFinal_clustered.pdf")
DimPlot(seurat_subFinal,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
DimPlot(seurat_subFinal,
        reduction = "umap",
      group.by = "subtype")
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_subFinal) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_subFinal <- NormalizeData(seurat_subFinal, verbose = FALSE)

# Save subsetted data
save(seurat_subFinal, file = "seurat_subFinal2.RData")