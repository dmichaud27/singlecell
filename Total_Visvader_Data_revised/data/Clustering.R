# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load seurat object
load("seurat_phase.RData")

# Printing out the most variable genes driving PCs
print(x = seurat_phase[["pca"]], 
      dims = 1:10, 
      nfeatures = 10)
use.pcs = 1:40

# Determine the K-nearest neighbor graph
seurat_phase <- FindNeighbors(seurat_phase, reduction = "pca",
                                   dims = use.pcs)

# Determine the clusters for various resolutions                                
seurat_clustered <- FindClusters(object = seurat_phase,
                                  resolution = seq(0.1, 2, 0.1))

sapply(grep("res",colnames(seurat_clustered@meta.data),value = TRUE),
       function(x) length(unique(seurat_clustered@meta.data[,x])))

save(seurat_clustered, file = "seurat_clustered.RData")

# Assign identity of clusters
Idents(object = seurat_clustered) <- "RNA_snn_res.1.9"

# Plot the tSNE
pdf("tSNE_clustered.pdf")
DimPlot(seurat_clustered, reduction = "tsne", label = TRUE, label.size = 6)
DimPlot(seurat_clustered, reduction = "tsne", group.by = "sample")
DimPlot(seurat_clustered, reduction = "tsne", group.by = "subtype")
DimPlot(seurat_clustered, reduction = "tsne", label = TRUE, label.size = 6, split.by = "mitoFr")
DimPlot(seurat_clustered, reduction = "tsne", label = TRUE, label.size = 6, split.by = "Phase")
dev.off()

# Plot the UMAP
pdf("UMAP_clustered.pdf")
DimPlot(seurat_clustered, reduction = "umap", label = TRUE, label.size = 6)
DimPlot(seurat_clustered, reduction = "umap", group.by = "sample")
DimPlot(seurat_clustered, reduction = "umap", group.by = "subtype")
DimPlot(seurat_clustered, reduction = "umap", label = TRUE, label.size = 6, split.by = "mitoFr")
DimPlot(seurat_clustered, reduction = "umap", label = TRUE, label.size = 6, split.by = "Phase")
dev.off()

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_clustered, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = "n_cells.csv")

# Defining the information in the seurat object of interest
columns_umap <- c(paste0("PC_", use.pcs),
             "ident",
             "UMAP_1", "UMAP_2")

# Defining the information in the seurat object of interest
columns_tsne <- c(paste0("PC_", use.pcs),
             "ident",
             "tSNE_1", "tSNE_2")

# Extracting this data from the seurat object
pc_data_umap <- FetchData(seurat_clustered, 
                          vars = columns_umap)

# Extracting this data from the seurat object
pc_data_tsne <- FetchData(seurat_clustered, 
                     vars = columns_tsne)

# Adding cluster label to center of cluster on tSNE
tsne_label <- FetchData(seurat_clustered, 
                        vars = c("ident", "tSNE_1", "tSNE_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(tSNE_1), y=mean(tSNE_2))

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_clustered, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_clustered) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_clustered <- NormalizeData(seurat_clustered, verbose = FALSE)

save(seurat_clustered, file = "seurat_clustered2.RData")