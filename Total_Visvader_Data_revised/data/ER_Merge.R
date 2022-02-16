# Single-cell RNA-seq analysis - ER Merge

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

# Create a Seurat object for each sample
for (file in c("AH0319_t", "MH0001_t", "MH0025_t", "MH0029_7c_t", "MH0029_9c_t", "MH0032_t", "MH0040_t", "MH0042_t", "MH0043_t", "MH0056_t", "MH0064_t", "MH0068_t", "MH0114_t", "MH0125_t", "MH0151_t", "MH0163_t", "MH0167_t", "MH0173_t", "PM0178_t", "PM0360_t")){
  seurat_data <- Read10X(data.dir = paste0("ER_samples/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}


# Create a merged Seurat object
merged_seuratHR <- merge(x = AH0319_t, 
                       y = c(MH0025_t, MH0001_t, MH0029_7c_t, MH0029_9c_t, MH0032_t, MH0040_t, MH0042_t, MH0043_t, MH0056_t, MH0064_t, MH0068_t, MH0114_t, MH0125_t, MH0151_t, MH0163_t, MH0167_t, MH0173_t, PM0178_t, PM0360_t),
                       add.cell.id = c("319_t", "025_t", "001_t", "029_7c_t", "029_9c_t", "032_t", "040_t", "042_t", "043_t", "056_t", "064_t", "068_t", "114_t", "125_t", "151_t", "163_t", "167_t", "173_t", "178_t", "360_t"))

# Add number of genes per UMI for each cell to metadata
merged_seuratHR$log10GenesPerUMI <- log10(merged_seuratHR$nFeature_RNA) / log10(merged_seuratHR$nCount_RNA)

# Compute percent mito ratio
merged_seuratHR$mitoRatio <- PercentageFeatureSet(object = merged_seuratHR, pattern = "^MT-")
merged_seuratHR$mitoRatio <- merged_seuratHR@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seuratHR@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "319_t"))] <- "319t"
metadata$sample[which(str_detect(metadata$cells, "001_t"))] <- "001t"
metadata$sample[which(str_detect(metadata$cells, "025_t"))] <- "025t"
metadata$sample[which(str_detect(metadata$cells, "029_7c_t"))] <- "029_7ct"
metadata$sample[which(str_detect(metadata$cells, "029_9c_t"))] <- "029_9ct"
metadata$sample[which(str_detect(metadata$cells, "032_t"))] <- "032t"
metadata$sample[which(str_detect(metadata$cells, "040_t"))] <- "040t"
metadata$sample[which(str_detect(metadata$cells, "042_t"))] <- "042t"
metadata$sample[which(str_detect(metadata$cells, "043_t"))] <- "043t"
metadata$sample[which(str_detect(metadata$cells, "056_t"))] <- "056t"
metadata$sample[which(str_detect(metadata$cells, "064_t"))] <- "064t"
metadata$sample[which(str_detect(metadata$cells, "068_t"))] <- "068t"
metadata$sample[which(str_detect(metadata$cells, "114_t"))] <- "114t"
metadata$sample[which(str_detect(metadata$cells, "125_t"))] <- "125t"
metadata$sample[which(str_detect(metadata$cells, "151_t"))] <- "151t"
metadata$sample[which(str_detect(metadata$cells, "163_t"))] <- "163t"
metadata$sample[which(str_detect(metadata$cells, "167_t"))] <- "167t"
metadata$sample[which(str_detect(metadata$cells, "173_t"))] <- "173t"
metadata$sample[which(str_detect(metadata$cells, "178_t"))] <- "178t"
metadata$sample[which(str_detect(metadata$cells, "360_t"))] <- "360t"

# Create subtype column
metadata$subtype <- NA
metadata$subtype[which(str_detect(metadata$cells, "_"))] <- "HR"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seuratHR@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seuratHR, file="merged_seuratHR.RData")
