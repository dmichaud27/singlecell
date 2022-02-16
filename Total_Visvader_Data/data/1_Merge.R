# Single-cell RNA-seq analysis - Merge

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

# Create a Seurat objects for HR Samples
for (file in c("AH0319", "MH0001", "MH0025", "MH0029_7c", "MH0029_9c", "MH0032", "MH0040", "MH0042", "MH0043", "MH0056", "MH0064", "MH0068", "MH0114_HR", "MH0125", "MH0151", "MH0163", "MH0167", "MH0173", "PM0178", "PM0360")){
  seurat_data <- Read10X(data.dir = paste0("/n/data1/bwh/medicine/guerriero/single-cell/dem11/visvader/GSE161529_RAW/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Create a merged Seurat object
merged_seuratHR <- merge(x = AH0319, 
                         y = c(MH0025, MH0001, MH0029_7c, MH0029_9c, MH0032, MH0040, MH0042, MH0043, MH0056, MH0064, MH0068, MH0114_HR, MH0125, MH0151, MH0163, MH0167, MH0173, PM0178, PM0360),
                         add.cell.id = c("319", "025", "001", "029_7c", "029_9c", "032", "040", "042", "043", "056", "064", "068", "114_HR", "125", "151", "163", "167", "173", "178", "360"))

# Add number of genes per UMI for each cell to metadata
merged_seuratHR$log10GenesPerUMI <- log10(merged_seuratHR$nFeature_RNA) / log10(merged_seuratHR$nCount_RNA)

# Compute percent mito ratio
merged_seuratHR$mitoRatio <- PercentageFeatureSet(object = merged_seuratHR, pattern = "^MT-")
merged_seuratHR$mitoRatio <- merged_seuratHR@meta.data$mitoRatio / 100

# Create metadata dataframe
metadataHR <- merged_seuratHR@meta.data

# Add cell IDs to metadata
metadataHR$cells <- rownames(metadataHR)

# Create sample column
metadataHR$sample <- NA
metadataHR$sample[which(str_detect(metadataHR$cells, "319_"))] <- "319"
metadataHR$sample[which(str_detect(metadataHR$cells, "001_"))] <- "001"
metadataHR$sample[which(str_detect(metadataHR$cells, "025_"))] <- "025"
metadataHR$sample[which(str_detect(metadataHR$cells, "029_7c_"))] <- "029_7c"
metadataHR$sample[which(str_detect(metadataHR$cells, "029_9c_"))] <- "029_9c"
metadataHR$sample[which(str_detect(metadataHR$cells, "032_"))] <- "032"
metadataHR$sample[which(str_detect(metadataHR$cells, "040_"))] <- "040"
metadataHR$sample[which(str_detect(metadataHR$cells, "042_"))] <- "042"
metadataHR$sample[which(str_detect(metadataHR$cells, "043_"))] <- "043"
metadataHR$sample[which(str_detect(metadataHR$cells, "056_"))] <- "056"
metadataHR$sample[which(str_detect(metadataHR$cells, "064_"))] <- "064"
metadataHR$sample[which(str_detect(metadataHR$cells, "068_"))] <- "068"
metadataHR$sample[which(str_detect(metadataHR$cells, "114_HR"))] <- "114_HR"
metadataHR$sample[which(str_detect(metadataHR$cells, "125_"))] <- "125"
metadataHR$sample[which(str_detect(metadataHR$cells, "151_"))] <- "151"
metadataHR$sample[which(str_detect(metadataHR$cells, "163_"))] <- "163"
metadataHR$sample[which(str_detect(metadataHR$cells, "167_"))] <- "167"
metadataHR$sample[which(str_detect(metadataHR$cells, "173_"))] <- "173"
metadataHR$sample[which(str_detect(metadataHR$cells, "178_"))] <- "178"
metadataHR$sample[which(str_detect(metadataHR$cells, "360_"))] <- "360"

# Create subtype column
metadataHR$subtype <- NA
metadataHR$subtype[which(str_detect(metadataHR$cells, "_"))] <- "HR"

# Rename columns
metadataHR <- metadataHR %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seuratHR@meta.data <- metadataHR

# Create .RData object to load at any time
save(merged_seuratHR, file="merged_seuratHR.RData")

# Create a Seurat object for HER2 samples
for (file in c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337")){
  seurat_data <- Read10X(data.dir = paste0("/n/data1/bwh/medicine/guerriero/single-cell/dem11/visvader/GSE161529_RAW", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Create a merged Seurat object
merged_seuratHER2 <- merge(x = AH0308, 
                       y = c(MH0031, MH0069, MH0161, MH0176, PM0337),
                       add.cell.id = c("AH0308", "MH0031", "MH0069", "MH0161", "MH0176", "PM0337"))

# Add number of genes per UMI for each cell to metadata
merged_seuratHER2$log10GenesPerUMI <- log10(merged_seuratHER2$nFeature_RNA) / log10(merged_seuratHER2$nCount_RNA)

# Compute percent mito ratio
merged_seuratHER2$mitoRatio <- PercentageFeatureSet(object = merged_seuratHER2, pattern = "^MT-")
merged_seuratHER2$mitoRatio <- merged_seuratHER2@meta.data$mitoRatio / 100

# Create metadata dataframe
metadataHER2 <- merged_seuratHER2@meta.data

# Add cell IDs to metadata
metadataHER2$cells <- rownames(metadataHER2)

# Create sample column
metadataHER2$sample <- NA
metadataHER2$sample[which(str_detect(metadataHER2$cells, "308_"))] <- "308"
metadataHER2$sample[which(str_detect(metadataHER2$cells, "031_"))] <- "031"
metadataHER2$sample[which(str_detect(metadataHER2$cells, "069_"))] <- "069"
metadataHER2$sample[which(str_detect(metadataHER2$cells, "161_"))] <- "161"
metadataHER2$sample[which(str_detect(metadataHER2$cells, "176_"))] <- "176"
metadataHER2$sample[which(str_detect(metadataHER2$cells, "337_"))] <- "337"

# Create subtype column
metadataHER2$subtype <- NA
metadataHER2$subtype[which(str_detect(metadataHER2$cells, "_"))] <- "HER2"

# Rename columns
metadataHER2 <- metadataHER2 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seuratHER2@meta.data <- metadataHER2

# Create .RData object to load at any time
save(merged_seuratHER2, file="merged_seuratHER2.RData")

# Create a Seurat object for each TNBC sample
for (file in c("MH0114", "MH0126", "MH0131", "MH0135", "MH0177", "MH4031", "SH0106", "Tum0554")){
  seurat_data <- Read10X(data.dir = paste0("/n/data1/bwh/medicine/guerriero/single-cell/dem11/visvader/GSE161529_RAW", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Create a merged Seurat object
merged_seuratTNBC <- merge(x = MH0114, 
                       y = c(MH0126, MH0131, MH0135, MH0177, MH4031, SH0106, Tum0554),
                       add.cell.id = c("MH0114", "MH0126", "MH0131", "MH0135", "MH0177", "MH4031", "SH0106", "Tum0554"))

# Add number of genes per UMI for each cell to metadata
merged_seuratTNBC$log10GenesPerUMI <- log10(merged_seuratTNBC$nFeature_RNA) / log10(merged_seuratTNBC$nCount_RNA)

# Compute percent mito ratio
merged_seuratTNBC$mitoRatio <- PercentageFeatureSet(object = merged_seuratTNBC, pattern = "^MT-")
merged_seuratTNBC$mitoRatio <- merged_seuratTNBC@meta.data$mitoRatio / 100

# Create metadata dataframe
metadataTNBC <- merged_seuratTNBC@meta.data

# Add cell IDs to metadata
metadataTNBC$cells <- rownames(metadataTNBC)

# Create sample column
metadataTNBC$sample <- NA
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "114_"))] <- "114"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "126_"))] <- "126"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "131_"))] <- "131"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "135_"))] <- "135"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "177_"))] <- "177"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "4031_"))] <- "4031"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "106_"))] <- "106"
metadataTNBC$sample[which(str_detect(metadataTNBC$cells, "554_"))] <- "554"

# Create subtype column
metadataTNBC$subtype <- NA
metadataTNBC$subtype[which(str_detect(metadataTNBC$cells, "_"))] <- "TNBC"

# Rename columns
metadataTNBC <- metadataTNBC %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seuratTNBC@meta.data <- metadataTNBC

# Create .RData object to load at any time
save(merged_seuratTNBC, file="merged_seuratTNBC.RData")

# Create a merged Seurat object
merged_seuratTotal <- merge(x = merged_seuratHR, 
                       y = c(merged_seuratHER2, merged_seuratTNBC))
                       
# Create metadata dataframe
metadataTotal <- merged_seuratTotal@meta.data

# Create .RData object to load at any time
save(merged_seuratTotal, file="merged_seuratTotal.RData")
