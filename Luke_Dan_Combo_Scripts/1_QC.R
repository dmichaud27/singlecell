# Single-cell RNA-seq analysis - QC

############################################################
# Directory for permanently saving data
# /n/data1/bwh/medicine/guerriero/single-cell/dem11/visvader

# Link to Seurat tutorial
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(DisHet)
library(SCINA)
library(pryr)

# Establish directories to load and save
# "exceldata" and "raw.dir" should remain the same
exceldata <- read.csv("/n/data1/bwh/medicine/guerriero/single-cell/kania/visvader/GSE161529_RAW/metadata.csv")
plot.dir <- "/n/scratch3/users/d/dem11/visvader/allsubtypes/plots"
data.dir <- "/n/scratch3/users/d/dem11/visvader/allsubtypes/data/rda"
raw.dir <- "/n/data1/bwh/medicine/guerriero/single-cell/kania/visvader/GSE161529_RAW"

# Selecting samples for analysis:
BC.samples <- exceldata[exceldata$type %in% c("er", "pr", "tnbc", "her2"), ]
BC.samples <- BC.samples[BC.samples$cells %in% "total", ]
sample.list <- c(BC.samples$Sample)
obj.list <- c()
obj.data.list <- c()

# Creating seurat objects from selected samples
for(n in 1:length(sample.list)) {
  if (sample.list[n] %in% BC.samples$Sample) {
    assign(paste0(sample.list[n],".data"), Read10X(data.dir = paste0(raw.dir, "/", sample.list[n])))
    assign(sample.list[n], CreateSeuratObject(counts = eval(parse(text = paste0(sample.list[n], ".data"))), project = sample.list[n], min.cells = 3, min.features = 100))
    obj.list = c(obj.list, eval(parse(text = sample.list[n])))
    obj.data.list <- c(obj.data.list, eval(parse(text = paste0(sample.list[n], ".data"))))
  }
}

gsms <- merge(obj.list[[1]], y = obj.list[2:length(obj.list)], add.cell.ids = sample.list, project = "gsms")
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

orig.table <- table(gsms@meta.data$orig)
gsms[["percent.mt"]] <- PercentageFeatureSet(gsms, pattern = "^MT-")

table(sample= gsms@meta.data$orig.ident,`Is live?`=gsms@meta.data$percent.mt<20) # 80% - 94%
table(sample= gsms@meta.data$orig.ident,`nFeature < 200` = gsms@meta.data$nFeature_RNA<200) # remove bottom 65%
table(sample= gsms@meta.data$orig.ident,`nFeature < 9000` = gsms@meta.data$nFeature_RNA<9000) # remove top five only
                       
# Create metadata dataframe
metadataTotal <- gsms@meta.data

# Create .RData object to load at any time
save(gsms, file="gsms_original.RData")

# Visualize the number of cell counts per sample
pdf("Ncells_before.pdf")
metadataTotal %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Visualize the number UMIs/transcripts per cell
pdf("NUMIs_before.pdf")
metadataTotal %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()

# Visualize the distribution of genes detected per cell via histogram
pdf("NGenes_histogram_before.pdf")
metadataTotal %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
pdf("NGenes_boxlpot_before.pdf")
metadataTotal %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
dev.off()

#Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("NgenesvsNUMIs_before.pdf")
metadataTotal %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("mitoRatio_before.pdf")
metadataTotal %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
pdf("NGenesperUMI_before.pdf")
metadataTotal %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

# Filter out low quality cells using selected thresholds - these will change with experiment
gsms <- subset(gsms, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 20)
filtered_seurat <- subset(x = merged_seuratTotal, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
save(filtered_seurat, file="seurat_filtered.RData")

# Create metadata dataframe
metadata_filtered <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
pdf("Ncells_after.pdf")
metadata_filtered %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Visualize the number UMIs/transcripts per cell
pdf("NUMIs_after.pdf")
metadata_filtered %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)
dev.off()

# Visualize the distribution of genes detected per cell via histogram
pdf("Ngenes_histogram_after.pdf")
metadata_filtered %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 800)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
pdf("Ngenes_boxplot_after.pdf")
metadata_filtered %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
dev.off()

#Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("NgenesvsNUMIs_after.pdf")
metadata_filtered %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
pdf("mitoRatio_after.pdf")
metadata_filtered %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
pdf("NgenesperUMI_after.pdf")
metadata_filtered %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()