# Single-cell RNA-seq - GSEA Analysis

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(limma)
library(topGO)
library(ggplot2)
library(msigdbr)
library(presto)
library(fgsea)

# Load seurat object
load("seurat_sub2.RData")

# Generate list of DE Genes
B.genes <- wilcoxauc(seurat_sub, group_by = "RNA_snn_res.0.5", seurat_assay = "RNA")
head(B.genes)
dplyr::count(B.genes, group)

# Set GSEA Collection
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Ranked DE gene list by cluster
B.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "2") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "3") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "4") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "5") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "6") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "7") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

B.genes %>%
  dplyr::filter(group == "8") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
  
B.genes %>%
  dplyr::filter(group == "9") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  
  
B.genes %>%
  dplyr::filter(group == "10") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  
  
B.genes %>%
  dplyr::filter(group == "11") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  
  
B.genes %>%
  dplyr::filter(group == "12") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  
  
B.genes %>%
  dplyr::filter(group == "13") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  
  
genes_0 <- B.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_0 <- deframe(genes_0)
head(ranks_0)

genes_1 <- B.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_1 <- deframe(genes_1)
head(ranks_1)

genes_2 <- B.genes %>%
  dplyr::filter(group == "2") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_2 <- deframe(genes_2)
head(ranks_2)

genes_3 <- B.genes %>%
  dplyr::filter(group == "3") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_3 <- deframe(genes_3)
head(ranks_3)

genes_4 <- B.genes %>%
  dplyr::filter(group == "4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_4 <- deframe(genes_4)
head(ranks_4)

genes_5 <- B.genes %>%
  dplyr::filter(group == "5") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_5 <- deframe(genes_5)
head(ranks_5)

genes_6 <- B.genes %>%
  dplyr::filter(group == "6") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_6 <- deframe(genes_6)
head(ranks_6)

genes_7 <- B.genes %>%
  dplyr::filter(group == "7") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_7 <- deframe(genes_7)
head(ranks_7)

genes_8 <- B.genes %>%
  dplyr::filter(group == "8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_8 <- deframe(genes_8)
head(ranks_8)

genes_9 <- B.genes %>%
  dplyr::filter(group == "9") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_9 <- deframe(genes_9)
head(ranks_9)

genes_10 <- B.genes %>%
  dplyr::filter(group == "10") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_10 <- deframe(genes_10)
head(ranks_10)

genes_11 <- B.genes %>%
  dplyr::filter(group == "11") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_11 <- deframe(genes_11)
head(ranks_11)

genes_12 <- B.genes %>%
  dplyr::filter(group == "12") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_12 <- deframe(genes_12)
head(ranks_12)

genes_13 <- B.genes %>%
  dplyr::filter(group == "13") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks_13 <- deframe(genes_13)
head(ranks_13)

# Run GSEA Analysis
fgseaRes_0 <- fgsea(fgsea_sets, stats = ranks_0, nperm = 10000)
fgseaResTidy_0 <- fgseaRes_0 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_0 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_1 <- fgsea(fgsea_sets, stats = ranks_1, nperm = 10000)
fgseaResTidy_1 <- fgseaRes_1 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_1 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_2 <- fgsea(fgsea_sets, stats = ranks_2, nperm = 10000)
fgseaResTidy_2 <- fgseaRes_2 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_2 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_3 <- fgsea(fgsea_sets, stats = ranks_3, nperm = 10000)
fgseaResTidy_3 <- fgseaRes_3 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_3 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_4 <- fgsea(fgsea_sets, stats = ranks_4, nperm = 10000)
fgseaResTidy_4 <- fgseaRes_4 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_4 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_5 <- fgsea(fgsea_sets, stats = ranks_5, nperm = 10000)
fgseaResTidy_5 <- fgseaRes_5 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_5 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_6 <- fgsea(fgsea_sets, stats = ranks_6, nperm = 10000)
fgseaResTidy_6 <- fgseaRes_6 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_6 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_7 <- fgsea(fgsea_sets, stats = ranks_7, nperm = 10000)
fgseaResTidy_7 <- fgseaRes_7 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_7 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_8 <- fgsea(fgsea_sets, stats = ranks_8, nperm = 10000)
fgseaResTidy_8 <- fgseaRes_8 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_8 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
  
fgseaRes_9 <- fgsea(fgsea_sets, stats = ranks_9, nperm = 10000)
fgseaResTidy_9 <- fgseaRes_9 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_9 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
  
fgseaRes_10 <- fgsea(fgsea_sets, stats = ranks_10, nperm = 10000)
fgseaResTidy_10 <- fgseaRes_10 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_10 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
  
fgseaRes_11 <- fgsea(fgsea_sets, stats = ranks_11, nperm = 10000)
fgseaResTidy_11 <- fgseaRes_11 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_11 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_12 <- fgsea(fgsea_sets, stats = ranks_12, nperm = 10000)
fgseaResTidy_12 <- fgseaRes_12 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_12 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

fgseaRes_13 <- fgsea(fgsea_sets, stats = ranks_13, nperm = 10000)
fgseaResTidy_13 <- fgseaRes_13 %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy_13 %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
  
# Plot Top GSEA Pathways
pdf("GSEA_Cluster0.pdf")
ggplot(fgseaResTidy_0 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster1.pdf")
ggplot(fgseaResTidy_1 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster2.pdf")
ggplot(fgseaResTidy_2 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster3.pdf")
ggplot(fgseaResTidy_3 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster4.pdf")
ggplot(fgseaResTidy_4 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster5.pdf")
ggplot(fgseaResTidy_5 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster6.pdf")
ggplot(fgseaResTidy_6 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster7.pdf")
ggplot(fgseaResTidy_7 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster8.pdf")
ggplot(fgseaResTidy_8 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster9.pdf")
ggplot(fgseaResTidy_9 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster10.pdf")
ggplot(fgseaResTidy_10 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster11.pdf")
ggplot(fgseaResTidy_11 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster12.pdf")
ggplot(fgseaResTidy_12 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()

pdf("GSEA_Cluster13.pdf")
ggplot(fgseaResTidy_13 %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
dev.off()