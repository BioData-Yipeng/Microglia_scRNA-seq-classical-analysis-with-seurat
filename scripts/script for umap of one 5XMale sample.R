getwd()
setwd("F:/data science/R/6_GSE 133433")
#library(dplyr)
library(Seurat)
#library(patchwork)
library(ggplot2)
library(presto)
#library(tidyverse)
library(data.table)
library(Matrix)
remotes::install_github("immunogenomics/presto")

# Creat seurat object
fivexmale.data <- Read10X(data.dir ="data/5xmale/")
fivexmale <- CreateSeuratObject(counts = fivexmale.data,project = "5xmale")
rownames(fivexmale)


#QCandFilter, check mitochondrial gene percentage, nFeature
fivexmale[["percent.mt"]] <- PercentageFeatureSet(fivexmale,pattern = "^GRCh38-MT-")
VlnPlot(fivexmale,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(fivexmale, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fivexmale, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dim(fivexmale)
fivexmale <- subset(fivexmale, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
dim(fivexmale) #a number of cells were filtered

#Analysis
fivexmale <- NormalizeData(fivexmale) #normalize data to align each cell
fivexmale <- FindVariableFeatures(fivexmale, nfeatures = 2000) #find the most varible features

top10 <- head(VariableFeatures(fivexmale),10)
plot3 <- VariableFeaturePlot(fivexmale) + theme(legend.position = "top")
plot4 <- LabelPoints(plot = plot3,points = top10,repel = TRUE) + theme(legend.position = "top")
plot3 + plot4

#output the most 1000 variable genes
top1000_genes <- head(VariableFeatures(fivexmale),1000)
top1000_gene_matrix <- GetAssayData(fivexmale, slot = "data")[top1000_genes,]

#calcualte and explore the frequency of each gene across all cells
top1000_nonzero <- Matrix::rowSums(top1000_gene_matrix>0)
top1000_frequency <- top1000_nonzero/ncol(top1000_gene_matrix)
top1000_gene_matrix <- as.data.frame(top1000_gene_matrix)
top1000_gene_matrix$gene_frequency <- top1000_frequency
top1000_gene_matrix$nozero <- top1000_nonzero
write.csv(top1000_gene_matrix,"fivexmale_top1000_variable_genes.csv")
#visualize the distribution of frequency
hist(top1000_frequency,breaks = 100)
plot(density(top1000_frequency))

#dimision reduction
fivexmale <- ScaleData(fivexmale)
fivexmale<- RunPCA(fivexmale, features = VariableFeatures(object = fivexmale))
View(fivexmale@reductions$pca@cell.embeddings) # view the default 50 PCs
View(fivexmale@reductions$pca@feature.loadings) # view how much each gene contribute to each PC
VizDimLoadings(fivexmale,dims = 1:2, nfeatures = 30) # view the most contribute genes to PC1 or PC2
print(fivexmale[["pca"]],dims=1:2,nfeatures=5) # print the top postive or negative gene contribute to PC1 or PC2
DimPlot(fivexmale, reduction = "pca") #PCA plot
DimHeatmap(fivexmale, dims = 1, cells = 500, balanced = TRUE) # heatmap show most contributed genes to PC1
ElbowPlot(fivexmale) #based on this pick a pc number for following analysis
fivexmale <- FindNeighbors(fivexmale,dims = 1:14)
fivexmale <- FindClusters(fivexmale, resolution = 0.2) # adjust the resolution to control culster numbers, low number mean few clusters
head(Idents(fivexmale),10)
fivexmale <- RunUMAP(fivexmale, dims = 1:14)
DimPlot(fivexmale, reduction = "umap", label = TRUE) + ggtitle("5X_Male UMAP")

#feature plot for some biomarker genes
FeaturePlot(fivexmale, features = "GRCh38-CX3CR1",label = TRUE, cols = c("gray", "red")) #CX3CR1 is the canonical biomarker for microglia
FeaturePlot(fivexmale, features = "GRCh38-CXCL10",label = TRUE, cols = c("gray", "red")) #CXCL10 is the inflammatory chemokines for microglia
FeaturePlot(fivexmale, features = "GRCh38-TREM2",label = TRUE, cols = c("gray", "red")) #TREM2 is the DAM disease associated microglia biomarker
FeaturePlot(fivexmale, features = "GRCh38-CD9",label = TRUE, cols = c("gray", "red")) #CD9 is the DAM disease associated microglia biomarker
FeaturePlot(fivexmale, features = "GRCh38-LPL",label = TRUE, cols = c("gray", "red")) #LPL is the DAM disease associated microglia biomarker
FeaturePlot(fivexmale, features = "GRCh38-ITGAX",label = TRUE, cols = c("gray", "red")) #ITGAX is DAM disease associated microglia biomarker

#vlnplot for some biomarker genes
VlnPlot(fivexmale, features = "GRCh38-CX3CR1")
VlnPlot(fivexmale, features = "GRCh38-CXCL10")
VlnPlot(fivexmale, features = "GRCh38-TREM2")
VlnPlot(fivexmale, features = "GRCh38-CD9")
VlnPlot(fivexmale, features = "GRCh38-LPL")
VlnPlot(fivexmale, features = "GRCh38-ITGAX")

#DE analysis marker genes for each cluster
de <- wilcoxauc(fivexmale, 'seurat_clusters') # by presto method
dim(de)

top_markers <- de %>%
  group_by(group) %>%
  filter(auc > 0.5, (logFC > 0.5 | logFC < -0.5), padj < 0.05) %>%
  slice_max(order_by = auc, n = 5)
#draw heatmap to show top different genes in each cluster
DoHeatmap(fivexmale, features = unique(top_markers$feature)) + NoLegend()

#commonly used
rm(list=ls())
gc()
