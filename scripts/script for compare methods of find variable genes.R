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
wtmale.data <- Read10X(data.dir ="data/WTMale/" )
wtmale <- CreateSeuratObject(counts = wtmale.data,project = "wtmale")
rownames(wtmale)


#QCandFilter, check mitochondrial gene percentage, nFeature
wtmale[["percent.mt"]] <- PercentageFeatureSet(wtmale,pattern = "^GRCh38-MT-")
VlnPlot(wtmale,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(wtmale, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wtmale, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dim(wtmale)
wtmale <- subset(wtmale, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
dim(wtmale) #623 cells were filtered

#Analysis
wtmale <- NormalizeData(wtmale) #normalize data to align each cell

#find most variable genes by various methods
wtmale_vst <- FindVariableFeatures(wtmale, selection.method = "vst",  nfeatures = 2000) #find the most varible features
wtmale_meanvarplot <- FindVariableFeatures(wtmale, selection.method = "mean.var.plot",  nfeatures = 2000) #find the most varible features
wtmale_dispersion <- FindVariableFeatures(wtmale, selection.method = "dispersion",  nfeatures = 2000) #find the most varible features

#visualization of variability features
top10_vst <- head(VariableFeatures(wtmale_vst),10)
top10_meanvarplot <- head(VariableFeatures(wtmale_meanvarplot),10)
top10_dispersion <- head(VariableFeatures(wtmale_dispersion),10)
top10_vst
top10_meanvarplot
top10_dispersion

plot3_vst <- VariableFeaturePlot(wtmale_vst) + theme(legend.position = "top")
plot4_vst <- LabelPoints(plot = plot3_vst,points = top10_vst,repel = TRUE) + theme(legend.position = "top")
plot3_vst + plot4_vst

plot3_meanvarplot <- VariableFeaturePlot(wtmale_meanvarplot) + theme(legend.position = "top")
plot4_meanvarplot <- LabelPoints(plot = plot3_meanvarplot,points = top10_meanvarplot,repel = TRUE) + theme(legend.position = "top")
plot3_meanvarplot + plot4_meanvarplot

plot3_dispersion <- VariableFeaturePlot(wtmale_dispersion) + theme(legend.position = "top")
plot4_dispersion <- LabelPoints(plot = plot3_dispersion,points = top10_dispersion,repel = TRUE) + theme(legend.position = "top")
plot3_dispersion + plot4_dispersion

#visualization of overlap of the top2000 most variable genes
install.packages("VennDiagram")
library(VennDiagram)
set_vst <- VariableFeatures(wtmale_vst)
set_meanvarplot <- VariableFeatures(wtmale_meanvarplot)
set_dispersion <- VariableFeatures(wtmale_dispersion)

venn.plot <- venn.diagram(x=list(Set1=set_vst, Set2=set_meanvarplot, Set3=set_dispersion),
category.names=c("vst","meanvarplot","dispersion"),
filename = NULL,
output = TRUE,
imagetype = "png",
col="black",
fill=c("skyblue","orange","green"),
alpha=0.5)
grid.draw(venn.plot)

#conclusion
#different methods of findvariablefeature produce different list of genes with about half overlap.
#commonly used
rm(list=ls())
gc()
