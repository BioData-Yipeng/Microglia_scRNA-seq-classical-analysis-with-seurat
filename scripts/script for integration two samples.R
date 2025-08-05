# This script is to integrate two samples to eliminate batch effects
getwd()
setwd("F:/data science/R/6_GSE 133433")
library(dplyr)
library(Seurat)

wtmale.data <- Read10X(data.dir ="data/WTMale/" )
fivexmale.data <- Read10X(data.dir ="data/5XMale/" )
wtmale <- CreateSeuratObject(counts = wtmale.data,project = "wtmale")
fivexmale <- CreateSeuratObject(counts = fivexmale.data,project = "5Xmale")
dim(wtmale)
dim(fivexmale)
merged_seurat <- merge(wtmale, y=list(fivexmale),
                     add.cell.ids=c("wtmale", "5Xmale"),
                     project = "merged")
table(merged_seurat$orig.ident)

View(merged_seurat@meta.data)

merged_seurat$condition <- ifelse(merged_seurat$orig.ident %in% c("wtmale"), "wt", "5x") #add metadata to distinguish groups

merged_seurat <- NormalizeData(merged_seurat)

merged_seurat@assays$RNA@layers$counts.wtmale@x[1:5] #counts before normalization
merged_seurat@assays$RNA@layers$data.wtmale@x[1:5] #counts after normalization

merged_seurat <- FindVariableFeatures(merged_seurat)

variable_features <- VariableFeatures(merged_seurat)
top10 <- head(variable_features,10)
p1 <- VariableFeaturePlot(merged_seurat)
p2 <- LabelPoints(p1,top10,repel = T)
p1+p2

merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat) #pca is run on the 2000 most variable features
View(merged_seurat.merged@reductions$pca@cell.embeddings)
View(merged_seurat.merged@reductions$pca@feature.loadings)
DimPlot(object = merged_seurat, reduction = "pca", group.by = "condition") #default is pc1 vs pc2
DimPlot(object = merged_seurat, reduction = "pca", dims=c(3,4), group.by = "condition") #check pc3 vs pc4
VizDimLoadings(object = merged_seurat, reduction = "pca", nfeatures=20, dims = 1:2) #check genes contribute most to each pc
DimHeatmap(object = merged_seurat, reduction = "pca", dims = 1, cells = 2000, nfeatures=50, balanced = TRUE) #visualize genes contribute to pc1
print(merged_seurat[["pca"]],dims=1,nfeatures=10)

ElbowPlot(merged_seurat,ndims=50) #check variance of each pc

merged_seurat <- FindNeighbors(merged_seurat,
                             dims = 1:20,
                             reduction = "pca")
merged_seurat <- FindClusters(merged_seurat,
                            resolution = 0.4,
                            cluster.name = "MergedCluster")
merged_seurat <- RunUMAP(merged_seurat,
                       dims=1:20,
                       reduction = "pca",
                       reduction.name = "MergedUMAP")
DimPlot(merged_seurat,
        reduction = "MergedUMAP",
        group.by = "condition") #this is used to visualize batch effects

Integrated_seurat <- IntegrateLayers(object = merged_seurat,
                                   method = HarmonyIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "harmony",
                                   verbose=TRUE) #harmony method of integration

Integrated_seurat <- FindNeighbors(Integrated_seurat,
                                 reduction = "harmony",
                                 dims = 1:20)
Integrated_seurat <- FindClusters(Integrated_seurat,
                                resultion=0.4,
                                cluster.name = "harmonycluster")
Integrated_seurat <- RunUMAP(Integrated_seurat,
                           reduction = "harmony",
                           dims = 1:20,
                           reduction.name = "harmonyUMAP")

#visualization clusters in wt vs 5X
DimPlot(Integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "harmonycluster",
        split.by ="condition",
        label = TRUE) + ggtitle("clusters wt vs 5X")

#visualization and compare effect of integration
beforeintegrate <- DimPlot(merged_seurat,
        reduction = "MergedUMAP",
        group.by = "condition") + ggtitle("before integration")
afterintegrate <- DimPlot(Integrated_seurat,
        reduction = "harmonyUMAP",
        group.by = "condition") + ggtitle("after integration")
beforeintegrate+afterintegrate

#by comparing wt vs 5x it can be seen group 4 numbers decreases and group 5 increase
#the following codes is to annotate and extract cells in group 5 and group 4
#different expression in gourp 5 vs other groups

# Find marker genes for "Cluster1"
Integrated_seurat <- JoinLayers(Integrated_seurat)
de_cluster4 <- FindMarkers(Integrated_seurat, 
                                 ident.1 = "4",        # Specify the cluster of interest
                                 test.use = "wilcox",         # Statistical test (e.g., Wilcoxon)
                                 min.pct = 0.25,              # Minimum percentage of cells expressing the gene
                                 logfc.threshold = 0.25)      # Log-fold change threshold

de_cluster5 <- FindMarkers(Integrated_seurat, 
                           ident.1 = "5",        # Specify the cluster of interest
                           test.use = "wilcox",         # Statistical test (e.g., Wilcoxon)
                           min.pct = 0.25,              # Minimum percentage of cells expressing the gene
                           logfc.threshold = 0.25)      # Log-fold change threshold

#cluster 5 highly express CD83, CCL3 which are biomarkers for activated microglia.
#cluster 5 represent those activated microglia which increases in 5X animal model.

#cluster 4 highly express those interferon response genes.
#the decrease number in 5x animals indicates lower interferon signal activation.

FeaturePlot(Integrated_seurat,
            reduction = "harmonyUMAP",
            features = "GRCh38-IFIT3",
            cols = c("gray","red"),
            split.by = "condition")
FeaturePlot(Integrated_seurat,
            reduction = "harmonyUMAP",
            features = "GRCh38-CD9",
            cols = c("gray","red"),
            split.by = "condition")


#Violin plot for the top 3 marker genes
VlnPlot(Integrated_seurat,
        split.by = "condition",
        features = "GRCh38-IFIT3",
        layer = "data")

VlnPlot(Integrated_seurat,
        split.by = "condition",
        features = "GRCh38-CD9",
        layer = "data")

DotPlot(Integrated_seurat,
        features = c("GRCh38-CD9", "GRCh38-CD83", "GRCh38-IFIT3", "GRCh38-IFIT1"),
        split.by = "condition") + RotatedAxis()


