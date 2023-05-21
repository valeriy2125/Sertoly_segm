install.packages('Seurat')
library(Seurat)

install.packages("tidyverse")
install.packages("devtools")

library(dplyr)
library(Seurat)
library(patchwork)



pbmc.data <- Read10X(data.dir = "../Downloads/scSeq2/")
pbmc <- CreateSeuratObject(counts = pbmc1.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mt < 1)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.4)

head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:20)

DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "../Downloads/scSeq/pbmc.rds")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv2(pbmc.markers, "../Downloads/scSeq/markers.csv")


cluster14.markers <- FindMarkers(pbmc, ident.1 = 14, ident.2 = c(2, 3, 11), min.pct = 0.25)
head(cluster14.markers, n = 5)

write.csv2(cluster14.markers, "../Downloads/scSeq/markersRTvsSC.csv")

new.cluster.ids <- c("Germ0", "Int1", "SC2", "SC3", "Int4", "Int5", "Germ6","Germ7", "Germ8", "Germ9", "Germ10", "SC11", "Int12", "Macro13", "RT", "Int15", "Int16", "Endot")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
