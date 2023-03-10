library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(Matrix)

#The chosen tissue is the heart
raw.data <- read.csv('FACS/Heart-counts.csv',row.names=1)
meta.data <- read.csv('metadata_FACS.csv')

plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]

rownames(meta.data) <- meta.data$plate.barcode
cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data)

# Find ERCC's and drop them from the raw data
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]

# Create the Seurat object with all the data
tiss <- CreateSeuratObject(counts = raw.data, min.cells = 5, min.genes = 5)

tiss <- AddMetaData(object = tiss, cell.meta.data)

# # Find and compute the percent of mitochondrial genes
tiss[["percent.mt"]] <- PercentageFeatureSet(tiss, pattern = "^MT-")

VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Mito genes were not found for some reason

#Filter out cells with few reads and few genes

tiss <- subset(tiss, subset = nFeature_RNA < 4500 & nCount_RNA < 2000000)

#Normalize the data
tiss <- NormalizeData(tiss)

#Identification of highly variable features
tiss <- FindVariableFeatures(tiss, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(tiss), 10)
plot1 <- VariableFeaturePlot(tiss)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
tiss <- ScaleData(tiss)

#Perform linear dimensional reduction
tiss <- RunPCA(tiss, features = VariableFeatures(object = tiss))
VizDimLoadings(tiss, dims = 1:2, reduction = "pca")
DimPlot(tiss, reduction = "pca")
DimHeatmap(tiss, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(tiss, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
ElbowPlot(tiss, 50)

#Cluster the cells
tiss <- FindNeighbors(tiss, dims = 1:30)
tiss <- FindClusters(tiss, resolution = 0.5)
head(Idents(tiss), 5)

#Run non-linear dimensional reduction UMAP
tiss <- RunUMAP(tiss, dims = 1:30)
DimPlot(tiss, reduction = "umap")

#Finding differentially expressed features (cluster biomarkers)
tiss.markers <- FindAllMarkers(tiss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tiss.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)






