################### SC3 ########################
library(SingleCellExperiment)
library(SC3)
library(mclust)
library(aricode)

result <- readRDS('./data_imp.rds') #### imputation results 
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(result),
    logcounts = log2(as.matrix(result) + 1)
  ), 
  colData = colnames(result)
)


rowData(sce)$feature_symbol <- rownames(sce)

label <- colnames(result)

K <- length(unique(label))


result <- sc3(sce, ks = K, gene_filter = FALSE)

ARI <- adjustedRandIndex(label, result@colData@listData$sc3_6_clusters) #### '6' is consistant, and equals to K

NMI <- NMI(result@colData@listData$sc3_6_clusters, label)


################# Seurat ########################
library(Seurat)
x.seurat <- CreateSeuratObject(result) #### imputation results 
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)

x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))  
x.seurat <- JackStraw(x.seurat, num.replicate = 100)
x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)
x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
x.seurat <- FindClusters(x.seurat, resolution = 0.5)

label <- colnames(result)

ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)

NMI <- NMI(x.seurat$seurat_clusters, as.factor(label))
