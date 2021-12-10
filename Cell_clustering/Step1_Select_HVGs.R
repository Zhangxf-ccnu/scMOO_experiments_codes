#################### select high variable genes #######################
library(Seurat)
x.seurat <- CreateSeuratObject(data)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)

filter_ID <- x.seurat@assays$RNA@var.features
data_hvgs <- data[filter_ID,]  #### preprocessed data with 2,000 HVGs before imputation