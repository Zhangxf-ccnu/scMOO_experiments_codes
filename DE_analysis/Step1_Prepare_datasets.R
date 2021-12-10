################### prepare datasets ########################
GSE75748_bulk_cell_type_ec <- read.csv("./GSE75748_bulk_cell_type_ec.csv", row.names = 1)
GSE75748_sc_cell_type_ec <- read.csv("./GSE75748_sc_cell_type_ec.csv", row.names = 1)

############## select 6 pairs of cell subpopulations containing DEC ################
## EC vs DEC  ############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(8:9, 10:12)]
bulk_group <- as.factor(c(rep("DEC",2),rep("EC",3)))
sc_data <- GSE75748_sc_cell_type_ec[,c(375:512,513:617)]
sc_group <- as.factor(c(rep("DEC",138),rep("EC",105)))

## H1 vs DEC ##############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(1:4, 8:9)]
bulk_group <- as.factor(c(rep("H1_ESC",4),rep("DEC",2)))
sc_data <- GSE75748_sc_cell_type_ec[,c(1:212,375:512)]
sc_group <- as.factor(c(rep("H1_ESC",212),rep("DEC",138)))

## H9 vs DEC ###############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(5:7, 8:9)]
bulk_group <- as.factor(c(rep("H9",3),rep("DEC",2)))
sc_data <- GSE75748_sc_cell_type_ec[,c(213:374,375:512)]
sc_group <- as.factor(c(rep("H9",162),rep("DEC",138)))

## HFF vs DEC ############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(8:9, 13:15)]
bulk_group <- as.factor(c(rep("DEC",2),rep("HFF",3)))
sc_data <- GSE75748_sc_cell_type_ec[,c(375:512,618:776)]
sc_group <- as.factor(c(rep("DEC",138),rep("HFF",159)))

## NPC vs DEC #############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(16:17, 8:9)]
bulk_group <- as.factor(c(rep("NPC",2),rep("DEC",2)))
sc_data <- GSE75748_sc_cell_type_ec[,c(777:949,375:512)]
sc_group <- as.factor(c(rep("NPC",173),rep("DEC",138)))

## TB vs DEC  ############
bulk_data <- GSE75748_bulk_cell_type_ec[,c(8:9, 18:19)]
bulk_group <- as.factor(c(rep("DEC",2),rep("TB",2)))
sc_data <- GSE75748_sc_cell_type_ec[,c(375:512,950:1000)]
sc_group <- as.factor(c(rep("DEC",138),rep("TB",51)))


################  select 2,000 HVGs before imputation  ################
library(Seurat)
x.seurat <- CreateSeuratObject(sc_data)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)

filter_ID <- x.seurat@assays$RNA@var.features

bulk_data <- as.matrix(bulk_data)[filter_ID ,]
sc_data <- as.matrix(sc_data)[filter_ID ,]