library(SingleCellExperiment)
library(TSCAN)
library(monocle)
#library(destiny)
# library(SLICER)
# library(ouija)
library(scater)
library(ggplot2)
library(igraph)


################# monocle 2 function #######################

my.monocle =  function(count, cellLabels){
  
  
  colnames(count) <- 1:ncol(count)
  geneNames <- rownames(count)
  rownames(count) <- 1:nrow(count)
  
  
  pd <- data.frame(timepoint = cellLabels)
  pd <- new("AnnotatedDataFrame", data=pd)
  fd <- data.frame(gene_short_name = geneNames)
  fd <- new("AnnotatedDataFrame", data=fd)
  
  dCellData <- newCellDataSet(count, phenoData = pd, featureData = fd, expressionFamily = uninormal())
  
  dCellData <- detectGenes(dCellData , min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(dCellData),
                                      num_cells_expressed >= 50))
  
  diff_test_res <- differentialGeneTest(dCellData[expressed_genes,],
                                        fullModelFormulaStr = "~timepoint",
                                        cores = 3)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  
  dCellData <- setOrderingFilter(dCellData, ordering_genes)
  
  dCellData <- reduceDimension(dCellData, max_components = 2,
                               method = 'DDRTree', norm_method='none')  
  
  dCellData <- orderCells(dCellData)
  
  cor.kendall = cor(dCellData@phenoData@data$Pseudotime, as.numeric(dCellData@phenoData@data$timepoint), 
                    method = "kendall", use = "complete.obs")
  
  lpsorder2 = data.frame(sample_name = colnames(count), State= dCellData@phenoData@data$State, 
                         Pseudotime = dCellData@phenoData@data$Pseudotime, rank = rank(dCellData@phenoData@data$Pseudotime))
  
  lpsorder_rank = dplyr::arrange(lpsorder2, rank)
  
  lpsorder_rank$Pseudotime = lpsorder_rank$rank
  lpsorder_rank = lpsorder_rank[-4]
  
  #row.names(lpsorder_rank) = lpsorder_rank$sample_name
  
  #i <- sapply(lpsorder_rank, is.factor)
  lpsorder_rank[1] <- lapply(lpsorder_rank[1], as.character)
  
  subpopulation <- data.frame(cell = colnames(count), sub = as.numeric(cellLabels)-1)
  
  POS <- TSCAN::orderscore(subpopulation, lpsorder_rank)[1]
  
  out = list(cor.kendall=cor.kendall, POS=POS)
  out
  
}