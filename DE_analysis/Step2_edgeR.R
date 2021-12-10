######################### edgeR function ###########################
library(edgeR)

my.edgeR <- function(count, group){
  
  # round data if there are non-integer columns
  intcol<-vector("logical")
  for(i in 1:ncol(count)){
    intcol<-c(intcol,is.integer(count[,i]))
  }
  if (!all(TRUE==intcol)) {
    warning("WARNING! Non-integer expression levels. Data rounded")
    count<-round(count)
  }
  
  # remove out zero count genes in all samples
  #count.sel = count[rowSums(count)>0,]
  count.sel = count
  # creat main edgeR object
  d = DGEList(count.sel, group=group)
  
  # do default normalisation for edgeR
  d = calcNormFactors(d)
  # estimate the common and tagwise dispersion
  d = estimateCommonDisp(d)
  d = estimateTagwiseDisp(d)
  
  # determine differentially expressed genes (using exact test)
  dest = exactTest(d)
  
  dest.sort = topTags(dest, n = dim(dest)[1])
  
  # need to extract table data from de.com object,
  # then select only required columns, plus calculate adjust p-value
  res.edgeR <- data.frame(rownames(dest$table),round(dest$table$logFC,digits=2),signif(p.adjust(dest$table$PValue,method='BH'),digits=3),signif(dest$table$PValue,digits=3))
  
  res.edgeR.sort <- data.frame(rownames(dest.sort$table),round(dest.sort$table$logFC,digits=2),signif(p.adjust(dest.sort$table$PValue,method='BH'),digits=3),signif(dest.sort$table$PValue,digits=3))
  
  
  # name
  names(res.edgeR) <- c('GeneID','log2FC','p.adjust','pvalue')
  names(res.edgeR.sort) <- c('GeneID','log2FC','p.adjust','pvalue')
  
  out = list(res.edgeR=res.edgeR, res.edgeR.sort = res.edgeR.sort)
  out
}
