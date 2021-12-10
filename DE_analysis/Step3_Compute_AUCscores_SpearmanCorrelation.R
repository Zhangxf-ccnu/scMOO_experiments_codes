################## compute AUC scores under different standards #################### 

bulk_data = as.matrix(bulk_data)[filter_ID ,]
bulk_deg = my.edgeR(bulk_data, bulk_group)$res.edgeR.sort

sc_data = as.matrix(sc_data)[filter_ID ,]
raw_deg = my.edgeR(sc_data, sc_group)$res.edgeR.sort

#### load imputation results of different methods ########
scMOO <- readRDS('./scmoo_imp.rds')
ALRA <- readRDS('./alra_imp.rds')
DrImpute <- readRDS('./drimpute_imp.rds')
SAVER <- readRDS('./saver_imp.rds')
VIPER <- readRDS('./viper_imp.rds')
MAGIC <- readRDS('./magic_imp.rds')
scVI <- readRDS('./scvi_imp.rds')
scImpute <- readRDS('./scimpute_imp.rds')
DCA <- readRDS('./dca_imp.rds')
SAVERX <- readRDS('./saverx_imp.rds')
WEDGE <- readRDS('./wedge_imp.rds')


scmoo_deg <- my.edgeR(scMOO, sc_group)$res.edgeR.sort  
alra_deg <- my.edgeR(ALRA, sc_group)$res.edgeR.sort
drimpute_deg <- my.edgeR(DrImpute, sc_group)$res.edgeR.sort
saver_deg <- my.edgeR(SAVER, sc_group)$res.edgeR.sort
viper_deg <- my.edgeR(VIPER, sc_group)$res.edgeR.sort
magic_deg <- my.edgeR(MAGIC, sc_group)$res.edgeR.sort
scvi_deg <- my.edgeR(scVI, sc_group)$res.edgeR.sort
scimpute_deg <- my.edgeR(scImpute, sc_group)$res.edgeR.sort
dca_deg <- my.edgeR(DCA, sc_group)$res.edgeR.sort
saverx_deg <- my.edgeR(SAVERX, sc_group)$res.edgeR.sort
wedge_deg <- my.edgeR(WEDGE, sc_group)$res.edgeR.sort



bulk_pvalue <- bulk_deg$p.adjust
raw_pvalue <- raw_deg$p.adjust
scmoo_pvalue <- scmoo_deg$p.adjust
dca_pvalue <- dca_deg$p.adjust
saverx_pvalue <- saverx_deg$p.adjust
alra_pvalue <- alra_deg$p.adjust
drimpute_pvalue <- drimpute_deg$p.adjust
saver_pvalue <- saver_deg$p.adjust
viper_pvalue <- viper_deg$p.adjust
magic_pvalue <- magic_deg$p.adjust
scvi_pvalue <- scvi_deg$p.adjust
scimpute_pvalue <- scimpute_deg$p.adjust
wedge_pvalue <- wedge_deg$p.adjust


n<-200  #### 400, 600, 800, 1000


bulk_tpr <- c()
raw_tpr <- c()
saver_tpr <- c()
scimpute_tpr <- c()
viper_tpr <- c()
alra_tpr <- c()
drimpute_tpr<-c()
magic_tpr<-c()
scvi_tpr<-c()
scmoo_tpr<-c()
dca_tpr<-c()
saverx_tpr<-c()
wedge_tpr<-c()

################

bulk_fpr <- c()
raw_fpr <- c()
saver_fpr <- c()
scimpute_fpr <- c()
viper_fpr <- c()
alra_fpr <- c()
drimpute_fpr<-c()
magic_fpr<-c()
scvi_fpr<-c()
scmoo_fpr<-c()
dca_fpr<-c()
saverx_fpr<-c()
wedge_fpr<-c()



degene = bulk_deg$GeneID[1:n]
nodegene = bulk_deg$GeneID[-c(1:n)]




for (i in 1:2000){
  saver_tpr[i] <- length(intersect(saver_deg[saver_deg$p.adjust<=saver_pvalue[i],]$GeneID,
                                  degene))/length(degene)
  scimpute_tpr[i] <- length(intersect(scimpute_deg[scimpute_deg$p.adjust<=scimpute_pvalue[i],]$GeneID,
                                degene))/length(degene)
  viper_tpr[i] <-length(intersect(viper_deg[viper_deg$p.adjust<=viper_pvalue[i],]$GeneID,
                                  degene))/length(degene)
  raw_tpr[i]<- length(intersect(raw_deg[raw_deg$p.adjust<=raw_pvalue[i],]$GeneID,
                                degene))/length(degene)
  alra_tpr[i] <-length(intersect(alra_deg[alra_deg$p.adjust<=alra_pvalue[i],]$GeneID,
                                 degene))/length(degene)
  magic_tpr[i] <- length(intersect(magic_deg[magic_deg$p.adjust<=magic_pvalue[i],]$GeneID,
                                  degene))/length(degene)
  drimpute_tpr[i] <- length(intersect(drimpute_deg[drimpute_deg$p.adjust<=drimpute_pvalue[i],]$GeneID,
                                degene))/length(degene)
  scvi_tpr[i] <- length(intersect(scvi_deg[scvi_deg$p.adjust<=scvi_pvalue[i],]$GeneID,
                                 degene))/length(degene)
  scmoo_tpr[i] <- length(intersect(scmoo_deg[scmoo_deg$p.adjust<=scmoo_pvalue[i],]$GeneID,
                                degene))/length(degene)
  dca_tpr[i] <- length(intersect(dca_deg[dca_deg$p.adjust<=dca_pvalue[i],]$GeneID,
                                degene))/length(degene)
  saverx_tpr[i] <- length(intersect(saverx_deg[saverx_deg$p.adjust<=saverx_pvalue[i],]$GeneID,
                                   degene))/length(degene)
  wedge_tpr[i] <- length(intersect(wedge_deg[wedge_deg$p.adjust<=wedge_pvalue[i],]$GeneID,
                                  degene))/length(degene)
  
  
  saver_fpr[i] <- length(intersect(saver_deg[saver_deg$p.adjust<=saver_pvalue[i],]$GeneID,
                                  nodegene))/length(nodegene)
  scimpute_fpr[i] <- length(intersect(scimpute_deg[scimpute_deg$p.adjust<=scimpute_pvalue[i],]$GeneID,
                                nodegene))/length(nodegene)
  viper_fpr[i] <- length(intersect(viper_deg[viper_deg$p.adjust<=viper_pvalue[i],]$GeneID,
                                  nodegene))/length(nodegene)
  raw_fpr[i] <- length(intersect(raw_deg[raw_deg$p.adjust<=raw_pvalue[i],]$GeneID,
                                nodegene))/length(nodegene)
  alra_fpr[i] <- length(intersect(alra_deg[alra_deg$p.adjust<=alra_pvalue[i],]$GeneID,
                                 nodegene))/length(nodegene)
  magic_fpr[i] <- length(intersect(magic_deg[magic_deg$p.adjust<=magic_pvalue[i],]$GeneID,
                                  nodegene))/length(nodegene)
  drimpute_fpr[i] <- length(intersect(drimpute_deg[drimpute_deg$p.adjust<=drimpute_pvalue[i],]$GeneID,
                                nodegene))/length(nodegene)
  scvi_fpr[i] <- length(intersect(scvi_deg[scvi_deg$p.adjust<=scvi_pvalue[i],]$GeneID,
                                 nodegene))/length(nodegene)
  scmoo_fpr[i] <- length(intersect(scmoo_deg[scmoo_deg$p.adjust<=scmoo_pvalue[i],]$GeneID,
                                nodegene))/length(nodegene)
  dca_fpr[i] <- length(intersect(dca_deg[dca_deg$p.adjust<=dca_pvalue[i],]$GeneID,
                                nodegene))/length(nodegene)
  saverx_fpr[i] <- length(intersect(saverx_deg[saverx_deg$p.adjust<=saverx_pvalue[i],]$GeneID,
                                   nodegene))/length(nodegene)
  wedge_fpr[i] <- length(intersect(wedge_deg[wedge_deg$p.adjust<=wedge_pvalue[i],]$GeneID,
                                  nodegene))/length(nodegene)
}



saver_auc <- 0
scimpute_auc <- 0
viper_auc <- 0
raw_auc <- 0
alra_auc <- 0
magic_auc <- 0
scvi_auc <- 0
drimpute_auc <- 0
scmoo_auc <- 0
dca_auc <- 0
saverx_auc <- 0
wedge_auc <- 0


for (i in 2:2000){
  saver_auc <- saver_auc + saver_tpr[i]*(saver_fpr[i]-saver_fpr[i-1])
  scimpute_auc <- scimpute_auc + scimpute_tpr[i]*(scimpute_fpr[i]-scimpute_fpr[i-1])
  viper_auc <- viper_auc + viper_tpr[i]*(viper_fpr[i]-viper_fpr[i-1])
  raw_auc <- raw_auc + raw_tpr[i]*(raw_fpr[i]-raw_fpr[i-1])
  alra_auc <- alra_auc + alra_tpr[i]*(alra_fpr[i]-alra_fpr[i-1])
  magic_auc <- magic_auc + magic_tpr[i]*(magic_fpr[i]-magic_fpr[i-1])
  scvi_auc <- scvi_auc + scvi_tpr[i]*(scvi_fpr[i]-scvi_fpr[i-1])
  drimpute_auc <- drimpute_auc + drimpute_tpr[i]*(drimpute_fpr[i]-drimpute_fpr[i-1])
  scmoo_auc <- scmoo_auc + scmoo_tpr[i]*(scmoo_fpr[i]-scmoo_fpr[i-1])
  dca_auc <- dca_auc + dca_tpr[i]*(dca_fpr[i]-dca_fpr[i-1])
  saverx_auc <- saverx_auc + saverx_tpr[i]*(saverx_fpr[i]-saverx_fpr[i-1])
  wedge_auc <- wedge_auc + wedge_tpr[i]*(wedge_fpr[i]-wedge_fpr[i-1])
  
}




auc_df_rank <- data.frame(Method = c('SAVER',"scImpute",  'VIPER', 'Observed', 'scMOO', 'ALRA',
                                       'MAGIC',  'DrImpute','DCA', 'SAVERX','scVI', 'WEDGE'),
                            AUC = c(saver_auc, scimpute_auc, viper_auc, raw_auc, scmoo_auc, 
                                    alra_auc, magic_auc, drimpute_auc, dca_auc, 
                                    saverx_auc, scvi_auc, wedge_auc))



################## compute Spearman correlation coefficient #################
saver_padj <- saver_deg$p.adjust[match(bulk_deg$GeneID, saver_deg$GeneID)]
raw_padj<- raw_deg$p.adjust[match(bulk_deg$GeneID, raw_deg$GeneID)]
scimpute_padj <- scimpute_deg$p.adjust[match(bulk_deg$GeneID, scimpute_deg$GeneID)]
viper_padj <- viper_deg$p.adjust[match(bulk_deg$GeneID, viper_deg$GeneID)]
alra_padj <- alra_deg$p.adjust[match(bulk_deg$GeneID, alra_deg$GeneID)]
drimpute_padj <- drimpute_deg$p.adjust[match(bulk_deg$GeneID, drimpute_deg$GeneID)]
scvi_padj <- scvi_deg$p.adjust[match(bulk_deg$GeneID, scvi_deg$GeneID)]
magic_padj <-magic_deg$p.adjust[match(bulk_deg$GeneID, magic_deg$GeneID)]
dca_padj <- dca_deg$p.adjust[match(bulk_deg$GeneID, dca_deg$GeneID)]
saverx_padj <- saverx_deg$p.adjust[match(bulk_deg$GeneID, saverx_deg$GeneID)]
wedge_padj <- wedge_deg$p.adjust[match(bulk_deg$GeneID, wedge_deg$GeneID)]
scmoo_padj = scmoo_deg$p.adjust[match(bulk_deg$GeneID, scmoo_deg$GeneID)]


cor.spearman = c()
bulk_padj = bulk_deg$p.adjust

cor.spearman[1] = cor(bulk_padj, saver_padj, method = c('spearman'))
cor.spearman[2] = cor(bulk_padj, raw_padj, method = c('spearman'))
cor.spearman[3] = cor(bulk_padj, scimpute_padj, method = c('spearman'))
cor.spearman[4] = cor(bulk_padj, viper_padj, method = c('spearman'))
cor.spearman[5] = cor(bulk_padj, alra_padj, method = c('spearman'))
cor.spearman[6] = cor(bulk_padj, drimpute_padj, method = c('spearman'))
cor.spearman[7] = cor(bulk_padj, magic_padj, method = c('spearman'))
cor.spearman[8] = cor(bulk_padj, scmoo_padj, method = c('spearman'))
cor.spearman[9] = cor(bulk_padj, dca_padj, method = c('spearman'))
cor.spearman[10] = cor(bulk_padj, saverx_padj, method = c('spearman'))
cor.spearman[11] = cor(bulk_padj, scvi_padj, method = c('spearman'))
cor.spearman[12] = cor(bulk_padj, wedge_padj, method = c('spearman'))










