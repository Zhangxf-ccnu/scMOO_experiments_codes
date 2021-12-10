############### run monocle 2 on Deng data ####################
### load imputation result of scMOO
deng <- readRDS('./deng_scmoo.rds')
deng_cellLabels = factor(colnames(deng),
                         levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                  '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                  'mid blastocyst', 'late blastocyst'))

deng.scmoo.monocle = my.monocle(deng, deng_cellLabels)

############### run monocle 2 on Petropoulos data ####################
### load imputation result of scMOO
petro <- readRDS('./petro_scmoo.rds')
petro_cellLabels = factor(colnames(petro),
                         levels=c('E3', 'E4', 'E5', 'E6', 'E7'))

petro.scmoo.monocle = my.monocle(petro, petro_cellLabels)

