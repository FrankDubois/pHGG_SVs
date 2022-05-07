library(mobster)
library(tidyr)
library(dplyr)
library(data.table)
library(parallel)
absmaf <- readRDS('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20200927mobster/20200927absmaf_mobst.rds')

x = unique(absmaf$TScode)[34]
lapply(unique(absmaf$TScode)[21:173], function(x){
  print(x)
  absmafX <- absmaf[TScode == x]
  mobsX <- mobster_fit(x = absmafX)
  saveRDS(mobsX, paste0('/xchip/beroukhimlab/Frank/DIPG/data/20200928mobster/20200928', x, 'mobst.rds'))
})

 
