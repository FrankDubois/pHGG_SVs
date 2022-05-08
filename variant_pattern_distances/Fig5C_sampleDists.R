library(data.table)
dists <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20201215sampleDists/20201215sampleDistMatrix_Jacc_Signa_ordered.txt')
colnames(dists)
table(dists$HistoneSVSignaGroup)
table(dists$HistoneGroup)
dists_K27M <- dists[HistoneGroup %in%  c('H3F3A_K27M', 'HIST1_K27M')]
idx <- colnames(dists_K27M) %in% c(dists_K27M$rn, 'rn', 'HistoneGroup', 'ComplxSVgr')
dists_K27M_s <- dists_K27M[,..idx]

dists_H3.1K27M_SimpleSV <- dists[HistoneSVSignaGroup == 'H3.1_K27M_1simpleSV']
MeanDists_H3.1K27M_SimpleSV <- as.data.table(colSums(dists_H3.1K27M_SimpleSV[,c(-1, -174, -175 )])/ nrow(dists_H3.1K27M_SimpleSV), keep.rownames ='ID')
colnames(MeanDists_H3.1K27M_SimpleSV) <- c('sid', 'meanDist')
dists_an <- dists[,c(1,174,175)]
MeanDists_H3.1K27M_SimpleSVan <- merge(x = MeanDists_H3.1K27M_SimpleSV, y = dists_an, by.x = 'sid', by.y = 'rn', all = TRUE)
MeanDists_H3.1K27M_SimpleSVan[,ComparisonGroup := 'H3.1_K27M_1simpleSV']

# g <- ggplot(pHGG_vstCbatmatXs, aes(x=cluster, y = MDM4, fill = Group))
# g + geom_boxplot( varwidth =TRUE, fill = 'grey') + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3) + theme_classic() + 
#   theme(text = element_text(size=28), axis.text.x = element_text(angle = 45, hjust = 1))

dists_H3.3_K27M_1simpleSV <- dists[HistoneSVSignaGroup == 'H3.3_K27M_1simpleSV']
MeanDists_H3.3_K27M_1simpleSV <- as.data.table(colSums(dists_H3.3_K27M_1simpleSV[,c(-1, -174, -175 )])/ nrow(dists_H3.3_K27M_1simpleSV), keep.rownames ='ID')
colnames(MeanDists_H3.3_K27M_1simpleSV) <- c('sid', 'meanDist')
#dists_an <- dists[,c(1,174,175)]
MeanDists_H3.3_K27M_1simpleSVan <- merge(x = MeanDists_H3.3_K27M_1simpleSV, y = dists_an, by.x = 'sid', by.y = 'rn', all = TRUE)
MeanDists_H3.3_K27M_1simpleSVan[,ComparisonGroup := 'H3.3_K27M_1simpleSV']

dists_H3.3_K27M_complexSV <- dists[HistoneSVSignaGroup == 'H3.3_K27M_complexSV']
MeanDists_H3.3_K27M_complexSV <- as.data.table(colSums(dists_H3.3_K27M_complexSV[,c(-1, -174, -175 )])/ nrow(dists_H3.3_K27M_complexSV), keep.rownames ='ID')
colnames(MeanDists_H3.3_K27M_complexSV) <- c('sid', 'meanDist')
#dists_an <- dists[,c(1,174,175)]
MeanDists_H3.3_K27M_complexSVan <- merge(x = MeanDists_H3.3_K27M_complexSV, y = dists_an, by.x = 'sid', by.y = 'rn', all = TRUE)
MeanDists_H3.3_K27M_complexSVan[,ComparisonGroup := 'H3.3_K27M_complexSV']

MeanDists_c <- rbind(MeanDists_H3.1K27M_SimpleSVan, MeanDists_H3.3_K27M_1simpleSVan, MeanDists_H3.3_K27M_complexSVan)
MeanDists_c[, Comparison := paste0(HistoneSVSignaGroup, '_vs_', ComparisonGroup)]
table(MeanDists_c$Comparison)
MeanDists_cs <- MeanDists_c[Comparison %in% c('H3.1_K27M_1simpleSV_vs_H3.1_K27M_1simpleSV', 'H3.1_K27M_1simpleSV_vs_H3.3_K27M_1simpleSV', 
                                              'H3.3_K27M_1simpleSV_vs_H3.3_K27M_1simpleSV', 'H3.3_K27M_1simpleSV_vs_H3.3_K27M_complexSV', 'H3.3_K27M_complexSV_vs_H3.3_K27M_complexSV', 
                                              'H3.1_K27M_1simpleSV_vs_H3.3_K27M_complexSV')]
require(ggplot2)
MeanDists_cs$Comparison <- gsub(pattern = '_', replacement = '', x = MeanDists_cs$Comparison )
g <- ggplot(MeanDists_cs, aes(x=Comparison, y = meanDist))
g + geom_boxplot( varwidth =TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3, fill = 'red') + 
  theme_classic() + theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))

MeanDists_cs2 <- MeanDists_c[Comparison %in% c('H3.3_K27M_1simpleSV_vs_H3.3_K27M_1simpleSV', 'H3.3_K27M_1simpleSV_vs_H3.3_K27M_complexSV', 'H3.3_K27M_complexSV_vs_H3.3_K27M_complexSV')]
g <- ggplot(MeanDists_cs2, aes(x=Comparison, y = meanDist))
g + geom_boxplot( varwidth =TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.3, fill = 'red') + 
  theme_classic() + theme(text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1))
write.table(MeanDists_cs2, '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220308figuretables/20220312Fig5C_20210106sampleDists.txt',quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)


Comparisons =unique(MeanDists_cs$Comparison)
xG = Comparisons[1]
WX_res <- lapply(Comparisons, function(xG){
  # Comparison
  print(xG)
  CL_memb <- MeanDists_cs[Comparison == xG,]
  x= MeanDists_cs$w[1]
  Comparisons_N <- Comparisons[Comparisons != xG]
  WX_resN <- rbindlist(lapply(Comparisons_N, function(xN){
    Wilcoxon_res <- wilcox.test(CL_memb[, meanDist], MeanDists_cs[Comparison == xN, meanDist], conf.int = TRUE)
    Wilcoxon_res_dt <- as.data.table(paste0( xG, '_VS_', xN))
    Wilcoxon_res_dt[,p_value := Wilcoxon_res$p.value]
    Wilcoxon_res_dt[,conf_int_low :=  Wilcoxon_res$conf.int[1]]
    Wilcoxon_res_dt[,conf_int_high :=  Wilcoxon_res$conf.int[2]]
    Wilcoxon_res_dt[,meanshift :=  Wilcoxon_res$estimate]
    Wilcoxon_res_dt[,ratio := mean(na.omit(CL_memb[, meanDist]))/mean(na.omit(MeanDists_cs[Comparison == xN, meanDist]))]
    return(Wilcoxon_res_dt)
  }))
  WX_resN[, sign_ratio := ifelse(ratio>1, "Enriched", "Exclusive")]
  return(WX_resN)
})
WX_res_dt <- rbindlist(WX_res)
require(qvalue)
qvalue(WX_res_dt$p_value, pi0 = 1)
WX_res_dt[p_value >= 1, p_value := 1]
WX_res_dt[p_value <= 0, p_value := 0]
WX_res_dts <- WX_res_dt[!duplicated(p_value)]
WX_res_dts[, qval := qvalue(p_value, pi0 = 1)$qvalues]
write.table(WX_res_dts, '/Users/fdubois/Dropbox/pHGG_data/data/20210106sampleDists/20210120sampleDists_byHistoneSVSignaGroup.txt',quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)


