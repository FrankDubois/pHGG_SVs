library(data.table)
library("RColorBrewer")
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
VarMat <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20201215comut_byComplexSVsigL20%/20201215comut_byComplexSVSigL20%_JaccClst.gct', skip = 'id')[-c(1:2),]
# take out Histone.snvs
VarMat_s <- VarMat[!(id %in% c('H3F3A:p.K28M', 'HIST1H3B:p.K28M', 'H3F3A:missense.snv', 'HIST1H3C:missense.snv', 'HIST1H3B:p.R27R'))]

TvarMat <- t(VarMat_s[,-1])
# sampleDists <- dist(x = TvarMat, upper = F,diag = T )
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix)
# pheatmap(sampleDistMatrix, cluster_rows =F , cluster_cols = F,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
library(ade4)
TvarMatnoNA <- na.omit(TvarMat)
df2 <- data.frame(apply(TvarMatnoNA, 2, function(x) as.numeric(as.character(x))))
rownames(df2) <- rownames(TvarMatnoNA)
df2[1:10, 1:10]
D2 <- dist.binary(df = na.omit(df2), method = 1, diag = T, upper = F)
sampleDistMatrix_Jacc <- as.matrix(D2)
rownames(sampleDistMatrix_Jacc)

pheatmap(sampleDistMatrix_Jacc, cluster_rows =F , cluster_cols = F,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

sampleDistMatrix_Jacc_dt <- as.data.table(sampleDistMatrix_Jacc, keep.rownames = TRUE)
metaD <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20201215sampleDists/20201215metaD_Sig3.txt')
sampleDistMatrix_Jacc_dt$rn[!(sampleDistMatrix_Jacc_dt$rn %in% metaD$Tumor_Sample_Barcode)]
sampleDistMatrix_Jacc_dt[, HistoneGroup := metaD[rn == Tumor_Sample_Barcode, Histone_group], by = rn]  
sampleDistMatrix_Jacc_dt[, HistoneSVSignaGroup := metaD[rn == Tumor_Sample_Barcode, Group], by = rn]  

sampleDistMatrix_Jacc_dt[, rowID :=  paste0(rn, ':', HistoneGroup)]
write.table(sampleDistMatrix_Jacc_dt,"/Users/fdubois/Dropbox/pHGG_data/data/20201215sampleDists/20201215sampleDistMatrix_Jacc_Signa_ordered.txt",quote=F,row.names=F,sep="\t")

rownames(sampleDistMatrix_Jacc) <- paste0(sampleDistMatrix_Jacc_dt$rn, ':', sampleDistMatrix_Jacc_dt$HistoneGroup)

pheatmap(sampleDistMatrix_Jacc, cluster_rows =F , cluster_cols = F,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
