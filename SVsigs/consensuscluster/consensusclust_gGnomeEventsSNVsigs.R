library(data.table)
load('/Users/fdubois/Dropbox/pHGG_Data/data/20191114SV_SNVsig_correlation/L1KL.lego96.DIPG.10.MAP14.WH.RData')
H_norm_SNV <- WH[[3]]
H_SNV <- WH[[2]]
rowSums(H_norm_SNV)
H_SNVnormalized <- apply(H_norm_SNV,2,function(x) x/sum(x))
colSums(H_SNVnormalized)

#gGnome SVsigs
load('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220203gGnomeSVsigs/L1KL.SV.pHGG.phi.1.7.RData')
WgGnome <- res[[1]]
HgGnome <- res[[2]]
indexgGnome <- colSums(WgGnome) > 1
KgGnome <- sum(indexgGnome)
WgGnome <- WgGnome[,indexgGnome]
HgGnome <- HgGnome[indexgGnome,]
colnames(WgGnome) <- paste("W",seq(1:ncol(WgGnome)),sep="")
rownames(HgGnome) <- colnames(WgGnome)
WgGnome.norm <- WgGnome
HgGnome <- HgGnome
for (i in 1:KgGnome) {
  WgGnome.norm[,i] <- WgGnome[,i]*rowSums(HgGnome)[i]
  HgGnome[i,] <- HgGnome[i,]*colSums(WgGnome)[i]
}
HgGnome <- apply(HgGnome,2,function(x) x/sum(x))
colSums(HgGnome)

H_SNVnormalized_mctred = sweep(H_SNVnormalized,1, apply(H_SNVnormalized,1,median,na.rm=T))
HgGnome_mctred = sweep(HgGnome,1, apply(HgGnome,1,median,na.rm=T))
library(ConsensusClusterPlus)
sc2 <- as.data.table(cbind(colnames(HgGnome), colnames(H_SNVnormalized)))
sc <- as.data.table(cbind(sort(colnames(HgGnome)), sort(colnames(H_SNVnormalized))))
colnames(sc) = c('gGnomeSVsigs', 'SNVsigs')

T_H_SNVnormalized_mctred <- as.data.table(t(H_SNVnormalized_mctred), keep.rownames = T)
T_H_SNVnormalized_mctred[, ID:= sc[rn==SNVsigs, gGnomeSVsigs], by = rn]
tT_H_SNVnormalized_mctred <- t(T_H_SNVnormalized_mctred[,-1])
colnames(tT_H_SNVnormalized_mctred) <- tT_H_SNVnormalized_mctred[15,]
tT_H_SNVnormalized_mctred <- tT_H_SNVnormalized_mctred[-15,]
sc2 <- as.data.table(cbind(colnames(HgGnome_mctred), colnames(tT_H_SNVnormalized_mctred)))
HgGnome_mctred_dt <- as.data.table(HgGnome_mctred, keep.rownames = T)
HgGnome_mctred_dt$rn<- gsub(pattern = 'W1', replacement = 'complex_unclassf.', x = HgGnome_mctred_dt$rn)
HgGnome_mctred_dt$rn<- gsub(pattern = 'W2', replacement = 'del', x = HgGnome_mctred_dt$rn)
HgGnome_mctred_dt$rn<- gsub(pattern = 'W3', replacement = 'dup', x = HgGnome_mctred_dt$rn)
HgGnome_mctred_dt$rn<- gsub(pattern = 'W4', replacement = 'inv', x = HgGnome_mctred_dt$rn)
HgGnome_mctred_dt$rn<- gsub(pattern = 'W5', replacement = 'tra', x = HgGnome_mctred_dt$rn)
tT_H_SNVnormalized_mctred <- as.data.table(tT_H_SNVnormalized_mctred, keep.rownames = T)
HgGnomeSNV_mctred <- rbind(HgGnome_mctred_dt, tT_H_SNVnormalized_mctred, use.names=TRUE)
write.table(HgGnomeSNV_mctred, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220115gGnomeEventSigSNVSigConsensusCluster/20220115_HgGnomeSNV_mctred.txt", sep ="\t", row.names = F, col.names = T, quote = F)

HgGnomeSNV_mctredMat <-  apply(HgGnomeSNV_mctred[,-1],2,FUN = as.numeric)
rownames(HgGnomeSNV_mctredMat) <- HgGnomeSNV_mctred[,1]

outputpths <- '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220115gGnomeEventSigSNVSigConsensusCluster/'
results <-  ConsensusClusterPlus(HgGnomeSNV_mctredMat,maxK=10,reps=1000,pItem=0.9,pFeature=0.9,
                                 title=outputpths,clusterAlg="hc",distance="spearman",seed=1262118388.71279, verbose = TRUE, plot = 'pdf')

cl = calcICL(results,title=outputpths,plot="pdf")
resultsCL2 <- results[[2]]
resCL2memb <- resultsCL2$consensusClass
resultsCL3 <- results[[3]]
resCL3memb <- resultsCL3$consensusClass
resultsCL4 <- results[[4]]
resCL4memb <- resultsCL4$consensusClass
resCl2_4 <- cbind(resCL2memb, resCL3memb, resCL4memb)

write.table(resCl2_4, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220115gGnomeEventSigSNVSigConsensusCluster/20220115_consensuscluster_Membership.txt", sep ="\t", row.names = T, col.names = T, quote = F)
resCl_ord <- fread('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220115gGnomeEventSigSNVSigConsensusCluster/20220115_consensuscluster_Membership_ord.txt')
resCl_ord <- resCl_ord[order(resCl_ord$resCL2memb), ] # sort membership samples by cluster

orderindx <- sapply(resCl_ord$sample, function(x){
  which(x == colnames(HgGnomeSNV_mctred))
})
#order original heatmap and safe
HgGnomeSNV_mctredOrdered <- HgGnomeSNV_mctred[,..orderindx]
HgGnomeSNV_mctredOrdered <- cbind(HgGnomeSNV_mctred[,1],HgGnomeSNV_mctredOrdered)
write.table(HgGnomeSNV_mctredOrdered, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220115gGnomeEventSigSNVSigConsensusCluster/20220115HgGnomeSNVmctredOrdClst2.txt", sep ="\t", row.names = F, col.names = T, quote = F)

