load('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191110SVsig/output/L1KL.SV.DIPG.phi.1.2.RData')
W <- res[[1]]
H <- res[[2]]
index <- colSums(W) > 1
K <- sum(index)
W <- W[,index]
H <- H[index,]
colnames(W) <- paste("W",seq(1:ncol(W)),sep="")
rownames(H) <- colnames(W)
pheatmap(H)
W.norm <- W
H.mid <- H
for (i in 1:K) {
  W.norm[,i] <- W[,i]*rowSums(H)[i]
  H.mid[i,] <- H[i,]*colSums(W)[i]
}
H.norm_SV <- apply(H.mid,2,function(x) x/sum(x))
num.sig <- dim(H.mid)[1]
require(qvalue)
require(data.table)

load('/Users/fdubois/Dropbox/pHGG_data/data/20191114SV_SNVsig_correlation/L1KL.lego96.DIPG.10.MAP14.WH.RData')
H_norm_SNV <- WH[[3]]
H_SNV <- WH[[2]]

colnames(H)[colnames(H) %in% colnames(H_SNV)]

sc <- cbind(sort(colnames(H)), sort(colnames(H_SNV)))

H_SNV_s <- H_SNV[,colnames(H_SNV)!='DIPG18T']
sc2 <- cbind(sort(colnames(H)), sort(colnames(H_SNV_s)))

orderindxSV <- sapply(sort(colnames(H)), function(x){
  which(x == colnames(H)) 
})
H_sorted <- H[,orderindxSV]

orderindxSNV <- sapply(sort(colnames(H_SNV_s)), function(x){
  which(x == colnames(H_SNV_s)) 
})
H_SNVs_sorted <- H_SNV_s[,orderindxSNV]

colnames(H_sorted) = colnames(H_SNVs_sorted)
H_sorted <- as.data.table(H_sorted, keep.rownames = TRUE)
H_SNVs_sorted <- as.data.table(H_SNVs_sorted, keep.rownames = TRUE)

H_sorted$rn <- gsub('W5', 'simple short dup', H_sorted$rn)
H_sorted$rn <- gsub('W4', 'complex large', H_sorted$rn)
H_sorted$rn <- gsub('W3', 'simple interchromosomal', H_sorted$rn)
H_sorted$rn <- gsub('W1', 'complex interchromosomal', H_sorted$rn)
H_sorted$rn <- gsub('W2', 'mixed simple', H_sorted$rn)

setkey(H_sorted, rn)
x= H_sorted$rn[1]
SVsigcorr_list <- lapply(H_sorted$rn, function(x){
  print(x)
  SVsigvals <- as.numeric(H_sorted[x,2:ncol(H_sorted)])
  y =1
  corrsSNV <- sapply(1:nrow(H_SNVs_sorted), function(y){
    SNV_sigvals <- as.numeric(H_SNVs_sorted[y,])
    corvals <- cor(x= SVsigvals, y= SNV_sigvals, method = 'spearman')
  })
  corSV <- cbind(x, as.data.table(t(corrsSNV)))
  colnames(corSV) <- c('SVsig', 1:14)
  return(corSV)
})
SVsnv_sigCordt <- rbindlist(SVsigcorr_list)
require(pheatmap)

SVsnv_sigCordt_mat <- as.matrix(SVsnv_sigCordt, rownames = 'SVsig')
max(SVsnv_sigCordt_mat)

pheatmap(SVsnv_sigCordt_mat, fontsize = 20)
# write.table(SVsnv_sigCordt_mat,"/Users/fdubois/Dropbox/pHGG_data/data/20200705SV_SNV_signaCorr/20200109SNVSVsigCorrMat.txt",quote=F,row.names=T,sep="\t")


SVsigcorrTest_list <- rbindlist(lapply(H_sorted$rn, function(x){
  print(x)
  SVsigvals <- as.numeric(H_sorted[x,2:ncol(H_sorted)])
  y =1
  corrsTestSNV <- rbindlist(lapply(1:nrow(H_SNVs_sorted), function(y){
    SNV_sigvals <- as.numeric(H_SNVs_sorted[y,])
    corTestvals <- cor.test(x= SVsigvals, y= SNV_sigvals, method = 'spearman', alternative = "two.sided")
    corTestvals_dt <- as.data.table(cbind(paste0('SNVsig:', y), corTestvals$p.value,corTestvals$estimate))
    colnames(corTestvals_dt) <- c('SNVsig', 'p_value', 'rho')
    return(corTestvals_dt)
    }))
  
  corTestSV <- cbind(x,corrsTestSNV)
  colnames(corTestSV)[1] <- c('SVsig')
  return(corTestSV)
}))
require(qvalue)
SVsigcorrTest_list$p_value <- as.numeric(SVsigcorrTest_list$p_value)
SVsigcorrTest_list [p_value > 1, p_value := 1]
SVsigcorrTest_list [p_value < 0, p_value := 0]
SVsigcorrTest_list[, qvalue := qvalue(p_value)$qvalue]
SVsigcorrTest_list[,log_q := -log10(qvalue)]
SVsigcorrTest_list[, sign_factor := ifelse(rho>0, 1, -1)]
SVsigcorrTest_list[, Log_qvalue_signed := sign_factor * log_q]


q_heatmap <- dcast(SVsigcorrTest_list,  SVsig ~ SNVsig, value.var = 'Log_qvalue_signed')
write.table(SVsigcorrTest_list,"...write_to_output_path...",quote=F,row.names=F,sep="\t")



