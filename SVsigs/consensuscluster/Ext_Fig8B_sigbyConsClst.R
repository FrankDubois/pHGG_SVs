library(data.table)
sv_snv_mat <- fread('/Volumes/xchip_beroukhimlab/Alex/pHGG/20220208_conssensuscluster/20220208_SVSNV_consencusclusterinputmat.txt')
resCl2_4 <- fread('/Volumes/xchip_beroukhimlab/Alex/pHGG/20220208_conssensuscluster/20220210_rescl2-4.txt')

i.sig = 1
dd <- lapply(seq(nrow(sv_snv_mat)), function(i.sig) {
  print(sv_snv_mat$V1[i.sig])
  idx1 <- colnames(sv_snv_mat) %in% resCl2_4[V2 == 1, V1]
  idx2 <- colnames(sv_snv_mat) %in% resCl2_4[V2 == 2, V1]
  # idx1 == !idx2
  wicxTest <-  wilcox.test(x = as.numeric(sv_snv_mat[i.sig, ..idx1]), y = as.numeric(sv_snv_mat[i.sig, ..idx2]), conf.int = T)#, conf.int =TRUE)
  ratio = mean(as.numeric(sv_snv_mat[i.sig, ..idx1]))/mean(as.numeric(sv_snv_mat[i.sig, ..idx2]))
  resW <- as.data.table(cbind(sv_snv_mat$V1[i.sig], wicxTest$p.value, ratio)) #,as.numeric(wicxTest$statistic)
  colnames(resW)[2] = 'p_value'
  return(resW) 
})
library(qvalue)
gene.list <- rbindlist(dd)
colnames(gene.list) = c('signature', 'p_value', 'Complex_vs_Simple_ratio')
gene.list$Complex_vs_Simple_ratio <- as.numeric(gene.list$Complex_vs_Simple_ratio)
gene.list[, sign.ratio := ifelse(Complex_vs_Simple_ratio>1, "High Activity", "Low Activity")]
gene.list[Complex_vs_Simple_ratio < 1, sign.ratio :=  "Low Activity"]
gene.list$p_value <- as.numeric(gene.list$p_value)
gene.list [p_value > 1, p_value := 1]
gene.list [p_value < 0, p_value := 0]
gene.list[, qvalue := qvalue(p_value)$qvalue]  

write.table(gene.list, "//Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220213SigbyConsClst/20220213SigbyConsClst.txt",
            sep ="\t", quote = F, col.names = T, row.names = F)

require(ggplot2)
require(ggrepel)
gene.list[, log2OR := log2(Complex_vs_Simple_ratio)]
gene.list[qvalue<0.2, labelf := signature]

ggplot(gene.list, aes(x=log2OR, y=-log10(qvalue), color = factor(sign.ratio)), size = incidence) + geom_point() +
  geom_text_repel(aes(label = labelf, vjust = 0.2, hjust = 0.5, size = 30)) + theme_classic() + theme(text = element_text(size = 25)) + xlab("log2OR")

ggplot(gene.list[p_value<1.1e-4], aes(x=SRSVsGroup, y=-log10(qvalue), color = factor(sign_ratio))) + geom_jitter(width = 0, height = 0.1) +
  geom_text_repel(aes(label = V1, vjust = 0.2, hjust = 0.5,  size = 30)) + theme_light() + theme(text = element_text(size = 16)) + xlab("SRSVsGroup")
