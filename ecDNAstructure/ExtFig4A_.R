require(data.table)
require(gUtils)
require(BSgenome.Hsapiens.UCSC.hg19)
require(gTrack)
gt.ge <- track.gencode()
#####
#subset magicList to samples with amps in PDGFRATAD
magiclist <- readRDS('/Users/fdubois/Dropbox/pHGG_data/data/20200228magiclistQC/20200228magiclistall.RDS')
colnames(magiclist)
#magiclist[,Start_Chr := as.numeric(Start_Chr)]
#Tad&genes
GZtad_gr_an <- readRDS('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20200127fhook/GZtad_gr_an.rds')
CGC <- readRDS("/Users/fdubois/Dropbox/pHGG_data/data/20190829SNVnhmmaf_dataprep/20190829cancergene_elements.rds")
colnames(CGC)
CGC_gr <- gr.sub(dt2gr(CGC[,c(1:3,6)]), 'chr', '')
# get number of samples with amp in each tile
GATKaCN <- fread('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20190914magicList/20190914GATKaCNseg.txt')
ICGCaCN <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20200222absCN_ICGCpGBM/20200222GATKaCNseg_aCNpladj.seg')
colnames(ICGCaCN) <- colnames(GATKaCN)[c(1:4,6)]
GATKaCN <- rbind(GATKaCN, ICGCaCN, fill = TRUE)
#match names
samplesMatch_all <- as.data.table(cbind(sort(unique(GATKaCN$track.name)), 
                                        sort(unique(magiclist$Tumor_Sample_Barcode))))#[!(Tumor_Sample_Barcode %in% c('GBM36_T', 'GBM84_T', 'GBM96_T', 'GBM97_T', 'GBM98_T')),Tumor_Sample_Barcode]))))

# samplesMatch <- as.data.table(cbind(sort(unique(GATKaCN$track.name)), 
#                                     sort(unique(magiclistCN_EGFR[!(Tumor_Sample_Barcode %in% c('GBM36_T', 'GBM84_T', 'GBM96_T', 'GBM97_T', 'GBM98_T')),Tumor_Sample_Barcode])))) 
colnames(samplesMatch_all) <- c('CN_pair', 'magiclist_Tcode')
GATKaCN[ ,Tumor_Sample_Barcode := samplesMatch_all[CN_pair == track.name, magiclist_Tcode], by = track.name]


###set CNcutoff####
CNcutoff = 2
#####
#class(magiclist$End_Chr)
magiclistCN_PDGFRA<- magiclist[Variant == 'CN' & (Start_Chr == 4 | End_Chr == 4)] ## SET Manual
#unique(magiclistCN_PDGFRA$End_Chr)
magiclistCN_PDGFRA[,Variant_Type := as.numeric(Variant_Type)]
unique(magiclistCN_PDGFRA$Tumor_Sample_Barcode)
#####
magiclistCN_PDGFRA[as.numeric(Variant_Type) >= CNcutoff, Amp := 1]
magiclistCN_PDGFRAamp <- magiclistCN_PDGFRA[Amp ==1]
unique(magiclistCN_PDGFRAamp$Tumor_Sample_Barcode)
sort(table(magiclistCN_PDGFRAamp$Hugo_Symbol), decreasing = TRUE)
min(magiclistCN_PDGFRAamp$Start_Pos)
max(magiclistCN_PDGFRAamp$End_Pos)
hist(magiclistCN_PDGFRAamp$Variant_Type, breaks = 'FD')
#magiclistCN_PDGFRAgain <- magiclistCN_PDGFRA[Variant_Type > 1.6]
unique(magiclistCN_PDGFRAamp[Hugo_Symbol == 'PDGFRA', Tumor_Sample_Barcode])
magiclistCN_PDGFRAamp_s <- magiclistCN_PDGFRAamp[Tumor_Sample_Barcode %in% unique(magiclistCN_PDGFRAamp[Hugo_Symbol == 'PDGFRA', Tumor_Sample_Barcode])]
min(magiclistCN_PDGFRAamp_s$Start_Pos)
max(magiclistCN_PDGFRAamp_s$End_Pos)
hist(magiclistCN_PDGFRAamp_s$Variant_Type, breaks = 'FD')
#####


GZtad_gr_an_dt <- as.data.table(GZtad_gr_an)
PDGFRAtad <- GZtad_gr_an_dt[c(1712:1719),] #set mamually
sum(PDGFRAtad$width)
PDGFRAtadgr <- dt2gr(PDGFRAtad)
bins = gr.sub(tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[4], 
                         tilewidth = 10000,
                         cut.last.tile.in.chrom = T),'chr', '')


PDGFRAtadbins <- bins %&% PDGFRAtadgr
PDGFRAtadbins <- bins %&% '4:52500000-58000000'
PDGFRAtadbins <- PDGFRAtadbins %$% CGC_gr

GATKaCN_PDGFRAamps <- GATKaCN[Tumor_Sample_Barcode %in% unique(magiclistCN_PDGFRAamp[Hugo_Symbol == 'PDGFRA', Tumor_Sample_Barcode]) & chrom == 4]
GATKaCN_PDGFRAamps_gr <- dt2gr(GATKaCN_PDGFRAamps)
GATKaCN_PDGFRAamps_an <- GATKaCN_PDGFRAamps_gr %$% GZtad_gr_an
GATKaCN_PDGFRAamps_andt <- as.data.table(GATKaCN_PDGFRAamps_an)

###PDGFRA####
# loop over samples and assign CN to bin
x = unique(GATKaCN_PDGFRAamps$Tumor_Sample_Barcode)[1]
CNlist <- lapply(unique(GATKaCN_PDGFRAamps$Tumor_Sample_Barcode), function(x){
  print(x)
  xCN <- GATKaCN_PDGFRAamps_andt[Tumor_Sample_Barcode == x]
  xCN_gr <- dt2gr(xCN)
  xBins <- PDGFRAtadbins %$% xCN_gr
  xBins_dt <- as.data.table(xBins)
  xBins_dt[, Tumor_Sample_Barcode := x]
  return(xBins_dt)
})
CNlist_dt <- rbindlist(CNlist)
# call amped bin
CNlist_dt[,amp := as.numeric(cn) >= CNcutoff]

# dcast bin ~ sample, fill = amp
CNlist_dt[,binID := paste0(start, ':', end)]
AmpMat <- dcast(CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'amp', fun.aggregate = mean)
CNMat <- dcast(CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'cn', fun.aggregate = mean)

#AmpMat <- AmpMat[,-10]
#CNMat <- CNMat[,-10]
#rowsum amp/bin -> plot
rAmpPercent <- sapply(1:nrow(AmpMat), function(y){
  mean(as.numeric(AmpMat[y,-1]))
})
#mean amp
rAmpMean <- sapply(1:nrow(CNMat), function(y){
  mean(as.numeric(CNMat[y,-1]))
})

#rowsum amp/bin -> plot
rAmpPercent_narm <- sapply(1:nrow(AmpMat), function(y){
  mean(na.omit(as.numeric(AmpMat[y,-1])))
})
#mean amp
rAmpMean_narm <- sapply(1:nrow(CNMat), function(y){
  mean(na.omit(as.numeric(CNMat[y,-1])))
})

PDGFRAtadbins_dt <- cbind(as.data.table(PDGFRAtadbins), rAmpPercent, rAmpMean, rAmpPercent_narm, rAmpMean_narm)
PDGFRAtadbins_dt[,binID := paste0(start, ':', end)]
PDGFRAtadbins_dts <- PDGFRAtadbins_dt[!is.na(rAmpPercent_narm)]
PDGFRAtad_amps_gr <- dt2gr(PDGFRAtadbins_dts)
write.table(PDGFRAtadbins_dts, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220316extDataFigs_sourceTbl/20220318extfiguretable4A_.txt", sep ="\t", row.names = F, col.names = T, quote = F)


gt_pHGGq <- gTrack('/Users/fdubois/Dropbox/pHGG_data/data/20200430CCDC26SVpHGGmergedH3K27ac/20200430MACStracks_pHGGmeanH3K27ac_qpois.bw', name = 'pHGG', yaxis.pretty = TRUE)
PDGFRAtad_Mean_amp_gt <- gTrack(data = PDGFRAtad_amps_gr, y.field = 'rAmpMean_narm', name = 'MeanCN', col = 'blue')
PDGFRAtad_per_amp_gt <- gTrack(data = PDGFRAtad_amps_gr, y.field = 'rAmpPercent_narm', name = '%amp', col = 'blue')
plot(c(gt.ge, gt_pHGGq, PDGFRAtad_Mean_amp_gt, PDGFRAtad_per_amp_gt), '4:52500000-57000000', cex.label = 0.3)
# plot(c(gt.ge, gt_cnmc710ATAC,gt_pHGGq, PDGFRAtad_Mean_amp_gt, PDGFRAtad_per_amp_gt, gg.jabbaSJ15$gt, gg.jabba004R$gt, gg.jabbaD28T$gt, gg.jabbaSJ11$gt, gg.jabba032P$gt, 
#        gg.jabbaD57T$gt, gg.jabbaSJ75$gt, gg.jabbaG17$gt, gg.jabbaD61T$gt, gg.jabbaSJ22$gt,gg.jabbaD06T$gt), '4:52500000-57000000', cex.label = 0.3)
# plot(c(gt.ge,gt_pHGGq, PDGFRAtad_per_amp_gt, PDGFRAtad_Mean_amp_gt, gg.jabbaSJ15$gt, gg.jabba004R$gt, gg.jabbaD28T$gt, gg.jabbaSJ11$gt, gg.jabba032P$gt, 
#        gg.jabbaD57T$gt, gg.jabbaSJ75$gt, gg.jabbaG17$gt, gg.jabbaD61T$gt, gg.jabbaSJ22$gt,gg.jabbaD06T$gt), '4:52500000-57000000', cex.label = 0.3)
# plot(c(gt.ge,gt_pHGGq, PDGFRAtad_per_amp_gt, PDGFRAtad_Mean_amp_gt, gg.jabbaSJ15$gt, gg.jabba004R$gt, gg.jabbaD28T$gt, gg.jabbaSJ11$gt, gg.jabba032P$gt, 
#        gg.jabbaD57T$gt, gg.jabbaSJ75$gt, gg.jabbaG17$gt, gg.jabbaD61T$gt, gg.jabbaSJ22$gt,gg.jabbaD06T$gt), '4:52500000-59000000', cex.label = 0.3)

AmpMat_dt <- merge(AmpMat, PDGFRAtadbins_dt, by = 'binID') 
CNMat_dt <- merge(CNMat, PDGFRAtadbins_dt, by = 'binID') 

colnames(AmpMat_dt)
AmpMat_dts <- AmpMat_dt[!is.na(rAmpPercent)]
nrow(AmpMat_dts)
which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA')
nrow(AmpMat_dts) - which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA')[1]
nrow(AmpMat_dts) - tail(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA'),1)
# let's go for 1.5 Mbp up and downstream

AmpMat_plmin15e5 <- AmpMat_dt[(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-200):(tail(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+200),]

AmpMat_min15e5 <- AmpMat_dt[(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-200):(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-1),]
AmpMat_pl15e5 <- AmpMat_dt[(tail(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+1):(tail(which(AmpMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+200),]
which(colnames(AmpMat_pl15e5) == 'seqnames')
Int_left <- colSums(AmpMat_min15e5[,2:(which(colnames(AmpMat_pl15e5) == 'seqnames')-1)], na.rm = T)
summary(Int_left)
Int_right <- colSums(AmpMat_pl15e5[,2:(which(colnames(AmpMat_pl15e5) == 'seqnames')-1)], na.rm = T)
summary(Int_right)

Int_left_dt <- as.data.table(Int_left, keep.rownames = TRUE)[,group := 'Int_left']
colnames(Int_left_dt)[2] = 'Integral'
Int_right_dt <- as.data.table(Int_right, keep.rownames = TRUE)[,group := 'Int_right']
colnames(Int_right_dt)[2] = 'Integral'
Int_gg <- rbind(Int_left_dt, Int_right_dt)

library(ggplot2)
g <- ggplot(Int_gg, aes(x=group, y = Integral))
g + geom_boxplot( varwidth =TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6) + 
  theme_classic() + theme(text = element_text(size=24))

wilcox.test(Int_left, Int_right, paired = T)
t.test(Int_left, Int_right, paired = T)

##CNmat####
CNMat_dts <- CNMat_dt[!is.na(rAmpPercent)]
nrow(CNMat_dts)

CNMat_plmin15e5 <- CNMat_dts[(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-200):(tail(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+200),]

CNMat_min15e5 <- CNMat_dts[(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-200):(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA')[1]-1),]
CNMat_pl15e5 <- CNMat_dts[(tail(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+1):(tail(which(CNMat_dts$geneName == 'FIP1L1, PDGFRA'),1)+200),]
colnames(CNMat_pl15e5)
Int_left <- colSums(CNMat_min15e5[,2:(which(colnames(AmpMat_pl15e5) == 'seqnames')-1)], na.rm = T)
summary(Int_left)
Int_right <- colSums(CNMat_pl15e5[,2:(which(colnames(AmpMat_pl15e5) == 'seqnames')-1)], na.rm = T)
summary(Int_right)

Int_left_dt <- as.data.table(Int_left, keep.rownames = TRUE)[,group := 'Int_left']
colnames(Int_left_dt)[2] = 'Integral'
Int_right_dt <- as.data.table(Int_right, keep.rownames = TRUE)[,group := 'Int_right']
colnames(Int_right_dt)[2] = 'Integral'
Int_gg <- rbind(Int_left_dt, Int_right_dt)

library(ggplot2)
g <- ggplot(Int_gg, aes(x=group, y = Integral))
g + geom_boxplot( varwidth =TRUE) + geom_dotplot(binaxis='y', stackdir='center', binwidth = 50, dotsize = 0.6) + 
  theme_classic() + theme(text = element_text(size=24))

wilcox.test(Int_left, Int_right, paired = T)
t.test(Int_left, Int_right, paired = T)

Int_gg2 <- cbind(Int_left_dt, Int_right_dt)
colnames(Int_gg2)[c(2,5)] = c('Int_left', 'Int_right')
Int_gg2 <- Int_gg2[,-c(3,4,6)]
p <- ggplot(Int_gg2) + geom_segment(aes(x=1, xend=2, y=Int_left, yend=Int_right), size=.75, show.legend=F) + 
  geom_vline(xintercept=1, linetype="dashed", size=.1) + 
  geom_vline(xintercept=2, linetype="dashed", size=.1) +
  xlim(.5, 2.5) +  labs(x="", y="∫Mean CN") + #+ ylim(0,(1.1*(max(df$`1952`, df$`1957`))))
  theme_classic() + theme(text = element_text(size=24)) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(1,2,1,2), "cm"))
p <- p + geom_text(label="Left", x=1, y = max(Int_gg$Integral), hjust=1.2, size=5)  # title
p <- p + geom_text(label="Right", x=2, y= max(Int_gg$Integral), hjust=-0.1, size=5)  # title
p

wilcox.test(Int_left, Int_right, paired = T)
t.test(Int_left, Int_right, paired = T)



