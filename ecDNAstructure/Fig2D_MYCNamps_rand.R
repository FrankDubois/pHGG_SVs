require(data.table)
require(gUtils)
require(BSgenome.Hsapiens.UCSC.hg19)
require(rtracklayer)

#####
#subset magicList to samples with amps in MYCN TAD
magiclist <- readRDS('/Users/fdubois/Dropbox/pHGG_data/data/20200228magiclistQC/20200228magiclistall.RDS')
colnames(magiclist)
#magiclist[,Start_Chr := as.numeric(Start_Chr)]
#class(magiclist$End_Chr)
magiclistCN_MYCN <- magiclist[Variant == 'CN' & (Start_Chr == 2 | End_Chr == 2)]
#unique(magiclistCN_MYCN$End_Chr)
magiclistCN_MYCN[,Variant_Type := as.numeric(Variant_Type)]
unique(magiclistCN_MYCN$Tumor_Sample_Barcode)
magiclistCN_MYCN[Variant_Type >= 2.1, Amp := 1]
magiclistCN_MYCNamp <- magiclistCN_MYCN[Amp ==1]
unique(magiclistCN_MYCNamp$Tumor_Sample_Barcode)
sort(table(magiclistCN_MYCNamp$Hugo_Symbol), decreasing = TRUE)
min(magiclistCN_MYCNamp$Start_Pos)
max(magiclistCN_MYCNamp$End_Pos)
hist(magiclistCN_MYCNamp$Variant_Type, breaks = 'FD')
#magiclistCN_MYCNgain <- magiclistCN_MYCN[Variant_Type > 1.6]
unique(magiclistCN_MYCNamp[Hugo_Symbol == 'MYCN', Tumor_Sample_Barcode])
magiclistCN_MYCNamp_s <- magiclistCN_MYCNamp[Tumor_Sample_Barcode %in% unique(magiclistCN_MYCNamp[Hugo_Symbol == 'MYCN', Tumor_Sample_Barcode])]
min(magiclistCN_MYCNamp_s$Start_Pos)
max(magiclistCN_MYCNamp_s$End_Pos)
hist(magiclistCN_MYCNamp_s$Variant_Type, breaks = 'FD')
#####

# tile MYCN TAD 10kb
genes = gr.sub(import('http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf')) # rtracklayer::import reads gtf and gr.sub replaces chr
GZtad_gr_an <- readRDS('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20200127fhook/GZtad_gr_an.rds')


GZtad_gr_an_dt <- as.data.table(GZtad_gr_an)
MYCN_ID2tad <- GZtad_gr_an_dt[c(1322:1334),]
MYCNtad <- GZtad_gr_an_dt[c(1327:1329),]
sum(MYCNtad$width)
MYCNtadgr <- dt2gr(MYCNtad)
bins = gr.sub(tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[2], 
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T),'chr', '')

MYCNtadbins <- bins %&% MYCNtadgr

# get number of samples with amp in each tile
GATKaCN <- fread('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20190914magicList/20190914GATKaCNseg.txt')
#match names
samplesMatch <- as.data.table(cbind(sort(unique(GATKaCN$track.name)), 
                                    sort(unique(magiclistCN_MYCN[!(Tumor_Sample_Barcode %in% c('GBM36_T', 'GBM86_T', 'GBM96_T', 'GBM97_T', 'GBM98_T')),Tumor_Sample_Barcode])))) # unique sample list ## of 207!! samples?? -> SJ doubled
colnames(samplesMatch) <- c('CN_pair', 'magiclist_Tcode')
GATKaCN[ ,Tumor_Sample_Barcode := samplesMatch[CN_pair == track.name, magiclist_Tcode], by = track.name]

GATKaCN_MYCNamps <- GATKaCN[Tumor_Sample_Barcode %in% unique(magiclistCN_MYCNamp[Hugo_Symbol == 'MYCN', Tumor_Sample_Barcode]) & chrom == 2]
GATKaCN_MYCNamps_gr <- dt2gr(GATKaCN_MYCNamps)
GATKaCN_MYCNamps_an <- GATKaCN_MYCNamps_gr %$% GZtad_gr_an
GATKaCN_MYCNamps_andt <- as.data.table(GATKaCN_MYCNamps_an)

###MYCN####
# loop over samples and assign CN to bin
x = unique(GATKaCN_MYCNamps$Tumor_Sample_Barcode)[1]
CNlist <- lapply(unique(GATKaCN_MYCNamps$Tumor_Sample_Barcode), function(x){
  xCN <- GATKaCN_MYCNamps_andt[Tumor_Sample_Barcode == x]
  xCN_gr <- dt2gr(xCN)
  xBins <- MYCNtadbins %$% xCN_gr
  xBins_dt <- as.data.table(xBins)
  xBins_dt[, Tumor_Sample_Barcode := x]
  return(xBins_dt)
})
CNlist_dt <- rbindlist(CNlist)
# call amped bin
CNlist_dt[,amp := as.numeric(cn) >= 2.1]

# dcast bin ~ sample, fill = amp
CNlist_dt[,binID := paste0(start, ':', end)]
AmpMat <- dcast(CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'amp', fun.aggregate = mean)
CNMat <- dcast(CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'cn', fun.aggregate = mean)

#rowsum amp/bin -> plot
rAmpPercent <- sapply(1:nrow(AmpMat), function(y){
  mean(na.omit(as.numeric(AmpMat[y,2:9])))
})

# MYCNtadbins_dt <- cbind(as.data.table(MYCNtadbins), rMean)
MYCNtad_amps_gr <- dt2gr(MYCNtadbins_dt)
require(gTrack)
gt.ge <- track.gencode()
gt_cnmc710ATAC <- gTrack('/Users/fdubois/Dropbox/pHGG_data/data/20191211ATACseq/GSM3596238_cnmc710t-atac1.rpm.bw', name ='710atac')
gt_434 <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191210MACStrsckbw/201911210_434_H3K27ac_FE.bdg.bw', name = '434')
gt_1840 <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191210MACStrsckbw/201911210_1840_H3K27ac_FE.bdg.bw', name = '1840')


MYCNtad_amps_gt <- gTrack(data = MYCNtad_amps_gr, y.field = 'rMean', name = '% amp', col = 'blue')

# plot(c(gt.ge, MYCNtad_amps_gt), '2:13160001-17880000')

#mean amp
rAmpMean <- sapply(1:nrow(CNMat), function(y){
  mean(na.omit(as.numeric(CNMat[y,2:9])))
})

MYCNtadbins_dt <- cbind(as.data.table(MYCNtadbins), rAmpPercent, rAmpMean)
MYCNtad_amps_gr <- dt2gr(MYCNtadbins_dt)
MYCNtad_Men_amp_gt <- gTrack(data = MYCNtad_amps_gr, y.field = 'rAmpMean', name = 'MeanCN', col = 'blue')

# plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNtad_Men_amp_gt), '2:13160001-17880000', cex.label = 0.2)

#Split into MYCN+ neighborhood / MYCN + ID2
#
#randomize amplicons


##### MYCN_ID2####

MYCN_ID2tadgr <- dt2gr(MYCN_ID2tad)
MYCN_ID2tadbins <- bins %&% MYCN_ID2tadgr

# loop over samples and assign CN to bin
x = unique(GATKaCN_MYCNamps$Tumor_Sample_Barcode)[1]
MYCN_ID2CNlist <- lapply(unique(GATKaCN_MYCNamps$Tumor_Sample_Barcode), function(x){
  print(x)
  xCN <- GATKaCN_MYCNamps_andt[Tumor_Sample_Barcode == x]
  xCN_gr <- dt2gr(xCN)
  xBins <- MYCN_ID2tadbins %$% xCN_gr
  xBins_dt <- as.data.table(xBins)
  xBins_dt[, Tumor_Sample_Barcode := x]
  return(xBins_dt)
})
MYCN_ID2CNlist_dt <- rbindlist(MYCN_ID2CNlist)
# call amped bin
MYCN_ID2CNlist_dt[,amp := as.numeric(cn) >= 2.1]
MYCN_ID2CNlist_dt[,binID := paste0(start, ':', end)]

MYCN_ID2CNlist_gr <- dt2gr(MYCN_ID2CNlist_dt)
MYCN_ID2CNlist_gr_an <- MYCN_ID2CNlist_gr %$% genes
MYCN_ID2CNlist_dt_an <- as.data.table(MYCN_ID2CNlist_gr_an)
  
# dcast bin ~ sample, fill = amp
MYCN_ID2AmpMat <- dcast(MYCN_ID2CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'amp', fun.aggregate = mean)
MYCN_ID2CNMat <- dcast(MYCN_ID2CNlist_dt, binID ~ Tumor_Sample_Barcode, value.var = 'cn', fun.aggregate = mean)
colnames(MYCN_ID2AmpMat)
#ID2MYCN : 
#rowsum amp/bin -> plot
MYCN_ID2rAmpPercent <- rbindlist(lapply(1:nrow(MYCN_ID2AmpMat), function(y){
  as.data.table(cbind(MYCN_ID2AmpMat[y,1],mean(na.omit(as.numeric(MYCN_ID2AmpMat[y,2:4])))))
}))
#mean amp
MYCN_ID2rAmpMean <- rbindlist(lapply(1:nrow(MYCN_ID2CNMat), function(y){
  as.data.table(cbind(MYCN_ID2CNMat[y,1],mean(na.omit(as.numeric(MYCN_ID2CNMat[y,2:4])))))
}))

MYCNID2bins_dt <- as.data.table(MYCN_ID2tadbins)
MYCNID2bins_dt[,BinID := paste0(start, ':', end)]
MYCNID2bins_dt_per <- merge(x = MYCNID2bins_dt, y = MYCN_ID2rAmpPercent, by.x = 'BinID', by.y = 'binID')
MYCNID2bins_dt_ <- merge(x = MYCNID2bins_dt_per, y = MYCN_ID2rAmpMean, by.x = 'BinID', by.y = 'binID')
MYCNID2_amps_gr <- dt2gr(MYCNID2bins_dt_)
# write.table(MYCNID2bins_dt_, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220308figuretables/20220308figuretable2G_MYCN_ID2.txt", sep ="\t", row.names = F, col.names = T, quote = F)

#MYCNID2_amps_gran <- MYCNID2_amps_gr %$% genes
#MYCNID2_amps_gran_dt <- as.data.table(MYCNID2_amps_gran)

MYCNID2_per_amp_gt <- gTrack(data = MYCNID2_amps_gr, y.field = 'V2.x', name = '%amp MYCN_ID2', col = 'blue')
MYCNID2_Men_amp_gt <- gTrack(data = MYCNID2_amps_gr, y.field = 'V2.y', name = 'MeanCN MYCN_ID2', col = 'blue')

# plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_per_amp_gt, MYCNID2_Men_amp_gt), '2:8000001-17880000', cex.label = 0.2)
# plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_per_amp_gt, MYCNID2_Men_amp_gt), '2:13160001-17880000', cex.label = 0.2)


#MYCN_neighborhood : 
#rowsum amp/bin -> plot
MYCN_NrAmpPercent <- rbindlist(lapply(1:nrow(MYCN_ID2AmpMat), function(y){
  as.data.table(cbind(MYCN_ID2AmpMat[y,1],mean(na.omit(as.numeric(MYCN_ID2AmpMat[y,5:9])))))
}))
#mean amp
MYCN_NrAmpMean <- rbindlist(lapply(1:nrow(MYCN_ID2CNMat), function(y){
  as.data.table(cbind(MYCN_ID2CNMat[y,1],mean(na.omit(as.numeric(MYCN_ID2CNMat[y,5:9])))))
}))

MYCN_Nbins_dt_per <- merge(x = MYCNID2bins_dt, y = MYCN_NrAmpPercent, by.x = 'BinID', by.y = 'binID')
MYCN_Nbins_dt_ <- merge(x = MYCN_Nbins_dt_per, y = MYCN_NrAmpMean, by.x = 'BinID', by.y = 'binID')
# write.table(MYCN_Nbins_dt_, "/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220308figuretables/20220308figuretable2G_MYCN_only.txt", sep ="\t", row.names = F, col.names = T, quote = F)

MYCN_N_amps_gr <- dt2gr(MYCN_Nbins_dt_)

MYCN_N_per_amp_gt <- gTrack(data = MYCN_N_amps_gr, y.field = 'V2.x', name = '%amp MYCN_N', col = 'blue')
MYCN_N_Mean_amp_gt <- gTrack(data = MYCN_N_amps_gr, y.field = 'V2.y', name = 'MeanCN MYCN_N', col = 'blue')

plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_per_amp_gt, MYCN_N_per_amp_gt), '2:8000001-17880000', cex.label = 0.2)
plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_per_amp_gt, MYCN_N_per_amp_gt), '2:13160001-17880000', cex.label = 0.2)

saveRDS(MYCNID2_Men_amp_gt, '/Users/fdubois/Dropbox/pHGG_data/data/20200615MYCNmeanamps/20200615MYCNID2_Mean_amp_gt.rds')
saveRDS(MYCN_N_Mean_amp_gt, '/Users/fdubois/Dropbox/pHGG_data/data/20200615MYCNmeanamps/20200615MYCN_N_Mean_amp_gt.rds')
saveRDS(MYCNID2_per_amp_gt, '/Users/fdubois/Dropbox/pHGG_data/data/20200615MYCNmeanamps/20200615MYCNID2_per_amp_gt.rds')
saveRDS(MYCN_N_per_amp_gt, '/Users/fdubois/Dropbox/pHGG_data/data/20200615MYCNmeanamps/20200615MYCN_N_per_amp_gt.rds')

gt_429 <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191210MACStrsckbw/201911210_429_H3K27ac_FE.bdg.bw', name = '429', y1 = 30)
plot(c(gt.ge, gt_429, MYCNID2_Men_amp_gt, MYCN_N_Mean_amp_gt), '2:8000001-17880000', cex.label = 0.2)
#pHGGmerged_q <- gTrack('/Users/fdubois/Dropbox/pHGG_data/data/20200430CCDC26SVpHGGmergedH3K27ac/20200430MACStracks_pHGGmeanH3K27ac_qpois.bw', name = 'pHGGmerg_q',  y1 = 400)
#plot(c(gt.ge, pHGGmerged_q, MYCNID2_Men_amp_gt, MYCN_N_Mean_amp_gt), '2:8000001-17880000', cex.label = 0.2)

plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_Men_amp_gt, MYCN_N_Mean_amp_gt), '2:8000001-17880000', cex.label = 0.2)
plot(c(gt.ge, gt_cnmc710ATAC,gt_434,gt_1840, MYCNID2_Men_amp_gt, MYCN_N_Mean_amp_gt), '2:13160001-17880000', cex.label = 0.2)


####H3K27acInputs####
gt_429i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_429inputNDe500.bw', name = '429')
gt_434i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_434inputNDe500.bw', name = '434')
gt_1840i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_1840inputNDe500.bw', name = '1840')
gt_2230i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_2230inputNDe500.bw', name = '2230')
gt_HSJ019i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_HSJ019inputNDe500.bw', name = 'HSJ019')
gt_JN27i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_JN27inputNDe500.bw', name = 'JN27')
gt_JN29i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_JN29inputNDe500.bw', name = 'JN29')
gt_JN57i <- gTrack('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20191209CHIPseqinputs/20191209_JN57inputNDe500.bw', name = 'JN57')

plot(c(gt.ge,gt_429i ,gt_434i,gt_1840i,gt_2230i, gt_HSJ019i, gt_JN27i,gt_JN29i, gt_JN57i,
       MYCNID2_Men_amp_gt, MYCN_N_Mean_amp_gt), '2:8000001-17880000', cex.label = 0.2)

