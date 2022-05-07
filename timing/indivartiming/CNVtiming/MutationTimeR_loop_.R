library("MutationTimeR")
require(data.table)
library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(parallel)
fls <- list.files('/xchip/beroukhimlab/Frank/DIPG/data/20201106vcf_reformating/')
fls2 <- gsub(pattern = '_mutTRform_filt.vcf', replacement = '', x = fls)
fpths <- paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201106vcf_reformating/', fls)
x = 49
#absSegs_sub <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200828palimsets/20200828absSegs_sub.rds')

absSegs_sub <- readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20200929mutationtimeR/20201025_absSegs_sub.rds')
metaD <- readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20200929mutationtimeR/20200929_metaD.rds')

genome_hg19msked_length <- seqlengths(BSgenome.Hsapiens.UCSC.hg19.masked)
names(genome_hg19msked_length) <- gsub(pattern = 'chr', replacement = '', x = names(genome_hg19msked_length))

clFLSf <- readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20201023mutationTimeR/20201023_clFLS_K2N3pi02.rds')
samplesMerge3 <- readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20201023mutationTimeR/20201023_samplesMerge3.rds')


mclapply(1:nrow(samplesMerge3), function(x){
  print(samplesMerge3[x,])
  print(samplesMerge3[x,1])
  vcf_X <- readVcf(fpths[x])
  bb_X_dt <- absSegs_sub[Sample == samplesMerge3[x, absSegTScode]]
  Wgd = mean(bb_X_dt$Ploidy) > 3
  bbX_gr <- dt2gr(dt = bb_X_dt, seqlengths=genome_hg19msked_length)
  clusterX <- clFLSf[TScode == samplesMerge3[x, mobstSamples]]
  if (nrow(clusterX) == 0) {mt_X <- mutationTime(vcf = vcf_X, cn = bbX_gr, n.boot=100, isWgd = wgd)
  } else {
    mt_X <- mutationTime(vcf = vcf_X, cn = bbX_gr, clusters = clusterX, n.boot=100, isWgd = wgd)
  }
  print(table(mt_X$V$CLS))
  info(header(vcf_X)) <- rbind(info(header(vcf_X)),mtHeader())
  info(vcf_X) <- cbind(info(vcf_X), mt_X$V)
  mcols(bbX_gr) <- cbind(mcols(bbX_gr),mt_X$T)
  plotSample(vcf_X,bbX_gr)
  bb_Xgr_dt_postMtimeR <- as.data.table(bbX_gr)
  write.table(x = bb_Xgr_dt_postMtimeR, file = 
                paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                       '_mutTimeR_CN.seg')
              ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)
  saveRDS(bb_Xgr_dt_postMtimeR, paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                                       '_mutTimeR_CN_dt.rds'))
  saveRDS(bbX_gr, paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                                       '_mutTimeR_CN_gr.rds'))
  writeVcf(vcf = vcf_X, paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                               '_mutTimeR_SNV.vcf'))
  saveRDS(vcf_X, paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                         '_mutTimeR_SNV_vcf.rds'))
  pdf(file = paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', as.character(samplesMerge3[x,1]),
                    '_samplePlot.pdf'), paper = 'USr')
  plotSample(vcf_X,bbX_gr)
  dev.off()
}, mc.preschedule = FALSE, mc.cores = 8)
  
  
  
  
  
  
  