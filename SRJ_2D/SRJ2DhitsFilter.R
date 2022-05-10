### this script filters the 2D hits
require(gUtils)
require(gTrack)
SRJs_A <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200227SRJ2D/20200227new2Dhitstbl_ql5_lt100.txt')
SRJs <- unique(SRJs_A$gene_i)

setkey(SRJs_A, gene_i)

x = SRJs[1]
SRJs_m <- lapply(SRJs, function(x){
  print(x)
  SRJs_A_s <- SRJs_A[x]
  min_p1 <- min(SRJs_A_s$pos_i)
  max_p1 <- max(SRJs_A_s$pos_i)
  min_p2 <- min(SRJs_A_s$pos_j)
  max_p2 <- max(SRJs_A_s$pos_j)
  if ((min_p1 == max_p1) | (min_p2 == max_p2)) { 
    print(paste('artefact:', x))
    return()}
  binspan_p1 <- max_p1-min_p1
  binspan_p2 <- max_p2-min_p2
  nsamples <- length(unique(SRJs_A_s$donor_unique_id))
  if (nsamples == 1) { 
    print(paste('only1sample:', x))
    return()}
  SRJbins <- as.data.table(cbind(x, SRJs_A_s$gene_j[1], SRJs_A_s$chr_i[1], min_p1,max_p1,binspan_p1, SRJs_A_s$chr_j[1],min_p2, max_p2, binspan_p2, 
                                 SRJs_A_s$pval[1], nsamples))
  colnames(SRJbins) <- c('SRJ_i','SRJ_j', 'chr1', 'start1','end1', 'span1', 'chr2', 'start2', 'end2', 'span2', 'pval', 'nsamples')
  return(SRJbins)
})
SRJs_mdt <- rbindlist(SRJs_m)
SRJs_mdt$pval <- as.numeric(SRJs_mdt$pval)
SRJs_mdt$span1 <- as.numeric(SRJs_mdt$span1)
SRJs_mdt$span2 <- as.numeric(SRJs_mdt$span2)
SRJs_mdt$nsamples <- as.numeric(SRJs_mdt$nsamples)
SRJs_mdt_s <- SRJs_mdt[nsamples >=2]
SRJs_mdt_s <- SRJs_mdt_s[(span1 < 1e6) & (span2 < 1e6)]

SRJs_A_ss <- SRJs_A[gene_i %in% SRJs_mdt_s$SRJ_i]

# remove hits with mean quality less than 20
SRJs_A_ss[, mean.quality := mean(QUALITY), by = hit_num]
SRJs_A_s2 <- SRJs_A_ss[mean.quality>25]

# remove hits with no ASDIS or ASSMB
SRJs_A_s2[, is.evdnce := any(EVDNC %in% c("ASSMB", "ASDIS")), by = hit_num]
SRJs_A_s3 <- SRJs_A_s2[is.evdnce == TRUE]

write.table(SRJs_A_s3, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200227SRJ2D/20200227_2Dhits_nle2_lt100.txt'
            ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)
