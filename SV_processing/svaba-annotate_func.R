#!/usr/bin/env Rscript
#This function reads in a svaba vcf and returns an annotated data.table 
make_nice_svaba_andcplot <- function(file) {
############### FUNCTIONS
load_svaba <- function(files, mc.cores=mc.cores, unlist=FALSE, bad.remove=TRUE) {

  snow4 <- lapply(files, function(x) {
    if (!file.exists(x)) {
      print(paste("File does not exist",x))
      if (unlist)
        return (GRanges())
      else
        return (GRangesList())
    }
    
    #print(paste(basename(x), "is", match(x, files), "of", length(files)))
    #aa <- skitools::ra_breaks(x)
    #aa2 <- gr.flipstrand(grl.unlist(aa))
    #dd <- mcols(aa)
    #grlix <- aa2$grl.ix
    #mcols(aa2) <- NULL
    #aa <- split(aa2, grlix)
    #mcols(aa) <- dd

    ## dont continue if empty
    if (!as.numeric(system(paste("grep -v '^#'", x, "| wc -l"), intern=TRUE))) {
      if (unlist)
        return (GRanges())
      else
        return (GRangesList())
    }
      
    ff <- data.table::fread(paste("grep -v '^#'", x),sep='\t')
    if (nrow(ff) == 0)
      return (ff)
    if (ncol(ff)==10) {
      setnames(ff, paste0("V",seq(1:10)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL"))
    } else {
      setnames(ff, paste0("V",seq(1:11)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"))
    }
    ff[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    ff$sample = gsub("(.*?)_.*","\\1",basename(x))
    ff[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    ff[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    ff[, NUMPARTS := as.integer(gsub(".*?NUMPARTS=([0-9]+).*", "\\1", INFO))]
    ff[, SCTG := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
    ff[, DISC_MAPQ := as.numeric(gsub(".*?DISC_MAPQ=([0-9]+).*", "\\1", INFO))]
    ff[, NM := as.integer(gsub(".*?;NM=([0-9]+).*", "\\1", INFO))]
    ff[, MATENM := as.integer(gsub(".*?;MATENM=([0-9]+).*", "\\1", INFO))]
    ff[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    print("...still formatting")
    ff[, REPSEQ := gsub(".*?;REPSEQ=([A-Z]+).*", "\\1", INFO)]
    ff[, REPSEQ := ifelse(grepl(";", REPSEQ), "", REPSEQ)] 
    ff[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
    ff[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    ff[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] ##{ xx <- gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO); if (!grepl(";",xx)) { xx } else { "" }}]
    ff[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    if ("TUMOR" %in% colnames(ff)) {
      ff[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
      ff[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
      ff[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
    }
    if ("NORMAL" %in% colnames(ff)) {
      ff[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
      ff[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
      ff[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
    }
    ff[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
    ff[, inv := strand[1] == strand[2], by=uid]
    ff[, altstrand := rev(strand), by=uid]
    ff[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
    ff[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
    ff[, end := start]

    ff[, c("ID","REF","ALT","QUAL", "INFO", "GENO") := NULL]
    
    bad.ix <- ff[grepl("^G|^M", seqnames), uid]
    ff <- ff[!uid %in% bad.ix]
    
    ## if (length(aa) == 0) {
    ##   if (unlist)
    ##     return (GRanges())
    ##   else
    ##     return (GRangesList())
    ## }
    ## bad.ix <- unique(grl.unlist(aa)[grepl("^G", as.character(seqnames(grl.unlist(aa))))]$grl.ix)
    ## #tab <- table(grl.unlist(aa)$grl.ix[grepl("ref", as.character(seqnames(grl.unlist(aa))))])
    ## #bad.ix <- unique(c(bad.ix, as.numeric(names(tab[tab == 2]))))
    ## if (length(bad.ix)) 
    ##   aa <- aa[-bad.ix]
    ## if (length(aa))
    ##   mcols(aa)$filename = basename(x);
    ## if (unlist)
    ##   return (grl.unlist(aa))
    ## else
    ##   return (aa)

    df <- as.data.frame(dt2gr(ff))[,setdiff(colnames(ff), c("seqnames","start","end","strand"))]
    df <- df[!duplicated(df$uid),]
    g <- dt2gr(ff[,.(seqnames, start, end, strand, uid)])
    g <- split(g, g$uid)
    stopifnot(all(sapply(g, length) == 2))
    mcols(g) <- df[match(names(g), df$uid),]
    return(g)
  })
  return (snow4)
}

#################

library(optparse)

option_list = list(
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input VCF file"),
  make_option(c("-o", "--output"), type = "character", default = "",  help = "Output annotation name"),
  make_option(c("-t", "--splitsupport"), type = "numeric", default = 4,  help = "Minimum number of tumor supporting reads for ASSMB"),
  make_option(c("-d", "--discsupport"), type = "numeric", default = 10,  help = "Minimum number of tumor support discordant for DSCRD"),
  make_option(c("-g", "--genes"), type = "numeric", default = 1,  help = "Add genes to the plot?"),
  make_option(c("-H", "--height"), type = "numeric", default = 10,  help = "Height"),
  make_option(c("-W", "--width"), type = "numeric", default = 10,  help = "Width"),
  make_option(c("-s", "--minsize"), type = "numeric", default = 0,  help = "Minimum size to plot for SV"),
  make_option(c("-S", "--svaba"), type = "numeric", default = 1,  help = "This is a svaba run. 0 for generic SV VCF") 
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

 if (is.null(file))
  stop(print_help(parseobj))

if (!nchar(opt$output)) {
  opt$output <- tools::file_path_sans_ext(opt$input)
}

print("...loading required packages")
library(gUtils) 
library(gplots)
library(plyr)
library(RCircos)

print("...importing tracks")
basedir <- "/xchip/beroukhimlab/Jeremiah/tracks"
genes <- readRDS("/Volumes/xchip_beroukhimlab/Jeremiah/tracks/gr.allgenes.rds")
genes <- genes[width(genes) < 2e6] 
#cgc <- readRDS(file.path(basedir, 'gr.allgenes.rds'))
#cgc = read.delim('/home/unix/marcin/DB/COSMIC/cancer_gene_census.tsv', strings = FALSE)
#cgc <- genes[genes$gene %in% cgc$Symbol]

## read in file
print("...reading VCF file")
opt$input <- file
if (opt$svaba) {
  suppressWarnings(bks <- load_svaba(opt$input)[[1]])
  mcols(bks)$tdisc  <- as.integer(gsub(".*?/.*?:.*?:.*?:.*?:.*?:.*?:(.*?):.*", "\\1", mcols(bks)$TUMOR))
  mcols(bks)$tsplit <- as.integer(gsub(".*?/.*?:.*?:.*?:.*?:.*?:(.*?):.*", "\\1", mcols(bks)$TUMOR))

  ## filter out discorant only with low support
  if (length(bks)) {
    print(paste("...filtering DSCRD breakpoints with < min support", opt$discsupport))
    bks <- bks[mcols(bks)$EVDNC != "DSCRD" | (mcols(bks)$EVDNC == "DSCRD" & mcols(bks)$tdisc >= opt$discsupport)]
  }
  
  ## filter out based on size
  if (length(bks)) {
    print(paste("...filtering breakpoints with < min span of", opt$minsize))
    bks <- bks[mcols(bks)$SPAN >= opt$minsize | mcols(bks)$SPAN < 0]
  }
  
  ## filter out assembly-only with low split support
  if (length(bks) && sum(bks$evidence == "ASSMB" & bks$tsplit < opt$splitsupport)) {
    print(paste("...filtering ASSMB with < min split support of", opt$minsize))
    bks <- bks[-(bks$evidence == "ASSMB" & bks$tsplit < opt$splitsupport)]
  }
  
} else {
  dt <- fread(paste("grep -v ^#", opt$input), sep="\t")
  setnames(dt, c("V1","V2","V3","V4","V5"), c("seqnames","start", "id","REF","ALT"))
  dt[, end := start]
  dt[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  dt[, uid := gsub("(,*?)_.*","\\1",id)]
  dt[, altstrand := rev(strand), by=uid]
  dt[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  dt[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9A-Z]+).*", "\\2", ALT)]

  dt[, V8 := NULL]
  gg <- dt2gr(dt)
  bks <- split(gg, gg$uid)
  mcols(bks)$INSERTION = ""
  mcols(bks)$HOMSEQ = ""
  mcols(bks)$sample = gsub("(.*?)\\..*", "\\1", basename(opt$input))
}

## bail early if nothing found
if (length(bks) == 0) {
  print("No rearrangements found")
  return()
} else {
  print(paste("Found", length(bks), "rearrangements"))
}

## annotate with gene overlaps
suppressWarnings(gru <- gr.val(grl.unlist(bks), genes, "gene"))
mcols(bks)$gene1 <- mcols(bks)$gene2 <- "Intergenic"
mcols(bks)$gene1 <- gru$gene[gru$grl.iix == 1]
mcols(bks)$gene2 <- gru$gene[gru$grl.iix == 2]                     

## annotate with genes within 100kb
suppressWarnings(gru <- gr.val(grl.unlist(bks), genes + 100e3, "gene"))
mcols(bks)$gene1_100kb <- gru$gene[gru$grl.iix == 1]
mcols(bks)$gene2_100kb <- gru$gene[gru$grl.iix == 2]

if (!"gene1" %in% colnames(mcols(bks)))
  mcols(bks)$gene1 <- ""
if (!"gene2" %in% colnames(mcols(bks)))
  mcols(bks)$gene2 <- ""
if (!"gene1_100kb" %in% colnames(mcols(bks)))
  mcols(bks)$gene1_100kb <- ""
if (!"gene2_100kb" %in% colnames(mcols(bks)))
  mcols(bks)$gene2_100kb <- ""

gr <- grl.unlist(bks)
llr <- gr2dt(gr)
llr$bk_msg <- llr$gstrand <- llr$gene <- llr$msg <- ""
llr$elem_num <- -1

gr.introns <- readRDS("/Volumes/xchip_beroukhimlab/Jeremiah/tracks/gr.introns.rds")
gr.exons <- readRDS("/Volumes/xchip_beroukhimlab/Jeremiah/tracks/gr.exons.rds")
gene.comp <- readRDS("/Volumes/xchip_beroukhimlab/Jeremiah/tracks/gr.intergenic.rds") 
gr.tub <- readRDS("/Volumes/xchip_beroukhimlab/Jeremiah/tracks/master_db_29062016_ranges.rds")

##
print("...annotating L1")
gr <- gr.val(gr, gr.tub, "type", ignore.strand=TRUE)
if ('type' %in% colnames(mcols(gr))) {
  llr$suspect_L1_retrotransposition <- nchar(gr$type) > 0
} else { 
  llr$suspect_L1_retrotransposition <- FALSE
}
llr[,frame := NA]
## overlap introns
print("...annotating introns")
fo <- gr2dt(gr.findoverlaps(gr, gr.introns))
if (nrow(fo)) {
  fo[, intron_num   := gr.introns$num[subject.id], by=subject.id]
  fo[, intron_str   := as.character(strand(gr.introns))[subject.id], by=subject.id]
  fo[, intron_frm   := gr.introns$frame[subject.id], by=subject.id]
  fo[, intron_gene  := gr.introns$gene[subject.id], by=subject.id]
  fo[, intron_start := start(gr.introns)[subject.id], by=subject.id]
  fo[, dist_after   := start - intron_start, by=subject.id]
  fo[, bk_msg := paste0(dist_after, " bp into intron ", intron_num, " of ", intron_gene, "(", intron_str, ")")]
  fo[, frame  := gr.introns$frame[subject.id], by=subject.id]
  llr$bk_msg[fo$query.id] <- fo$bk_msg
  llr$frame[fo$query.id] <- fo$frame
  llr$gene[fo$query.id] <- fo$intron_gene
  llr$gstrand[fo$query.id] <- fo$intron_str
  llr$elem_num[fo$query.id] <- fo$intron_num
}

## overlap exons
print("...annotating exons")
fo <- gr2dt(gr.findoverlaps(gr, gr.exons))
if (nrow(fo)) {
  fo[, exon_num := gr.exons$num[subject.id], by=subject.id]
  fo[, exon_str := as.character(strand(gr.exons))[subject.id], by=subject.id]
  fo[, exon_frm := gr.exons$frame[subject.id], by=subject.id]
  fo[, exon_gene := gr.exons$gene[subject.id], by=subject.id]
  fo[, exon_start := ifelse(strand(gr.exons)[subject.id] == '+', start(gr.exons)[subject.id], end(gr.exons)[subject.id]), by=subject.id]
  fo[, dist_after := abs(start - exon_start) + 1, by=subject.id]
  fo[, bk_msg := paste0(dist_after, " bp into exon ", exon_num, " of ", exon_gene, "(", exon_str, ")")]
  fo[, frame  := gr.introns$frame[subject.id], by=subject.id]
  llr$bk_msg[fo$query.id] <- fo$bk_msg
  llr$frame[fo$query.id] <- fo$frame
  llr$gene[fo$query.id] <- fo$exon_gene
  llr$gstrand[fo$query.id] <- fo$exon_str
  llr$elem_num[fo$query.id] <- fo$exon_num
}

## annotate fusions
print("...annotating fusions")
llr[, fusion := { f <- sum(grepl("intron|exon", bk_msg)) == 2 && (gene[1] != gene[2]); if (is.na(f)) { FALSE } else { f } } , by = uid]
llr[, sense  := all(gstrand != "") && ((strand[1] == strand[2] && gstrand[1] != gstrand[2]) || (strand[1] != strand[2] && gstrand[1] == gstrand[2])), by = uid]

## DECIDE WHICH PIECE IS "first" in the fusion gene (where Tx starts)
## RAR_STRAND        GENE_STRAND
## ++ or --          must be +- or -+. Starts on +
## +- or --          msut be -- or ++. Starts on side where gstrand == rarstrand
llr[, gorder := {

  mk <- c(NA,NA)
  ## inversion
  if (strand[1] == strand[2]) {
    if (gstrand[1] == '+')
      mk <- c(2,1)
    else
      mk <- c(1,2)
  }

  ## non type
  if (gstrand[1] == strand[1])
    mk <- c(1,2)
  else
    mk <- c(2,1)

  mk

  }, by=uid]

## decide if in frame
if ("frame" %in% colnames(llr)) {
  llr[, in_frame := {
    f = (frame[gorder[2]] - frame[gorder[1]]);
    xbases <- nchar(as.character(INSERTION[1])) - nchar(as.character(HOMSEQ[1]));
    if (is.na(f))
      FALSE
    else if (all(grepl('intron', bk_msg))) ## intron to intron
      f == 0
    else
      f + xbases %% 3 == 0
  }, by=uid]
  
  ## make the fusion messages
  llr[, in_frame_fusion := in_frame && fusion && sense, by=uid]
  llr[llr$in_frame_fusion, msg := paste("In-frame sense fusion from:", bk_msg[gorder[1]], "to", bk_msg[gorder[2]]), by=uid]
  llr[!llr$in_frame_fusion & llr$fusion & llr$sense, msg := paste("Out-of-frame sense fusion between:", bk_msg[1], "and", bk_msg[2]), by=uid]
  llr[!llr$in_frame_fusion & llr$fusion & !llr$sense, msg := paste("Anti-sense fusion between:", bk_msg[1], "and", bk_msg[2]), by=uid]

}
#### look for UTR overlap
## load refseq
                                        #gr.ref <- readRDS("/xchip/beroukhimlab/Jeremiah/tracks/gr.refseq.all.rds")
                                        #fo <- gr.findoverlaps(gr, gr.ref[gr.ref$type == "UTR"])
                                        #llr$bk_msg[intersect(fo$query.id, which(!nchar(llr$bk_msg)))] <- "UTR"


##llr$bk_msg[!grepl("intron|exon", llr$bk_msg)] <- ""
print("...annotating intergenic")
fo <- suppressWarnings(gr.findoverlaps(llr, gene.comp))
if (nrow(fo)) {
  fo[, left  := gene.comp$left[subject.id],  by=query.id]
  fo[, right := gene.comp$right[subject.id], by=query.id]
  fo[, right_start  := end(gene.comp)[subject.id], by=subject.id]
  fo[, left_end := start(gene.comp)[subject.id], by=subject.id]
  fo[, left_dist  := start - left_end, by=query.id]
  fo[, right_dist := right_start - end, by=query.id]
  fo[, mmm := paste(left_dist, "bp to right of", left, "and", right_dist, "bp to left of", right)]
  llr$msg[!nchar(llr$msg)] <- ""
  llr$bk_msg[fo$query.id] <- ifelse(!nchar(llr$bk_msg[fo$query.id]), fo$mmm, llr$bk_msg[fo$query.id])
}

## annotate genes withn 100kb
llr$gene100kb <- gr.val(dt2gr(llr), genes + 100e3, 'gene', sep="_")$gene
if (!"gene100kb" %in% colnames(llr))
  llr$gene100kb <- ""

if ("bk_msg" %in% colnames(llr)) {
  mcols(bks)$break1 <- llr[grl.iix==1, bk_msg]
  mcols(bks)$break2 <- llr[grl.iix==2, bk_msg]
} else {
  mcols(bks)$break1 <- mcols(bks)$break2 <- ""
}

if ("msg" %in% colnames(llr)) {
  mcols(bks)$fusion <- llr[grl.iix==1, msg]
} else {
  mcols(bks)$fusion <- ""
}

if ("suspect_L1_retrotranspotiion" %in% colnames(llr)) {
  mcols(bks)$suspect_L1_retrotransposition <- llr[grl.iix==1, suspect_L1_retrotransposition]
} else {
  mcols(bks)$suspect_L1_retrotransposition <- ""
}


## convert to data frame
print("...writing data")
dt.bks <- as.data.frame(gr2dt(grl.unlist(bks)))
dt.bks <- rename(dt.bks, c("seqnames"="chr1", "start"="pos1", "strand"="strand1","altchr"="chr2","altpos"="pos2","altstrand"="strand2"))
print(dt.bks)

if (opt$svaba) {
  dt.bks <- dt.bks[dt.bks$grl.iix==1, c("grl.ix", "chr1","pos1","strand1", "chr2","pos2","strand2", "FILTER","NORMAL","TUMOR","SPAN","sample","uid",
                     "EVDNC","NUMPARTS","SCTG","DISC_MAPQ","NM","MATENM","MAPQ","REPSEQ",
                     "HOMSEQ","INSERTION","TUMALT","TUMCOV","TUMLOD","NORMCOV","NORMALT",
                     "NORMLOD","gene1","gene2","gene1_100kb","gene2_100kb", "fusion","break1", "break2",
                     "suspect_L1_retrotransposition")]
} else {
  dt.bks <- dt.bks[dt.bks$grl.iix==1, c("grl.ix", "chr1","pos1","strand1", "chr2","pos2","strand2",
                     "sample","uid",
                     "gene1","gene2","gene1_100kb","gene2_100kb",
                     "fusion","break1", "break2",
                     "suspect_L1_retrotransposition")]
}


## write the output
write.table(dt.bks, file=paste0(unlist(strsplit(file,'[/]'))[8], ".txt"), quote=FALSE, sep='\t', row.names=FALSE)
##skitools::write.htab(dt.bks, file=paste0(opt$output, ".html"))

############## make the circos plot
print("...making Circos plot")
data(UCSC.HG19.Human.CytoBandIdeogram);
chr.exclude <- NULL;
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);

## get the gene label dat
#genes <- genes[width(genes) < 2e6]
#fo1 <- gr.findoverlaps(gr1+10e3, genes)
#fo2 <- gr.findoverlaps(gr2+10e3, genes)

## annoying bug with seqinfo clash on 'c'
#if (length(fo1) && length(fo2)) {
#  fo <- c(fo1,fo2)
#} else if (length(fo1)) { 
#  fo <- fo1
#} else {
#  fo <- fo2
#}

print("...constructing Circos plot")
try( {
  ## set the gene labels
 gene.dat <- data.frame(Chromosome=c(dt.bks$chr1, dt.bks$chr2), chromStart=c(dt.bks$pos1, dt.bks$pos2),
                       chromEnd=c(dt.bks$pos1,dt.bks$pos2), Gene=as.character(c(dt.bks$gene1, dt.bks$gene2)),
                       stringsAsFactors=FALSE)
  if (nrow(gene.dat))
  gene.dat <- gene.dat[nchar(gene.dat$Gene) > 0,]
  if (nrow(gene.dat))
   gene.dat <- gene.dat[!duplicated(gene.dat$Gene),]

  links = data.frame()
  if (length(bks))
    links = with(dt.bks, data.frame(Chromosome=chr1, chromStart=pos1, chromEnd=pos1, Chromosome.1=chr2, chromStart.1=pos2, chromEnd.1=pos2))

  gene.dat$Chromosome = ifelse(grepl("^chr", gene.dat$Chromosome), gene.dat$Chromosome, paste0("chr", gene.dat$Chromosome))
  links$Chromosome = ifelse(grepl("^chr", links$Chromosome), links$Chromosome, paste0("chr", links$Chromosome))
  links$Chromosome.1 = ifelse(grepl("^chr", links$Chromosome.1), links$Chromosome.1, paste0("chr", links$Chromosome.1))

  ## plot the PDF
  pdf(file=paste0(unlist(strsplit(file,'[/]'))[8],".pdf"), height=opt$height, width=opt$width, compress=TRUE, useDingbats=FALSE);
  RCircos.Set.Plot.Area();
  RCircos.Chromosome.Ideogram.Plot();
  if (opt$genes != 0 && nrow(gene.dat) > 0) {
  track.num <- 1
  RCircos.Gene.Connector.Plot(gene.dat, track.num, "in");
  track.num <- 2;
  name.col <- 4;
  RCircos.Gene.Name.Plot(gene.dat, name.col, track.num, "in");
}
if (nrow(links) > 0) {
  track.num = 1
  RCircos.Link.Plot(links, track.num, by.chromosome=TRUE) ## by.chromosome is for color
}
})
dev.off()
}
