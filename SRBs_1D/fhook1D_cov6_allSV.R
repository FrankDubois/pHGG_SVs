#########
#create a vcf file
library(data.table)
library(GenomicRanges)
library(gUtils)
library(fishHook)
DIPG.base <-  "/xchip/beroukhimlab/Frank/DIPG/data/20190709_svaba_calls" 

samplemaster <- fread('/xchip/beroukhimlab/Frank/DIPG/Rstuff/20200120fhook/20191125DIPGsampleMaster_cbio.txt')
#keep only the csvs
DIPG_csv <- fread('/xchip/beroukhimlab/Frank/DIPG/data/20200225fhook/20200225SVmafs.txt')

#load data and covariates here 
LOCAL  <- FALSE
if (LOCAL) {gtracks.base = "/xchip/beroukhimlab/Jeremiah/GenomeTracks"} else {gtracks.base="/xchip/beroukhimlab/Jeremiah/GenomeTracks" }

#set system environment
library(ggplot2)
library(plotly)
library(MASS)

library(gUtils)
library(roverlaps)
BINSIZE <- 5e4
OVERLAP <- 500
NEGATIVECONTROL <- FALSE
library(gUtils)

print("...loading tracks/covariates data")

oldFRAG <- fread(file.path(gtracks.base,"Fragile", "fragile_genes_smith.hg19fp.txt"))
setnames(oldFRAG, c("V3","V4","V5"),c("seqnames","start","end"))
oldFRAG <- dt2gr(oldFRAG)

  FRAG_o <- fread(file.path(gtracks.base,"Fragile", "fragile_genes_smith.hg19fp.txt")) ## Smith
  FRAG <- fread(file.path(gtracks.base,"Fragile", "fragile_genes_smith.hg19fp.txt")) ## Smith
  setnames(FRAG, c("V3","V4","V5","V8"),c("seqnames","start","end","name")) ## Smit
  fragfreq <- fread(file.path(gtracks.base,"Fragile", "fragile_breakability.csv"))
  setkey(FRAG, name)
  setkey(fragfreq, name)
  FRAG <- fragfreq[FRAG]
  FRAG <- FRAG[!is.na(frequency)]

  #FRAG <- gr2dt(dt2gr(FRAG[frequency > 1]) + 2.5e5)
  FRAG <-FRAG[frequency > 1] # w
  ff <- dt2gr(FRAG) ### combine in the zeros
  strand(ff) <- '+'
  ff <- gr.fix(ff, si)
  ff <- GenomicRanges::setdiff(si2gr(si), ff)
  ff$frequency <- 0
  FRAG <- rbind(fill=TRUE, FRAG, ff)  

FRAG <- dt2gr(FRAG)
FRAG$score <- FRAG$frequency 
## Hyper fragile
#### GENES
gr.genes = sort(gr.fix(gr.nochr(with(fread("/xchip/beroukhimlab/Jeremiah/tracks/genes.hg19.ucsc.txt", sep="\t"), GRanges(chr, IRanges(beg, end), gene=symbol))), si))
gr.genes <- gr.genes[!grepl("^ULK|^NBPF|^MIR|^LOC|^OR|^SNO|^FAM|^SMN|^NF1P2|^POTEB|^RGPD5|^RGPD2|^SNAR|^NBPF", gr.genes$gene) & width(gr.genes) < 3e6]
## add missing
gr.add <- c(GRanges(11, IRanges(118305205,118399539), gene="KMT2A"),
                        GRanges(14, IRanges(106032614,107288051), gene="IGH@"),
                        GRanges(22, IRanges(22380474,23265085), gene="IGL@"),
                        GRanges(2, IRanges(89890568,90274235), gene="IGK@"),
                        GRanges(4, IRanges(91048683,92523370), gene="CCSER1"))
gr.genes <- c(gr.genes, gr.add)

## make the gene size track
gr.genesize = sort(gr.fix(gr.nochr(with(fread("/xchip/beroukhimlab/Jeremiah/tracks/genes.hg19.ucsc.txt", sep="\t"), GRanges(chr, IRanges(beg, end), gene=symbol))), si))
gr.genesize <- gr.genesize[!grepl("^ULK|^NBPF|^MIR|^LOC|^OR|^SNO|^FAM|^SMN|^NF1P2|^POTEB|^RGPD5|^RGPD2|^SNAR|^NBPF", gr.genesize$gene) & width(gr.genesize) < 3e6]
gr.genesize$gene.size <- width(gr.genesize)
ff <- GenomicRanges::setdiff(si2gr(si), gr.genesize)
ff$gene = ""
ff$gene.size <- 0
gr.genesize$query.id <- gr.genesize$tile.id <- NULL
gr.genesize <- c(ff, gr.genesize) #hm?


## make the gene density track
grall <- si2gr(gUtils::si)
grall <- grall[seqnames(grall) %in% c(seq(22),"X")]
gr.density <- trim(gr.fix(gr.tile(grall, w=1e6) + 1e6, si))
fo <- gr2dt(gr.findoverlaps(gr.density, gr.genes))
fo[, gene.count := nrow(.SD), by=query.id]
fo[, gene.density := gene.count / width(gr.density)[query.id], by=query.id]
gr.density$gene.density <- 0
gr.density$gene.density[fo$query.id] <- fo$gene.density
gr.density <- gr.density[width(gr.density) > 2e6]
gene_density <- gr.density - 1e6
gene_density <- gr.fix(gene_density, si)
ff <- GenomicRanges::setdiff(si2gr(si), gene_density)
ff$gene.density <- 0
gene_density$query.id <- gene_density$tile.id <- NULL
gene_density <- c(ff, gene_density)
gene_density$score <- gene_density$gene.density
if(genesize){
  gene_density$score <- gene_density$score*width(gene_density)
}

#mappability
mappability <- readRDS("/xchip/beroukhimlab/Jeremiah/Projects/ICGC/TRACKS/marcin/DB/Tracks/wgEncodeCrgMapabilityAlign100mer.gr.rds")
SINE <- with(fread(file.path(gtracks.base, "repeat_masker/repeat_masker_hg19_SINE.bed")), GRanges(V1, IRanges(V2, V3)))
LINE <- with(fread(file.path(gtracks.base, "repeat_masker/repeat_masker_hg19_LINE.bed")), GRanges(V1, IRanges(V2, V3)))
LTR <-  with(fread(file.path(gtracks.base, "repeat_masker/repeat_masker_hg19_LTR.bed")), GRanges(V1, IRanges(V2, V3)))
LCR <-  with(fread(file.path(gtracks.base, "segmental_dups/segmental_duplications.merge.bed")), GRanges(gsub("chr", "", V1), IRanges(V2, V3)))
GC_1K    <- suppressWarnings(with(fread(file.path(gtracks.base,"GCcontent/hg19.1000.gc5.bed"), header=FALSE), GRanges(gsub("chr", "", V1), IRanges(V2, V3), score=ifelse(V4==".", NA, as.numeric(V4)))))

introns <- readRDS("/xchip/beroukhimlab/Jeremiah/tracks/gr.introns.rds")
exons <- readRDS("/xchip/beroukhimlab/Jeremiah/tracks/gr.exons.rds")

reptimedata = readRDS(gzcon(file('http://mskilab.com/fishHook/hg19/RT_NHEK_Keratinocytes_Int92817591_hg19.rds')))
reptime = Cov(data = reptimedata, field = 'score', name = 'ReplicationTiming')

eligible <- rtracklayer::import(file.path(gtracks.base,'HengLiMask/um75-hs37d5.covered.bed'))
#create GRanges object for tissue specific covariate
#set binsize and create bins
hypotheses <- gr.tile(si2gr(gUtils::si), BINSIZE) + OVERLAP
#remove Y chromomsome
#change name of eligible
gr.eligible <- eligible

#cancer genes
CFLPAD = 50e3
gr.genes = sort(gr.fix(gr.nochr(with(fread("/xchip/beroukhimlab/Jeremiah/tracks/genes.hg19.ucsc.txt", sep="\t"), GRanges(chr, IRanges(beg, end), gene=symbol))), si))
gr.genes <- gr.genes[!grepl("^ULK|^NBPF|^MIR|^LOC|^OR|^SNO|^FAM|^SMN|^NF1P2|^POTEB|^RGPD5|^RGPD2|^SNAR|^NBPF", gr.genes$gene) & width(gr.genes) < 3e6]
## add missing
gr.add <- c(GRanges(11, IRanges(118305205,118399539), gene="KMT2A"),
                        GRanges(14, IRanges(106032614,107288051), gene="IGH@"),
                        GRanges(22, IRanges(22380474,23265085), gene="IGL@"),
                        GRanges(2, IRanges(89890568,90274235), gene="IGK@"),
                        GRanges(4, IRanges(91048683,92523370), gene="CCSER1"))
gr.genes <- c(gr.genes, gr.add)
f <- fread("/xchip/beroukhimlab/Jeremiah/Projects/ICGC/data/Census_allFri_Jan_18_15_13_27_2019.csv")
cgc.genes <- gr.genes[gr.genes$gene %in% f[Tier ==1, `Gene Symbol`]]
cgc.genes <- c(cgc.genes, c(GRanges(11, IRanges(118305205,118399539), gene="KMT2A"),
                                        GRanges(14, IRanges(106032614,107288051), gene="IGH@"),
                                        GRanges(22, IRanges(22380474,23265085), gene="IGL@")))
covs <- list()
covSINE        <- Cov(SINE,        name = "SINE")
covLTR         <- Cov(LTR,         name = "LTR")
covRepTime        <- Cov(repl,        name = "ReplicationTiming", field = "score",    type="numeric")
covsMappability <- Cov(mappability, name = "Mappability",       field = "score",    type="numeric")
covFRAG        <- Cov(FRAG,        name = "HyperFragility",         field = "frequency",type="numeric")
covGC_1K       <- Cov(GC_1K,          name = "GC_1K",             field = "score",    type="numeric")

covsGenedens    <- Cov(gene_density, name = "gene_density", field = "score",type="numeric")
covs_intFRAG         <- Cov(oldFRAG, name = "IntervalFragility", type = "interval")
hetchromdata = readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20200127fhook/hetchromdata.rds')
hetchromCov = Cov(hetchromdata, name = 'Heterochromatin') ## instantiate interval covariate
##
####################overlapWithGenes############
#x is a GRanges object
overlapWithGenes <- function(x) {
  x <- gr2dt(x)
  genes <- dt2gr(fread("/xchip/beroukhimlab/ofer/tracks/gencode_v19.txt"))
  genes <- genes %Q% (gene_type =="protein_coding")
  overlap <-roverlaps::roverlaps(x, genes, index_only = TRUE)
  overlap[, gene:= genes[subject.id]$gene]
  overlap[, gene:=paste(gene, collapse = "_"), by = query.id]
  overlap[, tile.id := x$tile.id[query.id]]
  overlap <- overlap[!duplicated(query.id)]
  setkey(overlap, tile.id)
  setkey(x, tile.id)
  x <- overlap[, .(gene, tile.id)][x]
  return(x)
}
                                        #PIPELINE FOR THE 1D MODEL###
# A VCF FILE FROM SVABA AND RETURNS A LIST OF BREAKS####
library(rtracklayer)
genes_a = gr.sub(import('http://mskilab.com/fishHook/hg19/gencode.v19.genes.gtf')) # rtracklayer::import reads gtf and gr.sub replaces chr

 oneDModelDIPG  <- function(vcf.file, type = "SV", filepath = NULL, BINSIZE  = 50e3, r2 = FALSE, FDR = 0.1){
  OVERLAP <- 500
  hypotheses <- gr.tile(si2gr(gUtils::si), BINSIZE) + OVERLAP
  hypotheses <- hypotheses %$% GZtad_gr_an  
  if(type=="SV")  {
    raw.svs <- vcf.file     
    raw.svs1 <- data.table(seqnames = raw.svs$chr1, start = raw.svs$pos1, uid = raw.svs$uid, sample = raw.svs$sid, stringsAsFactors = FALSE)
    raw.svs2 <- data.table(seqnames = raw.svs$chr2, start = raw.svs$pos2, uid = raw.svs$uid, sample = raw.svs$sid, stringsAsFactors = FALSE)
    breakpoints.all <- rbind(raw.svs1, raw.svs2)
    breakpoints.all.gr <- dt2gr(breakpoints.all)
    print("...formatted data")
    print(paste("The number of breakpoints is", length(breakpoints.all.gr)))
    print(paste("The number of unique samples is", length(unique(breakpoints.all.gr$ID))))
    print(paste("The bin size is", BINSIZE))
    fishDIPG <- Fish(hypotheses = hypotheses[, 'geneName'], events = breakpoints.all.gr, eligible = gr.eligible, idcol = "ID")
    fishDIPG$mc.cores = 4 
    fishDIPG$covariates <- c(covSINE, reptime, covsMappability, covs_intFRAG, covGC_1K, hetchromCov)
    fishDIPG$score()
    print("...ran model")
    saveRDS(fishDIPG, '/xchip/beroukhimlab/..._.rds')
    sigBins.annot <- fishDIPG$res
    sigBins.annot <- sigBins.annot %$% genes_a
    #create maf that maps breaks to the sample they came from
    maf <- breakpoints.all.gr %**%  sigBins.annot
    sigBins.annot_an = sigBins.annot %$% GZtad_gr_an
    write.table(sigBins.annot_an, "/xchip/beroukhimlab/..._.txt", row.names = FALSE, sep = "\t")    
    write.table(maf, "/xchip/beroukhimlab/..._.txt", row.names = FALSE, sep = "\t")
  } else {
    print("type not specified or incorrectly specificed")
  } 
   if(!is.null(filepath)){
      write.table(sigBins.annot, file = filepath, row.names = FALSE, sep = "\t")
    } else if (r2) {
      if(type == "SV") { 
      return(r2 <-cor(fishDIPG$res$count.pred, fishDIPG$res$count, use = "complete.obs"))
    } else {
      return(r2 <- cor(fishIndels$res$count.pred, fishIndels$res$count, use = "complete.obs"))
    }}
     else {
      return(sigBins.annot)
    }
}
  
GZtad_gr_an = readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20200127fhook/GZtad_gr_an.rds')
colnames(DIPG_csv)
          
all.RSB <- oneDModelDIPG(DIPG_csv, type = "SV", BINSIZE = 5e4, FDR = .25)
#merge bins


