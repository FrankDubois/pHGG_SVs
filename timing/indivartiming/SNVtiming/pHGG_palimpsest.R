library(Palimpsest)
library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome of choice
#datadir <- "/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/R stuff/20200703palimpsest/LiC1162" 
#load(file.path(datadir,"vcf.RData"))
require(data.table)
pHGG_M2_maf <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200904palimsests/20200903M2mergeUniq.rds')
unique(pHGG_M2_maf$Tumor_Sample_Barcode)
pHGG_M2vcf <- pHGG_M2_maf
colnames(pHGG_M2vcf)
table(pHGG_M2vcf$Variant_Classification)
pHGG_M2vcf_ex <- pHGG_M2vcf[!(Variant_Classification %in% c('IGR', 'Intron', 'Silent', "3'UTR", "5'UTR"))]
colnames(pHGG_M2vcf_ex)
colnames(pHGG_M2vcf_ex)[c(2,3,10,7,9, 6)] <- c('CHROM', 'POS', 'Sample', 'REF', 'ALT', 'Type') 
pHGG_M2vcf_ex[Type == 'SNP', Type := 'SNV']
pHGG_M2vcf_ex[Type %in% c('DNP', 'TNP'), Type := 'INS']

table(pHGG_M2vcf_ex$Type)
colnames(pHGG_M2vcf_ex)
pHGG_M2vcf_ex[371475,]
pHGG_M2vcf_ex[is.na(POS), POS := Start_Position]
pHGG_M2vcf_ex[is.na(End_position), End_position := End_Position]
sum(pHGG_M2vcf_ex[is.na(POS), POS])
colnames(pHGG_M2vcf_ex)
pHGG_M2vcf_ex2 <- pHGG_M2vcf_ex[,-c(16,17)]
int_genes <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200513cbioList/20200528int_genes.txt')
drivers <- int_genes$gene
pHGG_M2vcf_ex2[Hugo_Symbol %in% int_genes$gene, Driver := paste0(Hugo_Symbol, ' ', Protein_Change)]

pHGG_M2vcfEX_annot <- annotate_VCF(vcf = pHGG_M2vcf_ex2, ref_genome = BSgenome.Hsapiens.UCSC.hg19)
saveRDS(pHGG_M2vcfEX_annot, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200904palimsests/20200904pHGG_M2vcfEX_annot.rds')
pHGG_M2vcfEX_annot <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200904palimsests/20200904pHGG_M2vcfEX_annot.rds')
#sigs?
SBS_input <- palimpsest_input(vcf = pHGG_M2vcfEX_annot, Type = "SBS")

# select desired COSMIC SBS reference signatures 
SBS_pHGG_names <- c("SBS1", "SBS3", "SBS5", "SBS6", "SBS13", "SBS4", "SBS12", "SBS16", "SBS17b", "SBS18", "SBS22",
                    "SBS23", "SBS24", 'SBS15', 'SBS8', 'SBS13')

for(new_name in c(SBS_pHGG_names[SBS_pHGG_names %!in% names(sig_cols)]))
sig_cols[new_name] <- signature_colour_generator(new_name)  ## generate colours for new signatures 

SBS_pHGG_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_pHGG_names,]

resdir <- '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/'
# calculate and plot the exposure of the signatures across the series
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_pHGG_sigs,
                                        threshold = 6, signature_colours = sig_cols, resdir = resdir, 
                                        input_vcf = pHGG_M2vcfEX_annot)
saveRDS(SBS_signatures_exp, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907SBS_signatures_exp.rds')
#SBS_signatures_exp <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200818palimsest/20200818SBS_signatures_exp.rds')
# This step is quite computationally-intensive for whole genome data. 
# For this example we restrict the analysis to coding mutations 
vcf.cod <- pHGG_M2vcfEX_annot[(!is.na(pHGG_M2vcfEX_annot$Driver) & pHGG_M2vcfEX_annot$Type=="SNV"),]
vcf.cod <- signature_origins(input = vcf.cod, Type = "SBS", input_signatures = SBS_pHGG_sigs,
                             signature_contribution = SBS_signatures_exp)

# Estimate and represent the cumulative contribution of signatures to each driver gene
int_genes <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200513cbioList/20200528int_genes.txt')
drivers <- int_genes$gene
matprob <- matrix(nrow=length(drivers),ncol=length(SBS_pHGG_names),dimnames=list(drivers, SBS_pHGG_names))
sig.cols <- paste0(rownames(SBS_pHGG_sigs),".prob")#grep("prob",colnames(vcf.cod))
for(i in 1:nrow(matprob)){
  g <- rownames(matprob)[i]
  ind <- which(vcf.cod$gene_name==g)
  matprob[i,] <- apply(vcf.cod[ind,sig.cols],2,sum,na.rm=T)
}
barplot(t(matprob),col = sig_cols,border = sig_cols,las=2)
legend("top",names(sig_cols)[names(sig_cols) %in% rownames(SBS_pHGG_sigs)],fill=sig_cols,ncol=5,
       cex=0.75,bty ="n",inset = c(0,-0.3),xpd = T)

# clonality CNA?
#datadir <- "/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/R stuff/20200703palimpsest/LiC1162" 

#####
lf <- list.files('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20190905absolute_segs/')
lfps <- paste0('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20190905absolute_segs/',lf)
x=1
absSegs <- rbindlist(lapply(lfps, fread))
colnames(absSegs)
colnames(absSegs)[1:4] <- c('Sample', "CHROM","POS_START","POS_END" )
absSegs_sub <- absSegs[,c(1:4, 9, 15, 37, 38)]
colnames(absSegs_sub)[6:8] <- c('ntot', "Nmin","Nmaj")
absSegs_sub[,LogR := log2(total_copy_ratio)]
metaD <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200823MetaDpp/20200823DIPGsampleMaster.txt')
absSegs_sub$Sample <- gsub(pattern = '-', replacement = '_', x = absSegs_sub$Sample)
metaD$Tumor_Sample_Barcode_L <- gsub(pattern = '-', replacement = '_', x = metaD$Tumor_Sample_Barcode_L)
metaD[175:179,Tumor_Sample_Barcode_L := paste0(Tumor_Sample_Barcode_L, '_T'), by = Tumor_Sample_Barcode_L]
#absSegs_sub$Sample <- gsub(pattern = 'Primary', replacement = 'P', x = absSegs_sub$Sample)
#absSegs_sub$Sample <- gsub(pattern = 'Recurrent', replacement = 'R', x = absSegs_sub$Sample)
absSegs_sub$Sample <- gsub(pattern = '048_pair', replacement = '048Primary', x = absSegs_sub$Sample)
absSegs_sub <- absSegs_sub[Sample != 'DIPG58_TOR_pair']
samplesMerge <- as.data.table(cbind(sort(unique(absSegs_sub$Sample)), sort(metaD$Tumor_Sample_Barcode_L)))
colnames(samplesMerge) <- c('absSeg', 'Meta')
samplesMerge$absSeg == samplesMerge$Meta

absSegs_sub[,Ploidy := metaD[Tumor_Sample_Barcode_L == Sample, ploidy], by = Sample]
#####

metaD_s <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200828palimsets/20200828metaD_s.txt')
absSegs_sub <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200828palimsets/20200828absSegs_sub.rds')


###### Calculate the Cancer Cell Fraction (CCF) of each mutation.#####
pHGG_M2vcfEX_annot2 <- pHGG_M2vcfEX_annot[!pHGG_M2vcfEX_annot$Sample %in% c('DIPG19T','DIPG58T'),]
absSegs_sub$Sample <- gsub(pattern = 'ICGCp', replacement = '', x = absSegs_sub$Sample)
metaD_s$Sample <- gsub(pattern = 'ICGCp', replacement = '', x = metaD_s$Sample)

samplesMerge2 <- as.data.table(cbind(sort(unique(pHGG_M2vcfEX_annot2$Sample)), sort(unique(absSegs_sub$Sample)), sort(metaD_s$Sample)))
colnames(samplesMerge2) <- c('M2vcf', 'absSeg', 'Meta')
#length(unique(pHGG_M2vcfEX_annot2$Sample))
#metaD_s[,TScode := Sample]
#metaD_s[,Sample := samplesMerge2[Meta == TScode, M2vcf], by= TScode]
#absSegs_sub[,TScode := Sample]
#absSegs_sub[,Sample := samplesMerge2[absSeg == TScode, M2vcf], by= TScode]
colnames(pHGG_M2vcfEX_annot2)[12] = 'Tumor_Varcount'
pHGG_M2vcfEX_annot2$Tumor_Depth = pHGG_M2vcfEX_annot2$Tumor_Varcount + pHGG_M2vcfEX_annot2$t_ref_count
max(pHGG_M2vcfEX_annot2$n_alt_count)
pHGG_M2vcfEX_annot2_notnorm <- pHGG_M2vcfEX_annot2[pHGG_M2vcfEX_annot2$n_alt_count>0,]
colnames(pHGG_M2vcfEX_annot2)[15] = 'Normal_Depth'

#pHGG_M2vcfEX_annot2$Normal_Depth = 80 # add column from original mafs


#####
vcf_cna <- cnaCCF_annot(vcf= pHGG_M2vcfEX_annot2, annot_data = metaD_s,cna_data = absSegs_sub, CCF_boundary = 0.95)
View(cnaCCF_annot)
saveRDS(vcf_cna, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907vcf_cna.rds')

# Estimate the contribution of each signature to clonal and subclonal mutations in each tumour
resdir <- file.path('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/',"Signatures_early_vs_late/");if(!file.exists(resdir)){dir.create(resdir)}

vcf.clonal <- vcf_cna[which(vcf_cna$Clonality=="clonal"),]
SBS_input_clonal <- palimpsest_input(vcf = vcf.clonal,Type = "SBS")
signatures_exp_clonal <- deconvolution_fit(input_matrices = SBS_input_clonal,
                                           input_signatures = SBS_pHGG_sigs, resdir =  resdir, input_vcf =vcf.clonal,
                                           save_signatures_exp = T)

vcf.subclonal <- vcf_cna[which(vcf_cna$Clonality=="subclonal"),]
SBS_input_subclonal <- palimpsest_input(vcf = vcf.subclonal,Type = "SBS")
signatures_exp_subclonal <- deconvolution_fit(input_matrices = SBS_input_subclonal,
                                              input_signatures = SBS_pHGG_sigs, resdir =  resdir, input_vcf =vcf.subclonal,
                                              save_signatures_exp = F)
#-------------------------------------------------------------------------------------------------
# 9] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir_parent = '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/'
resdir <- file.path(resdir_parent,"ChromosomeDups_timing/");if(!file.exists(resdir)){dir.create(resdir)}

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf_cna,cna_data = absSegs_sub,cyto=cytoband_hg19)
vcf_cna <- chrom_dup_time$vcf;point.mut.time <- chrom_dup_time$point.mut.time;cna_data <- chrom_dup_time$cna_data

# Visualising timing plots
chrTime_plot(vcf = vcf_cna, point.mut.time = point.mut.time, resdir = resdir,cyto = cytoband_hg19)


#-------------------------------------------------------------------------------------------------
# 10] Structural variant (SV) signature analysis:
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"SV_signatures/");if(!file.exists(resdir)){dir.create(resdir)}
library(bedr);
library(RCircos) # Loading dependencies necessary for SV annotation and CIRCOS plots
svabMaf <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200901palimsests/20200225SVmafs.txt')
colnames(svabMaf)
svabMaf[, svtype := ifelse(chr1 == chr2,ifelse(strand1 == strand2, ifelse(strand1 == '+', "h2hINV", "t2tINV"), ifelse(strand1 == "+", "DEL", "DUP")),"INTER")]

svabMaf_s <- svabMaf[,c(12,38,2:7, 24,25,27)]
unique(svabMaf_s$svtype)
intGenes <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200404cbio/20200414int_genes.txt')
colnames(svabMaf_s) <- c('Sample', "Type","CHROM_1","POS_1",'strand1', "CHROM_2","POS_2",'strand2',"Tumor_Varcount","Tumor_Depth","Normal_Depth")
svabMaf_s[,CHROM_1 := paste0('chr', CHROM_1), by = CHROM_1]
svabMaf_s[,CHROM_2 := paste0('chr', CHROM_2), by = CHROM_2]
class(svabMaf_s)
svabMaf_s <- as.data.frame(svabMaf_s)
# Preprocess SV inputs and annotate for further analysis:
SV.vcf <- preprocessInput_sv(input_data =  svabMaf_s,resdir = resdir)
SV.vcf_dt <- as.data.table(SV.vcf)
SV.vcf_dt[Associated.Gene.Name %in%intGenes$gene , Driver := paste0(Associated.Gene.Name, '_', Type), by = Associated.Gene.Name]

#-------------------------------------------------------------------------------------------------
# 11] Visualise the natural history of tumour samples:
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir_parent,"Natural_history/");if(!file.exists(resdir)){dir.create(resdir)}

samplesMergeSV <- as.data.table(cbind(sort(unique(vcf_cna$Sample)), sort(unique(SV.vcf_dt$Sample))))
colnames(samplesMergeSV) <- c('vcf_cna', 'SV_vcf_dt')

SV.vcf_dt[,TScode := Sample]
SV.vcf_dt[,Sample := samplesMergeSV[SV_vcf_dt == TScode, vcf_cna], by= TScode]

palimpsest_plotTumorHistories(vcf = vcf_d, sv.vcf = SV.vcf_dt,cna_data =  cna_data, 
                              point.mut.time = point.mut.time, 
                              clonsig = signatures_exp_clonal$sig_props, 
                              subsig = signatures_exp_subclonal$sig_props, 
                              msigcol = sig_cols,  msigcol.sv =sig_cols,  resdir = resdir)
sig_cols <- c(sig_cols,'black')
names(sig_cols)

#vcf_cna_s_dt[!is.na(Driver),SBS.Sig.max :='']
# Add signature probability columns to the original vcf
vcf_d <- merge(vcf_cna,vcf.cod,all=TRUE,sort=FALSE)
SV.vcf_dt[!is.na(Driver),SV.Sig.max :='']
sort(unique(vcf_d$Sample))
sort(unique(vcf.clonal$Sample))
sort(unique(vcf.subclonal$Sample))
vcf_ds <- vcf_d[vcf_d$Sample %in% vcf.subclonal$Sample,]

resdir <- file.path(resdir_parent,"Natural_history_sub/");if(!file.exists(resdir)){dir.create(resdir)}

palimpsest_plotTumorHistories(vcf = vcf_ds, sv.vcf = SV.vcf_dt,cna_data =  cna_data, 
                              point.mut.time = point.mut.time, 
                              clonsig = signatures_exp_clonal$sig_props, 
                              subsig = signatures_exp_subclonal$sig_props, 
                              msigcol = sig_cols,  msigcol.sv =sig_cols,  resdir = resdir)

saveRDS(vcf_ds, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907vcf_ds.rds')
saveRDS(vcf_d, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907vcf_d.rds')
saveRDS(SV.vcf_dt, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907SV.vcf_dt.rds')
saveRDS(cna_data, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907cna_data.rds')
saveRDS(point.mut.time, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200907palimsests/20200907point.mut.time.rds')

