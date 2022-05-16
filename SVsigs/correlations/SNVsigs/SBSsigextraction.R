library(devtools)
library(data.table)
library(readxl)
library(BiocManager)
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
devtools::install_github('raerose01/deconstructSigs')
library(deconstructSigs)

################################################################################
# read in and subset maf of snvs, subset to only those in paper
mut_master <- readRDS('/Users/Alex/Dropbox (Partners HealthCare)/pHGG_Data/data/20201119palimsest_indels/full/20201119pHGG_M2vcf_annot.rds')

unique(mut_master$Sample)

metadata <- read_xlsx('/Users/Alex/Dropbox (Partners HealthCare)/manuscript draft Rameen:Mimi edit/manuscript/old/resubmission/tables/supp. table 1_Metadata_bySample.xlsx', skip = 1)

sc <- cbind(sort(unique(mut_master$Sample)), sort(metadata$Tumor_Sample_Barcode))

mut_master <- as.data.table(mut_master)
mut_master_sub <- mut_master[!(Sample %in% c('DIPG19T', 'DIPG58T'))]

sc <- cbind(sort(unique(mut_master_sub$Sample)), sort(metadata$Tumor_Sample_Barcode)) # looks good


################################################################################
#can run this on default, does hg19  --> makes sbs matrix
sbs_matrix <- mut.to.sigs.input(mut.ref = mut_master_sub,
                                sample.id = 'Sample',
                                chr = 'CHROM',
                                pos = 'POS',
                                ref = 'REF',
                                alt = 'ALT',
                                sig.type = 'SBS')

#Warning:
#                       In mut.to.sigs.input(mut.ref = mut_master_sub, sample.id = "Sample",  :
#                       Some samples have fewer than 50 mutations:
#                       DIPG110.t, DIPG18T

cosmic_sigs <- fread('/Volumes/xchip_beroukhimlab/Alex/refs/COSMIC_v3_SBS_GRCh37.txt')

################################################################################
#### check counts of actual SBS events
sum(rowSums(sbs_matrix))
# [1] 6053980
mut_master_sub[, REF := gsub('-',NA, REF)]
mut_master_sub[, ALT := gsub('-',NA, ALT)]


mut_master_sub[, len := nchar(REF)]
mut_master_sub[, len_alt := nchar(ALT)]
mut_master_sub[, tot_len := len + len_alt]
length(which(mut_master_sub$tot_len <= 2))
# [1] 6053980 --> looks good to mee

################################################################################
# read in COSMICv3

cosmic_sigs <- fread('/Volumes/xchip_beroukhimlab/Alex/refs/COSMIC_v3_SBS_GRCh37.txt')
cosmic_sigs_matrix <- as.matrix(cosmic_sigs, rownames = 'Type')
cosmic_sigs_matrix <- t(cosmic_sigs_matrix)
cosmic_sigs_matrix <- as.data.frame(cosmic_sigs_matrix)
## closest sigs to de novo
# SBS21, SBS15, SBS3, SBS44, SBS13, SBS16, SBS6, SBS40, SBS26, SBS15, SBS1, SBS31, SBS18, SBS57
closest_sigs <- c('SBS21', 'SBS15', 'SBS3', 'SBS44', 'SBS13', 'SBS16', 'SBS6', 'SBS40', 'SBS26', 'SBS15', 'SBS1', 'SBS31', 'SBS18', 'SBS57')
cosmic_sigs_matrix_subset <- cosmic_sigs_matrix[rownames(cosmic_sigs_matrix) %in% closest_sigs,]

weights_mat <- NULL
product_mat <- NULL
for(i in 1: nrow(sbs_matrix)) {
  print(rownames(sbs_matrix)[i])
  sbs_sub <- sbs_matrix[i,]
  sample <- rownames(sbs_matrix)[i]
  
  extracted_sigs <- whichSignatures(tumor.ref = sbs_sub,
                                    signatures.ref = cosmic_sigs_matrix_subset,
                                    sample.id = sample,
                                    contexts.needed = TRUE,
                                    tri.counts.method = 'default',
                                    associated = closest_sigs)
  
  product <- extracted_sigs$product
  product_mat <- rbind(product_mat, product)
  
  unknown  <- extracted_sigs$unknown
  weights <- extracted_sigs$weights
  
  exposure <- cbind(weights, unknown)
  weights_mat <- rbind(weights_mat, exposure)
}

rowSums(weights_mat)
write.table(weights_mat, '/Volumes/xchip_beroukhimlab/Alex/pHGG/20220207_sbs_sigs/20220207_extracted_sbssigs_weightsmat.txt', row.names = T, col.names = T, sep= '\t', quote = F)
saveRDS(weights_mat, '/Volumes/xchip_beroukhimlab/Alex/pHGG/20220207_sbs_sigs/20220207_extracted_sbssigs_weightsmat.rds')
write.table(product_mat, '/Volumes/xchip_beroukhimlab/Alex/pHGG/20220207_sbs_sigs/20220207_extracted_sbssigs_productmat.txt', row.names = T, col.names = T, sep= '\t', quote = F)

######## 
# comparison of denovo sigs to cosmic v3 sigs

denovo<- fread('/Users/Alex/Dropbox (Partners HealthCare)/Alex&Frank/data/pHGG_CosineSimilarity/extdatafig5b_data.txt')
motif_adj <- NULL
for(i in 1:nrow(denovo)) {
  motif_adj <- rbind(motif_adj, paste0(unlist(strsplit(denovo$Motif[i], ""))[3], unlist(strsplit(denovo$Motif[i], ""))[1], unlist(strsplit(denovo$Motif[i], ""))[2], unlist(strsplit(denovo$Motif[i], ""))[4]))
}
denovo_mat <- as.matrix(denovo)
rownames(denovo_mat) <- motif_adj
denovo_mat <- denovo_mat[,-1]
cosmic_sigs_matrix_mat <-t(cosmic_sigs_matrix)
rownames(cosmic_sigs_matrix_mat) <- gsub("[^[:alnum:] ]", "", rownames(cosmic_sigs_matrix_mat))

rownames_in <- which(rownames(denovo_mat) %in% rownames(cosmic_sigs_matrix_mat))
sc <- cbind(rownames(cosmic_sigs_matrix_mat), rownames(denovo_mat))

denovo_mat_ord <- denovo_mat[order(rownames(denovo_mat)),]
cosmic_sigs_matrix_mat_ord <- cosmic_sigs_matrix_mat[order(rownames(cosmic_sigs_matrix_mat)),]
sc <- cbind(rownames(cosmic_sigs_matrix_mat_ord), rownames(denovo_mat_ord))

library(CoSim)
res <- cosine_similarity(denovo_mat_ord, cosmic_sigs_matrix_mat_ord)
write.table(res, '/Volumes/xchip_beroukhimlab/Alex/pHGG/20220207_sbs_sigs/20220505_denovoSBS_cosmic_cosinesimilarity.txt', row.names = T, col.names = T, sep = '\t', quote = F)


