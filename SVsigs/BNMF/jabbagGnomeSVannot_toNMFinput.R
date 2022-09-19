## this script reads in the JaBbA (https://github.com/mskilab/JaBbA) output, annotates and merges the JaBbA edges with the svaba SVs annotated by the Li et al 
## PCAWG method (https://github.com/cancerit/ClusterSV) to create the NMF input matrix for SV signature analysis
library(gGnome)
library(parallel)
fls <- list.files('/.../pHGG_Data/data/20220106jabba/out') # filepath to jabba.rds files (jabba output)
# get jabba alt edges ####
countmatrix_l <- mclapply(fls, function(xfl){
  print(xfl)
  ## load the graph for xfl
  jabbaX = gG(jabba = paste0('/Users/fdubois/Dropbox/pHGG_Data/data/20220106jabba/out/',xfl,'/jabba.rds'))
  ## Identify all supported SV event types
  jabbaX = events(jabbaX, verbose = T)
  jabbaX_edgesdt <- jabbaX$edgesdt # SVcalls
  jabbaX_ALTedgesdt <- jabbaX_edgesdt[type == 'ALT']
  jabbaX_ALTedgesdt[, sampleID := xfl]
  write.table(jabbaX_ALTedgesdt, paste0(
              '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220107gGnomeEventSigs/ALTedges/20220107gGnomeALTedges_',xfl,'.txt'), quote = FALSE, sep = "\t", row.names = F)
  if (nrow(jabbaX$meta$event) >=1) {
    retruntbl <- jabbaX$meta$event[, table(type)]
    retruntbl_dt <- as.data.table(retruntbl)
    retruntbl_dt[, sampleID := xfl]
    print(retruntbl_dt)
    return(retruntbl_dt)
  } else {
    retruntbl_dt <- as.data.table(t(c('no_events', 0)))
    colnames(retruntbl_dt) <- c('type', 'N')
    retruntbl_dt[, sampleID := xfl]
  }
}, mc.cores = 5, mc.silent = F)
# annotate alt edges ####
fls <- list.files('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220107gGnomeEventSigs/ALTedges')
fpth <- paste0('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220107gGnomeEventSigs/ALTedges/', fls)
ALTedges_l <- lapply(fpth, fread)
ALTedges <- rbindlist(ALTedges_l, fill = TRUE)

idx <- colnames(ALTedges) %in% c("cn", "type", "CHROM", "POS", 'REF', 'ALT', 'MATEID', 'SCTG', 'SPAN', 'HOMSEQ', 'INSERTION',
                                 'coord', 'mcoord', 'class', 'simple', 'cpxdm', 'tic', 'tyfonas', 'dm', 'bfb', 'fbi', 'tiny',
                                 "chromothripsis", "rigma", "del", "cn.max", "pyrgo", "dup", "chromoplexy", 'sampleID')
ALTedges_s <- ALTedges[,..idx]
length(unique(ALTedges_s$sampleID))

ALTedges_s[chromoplexy >0 , SVannot := 'chromoplexy']
ALTedges_s[!is.na(SVannot) & dup >0 , SVannot := paste(SVannot, 'dup')]
ALTedges_s[is.na(SVannot) & dup >0 , SVannot := 'dup']
ALTedges_s[!is.na(SVannot) & pyrgo >0 , SVannot := paste(SVannot, 'pyrgo')]
ALTedges_s[is.na(SVannot) & pyrgo >0 , SVannot := 'pyrgo']
ALTedges_s[!is.na(SVannot) & del >0 , SVannot := paste(SVannot, 'del')]
ALTedges_s[is.na(SVannot) & del >0 , SVannot := 'del']
ALTedges_s[!is.na(SVannot) & rigma >0 , SVannot := paste(SVannot, 'rigma')]
ALTedges_s[is.na(SVannot) & rigma >0 , SVannot := 'rigma']
ALTedges_s[!is.na(SVannot) & chromothripsis >0 , SVannot := paste(SVannot, 'chromothripsis')]
ALTedges_s[is.na(SVannot) & chromothripsis >0 , SVannot := 'chromothripsis']
ALTedges_s[!is.na(SVannot) & bfb >0 , SVannot := paste(SVannot, 'bfb')]
ALTedges_s[is.na(SVannot) & bfb >0 , SVannot := 'bfb']
ALTedges_s[!is.na(SVannot) & dm >0 , SVannot := paste(SVannot, 'dm')]
ALTedges_s[is.na(SVannot) & dm >0 , SVannot := 'dm']
ALTedges_s[!is.na(SVannot) & tyfonas >0 , SVannot := paste(SVannot, 'tyfonas')]
ALTedges_s[is.na(SVannot) & tyfonas >0 , SVannot := 'tyfonas']
ALTedges_s[!is.na(SVannot) & cpxdm >0 , SVannot := paste(SVannot, 'cpxdm')]
ALTedges_s[is.na(SVannot) & cpxdm >0 , SVannot := 'cpxdm']
ALTedges_s[!is.na(SVannot) & tic >0 , SVannot := paste(SVannot, 'tic')]
ALTedges_s[is.na(SVannot) & tic >0 , SVannot := 'tic']

table(ALTedges_s$simple)
table(ALTedges_s$SVannot)
ALTedges_s[!is.na(SVannot) & simple %in% c('INV1', 'INV2', 'INV3') , SVannot := paste(SVannot, 'inv')]
ALTedges_s[is.na(SVannot) &simple %in% c('INV1', 'INV2', 'INV3'), SVannot := 'inv']
ALTedges_s[!is.na(SVannot) & simple %in% c('TRA1', 'TRA2', 'TRA3', 'TRA4' ) , SVannot := paste(SVannot, 'tra')]
ALTedges_s[simple %in% c('TRA1', 'TRA2', 'TRA3', 'TRA4' ), SVannot := 'tra']
ALTedges_s[!is.na(SVannot) & simple == 'INVDUP1' , SVannot := paste(SVannot, 'invdup')]
ALTedges_s[is.na(SVannot) & simple == 'INVDUP1' , SVannot := 'invdup']
table(ALTedges_s$SVannot)

ALTedges_sRst <- ALTedges_s[is.na(SVannot)]
length(unique(ALTedges_sRst$sampleID))
# map to svaba SVs####
svmafs <- fread('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/draft/first_submission/tables/supp. table 2_202103210_SVmafs.txt')
colnames(svmafs)
svmafs[, svtype := ifelse(chr1 == chr2, ifelse(strand1 == strand2, ifelse(strand1 == '+', "h2hINV", "t2tINV"), 
                                               ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]

svmafs$sid <- gsub(pattern = '-', replacement = '_', x = svmafs$sid)
svmafs$sid <- gsub(pattern = '_STJ_WGS', replacement = '', x = svmafs$sid)
ALTedges_s$sampleID<- gsub(pattern = '_STJ_WGS', replacement = '', x = ALTedges_s$sampleID)
svmafs$sid <- gsub(pattern = 'rimary', replacement = '', x = svmafs$sid)
ALTedges_s$sampleID<- gsub(pattern = 'rimary', replacement = '', x = ALTedges_s$sampleID)
svmafs$sid <- gsub(pattern = 'ecurrent', replacement = '', x = svmafs$sid)
ALTedges_s$sampleID<- gsub(pattern = 'ecurrent', replacement = '', x = ALTedges_s$sampleID)
samplesmatch <- as.data.table(cbind(sort(unique(svmafs$sid)), sort(unique(ALTedges_s$sampleID))))
ALTedges_s$sampleID<- gsub(pattern = '-', replacement = '_', x = ALTedges_s$sampleID)
svmafs$sid <- gsub(pattern = '-pair', replacement = '', x = svmafs$sid)
svmafs$sid <- gsub(pattern = '_pair', replacement = '', x = svmafs$sid)
ALTedges_s$sampleID<- gsub(pattern = '_pair', replacement = '', x = ALTedges_s$sampleID)
ALTedges_s$sampleID<- gsub(pattern = '_T', replacement = '', x = ALTedges_s$sampleID)
ALTedges_s$sampleID<- gsub(pattern = 'OR', replacement = '_TOR', x = ALTedges_s$sampleID)
svmafs$sid <- gsub(pattern = 'ICGCp', replacement = '', x = svmafs$sid)

sort(unique(ALTedges_s$sampleID))[!(sort(unique(ALTedges_s$sampleID)) %in% sort(unique(svmafs$sid)))]
svmafs[, SVuid:= paste0(sid, ':', chr1, ':', pos1, ":", chr2, ':', pos2, ":", INSERTION)]
sum(as.numeric(duplicated(svmafs$SVuid)))
svmafs[duplicated(SVuid)]
ALTedges_s[is.na(INSERTION), INSERTION := '']
ALTedges_s[, SVuid:= paste0(sampleID, ':', coord, ":", mcoord, ":", INSERTION)]
sum(as.numeric(duplicated(ALTedges_s$SVuid)))
ALTedges_s$SVuid[!(ALTedges_s$SVuid %in% svmafs$SVuid)]
ALTedges_s_noMAP <- ALTedges_s[!(SVuid %in% svmafs$SVuid)]
SVsM <- merge(x = svmafs, y = ALTedges_s, by = 'SVuid', all = TRUE)
colnames(SVsM)
SVsMs <- SVsM[, c(1:12, 15:20, 22:23, 25:38, 46, 57, 62:63)]

write.table(SVsMs, 
            '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220202gGnomeSVs/20220202SVsAnnotationMerged.txt', quote = FALSE, sep = "\t", row.names = F)
table(SVsMs$SVannot)
SVsMs_ann <- as.data.table(SVsMs)
SVmaf <- readRDS('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/pHGG_SVmaf_clstSV.rds') # Li et al annotation of clustered SVs (PCAWG method)
SVmaf$sid <- gsub(pattern = '-', replacement = '_', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = '_STJ_WGS', replacement = '', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = 'rimary', replacement = '', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = 'ecurrent', replacement = '', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = '-pair', replacement = '', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = '_pair', replacement = '', x = SVmaf$sid)
SVmaf$sid <- gsub(pattern = 'ICGCp', replacement = '', x = SVmaf$sid)

SVsMs_ann[is.na(sid), sid := sampleID]
samplesmatch <- as.data.table(cbind(sort(unique(SVmaf$sid)), sort(unique(SVsMs_ann$sid))))
SVsMs_ann$sid[!(SVsMs_ann$sid %in% SVmaf$sid)]
SVmaf$sid[!(SVmaf$sid %in% SVsMs_ann$sid)]
SVmaf[, SVuid:= paste0(sid, ':', chr1, ':', pos1, ":", chr2, ':', pos2, ":", INSERTION)]
sum(duplicated(SVmaf$SVuid))
SVmafDups <- SVmaf[duplicated(SVuid)]
SVmafdeduped <- SVmaf[!duplicated(SVuid)]
colnames(SVmafdeduped)
SVmafdeduped_s <- SVmafdeduped[, c(61, 54)]

SVsM_ann <- merge(x= SVsMs_ann, y =  SVmafdeduped_s, by = "SVuid") # merge jabba and Li et al/ PCAWG annotations
table(SVsM_ann$SVannot)
SVsM_ann <- as.data.table(SVsM_ann)
SVsM_ann[is.na(SVannot) & clust_size <=2, SVannot:= tolower(svtype)]
SVsM_ann$SVannot <- gsub(pattern = 'h2hinv', replacement = 'inv', x = SVsM_ann$SVannot)
SVsM_ann$SVannot <- gsub(pattern = 't2tinv', replacement = 'inv', x = SVsM_ann$SVannot)
SVsM_ann$SVannot <- gsub(pattern = 'inter', replacement = 'tra', x = SVsM_ann$SVannot)
SVsM_ann[is.na(SVannot) & clust_size >2, SVannot:= 'complx_Unclassf']
sort(table(SVsM_ann$SVannot))
SVsM_ann[is.na(SVannot), SVannot:= tolower(gsub(pattern = '-like', replacement = '', x = class)), by = class]
sort(table(SVsM_ann$SVannot))

write.table(SVsM_ann, 
            '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220202gGnomeSVs/20220202SVsAnnCompNST_Merged.txt', quote = FALSE, sep = "\t", row.names = F)
# count combinations as 0.5 each
SVsM_ann_Homs <- SVsM_ann[!grepl(pattern = ' ', x = SVannot)]
sort(table(SVsM_ann_Homs$SVannot))

SVsM_ann_Divs1 <- SVsM_ann[grepl(pattern = ' ', x = SVannot)]
sort(table(SVsM_ann_Divs1$SVannot))
SVsM_ann_Divs1[, SVannot := unlist(strsplit(x = SVannot, split = " ",fixed = T))[1], by = SVannot]
sort(table(SVsM_ann_Divs1$SVannot))

SVsM_ann_Divs2 <- SVsM_ann[grepl(pattern = ' ', x = SVannot)]
# sort(table(SVsM_ann_Divs2$SVannot))
SVsM_ann_Divs2[, SVannot := unlist(strsplit(x = SVannot, split = " ",fixed = T))[2], by = SVannot]
sort(table(SVsM_ann_Divs2$SVannot))

SVsM_ann_BPs <- rbind(SVsM_ann_Homs, SVsM_ann_Homs, SVsM_ann_Divs1, SVsM_ann_Divs2, fill = T)
sort(table(SVsM_ann_BPs$SVannot))

countmatrix_M <- dcast(data = SVsM_ann_BPs, formula = sid ~ SVannot, fill = 0)

write.table(countmatrix_M, 
            '/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220202gGnomeSVs/20220202BPsSVAnnotCountmatrix.txt', quote = FALSE, sep = "\t", row.names = F)

countmatrix_Ms_mat <- as.matrix(countmatrix_M, rownames='sid')
rowSums(countmatrix_Ms_mat)
saveRDS(countmatrix_Ms_mat, '/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20220203gGnomeSVsigs/20220203gGnomeSVsigs_nmfInputMat.RDS') # input for NMF


