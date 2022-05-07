library(tidyverse)
require(data.table)
CNAwintable <- readRDS('/Users/fdubois/Dropbox/pHGG_data/data/20210118_SCNAwinningtables/20210118CNAwin_ann_sub.rds')
unique(CNAwintable$timing_state)
colnames(CNAwintable)
CNAwintable_s <- CNAwintable[,c(6, 29, 27)]
colnames(CNAwintable_s)
int_genesL <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20200513cbioList/20200528int_genes.txt')
VarCl4mat <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20200709VarCl4enrichment/20200709Varmatrix_NMForderedCL4_uniqArms.txt')
VarCl4mat_genes <- VarCl4mat$mut.name
VarCl4mat_genes <- sapply(VarCl4mat$mut.name, function(x){
  unlist(strsplit(x, split = '[:]'))[1]
})
int_genes <- unique(c(VarCl4mat_genes, int_genesL$gene))
#palimsets
#arm-level SCNAs
point.mut.time <- readRDS('/Users/fdubois/Dropbox/pHGG_data/data/20200907palimsests/20200907point.mut.time.rds')
unique(point.mut.time$Sample)
#SNVs
vcf_d <- readRDS('/Users/fdubois/Dropbox/pHGG_data/data/20200907palimsests/20200907vcf_d.rds')
vcf_dt <- as.data.table(vcf_d)
vcf_dt[, timing_state := paste0(Timing, '_', Clonality)]
unique(vcf_dt$timing_state)
vcf_dt[timing_state == 'NA_subclonal', timing_state2 := 'subclonal']
vcf_dt[timing_state == 'NA_clonal', timing_state2 := 'clonal_NA']
vcf_dt[timing_state == 'early_clonal', timing_state2 := 'early_clonal']
vcf_dt[timing_state == 'late_clonal', timing_state2 := 'late_clonal']
vcf_dt[timing_state == 'early_subclonal', timing_state2 := 'subclonal']
vcf_dt[timing_state == 'late_subclonal', timing_state2 := 'subclonal']
vcf_dt_s2 <- vcf_dt[!is.na(timing_state2)]
colnames(vcf_dt_s)
sort(unique(vcf_dt_s2$Sample))

vcf_dt_s <- vcf_dt_s2[Hugo_Symbol %in% int_genes]
vcf_dt_s[, Protein_Change := gsub(pattern = 'K28M', replacement = 'K27M', x = Protein_Change), by = Protein_Change]
vcf_dt_s[, Protein_Change := gsub(pattern = 'G35R', replacement = 'G34R', x = Protein_Change), by = Protein_Change]
vcf_dt_s[Hugo_Symbol %in% c('H3F3A', 'HIST1H3C', 'HIST1H3B') & Protein_Change %in% c('p.K27M', 'p.G34R'), Hugo_Symbol := paste0(Hugo_Symbol, Protein_Change)]
vcf_dt_s3 <- vcf_dt_s[!Variant_Classification %in% c('RNA', "5'Flank")]
table(vcf_dt_s3$Variant_Classification)
vcf_dt_sub <- vcf_dt_s3[,c(10, 65, 1)]
colnames(vcf_dt_sub) <- c("Sample","timing_state","Gene_Name")
colnames(CNAwintable_s) <- c("Sample","timing_state","Gene_Name")
unique(vcf_dt_sub$Sample)

unique(wintable_CN_plmstsSNV$Gene_Name)
#intGenes <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20200404cbio/20200414int_genes.txt')
vcf_dt_sub[, Gene_Name := paste0(Gene_Name, '_snv')]
wintable_CN_plmstsSNV <- rbind(CNAwintable_s, vcf_dt_sub)
sort(unique(wintable_CN_plmstsSNV$Sample))
unique(wintable_CN_plmstsSNV$timing_state)

#win_table_vcf = wintable_CN_plmstsSNV_deduped_X
## function to trim vcfs and make winning tables
mk_win_table <- function(win_table_vcf) {
  vcf_timing <- win_table_vcf %>%
    add_count(Gene_Name, timing_state) %>%
    mutate(result = case_when(grepl("early_clonal", timing_state) ~ 2,
                              grepl("clonal_NA", timing_state) ~ 1,
                              grepl("subclonal", timing_state) ~ -2,
                              grepl("late_clonal", timing_state) ~ -1,
                              TRUE ~ 0)) %>%
    dplyr::select(Gene_Name, result, n) %>%
    distinct() %>%
    dplyr::select(Gene_Name, result, n) %>%
    ## get weighted counts by multiplying events and occurrence
    mutate(weighted_result = result * n) %>%
    dplyr::rename(freq = n) %>%
    group_by(Gene_Name) %>%
    mutate(mean_result = sum(weighted_result)) %>%
    ungroup() %>%
    arrange(Gene_Name, mean_result) %>%
    distinct(Gene_Name, mean_result) %>%

    mutate(timing = case_when(mean_result < 0 ~ "LATE",
                              mean_result > 0 ~ "EARLY",
                              mean_result == 0 ~ "DRAW",
                              TRUE ~ "NONE"))
  ## let driver gene compete with one another: combination matrix
  wintbl_combn <- choose(nrow(vcf_timing), 2)
  snv_df_wintbl <- as.data.frame(matrix(apply(vcf_timing %>%
                                                dplyr::select(Gene_Name, mean_result),2, combn, m = 2),
                                        nrow = wintbl_combn),
                                 stringsAsFactors = F) %>%
    tbl_df() %>%
    ## rename V1 and V2 to P1 and P2, respectively
    dplyr::rename_at(vars(V1, V2), ~gsub("V", "P", .)) %>%
    ## where _vaf is clonality based estimate if event
    ## is ~clonal (early) or ~subclonal
    dplyr::rename(P1_vaf = V3,
                  P2_vaf = V4) %>%
    ## take a win/loss/draw call
    mutate(outcome = case_when(as.numeric(P1_vaf) > as.numeric(P2_vaf) ~ "P1",
                               as.numeric(P1_vaf) < as.numeric(P2_vaf) ~ "P2",
                               as.numeric(P1_vaf) == as.numeric(P2_vaf) ~ "D",
                               TRUE ~ NA_character_)) %>%
    dplyr::select(P1, P2, outcome)
  return(snv_df_wintbl)
}
VarCl4memb <- fread('/Users/fdubois/Dropbox/pHGG_data/data/20201227SRSVsVarcontext/20201227metaD_SRSVunique.txt')
table(VarCl4memb$SRSVsGroup)

wintable_CN_plmstsSNV$Sample <- gsub(pattern = '-', replacement = "_", x = wintable_CN_plmstsSNV$Sample)
wintable_CN_plmstsSNV$Sample <- gsub(pattern = '_WGS', replacement = "", x = wintable_CN_plmstsSNV$Sample)
samplesmatch <- as.data.table(cbind(sort(VarCl4memb$Tumor_Sample_Barcode), sort(unique(wintable_CN_plmstsSNV$Sample))))
colnames(samplesmatch) <- c('VarCl4', 'wintable_CN_plmstsSNV_Sample')

VarCl4memb_s <- VarCl4memb[SRSVsGroup %in% c('OncAmp', 'TADandNoncAmp'), c(4,11)]
#unique(wintable_CN_plmstsSNV$Sample)[!(unique(wintable_CN_plmstsSNV$Sample) %in% VarCl4memb$membership)]
VarCl4memb_s$Tumor_Sample_Barcode[!(VarCl4memb_s$Tumor_Sample_Barcode %in% unique(wintable_CN_plmstsSNV$Sample))]

wintable_CN_plmstsSNV <- wintable_CN_plmstsSNV[Gene_Name != 'CCDC26_snv']
wintable_CN_plmstsSNV[, ID := paste0(Sample, timing_state, Gene_Name)]
###
wintable_CN_plmstsSNV_deduped <- wintable_CN_plmstsSNV[!duplicated(ID)]
# Timing Plots ####
genes2plot_L = fread('/Users/fdubois/Dropbox/pHGG_data/data/20210118BTmodelHistoneGroups/20210118intgenes_pl_sub_s.txt')
#genes2plot_L <- genes2plot_L[1:29,]
testgenes <- unique(wintable_CN_plmstsSNV$Gene_Name)
sort(testgenes)
genes2plot_a <- unlist(lapply(genes2plot_L$genes, function(x){
  #print(x)
  unlist(grep(pattern = x, x = testgenes, value = TRUE))
}))
genes2plot <- genes2plot_a[!grepl(pattern = ',', x = genes2plot_a)]
genes2plot <- genes2plot[!grepl(pattern = 'WGD', x = genes2plot)]
genes2plot <- genes2plot[!genes2plot %in% c('ZNF133_snv', 'Mono-allelic Gain:7p15.1_Amp', 'Mono-allelic Gain:2p25.1_Amp', 
                                            'CN-LOH:11p')]
sort(unique(genes2plot))
library(BradleyTerryScalable)
library(igraph)
library(boot)
library(ggplot2)
library(ggthemes)

colnames(VarCl4memb_s)
x = unique(VarCl4memb_s$SRSVsGroup)[1]
lapply(unique(VarCl4memb_s$SRSVsGroup), function(x){
  print(x)
  VarCl4membX <- VarCl4memb_s[SRSVsGroup == x]
  wintable_CN_plmstsSNV_deduped_X <- wintable_CN_plmstsSNV_deduped[Sample %in% VarCl4membX$Tumor_Sample_Barcode]
  print(sort(table(wintable_CN_plmstsSNV_deduped_X$Gene_Name), decreasing = TRUE))
  # set plotgenes cutoff
  plotnames <- names(table(wintable_CN_plmstsSNV_deduped_X$Gene_Name)[table(wintable_CN_plmstsSNV_deduped_X$Gene_Name) >2]) ### set recurrence cutoff for plotting
  win_table_inputF <- mk_win_table(wintable_CN_plmstsSNV_deduped_X)
  table(c(win_table_inputF$P1, win_table_inputF$P2))
  testgenes <- unique(wintable_CN_plmstsSNV_deduped$Gene_Name)

  #### create btdata ####       
  data_4col <- codes_to_counts(win_table_inputF, c("P1", "P2", "D"))
  btdata <- btdata(data_4col, return_graph = TRUE)
    print(summary(btdata))
    fit_MAP <- btfit(btdata, 1.1)
    BTsummary2 <- summary(fit_MAP, SE = TRUE)
    BTsummary2_itemSummary <- BTsummary2$item_summary
    BTcoeffs2 <- coef(fit_MAP, as_df = TRUE )
    #?vcov
    #vcov(fit_MAP, ref = "home")
    #### get winning prob ####
    fitted(fit_MAP, as_df = TRUE)
    btprob(fit_MAP, as_df = TRUE)
    fit_MAP_prob <- btprob(fit_MAP)
    fit_map_df <- btprob(fit_MAP, as_df = TRUE)
    fit_map_mat <- as.data.frame(as.matrix(fit_MAP_prob))
    ## convert to rank-order df
    fit_map_mat_ranked <- fit_map_mat %>%
      rownames_to_column("driver") %>%
      mutate_at(vars(-driver), ~(1 - cume_dist(.)))
    ## rank order long table
    fit_map_mat_ranked_long_tbl <- fit_map_mat_ranked %>%
      gather(p2, rank, -driver) %>%
      distinct()
    
    # Timing Plots ####
    sort(genes2plot)
    fit_map_mat_long_tbl <- fit_map_mat %>%
      rownames_to_column() %>%
      ## select rows to plot matching genes found in
      ## known somatic driver genes
      filter(rowname %in% genes2plot) %>%
      gather(p2, prob, -rowname) %>%
      distinct()
  library(ggridges)
  library(scales)
  theme_set(theme_classic())

  bootBTm <- lapply(1:100, function(iter){
    sampI <- sample(1:nrow(wintable_CN_plmstsSNV_deduped_X), round(0.8*nrow(wintable_CN_plmstsSNV_deduped_X)))
    win_table_inputF_X <- mk_win_table(wintable_CN_plmstsSNV_deduped_X[sampI,])
    data4col_X <- codes_to_counts(win_table_inputF_X, c("P1", "P2", "D"))
    data4colX_bs <- data4col_X
    btdata <- btdata(data4colX_bs, return_graph = TRUE)
    #### fit model: btfit ####
    fit_MAP <- btfit(btdata, 1.1)
    BTcoeffs2 <- as.data.table(coef(fit_MAP, as_df = TRUE ))
    T_BTcoeffs2 <- transpose(BTcoeffs2, keep.names = 'rn', make.names = 'item')[-1,]
    return(T_BTcoeffs2)
  })
  bootBTm_dt <- rbindlist(bootBTm, use.names = T, fill = TRUE)

  idC = colnames(bootBTm_dt)[2]
  bootBTm_dtgg <- lapply(colnames(bootBTm_dt)[-1], function(idC){
    idC_X <- colnames(bootBTm_dt) == idC
    bootBTm_dt_X <- bootBTm_dt[,..idC_X]
    bootBTm_dt_X[, item := idC]
    colnames(bootBTm_dt_X)[1] = 'estimate'
    return(bootBTm_dt_X)
  })
  bootBTm_gg <- rbindlist(bootBTm_dtgg)
  bootBTm_gg$estimate <- as.numeric(bootBTm_gg$estimate)
  BTcoeffs3 <-  as.data.table(BTsummary2_itemSummary)[item %in% genes2plot]
  bootBTm_gg_nNA <- bootBTm_gg[!is.na(estimate)]
  bootBTm_gg_nNA_s <-  as.data.table(bootBTm_gg_nNA)[(item %in% genes2plot) & (item %in% plotnames)]
  write.table(bootBTm_gg_nNA_s, paste0('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220316extDataFigs_sourceTbl/20220319extFig7C', x, '.txt'),quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)
  
  pdf(file = paste0('/Volumes/GoogleDrive/My Drive/FD_Beroukhim Lab/pHGG/data/20220316extDataFigs_sourceTbl/20220319extFig7C', x, '.pdf'), width = 16, height = 20)
    plot(ggplot(bootBTm_gg_nNA_s, aes(x=estimate*-1, y=item, fill = estimate*-1)) + 
           geom_violin(fill = 'grey', scale = 'count') + geom_boxplot(width=0.2, fill = 'white', varwidth = TRUE) + 
         labs(title=paste("BT coef",x) , 
              subtitle=paste("early vs late :", length(unique(wintable_CN_plmstsSNV_deduped_X$Sample)), '/', nrow(VarCl4membX), 'Samples'), 
              caption="palimsestSNV, boots 100") +theme_ridges(font_size = 24)) # + theme_classic()
  dev.off()
})
