require(data.table)
fls <- list.files('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res')
fls_CNA_dts <- grep(pattern = '_CN_dt.rds', x = fls, value = TRUE, fixed = TRUE)
flsCNAs2 <- gsub(pattern = '_mutTimeR_CN_dt.rds', replacement = '', x = fls_CNA_dts)
x= 1
fpths <- paste0('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20201109mutationTimeR/20201109mutationTimeR_K2N3pi02_res/', fls_CNA_dts)

CNA_dts_Lsub <- lapply(1:length(fls_CNA_dts), function(x){
  print(flsCNAs2[x])
  CNA_dtX <- readRDS(fpths[x])
  CNA_dtX_fsub <- CNA_dtX[!is.na(type)]
  return(CNA_dtX_fsub)
})#, mc.cores = 6, mc.preschedule = FALSE)
CNA_dtsub <- rbindlist(CNA_dts_Lsub)
unique(CNA_dtsub$type)
table(CNA_dtsub$time.star)
mean(CNA_dtsub[time.star== '***', time])
mean(CNA_dtsub[time.star== '**', time])
mean(CNA_dtsub[time.star== '*', time])
hist(CNA_dtsub$time)
hist(CNA_dtsub$time, breaks = 'FD')
hist(CNA_dtsub$time, breaks = 'Scott')

library(gUtils)

CNA_dtsub_gr <- dt2gr(CNA_dtsub)
## arms####
armCoords <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210113canopy/20210113_sigarmDCoords.txt')
colnames(armCoords)[1] <- 'arm'
armCoords_gr <- dt2gr(armCoords)
CNAsub_arms__gr_ann <- CNA_dtsub_gr %$% armCoords_gr
CNAsub_arms_gr_anndt <- as.data.table(CNAsub_arms__gr_ann)
CNAsub_arms_gr_anndt[,eventName := paste0(type, ':', arm)]
CNAsub_arms_gr_anndt_sub <- CNAsub_arms_gr_anndt[arm !='']
CNAsub_arms_gr_anndt_sub[,meanTime := mean(na.omit(time)), by = c('Sample', 'eventName')]
CNAsub_arms_gr_anndt_sub$meanTime
CNAsub_arms_gr_anndt_subD <- CNAsub_arms_gr_anndt_sub[!is.na(arm)]
sort(unique(CNAsub_arms_gr_anndt_subD$eventName))
saveRDS(CNAsub_arms_gr_anndt_subD, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210114winningtablesCNAarms/20210114CNAwin_annarms.rds')
write.table(x = CNAsub_arms_gr_anndt_subD, 
            file = '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210114winningtablesCNAarms/20210114CNAwin_annarms.seg'
            ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)
CNAsub_arms_gr_anndt_subT <- CNAsub_arms_gr_anndt_subD[eventName %in% c('Bi-allelic Gain (WGD):1p', 'Bi-allelic Gain (WGD):1q', 'Mono-allelic Gain:1p',
                                                                        'Mono-allelic Gain:1q', 'Bi-allelic Gain (WGD):2p', 'Bi-allelic Gain (WGD):2q', 
                                                                        'Mono-allelic Gain:2q', 'Mono-allelic Gain:2p', 'Bi-allelic Gain (WGD):7p', 
                                                                        'Bi-allelic Gain (WGD):7q', 'Mono-allelic Gain:7p', 'Mono-allelic Gain:7q', 
                                                                        'CN-LOH:10p', 'CN-LOH:10q', 'CN-LOH:11p', 'CN-LOH:13q','CN-LOH:14q', 'CN-LOH:16p',
                                                                        'CN-LOH:16p, 16q', 'CN-LOH:17p', 'CN-LOH:18p', 'CN-LOH:18q')]


## GISTIC peaks####
GISTIC_gr <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210117canopy/20210117GISTIC_gr.RDS')
CNAsub_gr_ann <- CNA_dtsub_gr %$% GISTIC_gr
CNAsubgr_anndt <- as.data.table(CNAsub_gr_ann)
CNAsubgr_anndt[,eventName := paste0(type, ':', Descriptor)]
CNAsubgr_anndt_sub <- CNAsubgr_anndt[Descriptor !='']
CNAsubgr_anndt_sub[,meanTime := mean(na.omit(time)), by = c('Sample', 'eventName')]
CNAsubgr_anndt_sub$meanTime
CNAsubgr_anndt_subD <- CNAsubgr_anndt_sub[!is.na(Descriptor)]
sort(unique(CNAsubgr_anndt_subD$eventName))
saveRDS(CNAsubgr_anndt_subD, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20201110muttimeR/20201110CNAwin_ann.rds')
write.table(x = CNAsubgr_anndt_subD, 
            file = '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20201110muttimeR/20201110CNAwin_ann.seg'
            ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)


unique(grep(pattern = 'TP53', x = CNAsubgr_anndt_subD$eventName, value = TRUE))
CNAsubgr_anndt_subT <- CNAsubgr_anndt_subD[eventName %in% c('Mono-allelic Gain:4q12_Amp(PDGFRA)', 'Mono-allelic Gain:7q21.2_Amp(CDK6)', 
                                                       'Mono-allelic Gain:2p24.3_Amp(MYCN)', 'Mono-allelic Gain:2p25.1_Amp',
                                                       'Mono-allelic Gain:2p25.1_Amp(KIDINS220)', 
                                                       'Mono-allelic Gain:7p15.1_Amp', 'Mono-allelic Gain:7p11.2_Amp(EGFR)', 
                                                       'Mono-allelic Gain:7q31.2_Amp(MET)', 'Mono-allelic Gain:8q24.21_Amp(CCDC26)', 
                                                       'Mono-allelic Gain:2p25.1_Amp, 2p25.1_Amp(KIDINS220), 2p24.3_Amp(MYCN)',
                                                       'Mono-allelic Gain:7p15.1_Amp, 7p11.2_Amp(EGFR)', 'Mono-allelic Gain:8q24.21_Amp(MYC)', 
                                                       'CN-LOH:9p21.3(CDKN2A)_Del', 'CN-LOH:13q14.2(RB1)_Del', 'CN-LOH:17p13.3(TP53)_Del',
                                                       'CN-LOH:10q23.31(PTEN)_Del', 'CN-LOH:10q26.3(MGMT)_Del', 'CN-LOH:17q11.2(NF1)_Del',
                                                       'Mono-allelic Gain:8q24.21_Amp(MYC), 8q24.21_Amp(CCDC26), 8q24.3_Del', 
                                                       'Mono-allelic Gain:17q11.2(NF1)_Del, 17q25.3_Del'
                                                       
)]

CNAsubgr_anndt_subTcomb <- rbind(CNAsubgr_anndt_subT, CNAsub_arms_gr_anndt_subT,use.names=FALSE)
#CNAsubgr_anndt_subT[,timing_state := NA]
CNAsubgr_anndt_subTcomb[ meanTime < 0.25,timing_state := 'early_clonal']
CNAsubgr_anndt_subTcomb[meanTime >= 0.25,timing_state := 'clonal_NA']
CNAsubgr_anndt_subTcomb[ meanTime >= 0.5,timing_state := 'late_clonal']
CNAsubgr_anndt_subTcomb[ meanTime >= 0.75 ,timing_state := 'subclonal']
table(CNAsubgr_anndt_subTcomb$timing_state)
CNAsubgr_anndt_subTcombs <- CNAsubgr_anndt_subTcomb[!is.na(timing_state)]
saveRDS(CNAsubgr_anndt_subTcombs, '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210118_SCNAwinningtables/20210118CNAwin_ann_sub.rds')
write.table(x = CNAsubgr_anndt_subTcombs, 
            file = '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20210118_SCNAwinningtables/20210118CNAwin_ann_sub.seg'
            ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE)



unique(CNAsubgr_anndt_subTcomb$eventName)
unique(grep(pattern = 'TP53', x = CNAsubgr_anndt_subD$eventName, value = TRUE))

CNAsubgr_anndt_subTcomb_sub <- CNAsubgr_anndt_subTcomb[eventName %in% c('Mono-allelic Gain:4q12_Amp(PDGFRA)', 'Mono-allelic Gain:2p24.3_Amp(MYCN)', 
                                                                'Mono-allelic Gain:7p11.2_Amp(EGFR)', 'Mono-allelic Gain:8q24.21_Amp(CCDC26)', 
                                                                'Mono-allelic Gain:8q24.21_Amp(MYC)','CN-LOH:17p13.3(TP53)_Del'
)]
library(ggplot2)
g <- ggplot(CNAsubgr_anndt_subTcomb_sub, aes(x=Descriptor, y = meanTime))
g + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, fill="red", binwidth = 0.05) + 
  labs(title="GISTICpeak timing", subtitle="mutTimeR")+
  theme(text = element_text(size=20))




