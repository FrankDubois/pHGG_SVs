# This script plots the fishHook outputs 

require(data.table)
require(gTrack)

OVERLAP <- 0
BINSIZE  = 1e7
hypothesesTiles <- gr.tile(si2gr(gUtils::si), BINSIZE) + OVERLAP

fish50k_cov6 <- fread('/Volumes/xchip_beroukhimlab/Frank/DIPG/data/20200819fhook/20200330RSBs50k_cov6_allSV.txt')
colnames(fish50k_cov6)
fish50k_cov6[,logFDR := -log10(fdr)]
fish50k_cov6[,log_p := -log10(p)]
table(fish50k_cov6$seqnames)
fish50k_cov6_dt <- fish50k_cov6_dt[seqnames %in% c(1:22,'X','Y')]
fish50k_cov6_gr <- dt2gr(fish50k_cov6)

### This line orders the chromosomes correctly
fish50k_cov6_gr <- keepSeqlevels(fish50k_cov6_gr, as.character(c(1:22, 'X', 'Y')))

#gt_fish50k_cov62Pad= gTrack((fish50k_cov6_gr+1e4) , y.field = 'logFDR', col = 'blue', name = '-logFDR', lines = TRUE, y1 = 3, y0 = 0, 
#                            draw.paths =FALSE, sep.draw = FALSE, sep.bg.col = 'white')#, border = '6', y1 =  10)  #circles = TRUE,)
#gt_fish50k_cov62_cPad= gTrack((fish50k_cov6_gr +1e6) %Q% (!is.na(count)) , y.field = 'count', col = 'blue', name = 'SVcount', 
#                              lines = TRUE, y1 = 9,draw.paths =FALSE, sep.draw = FALSE, sep.bg.col = 'white')#, border = '6')  #circles = TRUE,)

### Create the gTrack object
gt_fish50k_cov6Pad= gTrack((fish50k_cov6_gr+1e4) , y.field = 'logFDR', col = 'blue', name = '-logFDR', bars = TRUE, y1 = 3, y0 = 0, 
                            draw.paths =FALSE, sep.draw = FALSE, sep.bg.col = 'white', sep.lwd = 0,  stack.gap = 0, yaxis.pretty = 3) #circles = TRUE,)

### Commented out the code below as not plotting the SV count graph currently
# gt_fish50k_cov6_cPad= gTrack((fish50k_cov6_gr +1e4) %Q% (!is.na(count)) , y.field = 'count', col = 'blue', name = 'SVcount'
#                               , bars = TRUE, y1 = 8.5,draw.paths =FALSE, sep.draw = FALSE, sep.bg.col = 'white',sep.lwd = 0, 
#                               stack.gap = 0, yaxis.pretty = 4)#, y.grid.lwd =0)

### only plot the FDR graph, not the SV count one
# plot(c(gt_fish50k_cov6_cPad, gt_fish50k_cov6Pad), plot.window = hypothesesTiles, sep.draw = FALSE, stack.gap = 0)#,y.grid.lwd =0, m.sep.lwd = FALSE)#,sep.lwd = 0)

### Only plotting the FDR graph 
plot(gt_fish50k_cov6Pad, plot.window = hypothesesTiles, sep.draw = FALSE, stack.gap = 0)#,y.grid.lwd =0, m.sep.lwd = FALSE)#,sep.lwd = 0)


nrow(as.data.table(hypothesesTiles))
# tilmA_mat <- matrix(data = 1:nrow(as.data.table(hypothesesTiles)), nrow = nrow(as.data.table(hypothesesTiles)), ncol = nrow(as.data.table(hypothesesTiles)))
# 
# plot(c(gt_fish50k_cov6_cPad, gt_fish50k_cov6Pad, gTrack(hypothesesTiles, mdata = tilmA_mat, cmap.max = 10, triangle = TRUE,
#                                                           sep.lwd = 0, sep.draw = FALSE,border = '0',
#                                                           colormaps = colorRampPalette(c("blue",'yellow', "red"))( 10))), 
#      plot.window = hypothesesTiles, sep.draw = FALSE, stack.gap = 0)#,y.grid.lwd =0, m.sep.lwd = FALSE)#,sep.lwd = 0)
