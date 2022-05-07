require(data.table)
#M2mafs <- readRDS('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200904palimsests/20200903M2mergeUniq.rds')
#colnames(M2mafs)
fls <- list.files('/xchip/beroukhimlab/Frank/DIPG/data/20200906M2vcfs')
fls2 <- gsub(pattern = '.cleaned-filtered.vcf', replacement = '', x = fls)
fls3 <- gsub(pattern = '-filtered.vcf', replacement = '', x = fls2)
fls4 <- gsub(pattern = '_tumor', replacement = '', x = fls3)
fls4 <- gsub(pattern = '.cleaned.rh', replacement = '', x = fls4)

fpths <- paste0('/xchip/beroukhimlab/Frank/DIPG/data/20200906M2vcfs/', fls)
x = 101
lapply(101:length(fls), function(x){
  print(fls4[x])
  vcfX_head <- readLines(fpths[x], n= 150)
  vcfX_head_ed <- c(vcfX_head,'##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="ref read count in tumor">', 
                    '##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="alt read count in tumor">')
  print(vcfX_head[150])
  vcfX <- fread(fpths[x], skip = 'CHROM')
  vcfxPass <- vcfX[FILTER == 'PASS']
  idx <- colnames(vcfxPass) == fls4[x]
  countsX <- vcfxPass[,..idx]
  #y = 2
  countsX_dt <- rbindlist(lapply(1:nrow(countsX), function(y){
    print(paste0(y,'/', nrow(countsX)))
    countsY <- unlist(strsplit(as.character(countsX[y]), "[:]"))[2]
    #print(countsY)
    t_ref_countY <- as.numeric(unlist(strsplit(countsY, "[,]"))[1])
    t_alt_countY <- as.numeric(unlist(strsplit(countsY, "[,]"))[2])
    countsY_dt <- as.data.table(cbind(paste0(vcfxPass$`#CHROM`[y],':',vcfxPass$POS[y]), t_ref_countY, t_alt_countY))
    return(countsY_dt)
  }))
  #countsX_dt$y[duplicated(countsX_dt$V1)]
  vcfxPass_c <- cbind(vcfxPass, countsX_dt) #hm
  #identical(vcfxPass_c$`00-001R`, vcfxPass_c$y)
  vcfxPass_c[,INFO := paste0(INFO,';', 't_ref_count=', t_ref_countY, ';', 't_alt_count=', t_alt_countY), by= INFO]
  #sum(duplicated(vcfxPass_c$INFO))
  writeLines(text = vcfX_head_ed, 
             con = paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201106vcf_reformating/', fls4[x], '_mutTRform_filt.vcf'))
  write.table(x = vcfxPass_c, file = paste0('/xchip/beroukhimlab/Frank/DIPG/data/20201106vcf_reformating/', fls4[x], '_mutTRform_filt.vcf')
              ,quote=FALSE, sep='\t', row.names=FALSE, col.names = TRUE, append = TRUE)
  return()
})
  
