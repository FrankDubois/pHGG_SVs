library(ggplot2)
library(reshape)
library(gplots)
library(reshape2)
library(gridExtra)

exposure_mat <- readRDS('/Volumes/xchip_beroukhimlab/Alex/pHGG/20220207_sbs_sigs/20220207_extracted_sbssigs_weightsmat.rds')
scale <- 0.8
.theme_ss <- theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12*scale, family="mono"),
        axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
        axis.title.x = element_text(face="bold",colour="black",size=14*scale),
        axis.title.y = element_text(face="bold",colour="black",size=14*scale),
        axis.text = element_text(size = 16*scale, family = "mono"),
        strip.text = element_text(lineheight=0.5),
        strip.text.x = element_text(size=10*scale,face='bold',angle=00),
        strip.text.y = element_text(size=10*scale,face="bold"),
        strip.background = element_rect(colour="black",fill="gray85"),
        panel.margin = unit(0.20,"lines"),
        plot.title=element_text(lineheight=1.0,face="bold",size=12*scale))
color.values = c("red","royalblue4","orange","cyan","darkgreen","black","magenta","darkviolet","gray50","green","darkred","darkblue","gray85","gold",rainbow(10)[5:20])
lego.colors <- c("cyan","black","red","gray50","green","magenta")
lego.colors <- c("cyan","red","yellow","purple","green","blue")

plot.signature.DNP <- function(W,title) {
  df1 <- data.frame(W)
  df1[df1 < 1.e-10] <- 0
  df1[,"feature"] <- rownames(W)
  df1 <- melt(df1,id.var="feature")
  colnames(df1) <- c("feature","signature","activity")
  x1 <- sapply(df1$feature,function(x) strsplit(x,"/")[[1]][1])
  x2 <- sapply(df1$feature,function(x) strsplit(x,"/")[[1]][2])
  x3 <- sapply(df1$feature,function(x) strsplit(x,"/")[[1]][3])
  x4 <- sapply(df1$feature,function(x) strsplit(x,"/")[[1]][4])
  x5 <- sapply(df1$feature,function(x) strsplit(x,"/")[[1]][5])
  df1[,"type1"] <- paste(x1,x2,sep="_")
  df1[,"type2"] <- paste(x4,x5,sep="_")
  #df1$feature <- factor(df1$feature,levels=rownames(W))
  K <- ncol(W)
  p = ggplot(df1)
  p = p + geom_bar(aes_string(x="feature",y="activity", fill="type1" ),stat="identity",position="identity")
  p = p + facet_grid(signature ~ ., scale = "free_y")
  p = p + .theme_ss
  p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
  p = p + ggtitle(title)
  p = p + xlab("Features") + ylab("Contributions")
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=10*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=10*scale))
  p = p + theme(strip.text.x = element_text(size = 10*scale, colour = "black", angle = 0))
  p = p + theme(strip.text.y = element_text(size = 10*scale, colour = "black", angle = 270))
  return(p)
}

plot.activity.barplot <- function(H.mid,H.norm,scale,tumor.type) {
  .theme_ss <- theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
          axis.text = element_text(size = 12*scale, family = "mono"))
  ordering <- order(colSums(H.mid),decreasing=T)
  H.mid <- H.mid[,ordering]
  #rownames(H.mid) <- paste("W",seq(1:nrow(H.mid)),sep="")
  H.norm <- H.norm[,ordering]
  #rownames(H.norm) <- paste("W",seq(1:nrow(H.norm)),sep="")
  sample.ordering <- colnames(H.mid)
  x1 <- melt(H.mid)
  x2 <- melt(H.norm)
  colnames(x1) <- c("Signature","Sample","Activity")
  colnames(x2) <- c("Signature","Sample","Activity")
  x1[,"class0"] <- c("Counts")
  x2[,"class0"] <- c("Fractions")
  df2 <- rbind(x1,x2)
  df2$class0 <- factor(df2$class0,c("Counts","Fractions"))
  df2$Sample <- factor(df2$Sample,sample.ordering)
  scale <- 1
  p = ggplot(df2,aes(x=factor(Sample),y=Activity,fill=factor(Signature)))
  p = p+geom_bar(stat="identity",position='stack',color='black',alpha=0.9)
  p = p + scale_fill_manual(values=c("red","cyan","yellow","blue","magenta","gray50","orange","darkgreen","brown","black","#CC79A7", "#0066CC", "#330013","#99FF00",rainbow(10)[1:10]))
  p = p + facet_grid(class0 ~ ., scale = "free_y")
  p = p + ggtitle(paste("Siganture Activities in",tumor.type,sep=" "))
  p = p + theme(plot.title=element_text(lineheight=1.0,face="bold",size=14*scale))
  p = p + xlab("Samples") + ylab("Signature Activities")
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.text.x = element_text(angle=90,vjust=0.5,size=8*scale,face="bold",colour="black"))
  p = p + theme(axis.text.y = element_text(size=10*scale,face="bold",colour="black"))
  p = p + theme(legend.title=element_blank())
  p = p + .theme_ss
  p = p + theme(legend.position="top")
  return(p)
}

p <- plot.signature.DNP(exposure_mat,'pHGG')
plot(p)
exposure_mat_t <- t(exposure_mat)

p1 <- plot.activity.barplot(exposure_mat_t,exposure_mat_t,1.0,'pHGG')
plot(p1)

