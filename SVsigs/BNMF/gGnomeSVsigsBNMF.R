### This script performs BNMF and generates plots
library(reshape)
library(gplots)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(parallel)

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

# as.matrix(lego.SV),200000,10,tol,Kcol,Kcol,1/phi
# V0 = lego.SV
# n.iter = 200000
# K = Kcol
# K0 = K
BayesNMF.L1.KL.phi.new <- function(V0,n.iter,a0,tol,K,K0,phi) {
  eps <- 1.e-50
  del <- 1.0
  #active_nodes <- colSums(V0) != 0
  #V0 <- V0[,active_nodes]
  V <- V0-min(V0) + eps
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  #W <- matrix(runif(N * K)*sqrt(Vmax),ncol=K)
  #H <- matrix(runif(M * K)*sqrt(Vmax),ncol=M)
  W <- matrix(runif(N * K)*sqrt(mean(V)),ncol=K)
  H <- matrix(runif(M * K)*sqrt(mean(V)),ncol=M)
  V.ap <- W %*% H + eps
  I <- array(1,dim=c(N,M))
  
  C <- N + M + a0 + 1
  b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K)
  #b0 <- sqrt(mean(V,na.rm=T)/K)
  lambda.bound <- b0/C
  lambda <- (colSums(W) + rowSums(H) + b0)/C
  lambda.cut <- 1.1 * lambda.bound
  
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+phi/lambda,M),ncol=M) + eps)
    V.ap <- W %*% H + eps
    W <- W * ((V/V.ap) %*% t(H))/t(matrix(rep(rowSums(H)+phi/lambda,N),ncol=N) + eps)
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
    n.lambda[[iter]] <- lambda
    if (iter %% 100 == 0) {
      like <- sum(V*log(V/V.ap)+V.ap-V)
     cat('like, phi, b0, C','\n', 
         like,phi, b0, C ,'\n')## db
      evid <- like/phi+sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
      if (exists('evid')) {
        error <- sum((V-V.ap)^2)
        #cat(iter,evid,like,error,del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
        cat('iter,evid,like,error,del,sum(colSums(W)!=0)','\n', iter,evid,like,error,del,sum(colSums(W)!=0),'\n')
        res <- list(W,H,like,evid,lambda,error) 
        save(res,file=paste(OUTPUT,paste("Bayes.temp.RData",sep="."),sep=""))
      }
      else{
        cat('no_evid')
      }
        
    }
    iter <- iter+1
  }
  cat("***********************************************************",'\n')
  if (exists('evid')){
  return(list(W,H,like,evid,lambda,error))
  }
  else{
    print('no_evid')
    return(0)
  }
}

######### Visualizing the activity of signatures
plot.activity.barplot <- function(H.mid,H.norm,scale,tumor.type) {
  .theme_ss <- theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
          axis.text = element_text(size = 12*scale, family = "mono"))
  ordering <- order(colSums(H.mid),decreasing=T)
  H.mid <- H.mid[,ordering]
  rownames(H.mid) <- paste("W",seq(1:nrow(H.mid)),sep="")
  H.norm <- H.norm[,ordering]
  rownames(H.norm) <- paste("W",seq(1:nrow(H.norm)),sep="")
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
  p = p + scale_fill_manual(values=c("red","cyan","yellow","blue","magenta","gray50","orange","darkgreen","brown","black",rainbow(10)[1:10]))
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


OUTPUT <- "/xchip/beroukhimlab/Frank/DIPG/data/20220203gGnomeSVsigs/out/"
system(paste("mkdir",OUTPUT,sep=" "))
tumor.type <- "pHGG"

require(data.table)
require(gUtils)
#basedir <- "/xchip/beroukhimlab/ofer/"


# run signatures analysis
reduce.mat <- readRDS('/xchip/beroukhimlab/Frank/DIPG/data/20220203gGnomeSVsigs/20220203gGnomeSVsigs_nmfInputMat.RDS')

# run signatures analysis
lego.SV <- t(reduce.mat)
Kcol <- dim(reduce.mat)[1]
tol <- 1.e-07
a0 <- 10
n.iter <- 20

phi <- 1
method <- paste("L1KL.SV",tumor.type,"phi",phi,sep=".")
mclapply(1:n.iter, function(xij){
# lapply(1:n.iter, function(xij){
  res <- BayesNMF.L1.KL.phi.new(as.matrix(lego.SV),200000,10,tol,Kcol,Kcol,1/phi)
  if (class(res) != 'list') {
    print('no_evid')
  } 
  else{
  save(res,file=paste(OUTPUT,paste(method,xij,"RData",sep="."),sep=""))
  W <- res[[1]]
  H <- res[[2]]
  index <- colSums(W) > 1
  K <- sum(index)
  W <- W[,index]
  H <- H[index,]
  colnames(W) <- paste("W",seq(1:ncol(W)),sep="")
  rownames(H) <- colnames(W)
  W1 <- W
  H1 <- H
  for (i in 1:K) {
    W1[,i] <- W[,i]*rowSums(H)[i]
    H1[i,] <- H[i,]*colSums(W)[i]
  }
  pdf(file=paste(OUTPUT,paste("signature",method,K,abs(round(res[[4]])),xij,"pdf",sep="."),sep=""),width=(15),height=(K*2))
  p <- plot.signature.DNP(W1,tumor.type)
  plot(p)
  dev.off()
  
  #mat.out <- list(W1,H1)
  #save(mat.out,file=paste(OUTPUT,paste(method,j,"RData",sep="."),sep=""))
  H.mid <- H1
  H.norm <- apply(H.mid,2,function(x) x/sum(x))
  p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
  pdf(file = paste0(OUTPUT,"activity.barplot.", xij, ".pdf"),width=15,height=12)
  plot(p1)
  dev.off()
  }
}, mc.cores = 6)


############# Activity plot
#load(paste0(OUTPUT, "L1KL.SV.DIPG.phi.1.7.RData"))
#H.mid <- res[[2]]
#H.norm <- apply(H.mid,2,function(x) x/sum(x))

#p1 <- plot.activity.barplot(H.mid,H.norm,1.0,"DIPG")
#pdf(file = paste0(OUTPUT,"activity.barplot.8.pdf"),width=15,height=12)
#plot(p1)
#dev.off()
