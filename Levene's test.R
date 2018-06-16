# CQR in UKBB for N=149 CREAM SNPs replicating in UKB

#dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
#
#install.packages("/home/sopjg2/R/cqr/SparseM_1.77.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/Matrix_1.2-12.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/MatrixModels_0.4-1.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/quantreg_5.35.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/iterators_1.0.9.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/foreach_1.4.4.tar.gz", Sys.getenv("R_LIBS_USER") )
#install.packages("/home/sopjg2/R/cqr/doParallel_1.0.11.tar.gz", Sys.getenv("R_LIBS_USER") )

rm(list=ls())

library(quantreg)
library(lawstat)
library(metafor)

#infolder  <- "/home/sopjg2/R/cqr/"
#outfolder <- "/home/sopjg2/R/cqr/"

infolder  <- "C:/Users/lopasy/Documents/"
outfolder <- "C:/Users/lopasy/Documents/"

shuffle <- function(x=x){
  x        <- as.numeric(x)
  ninds    <- length(x)
  d        <- as.data.frame(matrix(ncol=2,nrow=ninds))
  names(d) <- c("x","rands")
  d$x      <- x
  d$rands  <- runif(ninds,min=0,max=1)
  d        <- d[order(d$rands),]
  y        <- d$x
  return(y)
}

get_rntransform <- function(myvar){
  x         <- as.numeric(myvar)
  y         <- rnorm(length(x),mean=0,sd=1)
  z         <- as.data.frame(matrix(ncol=1,nrow=length(x)))
  z$x       <- x
  z$ord     <- 1:length(x)
  zz        <- z[order(z$x),]
  zz$y      <- sort(y)
  zzz       <- zz[order(zz$ord),]
  return(zzz$y)
}

#data1 <- read.table(file=paste(infolder,"CQR.raw",sep=""),header=TRUE)
data1 <- read.table(file="C:/Users/lopasy/Documents/height.raw",header=TRUE)

data1$FID <- NULL
data1$PAT <- NULL
data1$MAT <- NULL
data1$SEX <- NULL
data1$PHENOTYPE <- NULL

#data2 <- read.table(file=paste(infolder,"ukb_amblyopia_post-gwas_2017-10-26.phens",sep=""), header=TRUE)
data2 <- read.table(file="C:/Users/lopasy/Documents/height_pheno.txt", header=TRUE)
data2 <- data2[,c(2,3,14,16,17)]
data4 <- merge(data2,data1,by="IID")

data1 <- NULL
data2 <- NULL

qs    <- (1:19/20)
qss   <- qs^2

# need to rename 1 SNP
names(data4)[names(data4)=="rs9680365.lol_C"] <- "rs9680365_C"

# need to exclude 1 SNP
data4$rs9680365_C <- NULL

var_names   <- names(data4)
vars        <- as.matrix(var_names)
num_vars    <- nrow(vars)
num_snps    <- 0
snp_pos     <- matrix(nrow=num_vars,ncol=1)

for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  }
}



# LM: Loop through all SNPs
# ------------------------------------

results              <- as.data.frame(matrix(nrow=num_snps,ncol=8))
names(results)       <- c("Marker","Beta_OLS","LCI_OLS","UCI_OLS","P_OLS","MAF","P_LEVENE_OBS","P_LEVENE_PERM")

for (snp in 1:num_snps){
  
  results[snp,1]   <- vars[snp_pos[snp]]
  x                <- sum(data4[,vars[snp_pos[snp]]],na.rm=TRUE) / (2*sum(!is.na(data4[,vars[snp_pos[snp]]])))
  results[snp,6]   <- formatC(min(x, (1-x)),digits=2,format="f")
  
  m               <- min(table(data4[,vars[snp_pos[snp]]]))
  if(m > 50){
    
    myform          <- as.formula(paste("Height ~ ", vars[snp_pos[snp]], " + Sex + poly(Age,2) + Array", sep=""))
    
    mod_sumM <- summary(lm(formula=myform, data=data4))
    b        <- mod_sumM$coefficients[2,1]
    s        <- mod_sumM$coefficients[2,2]
    lci      <- b - (1.96*s)
    uci      <- b + (1.96*s)
    
    results[snp,2] <- formatC(b,digits=3,format="f")
    results[snp,3] <- formatC(lci,digits=3,format="f")
    results[snp,4] <- formatC(uci,digits=3,format="f")
    results[snp,5] <- formatC(mod_sumM$coefficients[2,4],digits=2,format="e")
  } 
}
results

# Levene's part 1: Loop through all SNPs to get observed p-value
# --------------------------------------------------------------

for (snp in 1:num_snps){
  m               <- min(table(data4[,vars[snp_pos[snp]]]))
  if(m > 50){
    z               <- data4[which(data4[,vars[snp_pos[snp]]]!="NA"),]
    z = z[!is.na(z$Height),]
    mod_lev         <- levene.test(z$Height, group=as.factor(z[,vars[snp_pos[snp]]]),location="mean")
    results[snp,7]  <- formatC(mod_lev$p.value,digits=3,format="e")
  } 
}
results

# Levene's part 2: Loop through all SNPs to get shuffled-phenotype p-values
# -------------------------------------------------------------------------
nperms         <- 100
myperms        <- as.data.frame(matrix(nrow=nperms,ncol=num_snps))
names(myperms) <- results[,1]

for (snp in 1:num_snps){
  
  cat("snp #",snp,":",vars[snp_pos[snp]])
  cat("\n")
  flush.console()
  
  m                 <- min(table(data4[,vars[snp_pos[snp]]]))
  if(m > 50){
    
    psig              <- 0
    pobs              <- as.numeric(results[snp,7])
    for (perm in 1:nperms){
      data4 = data4[!is.na(data4$Height),]
      data4$shuff_phen  <- shuffle(data4$Height)
      z                 <- data4[which(data4[,vars[snp_pos[snp]]]!="NA"),]
      mod_lev           <- levene.test(z$shuff_phen, group=as.factor(z[,vars[snp_pos[snp]]]),location="mean")
      myperms[perm,snp] <- mod_lev$p.value
      if(mod_lev$p.value < pobs){ psig = psig +1 }
    } # next perm
    
    results[snp,8]    <- formatC(as.numeric(psig/nperms),digits=3,format="e")
    
  } #end if
} #next snp

results$P_LEVENE_OBS_log = -log10(as.numeric(results$P_LEVENE_OBS))
results$P_LEVENE_PERM_log = -log10(as.numeric(results$P_LEVENE_PERM))
results = results[-which(results$P_LEVENE_PERM_log == "Inf"),]
m1 <- summary(lm(P_LEVENE_PERM_log ~ P_LEVENE_OBS_log,data=results)); m1; m1$coef[2,1]
results$P_LEVENE_OBS_corrected = 10^(m1$coef[2,1]*log10(as.numeric(results$P_LEVENE_OBS)))

d1 <- make_data1(results$P_LEVENE_p)
QQ1 <- ggplot()+ geom_ribbon(data=d1, fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(data=d1, aes(x=x1, y=x1), colour="red",size=1)+
  geom_point(data=d1, colour="black", size=2.5, aes(x=x1,y=y1))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN)
QQ1

file1=paste(outfolder,"cqr_2018-04-19_levene_",nperms,"perms.csv",sep=""); write.csv(results,          file=file1, row.names=TRUE)
file1=paste(outfolder,"cqr_2018-04-19_levene_perms.csv",sep=""); write.csv(myperms,          file=file1, row.names=TRUE)



#####

# QQ plot for Levene's test with 10k permutations

library(ggplot2)
library(grid)
library(gridExtra)


find_conf_intervals = function(row){
  i = row[1]
  len = row[2]
  if (i < 10000){
    return(c(-log10(qbeta(0.975,i,len-i+1)), -log10(qbeta(0.025,i,len-i+1))))
  } else { # Speed up
    return(c(NA,NA))
  }
}

theme_qq <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1.6), margin = margin(t = 15)),
      axis.title.y = element_text(colour = "black", size=rel(1.6), angle=90, margin = margin(r = 15)),
      axis.text.x = element_text(colour = "black", size=rel(1.4)),
      axis.text.y = element_text(colour = "black", size=rel(1.4)),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill="white"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    )   
}
theme_set(theme_qq())

#            plot.background = element_rect(fill="white"),

make_data1 <- function(data1){
  obs1 <- data1[order(data1,decreasing=F)]
  y1   <- -log10(obs1)
  x1   <- -log10( ppoints(length(obs1) ))
  d1   <- as.data.frame(cbind(x1,y1))
  N    <- dim(d1)[1]
  ci   <- 0.95
  d1$lci <- -log10(qbeta((1 - ci)/2, 1:N, N:1))
  d1$uci <- -log10(qbeta((1 + ci)/2, 1:N, N:1))
  return(d1)
}

scaleFUN <- function(x) sprintf("%.1f", x)

dataQQ=read.csv(file=paste(outfolder,"cqr_2018-04-19_levene_200perms.csv",sep=""),header=TRUE)


d1 <- make_data1(dataQQ$P_LEVENE_OBS)
QQ1 <- ggplot()+ geom_ribbon(data=d1, fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(data=d1, aes(x=x1, y=x1), colour="red",size=1)+
  geom_point(data=d1, colour="black", size=2.5, aes(x=x1,y=y1))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN)
QQ1

m1 <- summary(lm(y1~x1,data=d1)); m1; m1$coef[2,1]


gxescan = as.data.frame(matrix(nrow = ncol(myperms), ncol = 4))
names(gxescan) = c("SNP", "P_Levene_Observed", "P_PERM_JEZ", "P_Levene_PERM")
gxescan[, 1] = results$Marker
gxescan[, 2] = results$P_LEVENE_OBS
gxescan[, 3] = results$P_LEVENE_PERM


for(i in 1:ncol(myperms)){
  gxescan[i,4] = (length(which(as.numeric(myperms[,i]) < as.numeric(gxescan[i,2]))))/200
}
gxescan = na.omit(gxescan)
gxescan$P_Levene_PERM = formatC(gxescan$P_Levene_PERM, digits = 3, format = "e")

file1=paste(outfolder,"cqr_2018-04-19_levene_",nperms,"_gxescan.csv",sep=""); write.csv(gxescan,          file=file1, row.names=TRUE)
dataQQ=read.csv(file=paste(outfolder,"cqr_2018-04-19_levene_200_gxescan.csv",sep=""),header=TRUE)

d1 <- make_data1(dataQQ$P_Levene_Obs)
QQ1 <- ggplot()+ geom_ribbon(data=d1, fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(data=d1, aes(x=x1, y=x1), colour="red",size=1)+
  geom_point(data=d1, colour="black", size=2.5, aes(x=x1,y=y1))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN)
QQ1

m1 <- summary(lm(y1~x1,data=d1)); m1; m1$coef[2,1]







