##########################################
#                                        #
#            UKBB - PGRS x UNI           #
#                                        #
##########################################

rm(list=ls())

library(quantreg)
library(metafor)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)


get_ztransform <- function(myvar){
  x         <- as.numeric(myvar)
  m         <- mean(x,na.rm=TRUE)
  s         <- sd(x,na.rm=TRUE)
  y         <- (x - m)/s
  return(y)
}

get_score <- function(raw_file, score_file, skip=skip){ # raw_file=plink format file, score_file contains rsID and weight, skip=number of columns to skip in raw file before 1st SNP
  miss             <- as.integer(skip)
  raw              <- as.data.frame(raw_file)
  weights          <- as.data.frame(score_file)
  scores           <- as.data.frame(matrix(nrow=dim(raw)[1],ncol=(miss + dim(weights)[1])))
  nsnps            <- dim(weights)[1]
  scores[,1:miss]  <- raw[,1:miss]
  names(scores)[1:miss] <- names(raw)[1:miss]
  for (snp in 1:nsnps){
    raw_name         <- names(raw)[snp+miss]
    weight_name      <- weights[snp,1]
    if(raw_name==weight_name){
      weight                  <- weights[snp,2]
      scores[,snp+miss]       <- raw[,snp+miss] #*weight # for allele score, do not use weights
      names(scores)[snp+miss] <- raw_name
    } #end if
  } # next snp
  scores$score         <- rowSums(scores[,(miss+1):(miss+nsnps)],na.rm=TRUE)
  scores$nscore        <- get_ztransform(scores$score)
  return(scores$nscore)
}

theme_cqr <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1), margin = margin(t = 6)),
      axis.title.y = element_text(colour = "black", size=rel(1), angle=90, margin = margin(r = 6)),
      axis.text.x = element_text(colour = "black", size=rel(0.9)),
      axis.text.y = element_text(colour = "black", size=rel(0.9)),
      panel.grid.major.x = element_line(colour="#CCCCCC",linetype="solid",size=0.1), 
      panel.grid.minor.x = element_line(colour="#CCCCCC",linetype="dashed",size=0.1), 
      panel.grid.major.y = element_line(colour="#CCCCCC",linetype="solid",size=0.1),
      panel.grid.minor.y = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill="white"),
      legend.text = element_text(size=rel(0.7)),
      legend.title = element_text(size=rel(0.8)),
      legend.position=c(0.2,0.7),
      legend.key=element_blank(),
      plot.title = element_text(colour = "black", size=rel(1), margin = margin(b = 10),hjust = 0.5),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
    )   
}
theme_set(theme_cqr())

data1 <- read.csv("~/Documents/CQR_17_04_2018/AoO/AoO.raw", sep="")
withdrawn <- read.delim("~/w17351_2018-05-03.remove", header=FALSE)
Euro <- read.table("~/Documents/CQR_17_04_2018/AoO/ukb_v2_10xSD_PC-derived_Europeans.keep", quote="\"", comment.char="")
data2 <- read.csv("~/Documents/CQR_17_04_2018/AoO/cnv_all_details_included.txt", sep="")
data2 = data2[which(!data2$FID %in% withdrawn$V1),]
data2 = data2[which(Euro$V1 %in% data2$FID),]
data2 = data2[which(data2$het_miss_out == 0),]
data2 = data2[which(data2$in_kinship == 0),]
data2 = data2[!is.na(data2$AgeSpexWear),]


data1$FID <- NULL
data1$PAT <- NULL
data1$MAT <- NULL
data1$SEX <- NULL
data1$PHENOTYPE <- NULL

data2 <- read.table(file="D:/cygwin64/home/raven/UKB/qc_full/ukb_geno_mse2017-02-16.txt", header=TRUE)
data3 <- data2[,c("IID","UniEdu")]
data4 <- merge(data3,data1,by="IID")
data1 <- data4

#data2 <- read.table(file="D:/cygwin64/home/raven/UKB/qc_full/ukb_amblyopia_post-gwas_2017-10-26.phens", header=TRUE)

data3 <- data2[,c("IID","Age","Sex_matched","avMSE","Geno_array","Heterozygosity",
                  "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]

data4 <- merge(data2,data1,by="IID")

data1 <- NULL
data2 <- NULL
data3 <- NULL

qs    <- (1:19/20)
qss   <- qs^2

# need to rename 1 SNP
names(data4)[names(data4)=="rs9680365.a_C"] <- "rs9680365_C"

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

data1 <- read.csv("~/Documents/CQR_17_04_2018/cream2017_ukbb_replicated.csv")
data2        <- data1[,c("ukbSNP","ukbCHR","ukbPOS","GENE.1","ukbA1","ukbA2","ukbBETA","ukbSE","ukbP_bolt")]
data2$EA     <- ifelse(data2$ukbBETA < 0, as.character(data2$ukbA1), as.character(data2$ukbA2))

data3         <- as.data.frame(matrix(nrow=149,ncol=1))
data3$start   <- var_names[snp_pos[1]:snp_pos[149]]
pieces        <- strsplit(data3$start,"_")
data3$ukbSNP  <- as.list(sapply(pieces, "[", 1))
data3$plinkA1 <- as.list(sapply(pieces, "[", 2))
data3$V1      <- NULL
data3$start   <- NULL
data6         <- merge(data2,data3,by="ukbSNP")
data6$switch  <- ifelse(data6$EA==data6$plinkA1,0,1)

# ensure SNP is coded so that test allele produces a negative beta
for (n in 1:num_snps){
  if(data6$switch[n]==1){
    old_snp  <- paste(data6$ukbSNP[n],"_",data6$plinkA1[n],sep="")
    new_snp  <- paste(data6$ukbSNP[n],"_",data6$EA[n],sep="")
    curr_col<-which(colnames(data4)==old_snp)
    data4[,curr_col] <- 2 - data4[,curr_col]
    names(data4)[curr_col] <- new_snp
  }
}

# collect SNP names again after re-coding
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

# Create Allele Score for avMSE (i.e. simply count number of risk alleles)

data5 <- read.csv("~/Documents/jez_ukb_summary.csv")
names(data5)[1] = "Marker"
data4$PGRS   <- get_score (data4,data5[,c("Marker", "Beta_OLS")],skip=118)
data4$Age2 <- (data4$Age)^2
m1 <- summary(lm(avMSE  ~ Sex + Age + Age2       ,data=data4))
m2 <- summary(lm(avMSE  ~ Sex + Age + Age2 + PGRS,data=data4))
m1
m2
m2$adj.r.squared - m1$adj.r.squared

# PGRS x UniEdu plot (10 lines)
# ----------------------------
pdata              <- expand.grid(avMSE=1, Sex=1, Age=58, Age2=3364, Geno_array=1, PGRS=c(-2,-1,0,1,2), UniEdu.x=c(0,1))
preds_est          <- as.data.frame(matrix(ncol=19,nrow=10))
preds_lci          <- as.data.frame(matrix(ncol=19,nrow=10))
preds_uci          <- as.data.frame(matrix(ncol=19,nrow=10))

names(preds_est)   <- c(paste("Q",1:19,sep=""))
names(preds_lci)   <- c(paste("Q",1:19,sep=""))
names(preds_uci)   <- c(paste("Q",1:19,sep=""))

preds_est$group    <- c("N1","N2","N3","N4","N5","U1","U2","U3","U4","U5")
preds_lci$group    <- c("N1","N2","N3","N4","N5","U1","U2","U3","U4","U5")
preds_uci$group    <- c("N1","N2","N3","N4","N5","U1","U2","U3","U4","U5")

myform            <- as.formula("AgeSpexWear ~ Sex + Age + Age2 + Geno_array + PGRS + UniEdu.x + PGRS:UniEdu.x")

for (i in 1:19){
  q                            <- i/20
  mod_qr                       <- try(rq(formula=myform, data=data4, tau=q))
  preds                        <- as.data.frame(predict.rq(mod_qr,pdata,interval = "confidence", level=0.95))
  preds_est[,i]                <- preds$fit
  preds_lci[,i]                <- preds$lower
  preds_uci[,i]                <- preds$higher
} # next quantile

est_preds                    <- melt(preds_est,id.vars="group")
lci_preds                    <- melt(preds_lci,id.vars="group")
uci_preds                    <- melt(preds_uci,id.vars="group")
names(est_preds)             <- c("group","quantile","avMSE")
names(lci_preds)             <- c("group","quantile","lci")
names(uci_preds)             <- c("group","quantile","uci")
all_preds                    <- cbind(est_preds,lci_preds$lci,uci_preds$uci)
names(all_preds)             <- c("group","quantile","avMSE","lci","uci")
all_preds$group              <- as.factor(all_preds$group)
all_preds$quantile           <- sort(rep(qs,10))
all_preds$edu_group          <- as.factor(c(rep(c(1,1,1,1,1,2,2,2,2,2),19)))
all_preds$pgrs_group         <- as.factor(c(rep(c(1,2,3,4,5),38)))

all_preds          <- all_preds[order(all_preds$group),]
all_preds$baseline <- rep(as.numeric(all_preds[which(all_preds$group=="N3"),]$avMSE),10)
all_preds$avMSEx   <- all_preds$avMSE - all_preds$baseline
all_preds$lcix     <- all_preds$lci - all_preds$baseline
all_preds$ucix     <- all_preds$uci - all_preds$baseline

plot1 <- ggplot(all_preds, aes(x=quantile, y=avMSE))+
  geom_ribbon(aes(x=quantile, ymin=lci, ymax=uci, fill=group),alpha=0.2)+
  geom_line(aes(colour=pgrs_group, linetype=edu_group),size=1)+
  scale_x_continuous(limits=c(0.05,0.95), breaks=seq(0.1,0.9,by=0.2), labels=seq(0.1,0.9,by=0.2)) +
  scale_y_continuous(limits=c(5,50), breaks=seq(5,50,by=5)) +
  #scale_fill_manual(guide=FALSE, values = c("blue","dark blue","black","dark red","red","blue","dark blue","black","dark red","red"))+
  scale_fill_manual(guide=FALSE, values = c(rep("grey",10)))+
  scale_colour_manual(name=NULL, values = c("blue","turquoise","green","orange","red"), labels = c("PGRS -2 SD","PGRS -1 SD","PGRS Average","PGRS +1 SD","PGRS +2 SD"))+
  scale_linetype_manual(name=NULL, values=c("solid","dotted"),labels=c("No degree","University degree"))+
  labs(x="Quantile",y="Age of Onset (Years)") +
  theme(legend.key.width = unit(1,"cm"))+
  theme(legend.spacing = unit(0.01,"cm"))+
  theme(legend.justification=c(1,0), legend.position=c(0.99,0.01))
plot1


plot2 <- ggplot(all_preds, aes(x=quantile, y=avMSEx))+
  geom_ribbon(aes(x=quantile, ymin=lcix, ymax=ucix, fill=group),alpha=0.06)+
  geom_line(aes(colour=pgrs_group, linetype=edu_group),size=1)+
  scale_x_continuous(limits=c(0.05,0.95), breaks=seq(0.1,0.9,by=0.2), labels=seq(0.1,0.9,by=0.2)) +
  scale_y_continuous(limits=c(-22.5,10), breaks=seq(-25,15,by=5)) +
  scale_fill_manual(guide=FALSE, values = c("blue","turquoise","green","orange","red","blue","turquoise","green","orange","red"))+
  #scale_fill_manual(guide=FALSE, values = c(rep("grey",10)))+
  scale_colour_manual(guide=FALSE, values = c("blue","turquoise","green","orange","red"), labels = c("PGRS -2 SD","PGRS -1SD","PGRS Average","PGRS +1SD","PGR +2SD"))+
  scale_linetype_manual(guide=FALSE, values=c("solid","dotted"),labels=c("No degree","University degree"))+
  labs(x="Quantile",y="Age of Onset (Years)") +
  theme(legend.key.width = unit(1,"cm"))+
  theme(legend.justification=c(1,0), legend.position=c(0.9,0.1))
plot2

tiff(filename = "C:/Users/Alfred/UKB_PGRSxEdu_AgeOfOnset.tif",
     width = 12, height = 20, units = "cm", pointsize = 12, compression = "lzw",bg = "white", res = 300, type = "windows")
grid.arrange(plot1,plot2,#plot4,plot5,plot6,
             #left=textGrob("Polygenic risk score effect size (D)", rot=90, vjust=1, gp=gpar(cex=1.1)),
             #bottom=textGrob ("Quantile", vjust=0, gp=gpar(cex=1.1)),
             ncol=1,clip=FALSE, heights=c(2,1))
dev.off()





