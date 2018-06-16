library(data.table); library(ggplot2); library(reshape)

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
      axis.title.x = element_text(colour = "black", size=rel(1), margin = margin(t = 6)),
      axis.title.y = element_text(colour = "black", size=rel(1), angle=90, margin = margin(r = 6)),
      axis.text.x = element_text(colour = "black", size=rel(0.9)),
      axis.text.y = element_text(colour = "black", size=rel(0.9)),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill="white"),
      legend.text = element_text(size=rel(0.7)),
      legend.title = element_text(size=rel(0.8)),
      legend.position=c(0.2,0.7),
      legend.key=element_blank(),
      plot.title = element_text(colour = "black", size=rel(1), margin = margin(b = 10),hjust = 0.5),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"))   
}
theme_set(theme_qq())

make_data2 <- function(data1){
  my1        <- as.data.frame(data1)
  names(my1) <- c("obs","bin")
  my2        <- my1[order(my1$obs),]
  my2$y1     <- -log10(my2$obs)
  x1         <- -log10( ppoints(length(my2$obs) ))
  my2$x1     <- x1
  N          <- dim(my2)[1]
  ci         <- 0.95
  my2$lci    <- -log10(qbeta((1 - ci)/2, 1:N, N:1))
  my2$uci    <- -log10(qbeta((1 + ci)/2, 1:N, N:1))
  return(my2)
}

scaleFUN <- function(x) sprintf("%.1f", x)


###########################################################################################################################
################################################# Intercept correction ####################################################
###########################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms) == "Int_only_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms = perms[,perm_num, with=F]; names(perms) = chr
stopwords = c("./mr_perms_"); colnames(perms) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms))
stopwords = c(".csv"); colnames(perms) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms))



dataMAF <- read.csv("~/Documents/jez_ukb_summary.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)


res            <- as.data.frame(t(perms))
res$SNP     <- rownames(res)
res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQint <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQint

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms = do.call(rbind.data.frame, chr) 
perms$pcor_int = -log10(perms$P_MR); perms$pcor_int = perms$pcor_int/m1$coef[2,1]
perms$pcor_int= 10^-perms$pcor_int; perms$pcor_int_form = formatC(perms$pcor_int, digits = 2, format = "e")

length(which(perms$pcor_int < 0.05/146))
###########################################################################################################################
###########################################################################################################################

###########################################################################################################################
################################################### Beta1 correction ######################################################
###########################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms2 = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms2) == "Int_qs_vs_Int_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms2 = perms2[,perm_num, with=F]; names(perms2) = chr
stopwords = c("./mr_perms_"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))
stopwords = c(".csv"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))



dataMAF <- read.csv("~/Documents/jez_ukb_summary.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)




res            <- as.data.frame(t(perms2))
res$SNP     <- rownames(res)
res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQqs <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQqs

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms2 = do.call(rbind.data.frame, chr) 
perms2$pcor_beta = -log10(perms2$P_CQR1); perms2$pcor_beta = perms2$pcor_beta/m1$coef[2,1]
perms2$pcor_beta = 10^-perms2$pcor_beta; perms2$pcor_beta_form = formatC(perms2$pcor_beta, digits = 2, format = "e")
perms2 = perms2[,c(1,14,15)]

length(which(perms2$pcor_beta < 0.05/146))
perms = merge(perms,perms2,"SNP")
###########################################################################################################################
###########################################################################################################################

###########################################################################################################################
################################################### Beta2 correction ######################################################
###########################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms2 = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms2) == "Int_qs_qss_vs_Int_qs_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms2 = perms2[,perm_num, with=F]; names(perms2) = chr
stopwords = c("./mr_perms_"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))
stopwords = c(".csv"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))



dataMAF <- read.csv("~/Documents/jez_ukb_summary.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)




res            <- as.data.frame(t(perms2))
res$SNP     <- rownames(res)
res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQqss <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQqss

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms2 = do.call(rbind.data.frame, chr) 
perms2$pcor_beta2 = -log10(perms2$P_CQR2); perms2$pcor_beta2 = perms2$pcor_beta2/m1$coef[2,1]
perms2$pcor_beta2 = 10^-perms2$pcor_beta2; perms2$pcor_beta2_form = formatC(perms2$pcor_beta2, digits = 2, format = "e")
perms2 = perms2[,c(1,14,15)]

length(which(perms2$pcor_beta2 < 0.05/146))
perms = merge(perms,perms2,"SNP")



perms2 = perms[order(as.numeric(perms$pcor_int)),]
cor.test(perms2$pcor_int,perms2$pcor_beta, method = "spearman")
cor.test(perms2$pcor_int,perms2$pcor_beta2, method = "spearman")
cor.test(perms2$pcor_beta,perms2$pcor_beta2, method = "spearman")







for(i in 1:nrow(perms2)){
  threshold = 0.05/(nrow(perms2)+1)
  if(as.numeric(perms2[i,15]) < threshold){
    perms2[i,15] = paste(perms2[i,15], "*", sep="")
  } else{
    perms2[i,15] = paste(perms2[i,15])
  }
}

for(i in 1:nrow(perms2)){
  threshold = 0.05/(nrow(perms2)+1)
  if(as.numeric(perms2[i,17]) < threshold){
    perms2[i,17] = paste(perms2[i,17], "*", sep="")
  } else{
    perms2[i,17] = paste(as.numeric(perms2[i,17]))
  }
}

for(i in 1:nrow(perms2)){
  threshold = 0.05/(nrow(perms2)+1)
  if(as.numeric(perms2[i,19]) < threshold){
    perms2[i,19] = paste(perms2[i,19], "*", sep="")
  } else{
    perms2[i,19] = paste(perms2[i,19])
  }
}


mr$uno = paste(mr$Beta_MR, "[", mr$LCI_MR, ";", mr$UCI_MR, "]", sep = "")
mr$dos = paste(mr$BETA_CQR1, "[", mr$LCI_CQR1, ";", mr$UCI_CQR1, "]", sep = "")
mr$tres = paste(mr$BETA_CQR2, "[", mr$LCI_CQR2, ";", mr$UCI_CQR2, "]", sep = "")



###################################################################################################################

################################################### Height ########################################################

###################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms) == "Int_only_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms = perms[,perm_num, with=F]; names(perms) = chr
stopwords = c("./mr_perms_"); colnames(perms) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms))
stopwords = c(".csv"); colnames(perms) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms))



dataMAF <- read.csv("~/Documents/CQR_17_04_2018/height_ols.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)


res            <- as.data.frame(t(perms))
res$SNP     <- rownames(res)
#res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQint <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQint

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms = do.call(rbind.data.frame, chr) 
perms$pcor_int = -log10(perms$P_MR); perms$pcor_int = perms$pcor_int/m1$coef[2,1]
perms$pcor_int= 10^-perms$pcor_int; perms$pcor_int_form = formatC(perms$pcor_int, digits = 2, format = "e")

length(which(perms$pcor_int < 0.05/146))
###########################################################################################################################
###########################################################################################################################

###########################################################################################################################
################################################### Beta1 correction ######################################################
###########################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms2 = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms2) == "Int_qs_vs_Int_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms2 = perms2[,perm_num, with=F]; names(perms2) = chr
stopwords = c("./mr_perms_"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))
stopwords = c(".csv"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))



dataMAF <- read.csv("~/Documents/CQR_17_04_2018/height_ols.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)




res            <- as.data.frame(t(perms2))
res$SNP     <- rownames(res)
#res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQqs <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQqs

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms2 = do.call(rbind.data.frame, chr) 
perms2$pcor_beta = -log10(perms2$P_CQR1); perms2$pcor_beta = perms2$pcor_beta/m1$coef[2,1]
perms2$pcor_beta = 10^-perms2$pcor_beta; perms2$pcor_beta_form = formatC(perms2$pcor_beta, digits = 2, format = "e")
perms2 = perms2[,c(1,14,15)]

length(which(perms2$pcor_beta < 0.05/146))
perms = merge(perms,perms2,"SNP")
###########################################################################################################################
###########################################################################################################################

###########################################################################################################################
################################################### Beta2 correction ######################################################
###########################################################################################################################

chr = list.files(pattern = "mr_perms.", full.names = TRUE)
chr = lapply(chr, fread)
perms2 = do.call(cbind.data.frame, chr); perm_num = perms[,which(colnames(perms2) == "Int_qs_qss_vs_Int_qs_p")]
chr = list.files(pattern = "mr_perms.", full.names = TRUE)
perms2 = perms2[,perm_num, with=F]; names(perms2) = chr
stopwords = c("./mr_perms_"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))
stopwords = c(".csv"); colnames(perms2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(perms2))



dataMAF <- read.csv("~/Documents/CQR_17_04_2018/height_ols.csv")
dataMAF        <- dataMAF[,c("SNP","MAF")]
stopwords = c("_A","_C","_T","_G"); dataMAF$SNP = gsub(paste0(stopwords, collapse = "|"), "", dataMAF$SNP)
rm(chr, stopwords,perm_num)




res            <- as.data.frame(t(perms2))
res$SNP     <- rownames(res)
#res = res[-which(res$SNP == "rs74764079"),]; res = res[-which(res$SNP == "rs73730144"),]; res = res[-which(res$SNP == "rs17837871"),]

dataM          <- merge(dataMAF,res,by="SNP")
dataMM         <- melt(dataM, id.vars = c("SNP","MAF"))
names(dataMM)  <- c("SNP","MAF","Permutation","Pval")
dataMM$logP    <- log10(dataMM$Pval)
dataMM$MAFbin  <- as.factor(ceiling(dataMM$MAF*10))
d2_1           <- make_data2(dataMM[,c("Pval","MAFbin")])
QQqss <- ggplot(data=d2_1)+ geom_ribbon( fill="light grey", aes(x=x1,ymin=lci, ymax=uci))+
  geom_line(aes(x=x1, y=x1), colour="red",size=1)+
  geom_line(size=1.5, aes(x=x1,y=y1,colour=bin))+
  labs(x=expression(paste("Expected -",log[10],"(P)")),y=expression(paste("Observed -",log[10],"(P)")))+
  scale_colour_manual(values=c("blue", "green", "orange", "purple", "pink"), name="MAF bin", labels=c("0.0 - 0.1","0.1 - 0.2","0.2 - 0.3","0.3 - 0.4","0.4 - 0.5"))+
  scale_x_continuous(labels=scaleFUN)+
  scale_y_continuous(labels=scaleFUN,limits=c(0,10)) +
  ggtitle("CQR Intercept (quantiles 0.05 - 0.95)")
QQqss

m1 <- summary(lm(y1~x1,data=d2_1)); m1; m1$coef[2,1]
chr = list.files(pattern = "mr_rs.", full.names = TRUE)
chr = lapply(chr, fread); perms2 = do.call(rbind.data.frame, chr) 
perms2$pcor_beta2 = -log10(perms2$P_CQR2); perms2$pcor_beta2 = perms2$pcor_beta2/m1$coef[2,1]
perms2$pcor_beta2 = 10^-perms2$pcor_beta2; perms2$pcor_beta2_form = formatC(perms2$pcor_beta2, digits = 2, format = "e")
perms2 = perms2[,c(1,14,15)]

length(which(perms2$pcor_beta2 < 0.05/146))
perms = merge(perms,perms2,"SNP")







