plain_eu=ped[which(!ped$FID %in% z$V1),]



write.table(unrelated[,1:2], "remove_unrelated.txt", row.names = F, quote = F)



china_cnvs = as.data.frame(matrix(nrow = 54, ncol = 1))
for(i in 47:100){
  china_cnvs[i-46,1] = table(z[,i])[2]
}

china_cnvs = which(!is.na(china_cnvs$V1))


for(i in china_cnvs){
  print(summary(lm(avMSE~Sex+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+z[,i],z)))
}


ukb5 = ukb2[which(!ukb2$FID %in% withdrawn$V1),]
ukb5 = ukb5[which(ukb5$outlierMSE == 0),]


Europeans = read.table("~/Documents/CQR_17_04_2018/GWAS/ukb_v2_10xSD_PC-derived_Europeans.keep", quote="\"", comment.char="")

ukb5 = cnv_covar[which(!cnv_covar$FID %in% withdrawn$V1),]
#ukb5=merge(ukb5,outlier,1)
ukb5 = ukb5[which(ukb5$outlierMSE == 0),]
ukb5 = ukb5[which(ukb5$White_British == 1),]
ukb5 = merge(ukb5,het_out,1)
ukb5 = ukb[which(ukb$het_miss_out == 0),]
ukb5 = ukb5[which(ukb5$in_kinship == 0),]

ukb5 = ukb5[which(ukb5$FID %in% Europeans$V1),]
ukb5 = merge(ukb5,ukb,1)

uk = merge(ukb5,spex_all,1)
uk = uk[!is.na(uk$AgeSpexWear.y),]


#Kolmogorov-Smirnov (procedure 1)
ukb5 = ukb[which(ukb$het_miss_out == 0),]
ukb5 = ukb5[which(ukb$FID %in% Europeans$V1),]
ukb5 = ukb5[which(ukb5$in_kinship == 0),]
ukb5 = ukb5[which(ukb5$FID %in% true_europeans_remove$V1),]

cnvs = 47:100
pvals = as.data.frame(matrix(nrow = length(cnvs), ncol = 2))
names(pvals) = c("CNV", "P")
for(i in cnvs){
  ukb_cnv = ukb5[which(ukb5[,i] == 1),]
  ukb_no_cnv = ukb5[which(ukb5[,i] == 0),]
  a = tryCatch(ks.test(ukb_cnv$avMSE, ukb_no_cnv$avMSE)$p.value, error=function(err) NA)

  pvals[i-46,2] = a
  pvals[i-46,1] = names(ukb5)[i]
}

#Kolmogorov-Smirnov (procedure 2)
ukb = ukb[which(ukb$het_miss_out == 0),]
ukb2 = ukb[which(ukb$FID %in% Europeans$V1),]
ukb2 = ukb2[which(ukb2$FID %in% true_europeans_remove$V1),]
cnvs = c(47:100)
num_shuffles = 1000
perms = as.data.frame(matrix(nrow = num_shuffles, ncol = 1))
pvals = as.data.frame(matrix(nrow = length(cnvs), ncol = 2))
names(pvals) = c("CNV", "P")
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
for(i in cnvs){
  for(j in 1:num_shuffles){
    ukb2$shuffle = shuffle(ukb2[,i])
    ukb_cnv = ukb2[which(ukb2$shuffle == 1),]
    ukb_no_cnv = ukb2[which(ukb2$shuffle == 0),]
    perms[j,1] = tryCatch(ks.test(ukb_cnv$avMSE, ukb_no_cnv$avMSE)$p.value, error=function(err) NA)
  }
  
  ukb_cnv = ukb2[which(ukb2[,i] == 1),]
  ukb_no_cnv = ukb2[which(ukb2[,i] == 0),]
  a = tryCatch(ks.test(ukb_cnv$avMSE, ukb_no_cnv$avMSE)$p.value, error=function(err) NA)
  if (is.numeric(a)) pvals[i-46,2] = length(which(perms[,1] < a))/num_shuffles
  pvals[i-46,1] = names(ukb2)[i]
}


a = lm(AgeSpexWear~CNV_TAR_del+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+UniEdu,ukb5)


z = cnv_covar[which(!cnv_covar$FID %in% remove_spex_true_europeans_excluded$V1),]
z = z[!is.na(z$CNV_22q11.2dup),]
z$CNV_22q11.2dup[z$CNV_22q11.2dup == 1] = "Cases"
z$CNV_22q11.2dup[z$CNV_22q11.2dup == 0] = "Controls"
table(z$CNV_22q11.2dup)

ppi = 70
p = ggplot(z, aes(factor(CNV_22q11.2dup), AgeSpexWear, fill=factor(CNV_22q11.2dup))) +
  geom_violin(trim = FALSE, width = 1.1) + 
  geom_boxplot(width = 0.08,outlier.shape=NA, fill = "white") +
  scale_fill_brewer(palette="Blues") +
  guides(fill=FALSE) +
  labs(x = "CNV_22q11.2dup", y = "Years") +
  theme_classic() +
  theme(axis.text= element_text(size=10,colour = "black"), axis.title=element_text(size=11,face="bold"))

file_out = paste("F:/Full_ukb/May/", "CNV_22q11.2dup_spex", ".png", sep = "")
png(file_out, width = 4*ppi, height = 4.5*ppi, res = ppi)
print(p)
dev.off()





z = z[!is.na(z$AgeSpexWear),]
zz1 = z[which(z$CNV_22q11.2dup == 0),]
zz2 = z[which(z$CNV_22q11.2dup == 1),]
vioplot(zz1$AgeSpexWear, zz2$AgeSpexWear)






