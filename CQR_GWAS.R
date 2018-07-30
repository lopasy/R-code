library(data.table);library(qqman)

bim = read.table("D:/trues.bim", quote="\"", comment.char="")
names(bim)[7] = "SNP"

chr = list.files(pattern = "meta.", full.names = TRUE)
chr = lapply(chr, fread)
mr10 = do.call(rbind.data.frame, chr)

mr10 = mr10[-(which(mr10$SNP == "SNP")),]

mr10 = mr10[!duplicated(mr10$SNP),]

stopwords = c("_A","_C","_T","_G"); stopwords2 = c("A","C","T","G")
mr10$SNP = gsub(paste0(stopwords, collapse = "|"), "", mr10$SNP)
mr10$SNP = gsub(paste0(stopwords2, collapse = "|"), "", mr10$SNP)
mr = merge(mr10, bim, "SNP")

mr$P_MR = as.numeric(as.character(mr$P_MR))
mr$P_CQR1 = as.numeric(as.character(mr$P_CQR1))
mr$P_CQR2 = as.numeric(as.character(mr$P_CQR2))

mr_m = mr[which(mr$P_MR < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "P_MR", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))
mr_m = mr[which(mr$P_CQR1 < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "P_CQR1", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))
mr_m = mr[which(mr$P_CQR2 < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "P_CQR2", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))

qq(mr$P_MR); qq(mr$P_CQR1); qq(mr$P_CQR2)

test = mr[which(mr$P_MR < 0.05/(3000000)),]
test = test[which(test$V1 == 14),]

test = mr[which(mr$P_CQR1 < 0.05/(3000000)),]
test = mr[which(mr$P_CQR2 < 0.05/(3000000)),]



# Calculate lambda gc (??gc)
chisq = qchisq(mr$P_MR,1, lower.tail = F);int_c = median(chisq)/qchisq(0.5,1)
newchisq <- chisq/int_c; mr$int_cc <- pchisq(newchisq, df=1,lower.tail=FALSE)
chisq = qchisq(mr$P_CQR1,1, lower.tail = F);int_l = median(chisq)/qchisq(0.5,1)
newchisq <- chisq/int_l; mr$int_ll <- pchisq(newchisq, df=1,lower.tail=FALSE)
chisq = qchisq(mr$P_CQR2,1, lower.tail = F);int_q = median(chisq)/qchisq(0.5,1)
newchisq <- chisq/int_q; mr$int_qq <- pchisq(newchisq, df=1,lower.tail=FALSE)

mr_m = mr[which(mr$int_cc < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "int_cc", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))
mr_m = mr[which(mr$int_ll < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "int_ll", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))
mr_m = mr[which(mr$int_qq < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "int_qq", suggestiveline = -log10(0.000001/3),
          genomewideline = -log10(0.05/(3000000)), cex = 0.1, col = c("goldenrod3", "darkblue"))

chisq = qchisq(mr$int_cc,1, lower.tail = F); median(chisq)/qchisq(0.5,1)
chisq = qchisq(mr$int_ll,1, lower.tail = F); median(chisq)/qchisq(0.5,1)
chisq = qchisq(mr$int_qq,1, lower.tail = F); median(chisq)/qchisq(0.5,1)

qq(mr$int_cc); qq(mr$int_ll); qq(mr$int_qq)

test = mr[which(mr$int_cc < 0.05/(3000000)),]
test = mr[which(mr$int_ll < 0.05/(3000000)),]
test = mr[which(mr$int_qq < 0.05/(3000000)),]


# Correction for 90%
data3<-subset(mr,select=c("P_MR"))
plen <- nrow(data3)
data3<- data3[order(data3[,1]), ]
pran <- rchisq(plen,df=1,ncp = 0)
pran <- as.data.frame(pran)
pran <- pran[order(pran[,1]), ]
data4 <- data.frame(cbind(data3,pran))
plen90<-0.9*plen
data5 <-data4[1:plen90,]
lambda<-mean(data5$P_MR)/mean(data5$pran)

chisq = qchisq(1-mr$P_MR,1)
newchisq <- chisq/lambda

mr$cc_90 <- pchisq(newchisq, df=1,lower.tail=FALSE)
mr_m = mr[which(mr$cc_90 < 0.0005),]
manhattan(mr_m, chr = "V1", bp = "V4", p = "cc_90", genomewideline = -log10(0.05/(3000000)), cex = 0.1)





