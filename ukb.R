library(readr)
unrelated <- read_delim("~/Documents/ukb_euro_het_unrelated.keep", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) # For ukb_euro_het_unrelated.keep
names(unrelated) <- c("FID","IID")
ukb <- read_delim("~/Documents/ukb_geno_mse2017-02-16.txt",  " ", escape_double = FALSE)

data3 <- ukb[,c(1,3,5,8,9,11,13,15:29)]
data4 <- merge(unrelated,data3,by="FID")
#data5 <- ukb[which(ukb$IID %in% unrelated$X1),]
table(data4$Ethnicity)


ukb = data4[which(data4$Ethnicity == "White"),]
ukb = na.omit(data4) # Remove those with missing covariates

' My incorrect way
ukb = ukb[which(unrelated$X1 %in% ukb$IID),] # Select IDs in ukb pheno file from unrelated file
table(ukb$Ethnicity)
ukb = na.omit(ukb) # Remove those with missing covariates'
dim(ukb)
# 6965   18


# Mke a list of unrelated Europeans with no missing covariates
a = ukb[,c(1,2)]

write.table(a, file = "unrelated.txt", row.names = F, col.names = F, quote = F)
write.table(ukb, file = "cov.txt", row.names = F,  quote = F)



a = lm(avMSE ~ Sex+Age+UniEdu+GenoArray+GenoPCA1+GenoPCA2+GenoPCA3+GenoPCA4+GenoPCA5+GenoPCA6+GenoPCA7+GenoPCA8+GenoPCA9+GenoPCA10+GenoPCA11+GenoPCA12+GenoPCA13+GenoPCA14+GenoPCA15, data = ukb)
ukb$res = a$residuals

write.table(ukb, file = "residuals_gxg_15PCA.txt", row.names = F,  quote = F)


