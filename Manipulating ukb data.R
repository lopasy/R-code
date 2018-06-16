library(readr)
unrelated <- read_delim("~/Documents/ukb_euro_het_unrelated.keep", " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) # For ukb_euro_het_unrelated.keep
names(unrelated) <- c("FID","IID")
ukb <- read_delim("~/Documents/ukb_geno_mse2017-02-16.txt",  " ", escape_double = FALSE)

data3 <- ukb[,c(1,3,5,8,9,11,13,15:24)]
data4 <- merge(unrelated,data3,by="FID")
#data5 <- ukb[which(ukb$IID %in% unrelated$X1),]
table(data4$Ethnicity)

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
