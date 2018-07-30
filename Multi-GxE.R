# Multi-GxE

master_true = read.csv("D:/master_true.txt", sep="")
ukb9448 = read.csv("D:/CQR GWAS/ukb9448.csv")

test = ukb9448[,c(1,2,5)]
names(test) = c("FID", "North", "East")

test = merge(master_true, test, "FID")
test = test[!is.na(test$North),]
test = test[-which(test$North == -1),]


write.table(test[,1:2], "multi_sampleID.txt", row.names = F, col.names = F, quote = F)

write.table(test[,c("UniEdu", "North", "East")], "multi_gxe_env.txt", row.names = F, col.names = F, quote = F)


test2 = as.data.frame(t(test))
colnames(test2) = test2["FID",]
test2 = test2["avMSE",]
test2[2,] = 0
write.table(test2, "multi_gxe_pheno.txt", quote = F)


test$intercept = 1
write.table(test[,c(52,3:13,22,19)], "multi_gxe_cov.txt", row.names = F, col.names = F , quote = F)
