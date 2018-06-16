a = lm(predicted_avMSE ~ PC1 + PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10 + Age + Sex + UniEdu + Array, predicted)
predicted_uniedu = as.data.frame(a$residuals)
a = lm(predicted_avMSE ~ PC1 + PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10 + Age + Sex + EduYearsHigh + Array, predicted)
predicted_eduyearshigh = as.data.frame(a$residuals)


a = as.data.frame(sample(predicted$FID, 50000, F))
a$IID = a$`sample(predicted$FID, 50000, F)`
write.table(a, "predicted_keep.txt", row.names = F, col.names = F, quote = F)



a = lm(avMSE ~ PC1 + PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10 + Age + Sex + UniEdu + Array, true)
uniedu = as.data.frame(a$residuals)
a = lm(avMSE ~ PC1 + PC2+ PC3+ PC4+ PC5+ PC6+ PC7+ PC8+ PC9+ PC10 + Age + Sex + EduYearsHigh + Array, true)
eduyearshigh = as.data.frame(a$residuals)
a = as.data.frame(sample(true$FID, 50000, F))


a$IID = a$`sample(true$FID, 50000, F)`
write.table(a, "true_keep.txt", row.names = F, col.names = F, quote = F)



true_pheno = true[,c(1,2,25:49)]
write.table(true_pheno, "true_pheno.txt", row.names = F, col.names = F, quote = F)
predicted_pheno = predicted[,c(1,2,24:36)]
write.table(predicted_pheno, "predicted_pheno.txt", row.names = F, col.names = F, quote = F)



predicted_uni = predicted[,c(1,2,21)]
predicted_edu = predicted[,c(1,2,23)]
write.table(predicted_uni, "predicted_uni.txt", row.names = F, col.names = F, quote = F)
write.table(predicted_edu, "predicted_edu.txt", row.names = F, col.names = F, quote = F)
true_uni = true[,c(1,2,21)]
true_edu = true[,c(1,2,24)]
write.table(true_uni, "true_uni.txt", row.names = F, col.names = F, quote = F)
write.table(true_edu, "true_edu.txt", row.names = F, col.names = F, quote = F)


true_qcovar = true[,c(1,2,4:13,19)]
true_covar = true[,c(1,2,3,22)]
write.table(true_qcovar, "true_qcovar.txt", row.names = F, col.names = F, quote = F)
write.table(true_covar, "true_covar.txt", row.names = F, col.names = F, quote = F)
predicted_qcovar = predicted[,c(1,2,4:13,19)]
predicted_covar = predicted[,c(1,2,3,22)]
write.table(predicted_qcovar, "predicted_qcovar.txt", row.names = F, col.names = F, quote = F)
write.table(predicted_covar, "predicted_covar.txt", row.names = F, col.names = F, quote = F)








