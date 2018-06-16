ukb = ukb[,c(1,2,3,5,8,9,11,13,15:24)]
ukb = na.omit(ukb)
table(ukb$Ethnicity)
ukb2 = ukb[which(ukb$Ethnicity == "White"),]

a = lm(avMSE ~ Sex + Age + UniEdu + GenoArray + GenoPCA1 + GenoPCA2 + GenoPCA3 + GenoPCA4 + 
         GenoPCA5 + GenoPCA6 + GenoPCA7 + GenoPCA8 + GenoPCA9 + GenoPCA10, data = ukb2)

a = as.data.frame(a$residuals)
ukb2$residuals = a$`a$residuals`

write.table(ukb2, file = "residuals_gxg_revised.txt", row.names = F, quote = F)

a = a[which(unrelated$X1 %in% a$FID),]
table(a$Ethnicity)


ukb = ukb[,c(1,2,3,5,8,9,11,13,15:29)]
ukb = na.omit(ukb)
table(ukb$Ethnicity)
ukb2 = ukb[which(ukb$Ethnicity == "White"),]

a = lm(avMSE ~ Sex + Age + UniEdu + GenoArray + GenoPCA1 + GenoPCA2 + GenoPCA3 + GenoPCA4 + 
         GenoPCA5 + GenoPCA6 + GenoPCA7 + GenoPCA8 + GenoPCA9 + GenoPCA10 + GenoPCA11 +
         GenoPCA12 + GenoPCA13 + GenoPCA14 + GenoPCA15, data = ukb2)

a = as.data.frame(a$residuals)
ukb2$residuals = a$`a$residuals`

write.table(ukb2, file = "residuals_gxg_15PCA.txt", row.names = F, quote = F)
