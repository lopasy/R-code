# Sample size creation for Age of Onset replication

test = cnv[which(cnv$Ethnicity.x == "Chinese"),] # Restrict to Chinese only (N_remaining = 1504)
test = test[which(test$het_miss_out == 0),] # Remove heterozygosity outliers (N_remaining = 150)
test = test[which(test$in_kinship == 0),] # Keep only unrelated as provided by UKB (N_remaining = 1353)
test = test[!is.na(test$AgeSpexWear),] # Remove people with missing phenotype (N_remaining = 1091)
test = test[!is.na(test$UniEdu),] # Remove people without education information (N_remaining = 1045)

write.table(test[,1:2], "Chinese_AoO_ID.txt", row.names = F, col.names = F, quote = F)
write.table(test[,c(1:4,6,12,14,18:27)], "Chinese_AoO_full_cov.txt", row.names = F, quote = F)


# Sample size creation for avMSE replication

test = cnv[which(cnv$Ethnicity.x == "Chinese"),] # Restrict to Chinese only (N_remaining = 1504)
test = test[which(test$het_miss_out == 0),] # Remove heterozygosity outliers (N_remaining = 150)
test = test[which(test$in_kinship == 0),] # Keep only unrelated as provided by UKB (N_remaining = 1353)
test = test[which(test$outlierMSE.x == 0),] # Remove individuals with ocular comorbidities (N_remaining = 1203)
test = test[!is.na(test$avMSE),] # Remove people with missing phenotype (N_remaining = 416)
test = test[!is.na(test$UniEdu),] # Remove people without education information (N_remaining = 401)

write.table(test[,1:2], "Chinese_avMSE_ID.txt", row.names = F, col.names = F, quote = F)
write.table(test[,c(1:4,6,12,14,18:27)], "Chinese_avMSE_full_cov.txt", row.names = F, quote = F)
