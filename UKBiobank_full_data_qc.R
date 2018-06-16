# UK Biobank quality control. There are 3 steps: data preparation for clumping, quantitative avMSE as binary and education as quantitative.
# 1. Prepare two files for normal avMSE and predicted_avMSE

#### Normal avMSE ####

# Combine FID with full UKB data
{ukb = cbind(fam,data)
colnames(ukb)[1] <- "FID"; colnames(ukb)[2] <- "IID"

# Remove individuals who have withdrawn
colnames(withdrawn)[1] <- "FID"
ukb = ukb[which(!(ukb$FID %in% withdrawn$FID)),]   # 8 people removed, 488369 remaining

# Merge with the file that has avMSE, predMSE, Eye disorder Outliers and UniEdu information
ukb2 = merge(ukb, predicted, by = 1) # 3 individuals removed (don't know why),   488366 remaining

# Remove individuals with other eye disorders
ukb2 = ukb2[which(ukb2$outlierMSE == 0),]          # 48142 people removed, 440224 remaining

# Remove non-British
ukb2 = ukb2[which(ukb2$in_white_British_ancestry_subset == 1),] # 69938 people removed, 370286 remaining

# Remove heterozigosity outliers
ukb2 = ukb2[which(ukb2$het_missing_outliers == 0),] # 648 people removed, 369638 remaining

# Remove related individuals
ukb2 = merge(ukb2,unrelated, by = 1) # 119722 people removed, 249916 remaining

# Remove gender mismatches
ukb2$conc_sex = paste(ukb2$Submitted_Gender, ukb2$Inferred_Gender, sep = "") # No one removed

# Select columns of interest
ukb3 = ukb2[,c(1,2,9,32:46,77,81,86,128)]

# Clean up
ukb3$conc_sex[ukb3$conc_sex == "MM"] <- 0; ukb3$conc_sex[ukb3$conc_sex == "FF"] <- 1
ukb3$genotyping_array[ukb3$genotyping_array == "UKBB"] <- 0; ukb3$genotyping_array[ukb3$genotyping_array == "UKBL"] <- 1
colnames(ukb3)[1] <- "FID"; colnames(ukb3)[2] <- "IID" 
colnames(ukb3)[22] <- "Sex"; colnames(ukb3)[3] <- "Array"

# Remove individuals without avMSE
ukb4 = na.omit(ukb3)} # 61426 individuals remaining with complete data

table(ukb4$UniEdu)[2]/table(ukb4$UniEdu)[1] # Education prevalence, 0.57
mean(ukb4$avMSE) # -0.26
sd(ukb4$avMSE) # 2.67


write.table(ukb4, file = "avMSE.txt", quote = F, row.names = F)




#### BOLT ####

# Combine FID with full UKB data
{ukb = cbind(fam,data)
colnames(ukb)[1] <- "FID"; colnames(ukb)[2] <- "IID"

# Remove individuals who have withdrawn
colnames(withdrawn)[1] <- "FID"
ukb = ukb[which(!(ukb$FID %in% withdrawn$FID)),]   # 8 people removed, 488369 remaining

# Merge with the file that has avMSE, predMSE, Eye disorder Outliers and UniEdu information
ukb2 = merge(ukb, predicted, by = 1) # 3 individuals removed (don't know why),   488366 remaining

# Remove individuals with other eye disorders
ukb2 = ukb2[which(ukb2$outlierMSE == 0),]          # 48142 people removed, 440224 remaining

# Remove non-British
ukb2 = ukb2[which(ukb2$in_white_British_ancestry_subset == 1),] # 69938 people removed, 370286 remaining

# Remove heterozigosity outliers
ukb2 = ukb2[which(ukb2$het_missing_outliers == 0),] # 648 people removed, 369638 remaining

# Remove gender mismatches
ukb2$conc_sex = paste(ukb2$Submitted_Gender, ukb2$Inferred_Gender, sep = "") # No one removed
toMatch <- c("MM", "FF")
ukb2 <- ukb2[grepl(paste(toMatch, collapse="|"), ukb2$conc_sex), ]

# Select columns of interest
ukb3 = ukb2[,c(1,2,9,32:46,77,81,86,127)]

# Clean up
ukb3$conc_sex[ukb3$conc_sex == "MM"] <- 0; ukb3$conc_sex[ukb3$conc_sex == "FF"] <- 1
ukb3$genotyping_array[ukb3$genotyping_array == "UKBB"] <- 0; ukb3$genotyping_array[ukb3$genotyping_array == "UKBL"] <- 1
ukb3$conc_sex[ukb3$conc_sex == "MM"] <- 0; ukb3$conc_sex[ukb3$conc_sex == "FF"] <- 1
colnames(ukb3)[1] <- "FID"; colnames(ukb3)[2] <- "IID" 
colnames(ukb3)[22] <- "Sex"; colnames(ukb3)[3] <- "Array"

# Remove individuals without avMSE
ukb4 = na.omit(ukb3)} # 86286 related individuals remaining with complete data

table(ukb4$UniEdu)[2]/table(ukb4$UniEdu)[1] # Education prevalence, 0.54
mean(ukb4$avMSE) # -0.229
sd(ukb4$avMSE) # 2.66


write.table(ukb4, file = "avMSE.txt", quote = F, row.names = F)
write.table(ukb4[,c(1,2)], file = "avMSE_ID.txt", quote = F, row.names = F)



#### Predicted avMSE ####

# Combine FID with full UKB data
{ukb = cbind(fam,data)
colnames(ukb)[1] <- "FID"; colnames(ukb)[2] <- "IID"

# Remove individuals who have withdrawn
colnames(withdrawn)[1] <- "FID"
ukb = ukb[which(!(ukb$FID %in% withdrawn$FID)),]   # 8 people removed, 488369 remaining

# Merge with the file that has avMSE, predMSE, Eye disorder Outliers and UniEdu information
ukb2 = merge(ukb, predicted, by = 1) # 3 individuals removed (don't know why),   488366 remaining

# Remove individuals with other eye disorders
ukb2 = ukb2[which(ukb2$outlierMSE == 0),]          # 48142 people removed, 440224 remaining

# Remove non-British
ukb2 = ukb2[which(ukb2$in_white_British_ancestry_subset == 1),] # 69938 people removed, 370286 remaining

# Remove heterozigosity outliers
ukb2 = ukb2[which(ukb2$het_missing_outliers == 0),] # 648 people removed, 369638 remaining

# Remove gender mismatches
ukb2$conc_sex = paste(ukb2$Submitted_Gender, ukb2$Inferred_Gender, sep = "") # No one removed
toMatch <- c("MM", "FF")
ukb2 <- ukb2[grepl(paste(toMatch, collapse="|"), ukb2$conc_sex), ]

# Select columns of interest
ukb3 = ukb2[,c(1,2,9,32:46,77,82,86,127)]

# Clean up
ukb3$conc_sex[ukb3$conc_sex == "MM"] <- 0; ukb3$conc_sex[ukb3$conc_sex == "FF"] <- 1
ukb3$genotyping_array[ukb3$genotyping_array == "UKBB"] <- 0; ukb3$genotyping_array[ukb3$genotyping_array == "UKBL"] <- 1
ukb3$conc_sex[ukb3$conc_sex == "MM"] <- 0; ukb3$conc_sex[ukb3$conc_sex == "FF"] <- 1
colnames(ukb3)[1] <- "FID"; colnames(ukb3)[2] <- "IID" 
colnames(ukb3)[22] <- "Sex"; colnames(ukb3)[3] <- "Array"

# Remove individuals without avMSE
ukb4 = na.omit(ukb3)} # 240125 related individuals remaining with complete data

table(ukb4$UniEdu)[2]/table(ukb4$UniEdu)[1] # Education prevalence, 0.44
mean(ukb4$predicted_avMSE) # Mean is -0.32
sd(ukb4$predicted_avMSE) # SD is 1.54

write.table(ukb4, file = "predictedMSE.txt", quote = F, row.names = F)
write.table(ukb4[,c(1,2)], file = "predictedMSE_ID.txt", quote = F, row.names = F)


#### Edu as quantitative ####

# Merge with the file containing UniEdu as quantitative
edu = merge(ukb4, quant, by = 1)
edu = edu[,c(1:22,35)]
colnames(edu)[19] <- "Age"; colnames(edu)[20] <- "avMSE"; colnames(edu)[22] <- "Sex"

# Split the file to separate those with known EducationAgeOLD and missing
split1 = na.omit(edu) # 55922 individuals with no missing UniEdu data

'table(split1$UniEdu, split1$EducationAgeOLD)
       13    14    15    16    17    18    19    20    21    22    23    24    25    26
  0   513   500 14984 18757  6545  7126  1666  1037  2505   762   355   154   179   341
  1     1     1    40    83    38    83    26    23   137    38    15     3     3     7' # This how many individuals have both types of information

library(ggplot2); library(gridExtra)
a = as.data.frame(table(split1$UniEdu, split1$EducationAgeOLD))
colnames(a)[1] <- "Education"; colnames(a)[2] <- "Years of Education"; colnames(a)[3] <- "Number of Individuals"
b = table(split1$UniEdu, split1$EducationAgeOLD)
ggplot(data=a, aes_string(x = "`Years of Education`", y= "`Number of Individuals`", fill = "Education")) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black")+
  annotation_custom(tableGrob(t(b)), xmin=10, xmax=16, ymin=11000, ymax=16000)+
  theme_bw()+
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black"))+
  ggtitle("People with Known UniEdu and EducationAgeOLD information")

split2 = edu[which(!(edu$FID %in% split1$FID)),] # 30364 individuals who don't have EducationAgeOLD
'table(split2$UniEdu)
0     1 
551 29813'  # This is how many individuals don't have EducationAgeOLD but have UniEdu information 

# Those who have UniEdu = 1, should at least have EducationAgeOLD = 21, so change those to 21
split2$EducationAgeOLD[split2$UniEdu == 1] <- 21

edu = rbind(split1, split2)
a = as.data.frame(table(edu$UniEdu, edu$EducationAgeOLD))
colnames(a)[1] <- "Education"; colnames(a)[2] <- "Years of Education"; colnames(a)[3] <- "Number of Individuals"
b = table(edu$UniEdu, edu$EducationAgeOLD)
ggplot(data=a, aes_string(x = "`Years of Education`", y= "`Number of Individuals`", fill = "Education")) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black")+
  annotation_custom(tableGrob(t(b)), xmin=9, xmax=15, ymin=15000, ymax=27000)+
  theme_bw()+
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black"))+
  ggtitle("People with Known UniEdu and EducationAgeOLD information")

a = na.omit(edu$EducationAgeOLD)
sd(a) # Continious environmental variable SD is 2.626

#### avMSE as binary ####

edu$cc[edu$avMSE <= -0.5] <- 1; edu$cc[edu$avMSE > -0.5] <- 0
table(edu$cc)[2]/table(edu$cc)[1] # Disease prevalence 0.50
table(edu$cc)[1]/table(edu$cc)[2] # ~2 controls per case
write.table(edu, file = "avMSE.txt", quote = F, row.names = F)





