# CNV analysis
library(data.table)

# Arrange UKB dataset
{cnv = fread("F:/CNV Dataset 22Feb2018_Optometry.dat")
names = fread("F:/Full_ukb/March/ukb_sqc_descriptions.txt")
ukb = fread("F:/Full_ukb/March/ukb_sqc_v2.txt")
names(ukb) = colnames(names)
fam = read.table("F:/Full_ukb/March/ukb_snp_v2.fam")
fam = fam[,1:2]; names(fam) = c("FID","IID")

# Merge with CNV data
ukb = cbind(fam,ukb)
names(ukb)[3] = "array_ID"
ukb = merge(cnv,ukb,"array_ID")
ukb = ukb[,c(1,2,5:71)]

# Obtain relevant samples and remove individuals not passing quality control
true = fread("~/Documents/Useless/master_true.txt")
true = merge(true, ukb, "FID")
true = true[which(true$Filter_CNVQC == 1),]
true = as.data.frame(true)

rm(cnv, fam, names, ukb)}

{cnvs = c(47:60, 63:92, 94:100)
results = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 5))
names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV")

# If only main effects desired
for (i in cnvs){
  
  a = summary(lm(avMSE ~ Geno_array + UniEdu + Sex + Age + 
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i], true))
  
  results[i-46,1] = names(true)[i]
  results[i-46,2] = a$coefficients[16,1]
  results[i-46,3] = results[i-46,2] - 1.96*a$coefficients[16,2]
  results[i-46,4] = results[i-46,2] + 1.96*a$coefficients[16,2]
  results[i-46,5] = a$coefficients[16,4]
}}

# To test for interaction
{cnvs = c(47:51, 53:59, 63, 65, 66, 68:70, 72, 73, 75, 76, 78, 79, 81, 83, 84, 86:92, 94, 96, 98)
results2 = as.data.frame(matrix(nrow = length(cnvs)+17, ncol = 9))
names(results2) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "BETA_GxE", "LCI_GxE", "UCI_GxE", "P_GxE")

for (i in cnvs){
  
  a = summary(lm(avMSE ~ Geno_array + UniEdu + Sex + Age + 
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i] + UniEdu*true[,i], true))
  
  results2[i-46,1] = names(true)[i]
  results2[i-46,2] = a$coefficients[16,1]
  results2[i-46,3] = results2[i-46,2] - 1.96*a$coefficients[16,2]
  results2[i-46,4] = results2[i-46,2] + 1.96*a$coefficients[16,2]
  results2[i-46,5] = a$coefficients[16,4]
  
  results2[i-46,6] = a$coefficients[17,1]
  results2[i-46,7] = results2[i-46,6] - 1.96*a$coefficients[17,2]
  results2[i-46,8] = results2[i-46,6] + 1.96*a$coefficients[17,2]
  results2[i-46,9] = a$coefficients[17,4]
}}

{counts = as.data.frame(matrix(nrow = 54, ncol = 3))
names(counts) = c("Zero", "One", "Two")
for (i in 3:56){
  counts[i-2,1] = table(true[,i])[1]
  counts[i-2,2] = table(true[,i])[2]
  counts[i-2,3] = table(true[,i])[3]
}}

b = true[which(true$Medical_CNV == 1),]

{counts = as.data.frame(matrix(nrow = 54, ncol = 3))
  names(counts) = c("Zero", "One", "Two")
  for (i in 3:56){
    counts[i-2,1] = table(a[,i])[1]
    counts[i-2,2] = table(a[,i])[2]
    counts[i-2,3] = table(a[,i])[3]
  }}

{cnvs = c(47)
  a = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 5))
  names(a) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + b[,i], true))
    
    a[i-5,1] = names(true)[i]
    a[i-5,2] = a$coefficients[16,1]
    a[i-5,3] = a[i-5,2] - 1.96*a$coefficients[16,2]
    a[i-5,4] = a[i-5,2] + 1.96*a$coefficients[16,2]
    a[i-5,5] = a$coefficients[16,4]
  }}


### Predicted
names(ukb2)[1]="array_ID"
fam = read.table("F:/Full_ukb/March/ukb_snp_v2.fam", quote="\"", comment.char="")
fam = fam[,1:2]; names(fam) = c("FID","IID")

ukb2 = cbind(fam,ukb2)
z = merge(cnv,ukb2,"array_ID")
z = z[,1:71]
zz = merge(predicted, z, "FID")
zz = zz[which(zz$Filter_CNVQC == 1),]


{cnvs = c(47:100)
  results = as.data.frame(matrix(nrow = length(cnvs), ncol = 5))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(predicted_avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + zz[,i], zz))
    
    results[i-41,1] = names(zz)[i]
    results[i-41,2] = a$coefficients[16,1]
    results[i-41,3] = results[i-41,2] - 1.96*a$coefficients[16,2]
    results[i-41,4] = results[i-41,2] + 1.96*a$coefficients[16,2]
    results[i-41,5] = a$coefficients[16,4]
  }}

{cnvs = c(47:60, 63, 65, 66, 68:73, 75, 76, 78:81, 83:92, 94, 96, 98, 100)
  results2 = as.data.frame(matrix(nrow = length(cnvs)+17, ncol = 9))
  names(results2) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "BETA_GxE", "LCI_GxE", "UCI_GxE", "P_GxE")
  
  for (i in cnvs){
    
    a = summary(lm(predicted_avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + zz[,i] + UniEdu*zz[,i], zz))
    
    results2[i-46,1] = names(zz)[i]
    results2[i-46,2] = a$coefficients[16,1]
    results2[i-46,3] = results2[i-46,2] - 1.96*a$coefficients[16,2]
    results2[i-46,4] = results2[i-46,2] + 1.96*a$coefficients[16,2]
    results2[i-46,5] = a$coefficients[16,4]
    
    results2[i-46,6] = a$coefficients[17,1]
    results2[i-46,7] = results2[i-46,6] - 1.96*a$coefficients[17,2]
    results2[i-46,8] = results2[i-46,6] + 1.96*a$coefficients[17,2]
    results2[i-46,9] = a$coefficients[17,4]
  }}









# CNV analysis
library(data.table)

# Arrange UKB dataset
{cnv = fread("C:/Users/Alfred/Documents/ukb_geno_mse2017-08-29_cnv.csv")
  
  ###############################
  ### No filter, no education ###
  ###############################
  
  true = fread("D:/master_true.txt")
  true = true[,1]
  true = merge(true, cnv, "FID")
  true = as.data.frame(true)}

{cnvs = c(47:60,63:92,94:100)
  results = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 6))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "Freq")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Geno_array + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i], true))
    
    results[i-46,1] = names(true)[i]
    results[i-46,2] = a$coefficients[15,1]
    results[i-46,3] = results[i-46,2] - 1.96*a$coefficients[15,2]
    results[i-46,4] = results[i-46,2] + 1.96*a$coefficients[15,2]
    results[i-46,5] = a$coefficients[15,4]
    results[i-46,6] = formatC(table(true[,i])[2]/nrow(true), digits = 5, format = "f")
  }}

write.table(results, "my_unrelated_nofilter_nouniedu.txt", row.names = F, quote = F)

############################
### Filter, no education ###
############################

true = fread("D:/master_true.txt")
true = true[,1]
true = merge(true, cnv, "FID")
true = as.data.frame(true)
true = as.data.frame(true[which(true$Filter_CNVQC == 1),])

{cnvs = c(47:60,63:92,94:100)
  results = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 6))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "Freq")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Geno_array + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i], true))
    
    results[i-46,1] = names(true)[i]
    results[i-46,2] = a$coefficients[15,1]
    results[i-46,3] = results[i-46,2] - 1.96*a$coefficients[15,2]
    results[i-46,4] = results[i-46,2] + 1.96*a$coefficients[15,2]
    results[i-46,5] = a$coefficients[15,4]
    results[i-46,6] = formatC(table(true[,i])[2]/nrow(true), digits = 5, format = "f")
  }}

write.table(results, "my_unrelated_filter_nouniedu.txt", row.names = F, quote = F)


############################
### No filter, education ###
############################

true = fread("D:/master_true.txt")
true = true[,1]
true = merge(true, cnv, "FID")
true = as.data.frame(true)}

{cnvs = c(47:60,63:92,94:100)
  results = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 6))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "Freq")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i], true))
    
    results[i-46,1] = names(true)[i]
    results[i-46,2] = a$coefficients[16,1]
    results[i-46,3] = results[i-46,2] - 1.96*a$coefficients[16,2]
    results[i-46,4] = results[i-46,2] + 1.96*a$coefficients[16,2]
    results[i-46,5] = a$coefficients[16,4]
    results[i-46,6] = formatC(table(true[,i])[2]/nrow(true), digits = 5, format = "f")
  }}

write.table(results, "my_unrelated_nofilter_nouniedu.txt", row.names = F, quote = F)









true = fread("D:/master_true.txt")
true = true[,1]
true = merge(true, cnv, "FID")
true = as.data.frame(true)
true = as.data.frame(true[which(true$Filter_CNVQC == 1),])

{cnvs = c(47:60,63:92,94:100)
  results = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 6))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "Freq")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Geno_array + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i], true))
    
    results[i-46,1] = names(true)[i]
    results[i-46,2] = a$coefficients[15,1]
    results[i-46,3] = results[i-46,2] - 1.96*a$coefficients[15,2]
    results[i-46,4] = results[i-46,2] + 1.96*a$coefficients[15,2]
    results[i-46,5] = a$coefficients[15,4]
    results[i-46,6] = formatC(table(true[,i])[2]/nrow(true), digits = 5, format = "f")
  }}

write.table(results, "my_unrelated_filter_nouniedu.txt", row.names = F, quote = F)






























# To test for interaction
{cnvs = c(47:51, 53:59, 63, 65, 66, 68:70, 72, 73, 75, 76, 78, 79, 81, 83, 84, 86:92, 94, 96, 98)
  results2 = as.data.frame(matrix(nrow = length(cnvs)+17, ncol = 9))
  names(results2) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "BETA_GxE", "LCI_GxE", "UCI_GxE", "P_GxE")
  
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + true[,i] + UniEdu*true[,i], true))
    
    results2[i-46,1] = names(true)[i]
    results2[i-46,2] = a$coefficients[16,1]
    results2[i-46,3] = results2[i-46,2] - 1.96*a$coefficients[16,2]
    results2[i-46,4] = results2[i-46,2] + 1.96*a$coefficients[16,2]
    results2[i-46,5] = a$coefficients[16,4]
    
    results2[i-46,6] = a$coefficients[17,1]
    results2[i-46,7] = results2[i-46,6] - 1.96*a$coefficients[17,2]
    results2[i-46,8] = results2[i-46,6] + 1.96*a$coefficients[17,2]
    results2[i-46,9] = a$coefficients[17,4]
  }}

{counts = as.data.frame(matrix(nrow = 54, ncol = 3))
  names(counts) = c("Zero", "One", "Two")
  for (i in 3:56){
    counts[i-2,1] = table(true[,i])[1]
    counts[i-2,2] = table(true[,i])[2]
    counts[i-2,3] = table(true[,i])[3]
  }}

b = true[which(true$Medical_CNV == 1),]

{counts = as.data.frame(matrix(nrow = 54, ncol = 3))
  names(counts) = c("Zero", "One", "Two")
  for (i in 3:56){
    counts[i-2,1] = table(a[,i])[1]
    counts[i-2,2] = table(a[,i])[2]
    counts[i-2,3] = table(a[,i])[3]
  }}

{cnvs = c(47)
  a = as.data.frame(matrix(nrow = length(cnvs)+3, ncol = 5))
  names(a) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(avMSE ~ Array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + b[,i], true))
    
    a[i-5,1] = names(true)[i]
    a[i-5,2] = a$coefficients[16,1]
    a[i-5,3] = a[i-5,2] - 1.96*a$coefficients[16,2]
    a[i-5,4] = a[i-5,2] + 1.96*a$coefficients[16,2]
    a[i-5,5] = a$coefficients[16,4]
  }}


### Predicted
names(ukb2)[1]="array_ID"
fam = read.table("F:/Full_ukb/March/ukb_snp_v2.fam", quote="\"", comment.char="")
fam = fam[,1:2]; names(fam) = c("FID","IID")

ukb2 = cbind(fam,ukb2)
z = merge(cnv,ukb2,"array_ID")
z = z[,1:71]
zz = merge(predicted, z, "FID")
zz = zz[which(zz$Filter_CNVQC == 1),]


{cnvs = c(47:100)
  results = as.data.frame(matrix(nrow = length(cnvs), ncol = 5))
  names(results) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV")
  
  # If only main effects desired
  for (i in cnvs){
    
    a = summary(lm(predicted_avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + zz[,i], zz))
    
    results[i-41,1] = names(zz)[i]
    results[i-41,2] = a$coefficients[16,1]
    results[i-41,3] = results[i-41,2] - 1.96*a$coefficients[16,2]
    results[i-41,4] = results[i-41,2] + 1.96*a$coefficients[16,2]
    results[i-41,5] = a$coefficients[16,4]
  }}

{cnvs = c(47:60, 63, 65, 66, 68:73, 75, 76, 78:81, 83:92, 94, 96, 98, 100)
  results2 = as.data.frame(matrix(nrow = length(cnvs)+17, ncol = 9))
  names(results2) = c("CNV_NAME", "BETA_CNV", "LCI_CNV", "UCI_CNV", "P_CNV", "BETA_GxE", "LCI_GxE", "UCI_GxE", "P_GxE")
  
  for (i in cnvs){
    
    a = summary(lm(predicted_avMSE ~ Geno_array + UniEdu + Sex + Age + 
                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + zz[,i] + UniEdu*zz[,i], zz))
    
    results2[i-46,1] = names(zz)[i]
    results2[i-46,2] = a$coefficients[16,1]
    results2[i-46,3] = results2[i-46,2] - 1.96*a$coefficients[16,2]
    results2[i-46,4] = results2[i-46,2] + 1.96*a$coefficients[16,2]
    results2[i-46,5] = a$coefficients[16,4]
    
    results2[i-46,6] = a$coefficients[17,1]
    results2[i-46,7] = results2[i-46,6] - 1.96*a$coefficients[17,2]
    results2[i-46,8] = results2[i-46,6] + 1.96*a$coefficients[17,2]
    results2[i-46,9] = a$coefficients[17,4]
  }}











# Create ped file
cnv = fread("C:/Users/Alfred/Documents/ukb_geno_mse2017-08-29_cnv.csv")

# Obtain relevant samples and remove individuals not passing quality control
true = fread("D:/master_true.txt")
true = true[,1]
true = merge(true, cnv, "FID")
true = as.data.frame(true[which(true$Filter_CNVQC == 1),])

true = true[,c(1:2,6,4,47:100)]
true$PAT = 0
true$MAT = 0
true$Sex[true$Sex == 1] = 2; true$Sex[true$Sex == 0] = 1
true = true[,c(1,2,59,60,4,3,5:58)]
true[is.na(true)] = 0

cnvs = true[,7:60]

ped = as.data.frame(matrix(nrow = nrow(true), ncol = 2*ncol(cnvs)))
j = 1
for(i in 1:ncol(cnvs)){
  ped[,j] = cnvs[,i]
  ped[,j + 1] = 0
  j = j + 2
}
ped[ped == 1] = 2
ped[ped == 0] = 1

true = true[,1:6]
ped = cbind(true,ped)

write.table(ped, "cnv_europeans.ped", quote = F, row.names = F, col.names = F)



a = as.data.frame(names(cnvs[,1:54]))
a$CHR = 0
a$BP = seq(1,nrow(a),1)
a$CM = 0

write.table(a[,c(2,1,4,3)], "cnv_related.map", row.names = F, col.names = F, quote = F)
