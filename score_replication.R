library(tidyr)
ma = read.csv("~/UKB/March/CQR/ma_cream2017_ukbb_rep_hits_allchr_pheno_alfredID.csv")
cream = read.csv("~/UKB/March/CQR/cream2017_ukbb_replicated.csv")


# Create first 6 ped columns
ped1 = as.data.frame(ma$alfred_ID1)
ped1[,2] = ped1$`ma$alfred_ID1`
ped1$MAT = 0
ped1$PAT = 0
ped1$Sex = 2
ped1$Phenotype = ma$avMSE

# Create ped genotypes
ped2 = as.data.frame(ma[,4:152])
ped3 = as.data.frame(matrix(nrow = nrow(ped2), ncol = ncol(ped2)*2))

j = 1
for (i in 1:ncol(ped2)){
  ped3[,j] = ped2[,i]
  names(ped3)[j] = names(ped2)[i]
  ped3[,j+1] = ped2[,i]
  names(ped3)[j+1] = names(ped2)[i]
  ped3[,j+1] = ifelse(ped3[,j+1] == 2,1, ifelse(ped3[,j+1] == 1,0,0))
  j = j + 2
}
ped3[ped3 == 2] = 1 # Separates alleles, e.g. people with two risk  alleles become homozygous
ped3[ped3 == 1] = 2 # Recoding risk allele
ped3[ped3 == 0] = 1 # Recoding other allele

ped = cbind(ped1,ped3)
write.table(ped, "ma.ped", row.names = F, col.names = F, quote = F)

# Create map file
map = as.data.frame(names(ma[,4:152]))
map$chr = 0
map$cm = 0
map$bp = 0

write.table(map[,c(2,1,3,4)], "ma.map", row.names = F, col.names = F, quote = F)



# Recode bim file with effect allele
ma.bim = read.table("~/ma.bim", quote="\"", comment.char="")
ma.bim = separate(ma.bim, V2, into = c("SNP","EA"), sep = "_")
write.table(ma.bim[,c(1,2,4,5,3,7)], "ma.bim", row.names = F, col.names = F, quote = F)



# Define variables 
{var_names   = names(ma)
vars        = as.matrix(var_names)
num_vars    = nrow(vars)
num_snps    = 0
snp_pos     = matrix(nrow = num_vars, ncol = 1)

for (n in 1:num_vars){
  if (substr(vars[n],1,2) == "rs"){
    num_snps = num_snps + 1
    snp_pos[num_snps] = n
  }
}

{ data2        = cream[, c("ukbSNP", "ukbCHR", "ukbPOS", "GENE.1", "ukbA1", "ukbA2", "ukbBETA", "ukbSE", "ukbP_bolt")]
  data2$EA     = ifelse(data2$ukbBETA < 0, as.character(data2$ukbA1), as.character(data2$ukbA2))
  
  data3        = as.data.frame(matrix(nrow = nrow(data2), ncol = 1))
  data3$start  = var_names[snp_pos[1]:snp_pos[nrow(data2)]]; data3$V1 = NULL
  data5        = separate(data = data3, col = start, into = c("ukbSNP", "plinkA1"), sep = "_")
  data6        = merge(data2, data5, "ukbSNP")
  data6$switch = ifelse(data6$EA == data6$plinkA1,0,1)}

for (n in 1:num_snps){
  if(data6$switch[n] == 1){
    old_snp  = paste(data6$ukbSNP[n],"_", data6$plinkA1[n], sep="")
    new_snp  = paste(data6$ukbSNP[n],"_", data6$EA[n], sep="")
    curr_col = which(colnames(ma) == old_snp)
    ma[, curr_col] = 2 - ma[, curr_col]
    names(ma)[curr_col] = new_snp
  }
}

# collect SNP names again after re-coding
var_names   = names(ma)
vars        = as.matrix(var_names)
num_vars    = nrow(vars)
num_snps    = 0
snp_pos     = matrix(nrow = num_vars, ncol = 1)

for (n in 1:num_vars){
  if (substr(vars[n],1,2) == "rs"){
    num_snps = num_snps + 1
    snp_pos[num_snps] = n
  }}}

ma2 = ma[,4:152]
cream$names = names(ma2)
cream = separate(cream, names, into = c("test","ea"), sep = "_")
cream$test = ma.bim$EA
# Create score file
write.table(cream[,c(1,22,15)], "ols_score.txt", row.names = F, col.names = F, quote = F)



# Variance explained
ma.profile = read.csv("~/ma.profile", sep="")
summary(lm(PHENO~SCORE,ma.profile))















yp = read.csv("~/UKB/March/CQR/yp_cream2017_ukbb_rep_hits_allchr_pheno_alfredID.csv")
cream = read.csv("~/UKB/March/CQR/cream2017_ukbb_replicated.csv")


# Create first 6 ped columns
ped1 = as.data.frame(yp$alfred_ID1)
ped1[,2] = ped1$`yp$alfred_ID1`
ped1$MAT = 0
ped1$PAT = 0
ped1$Sex = yp$sex
ped1$Phenotype = -9

# Create ped genotypes
ped2 = as.data.frame(yp[,13:161])
ped3 = as.data.frame(matrix(nrow = nrow(ped2), ncol = ncol(ped2)*2))

j = 1
for (i in 1:ncol(ped2)){
  ped3[,j] = ped2[,i]
  names(ped3)[j] = names(ped2)[i]
  ped3[,j+1] = ped2[,i]
  names(ped3)[j+1] = names(ped2)[i]
  ped3[,j+1] = ifelse(ped3[,j+1] == 2,1, ifelse(ped3[,j+1] == 1,0,0))
  j = j + 2
}
ped3[ped3 == 2] = 1 # Separates alleles, e.g. people with two risk  alleles become homozygous
ped3[ped3 == 1] = 2 # Recoding risk allele
ped3[ped3 == 0] = 1 # Recoding other allele

ped = cbind(ped1,ped3)
write.table(ped, "yp.ped", row.names = F, col.names = F, quote = F)

# Create map file
map = as.data.frame(names(yp[,13:161]))
map$chr = 0
map$cm = 0
map$bp = 0

write.table(map[,c(2,1,3,4)], "yp.map", row.names = F, col.names = F, quote = F)



# Recode bim file with effect allele
yp.bim = read.table("~/yp.bim", quote="\"", comment.char="")
yp.bim = separate(yp.bim, V2, into = c("SNP","EA"), sep = "_")
write.table(yp.bim[,c(1,2,4,5,3,7)], "yp.bim", row.names = F, col.names = F, quote = F)



# Define variables 
{var_names   = names(yp)
  vars        = as.matrix(var_names)
  num_vars    = nrow(vars)
  num_snps    = 0
  snp_pos     = matrix(nrow = num_vars, ncol = 1)
  
  for (n in 1:num_vars){
    if (substr(vars[n],1,2) == "rs"){
      num_snps = num_snps + 1
      snp_pos[num_snps] = n
    }
  }
  
  { data2        = cream[, c("ukbSNP", "ukbCHR", "ukbPOS", "GENE.1", "ukbA1", "ukbA2", "ukbBETA", "ukbSE", "ukbP_bolt")]
    data2$EA     = ifelse(data2$ukbBETA < 0, as.character(data2$ukbA1), as.character(data2$ukbA2))
    
    data3        = as.data.frame(matrix(nrow = nrow(data2), ncol = 1))
    data3$start  = var_names[snp_pos[1]:snp_pos[nrow(data2)]]; data3$V1 = NULL
    data5        = separate(data = data3, col = start, into = c("ukbSNP", "plinkA1"), sep = "_")
    data6        = merge(data2, data5, "ukbSNP")
    data6$switch = ifelse(data6$EA == data6$plinkA1,0,1)}
  
  for (n in 1:num_snps){
    if(data6$switch[n] == 1){
      old_snp  = paste(data6$ukbSNP[n],"_", data6$plinkA1[n], sep="")
      new_snp  = paste(data6$ukbSNP[n],"_", data6$EA[n], sep="")
      curr_col = which(colnames(yp) == old_snp)
      yp[, curr_col] = 2 - yp[, curr_col]
      names(yp)[curr_col] = new_snp
    }
  }
  
  # collect SNP names again after re-coding
  var_names   = names(yp)
  vars        = as.matrix(var_names)
  num_vars    = nrow(vars)
  num_snps    = 0
  snp_pos     = matrix(nrow = num_vars, ncol = 1)
  
  for (n in 1:num_vars){
    if (substr(vars[n],1,2) == "rs"){
      num_snps = num_snps + 1
      snp_pos[num_snps] = n
    }}}

yp2 = yp[,13:161]
cream$names = names(yp2)
cream = separate(cream, names, into = c("test","ea"), sep = "_")
#cream$test = ma.bim$EA
# Create score file
write.table(cream[,c(1,22,15)], "ols_score.txt", row.names = F, col.names = F, quote = F)



# Variance explained
yp.profile = read.csv("~/yp.profile", sep="")
names(yp)[1] = "FID"; yp.profile = merge(yp,yp.profile,"FID")
summary(lm(Yr7_avMSE~sex + age7 + SCORE,yp.profile))
summary(lm(Yr10_avMSE~sex + age10 + SCORE,yp.profile))
summary(lm(Yr11_avMSE~sex + age11 + SCORE,yp.profile))
summary(lm(Yr12_avMSE~sex + age12 + SCORE,yp.profile))
summary(lm(Yr15_avMSE~sex + age15 + SCORE,yp.profile))















