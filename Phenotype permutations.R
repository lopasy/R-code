##########################
### Permute phenotypes ###-------------------------### Not quite correct
##########################

## Data prep ##
EduYearsHigh3 = EduYearsHigh; gg = gg[,7:46] # gg is transposed bed file
{EduYearsHigh3 = cbind(EduYearsHigh3,gg)
##

start.time <- Sys.time()

permutation = 10000
snps = 27:66

## Initialize vector ##
y = matrix(0, permutation, max(snps)-min(snps)+1)
rep = 4
matList = lapply(1:ncol(y), function(x) replicate(rep, y[,x])) ## This is an array

## Linear regression for each SNP
out = matrix(nrow = permutation, ncol = rep)
for (j in snps) {
  for (i in 1:permutation) {
  EduYearsHigh3$avMSE = sample(EduYearsHigh3$avMSE, nrow(EduYearsHigh3), F) # label swapping
  
  model = summary(lm(as.formula(paste(colnames(EduYearsHigh3)[20], "~",
                      paste(colnames(EduYearsHigh3)[c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 22, 26, j)], collapse = "+"),
                      sep = "")), data=EduYearsHigh3))
  
  out[i, 1] = names(EduYearsHigh3)[j]
  out[i, 2] = model$coefficients[16,1]   # snp beta
  out[i, 3] = model$coefficients[16,2]   # snp se
  out[i, 4] = model$coefficients[16,4]   # snp p
  
  matList[[j-26]] = out
  }
}

matList = lapply(matList, function(x) {colnames(x) = c("SNP", "BETA", "SE", "P");x}) # Assign colnames

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken}




############################
### Find permuted p-vals ###
############################

## Create plink.linear output dataset
plink = plink[which(plink$TEST == "ADD"),]
plink$conc = paste0(plink$SNP, sep = "_", plink$A1)
plink = plink[,c(10,7,8,9)]

## Extracts p-vals for each SNP
permuted = matrix(0, permutation, max(snps)-min(snps)+1)
for (i in 1:length(matList)){
  permuted[,i] = matList[[i]][,4]
}

colnames(permuted) = plink[,1]
permuted = rbind(permuted, plink[which(plink$conc == colnames(permuted)), 4]) # Add original p-vals to permuted p-vals
ordered = apply(permuted, 2, sort,decreasing=F) # Sort from lowest to highest p-val


pvals = matrix(0, ncol(ordered), 3)
for (i in 1:ncol(ordered)){
  pvals[i,2] = which(plink[i,4] == ordered[,i]) # Obtains the rank of the SNP in empirical p-value distribution (under no association)
}

## Final data clean-up
pvals[,3] = pvals[,2]/nrow(ordered) # Permutation p-vals depend on the rank assigned
pvals[,1] = colnames(permuted)
pvals = cbind(pvals, plink[,4])
colnames(pvals) = c("SNP", "N_perm", "Perm_P", "P_linear")
pvals = as.data.frame(pvals)
pvals$Perm_P = as.numeric(paste(pvals$Perm_P))





