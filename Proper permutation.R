

gg = gg[,7:46] # gg is transposed bed file
{EduYearsHigh3 = EduYearsHigh
 EduYearsHigh3 = cbind(EduYearsHigh3,gg)

  
  permutation = 10000
  snps = 27:66
  
  ## Initialize vector ##
  sig = matrix(0, permutation, 1)
  
  i = 1
  j = 27
  
  for (i in 1:permutation) {
    out = matrix(nrow = max(snps)-min(snps)+1, ncol = 1)
    for (j in snps) {
      EduYearsHigh3$avMSE = sample(EduYearsHigh3$avMSE, nrow(EduYearsHigh3), F) # label swapping
      
      model = summary(lm(as.formula(paste(colnames(EduYearsHigh3)[20], "~",
                                          paste(colnames(EduYearsHigh3)[c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 22, 26, j)], collapse = "+"),
                                          sep = "")), data=EduYearsHigh3))
      
      out[j-min(snps)+1, 1] = model$coefficients[16,4]   # snp p

    }
    sig.pvals = length(which(out[, 1] < 0.00125))      # how many significant
    sig[i,] = sig.pvals
    out = NULL
  }
  
  EduYearsHigh3 = EduYearsHigh
  EduYearsHigh3 = cbind(EduYearsHigh3,gg)
  
  
  permutation = 10000
  snps = 27:66
  
  ## Initialize vector ##
  sig_0.05 = matrix(0, permutation, 1)
  
  i = 1
  j = 27
  
  for (i in 1:permutation) {
    out = matrix(nrow = max(snps)-min(snps)+1, ncol = 1)
    for (j in snps) {
      EduYearsHigh3$avMSE = sample(EduYearsHigh3$avMSE, nrow(EduYearsHigh3), F) # label swapping
      
      model = summary(lm(as.formula(paste(colnames(EduYearsHigh3)[20], "~",
                                          paste(colnames(EduYearsHigh3)[c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 19, 22, 26, j)], collapse = "+"),
                                          sep = "")), data=EduYearsHigh3))
      
      out[j-min(snps)+1, 1] = model$coefficients[16,4]   # snp p
      
    }
    sig.pvals = length(which(out[, 1] < 0.05))      # how many significant
    sig_0.05[i,] = sig.pvals
    out = NULL
  }
  
}

#########################################
### If needed to produce permutations ###
#########################################
summary(lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,true))$fstat[1]

fstats = numeric(10000)
for (i in 1:10000){
  ge = lm(sample(avMSE)~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,true)
  fstats [i] = summary(ge)$fstat[1]
}
length(fstats[fstats > 306.9984])/10000
mean(fstats)
hist(fstats)




