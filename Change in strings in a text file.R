library(stringr)

fileName <- 'avMSE.R'

for(i in 1:149){
  
a = readChar(fileName, file.info(fileName)$size)

a = str_replace(a, "library(quantreg)", "library(quantreg)")
a = str_replace(a, "22,49,50", paste("22,49,", i + 49, sep = ""))
#a = gsub(paste(",c(3,19,22,49,50)", sep = ""), paste(",c(3,19,22,49,", i + 49, sep = ""), a)

a = gsub("06_cqr_mr_perms_int_p", paste("06_cqr_mr_perms_int_p", i, sep = ""), a)
a = gsub("06_cqr_mr_perms_beta1_p", paste("06_cqr_mr_perms_beta1_p", i, sep = ""), a)
a = gsub("06_cqr_mr_perms_beta2_p", paste("06_cqr_mr_perms_beta2_p", i, sep = ""), a)

a = gsub("06_cqr_mr15_85_perms_int_p", paste("06_cqr_mr15_85_perms_int_p", i, sep = ""), a)
a = gsub("06_cqr_mr15_85_perms_beta1_p", paste("06_cqr_mr15_85_perms_beta1_p", i, sep = ""), a)
a = gsub("06_cqr_mr15_85_perms_beta2_p", paste("06_cqr_mr15_85_perms_beta2_p", i, sep = ""), a)


writeLines(a, paste("avMSE", i, ".R", sep = ""))
}


fileName <- 'G:/Work/Rcode/BEDMatrix2.R'

snp = 10000
for(i in 1:45){
  
  a = readChar(fileName, file.info(fileName)$size)
  
  #a = str_replace(a, "ref = 1", paste("ref = ", i, sep = ""))
  a = str_replace(a, "num_snps = missing", paste("num_snps = ", "missing[", snp-9999, ":", snp, ",]", sep = ""))
  snp = snp + 10000
  files = file(paste("avMSE_GWAS", i, ".R", sep = ""))
  writeLines(a, files)
  close(files)
}


