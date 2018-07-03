eye_color <- read.csv("D:/eye_color.raw", sep=""); eye_color = eye_color[,c(1,7:27)]
master_true <- read.csv("D:/master_true.txt", sep="")

final = merge(master_true,eye_color,"FID"); final = final[,c(1:22,48:70)]

name = names(final)[25:45]


num_snps = 21
snp = 1:num_snps
results_matrix=as.data.frame(matrix(nrow=(num_snps*(num_snps-1))/2,ncol=7))
colnames(results_matrix)=c("SNP1check","SNP2check", "beta_SNP1:SNP2","se_SNP1:SNP2","P_SNP1:SNP2","P_SNP1","P_SNP2")

for (snp1 in 1:(num_snps-1)){
  for(snp2 in (snp1+1):num_snps) {
    myform = as.formula(paste("avMSE ~ ", name[snp1], "*", name[snp2], "+ UniEdu + Sex + Age + Array + PC1 + PC2 + 
                          PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", sep = ""))
    
    a = summary(lm(myform, data= final))
    results_matrix[snp,1] <- name[snp1] #SNP1 name check
    results_matrix[snp,2] <- name[snp2] #SNP2 name check
    results_matrix[snp,3] <- round(ifelse(dim(a$coefficients)[1] == 17,NA,a$coefficients[18,1]),3)
    results_matrix[snp,4] <- round(ifelse(dim(a$coefficients)[1] == 17,NA,a$coefficients[18,2]),3)
    results_matrix[snp,5] <- ifelse(dim(a$coefficients)[1] == 17,NA,a$coefficients[18,4])
    results_matrix[snp,6] <- ifelse(dim(a$coefficients)[1] == 17,NA,a$coefficients[2,4])
    results_matrix[snp,7] <- ifelse(dim(a$coefficients)[1] == 17,NA,a$coefficients[3,4])
    
    snp = snp + 1
  }
}    
results_matrix = results_matrix[1:(num_snps*(num_snps-1)/2),]
