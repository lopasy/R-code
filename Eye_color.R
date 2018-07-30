eye_color = read.csv("D:/eye_color.raw", sep=""); eye_color = eye_color[,c(1,7:27)]
master_true = read.csv("D:/master_true.txt", sep="")

final = merge(master_true,eye_color,"FID"); final = final[,c(1:22,50:70)]

name = names(final)[23:43]

ukb9448 = read.csv("D:/ukb9448.csv")
ukb9448 = ukb9448[,c(1,2,5)]
names(ukb9448) = c("FID", "North", "East")

final = merge(final, ukb9448, "FID")
final = final[!is.na(final$North),]
final = final[-which(final$North == -1),]



num_snps = 21
snp = 1:num_snps
results_matrix=as.data.frame(matrix(nrow=num_snps,ncol=8))
colnames(results_matrix)=c("SNP","beta","LCI","UCI","P","MAF","North","East")

for (snp1 in snp){
  myform = as.formula(paste("avMSE ~ ", name[snp1], "+ UniEdu + Sex + Age + Array + PC1 + PC2 + 
                            PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + North + East", sep = ""))
  
  a = summary(lm(myform, data= final))
  results_matrix[snp1,1] = name[snp1] #SNP1 name check
  results_matrix[snp1,2] = a$coefficients[2,1]
  results_matrix[snp1,3] = a$coefficients[2,1] - 1.96*a$coefficients[2,2]
  results_matrix[snp1,4] = a$coefficients[2,1] + 1.96*a$coefficients[2,2]
  results_matrix[snp1,5] = a$coefficients[2,4]
  results_matrix[snp1,7] = a$coefficients[17,4]
  results_matrix[snp1,8] = round(a$coefficients[18,4],4)
  
  x = sum(final[, name[snp1]], na.rm = TRUE) / (2*sum(!is.na(final[,name[snp1]])))

  results_matrix[snp1,6] = as.numeric(formatC(min(x, (1-x)), digits = 5, format = "f"))
  
}   

length(which(results_matrix$P < 0.05/num_snps))



