library(data.table)
CQR = fread("D:/CQR.raw")
master = fread("D:/master_true.txt")


final = merge(master,CQR,"FID")
final = final[,c(3:13,19:22,55:200,202:204)]


results_matrix=as.data.frame(matrix(nrow=12000,ncol=5))
colnames(results_matrix)=c("SNP1check","SNP2check", "beta_SNP1:SNP2","se_SNP1:SNP2","P_SNP1:SNP2")

num_snps = length(16:ncol(final))
snp = 1
for (snp1 in 16:(num_snps-1)){
  for(snp2 in (snp1+1):num_snps) {
    
    a = names(final)[snp1]; b = names(final)[snp2]
    myform<-as.formula(paste("avMSE ~ Age+Array+Sex+UniEdu+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+",a, "*", b,  sep=""))
    res=summary(lm(myform,final))
    
    results_matrix[snp,1] <-a #SNP1 name check
    results_matrix[snp,2] <-b #SNP2 name check
    results_matrix[snp,3] <-ifelse(is.null(res),NA,res$coefficients[18,1])       #beta_SNP1
    results_matrix[snp,4] <-ifelse(is.null(res),NA,res$coefficients[18,2])      #se_SNP1
    results_matrix[snp,5]<-ifelse(is.null(res),NA,res$coefficients[18,4])       #P_SNP1
    
    cat("snp",b,":", "\n")
    
    #print(a$tTable[8,5])
    snp = snp + 1
  }
}

#write.csv(results_matrix, file="E:/CSV/SNP_SNP_AGE_incomplete.csv", row.names=FALSE)
