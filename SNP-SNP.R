library (nlme)
{var_names<-names(all[,c(11:161)])
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
testFunction <- function (formula_in) {
  myformula <- as.formula(formula_in)
  return(tryCatch(summary(do.call("lme", args = list(myformula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                                     correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                                     na.action = na.omit,  method="ML", data=all))), error=function(e) NULL))
}
results_matrix=matrix(nrow=11175,ncol=5)
colnames(results_matrix)=c("SNP1check","SNP2check", "beta_SNP1:SNP2","se_SNP1:SNP2","P_SNP1:SNP2")
results_matrix<-as.data.frame(results_matrix)

snp = 1:num_snps
for (snp1 in 1:(num_snps-1)){
  for(snp2 in (snp1+1):num_snps) {
    
    b = paste0("+", vars[snp_pos[snp2]])
    formula<-as.formula(paste("mse_at_visit ~ ", 
                              factor(vars[snp_pos[snp1]]), b,
                              " + poly(I(age_at_visit - 7.5),4) + ",
                              factor(vars[snp_pos[snp1]]): factor(vars[snp_pos[snp2]]),
                              " + I(age_at_visit - 7.5):",vars[snp_pos[snp1]],
                              " + I(age_at_visit - 7.5):",vars[snp_pos[snp2]],  sep=""))
    a<-testFunction(formula)
    
    results_matrix[snp,1] <-vars[snp_pos[snp1]] #SNP1 name check
    results_matrix[snp,2] <-vars[snp_pos[snp2]] #SNP2 name check
    results_matrix[snp,3] <-ifelse(is.null(a),NA,a$tTable[8,1])       #beta_SNP1
    results_matrix[snp,4] <-ifelse(is.null(a),NA,a$tTable[8,2])      #se_SNP1
    results_matrix[snp,5]<-ifelse(is.null(a),NA,a$tTable[8,5])       #P_SNP1
    
    print(a$tTable[8,5])
    snp = snp + 1
  }
}}

write.csv(results_matrix, file="E:/CSV/SNP_SNP_AGE_incomplete.csv", row.names=FALSE)

