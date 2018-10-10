#######################################
### OLS, CQR and MR analysis in UKB ###
#######################################
{library(data.table); library(caret)
  
  #Load and clean data
  cream =fread("D:/CQR.raw")
  master = fread("D:/master_true.txt")
  final = merge(master, cream, "FID")
  final = final[,c(3:13, 19:22, 55:204)]; final = final[, -162]
  
  # This will be used for plotting SNP names
  stopwords = c("_A","_C","_T","_G")
  names(final) = gsub(paste0(stopwords, collapse = "|"), "", colnames(final))
  
  flds = createFolds(final$avMSE, k = 10, list = TRUE, returnTrain = FALSE)
  training = final[-flds[[1]],]; testing = final[flds[[1]],]
  
  # Define variables 
  var_names   = names(training)
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
  
  results = as.data.frame(matrix(nrow = num_snps, ncol = 3))
  names(results) = c("Beta_SNP", "Beta_GxE", "SNP")
  
  #####################
  ### Main analysis ###
  #####################
  for (snp in 1:num_snps){
    name1 = names(training)[snp+15]
    
    # Quantile regression
    myform = as.formula(paste("avMSE ~ ", vars[snp_pos[snp]], "*UniEdu+Sex+Age+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC9+PC10", sep = ""))
    mod_ols = summary(lm(myform, training))
    
    results[snp,1] = mod_ols$coefficients[2,1]
    results[snp,2] = mod_ols$coefficients[17,1]
    results[snp,3] = name1
    
  }
}

library(tidyverse)

  score_marginal = as.data.frame(matrix(nrow = nrow(testing), ncol = 1))
  names(score_marginal) = "PRS"
  
  for(i in 1:nrow(testing)){
    PRS = 0
    for(j in 16:ncol(testing)){
      
      name = names(testing)[j]
      beta = results$Beta_SNP[which(results$SNP == name)]
      weighted = beta*testing[i,name,with=F]
      if(is.na(weighted) == F){
        PRS = PRS + weighted
      }
    }
    score_marginal[i,1] = PRS
  }
  
  testing$PRS_marginal = score_marginal$PRS
  
  

  score_GxE = as.data.frame(matrix(nrow = nrow(testing), ncol = 1))
  names(score_GxE) = "PRS"
  
  for(i in 1:nrow(testing)){
    PRS = 0
    for(j in 16:(ncol(testing)-1)){
      name = names(testing)[j]
      beta = results$Beta_GxE[which(results$SNP == name)]
      weighted = beta*testing[i,name,with=F]
      if(is.na(weighted) == F){
        PRS = PRS + weighted
      }
    }
    score_GxE[i,1] = PRS
  }
  
  testing$PRS_GxE = score_GxE$PRS
  
  cor.test(testing$PRS_marginal, testing$avMSE)
  simple = lm(avMSE~PRS_marginal,testing)
  summary(simple)
  mean(testing$PRS_marginal); sd(testing$PRS_marginal)
  
  par(mfrow=c(1,3))
  hist(testing$PRS_marginal, breaks = 100)
  hist(testing$PRS_GxE, breaks = 100)
  
  cor.test(testing$PRS_GxE, testing$avMSE)
  mean(testing$PRS_GxE); sd(testing$PRS_GxE)
  simple = lm(avMSE~PRS_GxE,testing)
  summary(simple)

  testing$PRS_combined = testing$PRS_marginal + testing$PRS_GxE
  hist(testing$PRS_combined, breaks = 100)
  cor.test(testing$PRS_combined, testing$avMSE)
  mean(testing$PRS_combined); sd(testing$PRS_combined)
  simple = lm(avMSE~PRS_combined,testing)
  summary(simple)
  














