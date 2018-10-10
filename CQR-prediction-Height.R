#######################################
### OLS, CQR and MR analysis in UKB ###
#######################################
{
  library(quantreg); library(data.table); library(caret)
  
  #Load and clean data
  cream =fread("D:/height.raw")
  master = fread("D:/height_pheno.txt")
  final = merge(master, cream, "FID")
  final = final[,c(3:17, 23:173)]; final = final[, -146]
  final = final[!is.na(final$Height),]
  
  # This will be used for plotting SNP names
  stopwords = c("_A","_C","_T","_G")
  names(final) = gsub(paste0(stopwords, collapse = "|"), "", colnames(final))
  
  flds = createFolds(final$Height, k = 10, list = TRUE, returnTrain = FALSE)
  training = final[-flds[[7]],]; testing = final[flds[[7]],]
  
  
  # Miscellaneous
  qs = seq(10,90,10)/100; qss = qs^2
  c = 1:length(qs); j = 1; rm(cream, master)
  
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
  
  # Create variable were results will be stored
  df2 = as.data.frame(matrix(nrow = num_snps*length(qs), ncol = 5))
  names(df2) = c("x","intercept","slope","se", "SNP")
  df2$x = rep(qs)
  
  results = as.data.frame(matrix(nrow = num_snps, ncol = 2))
  names(results) = c("Beta_OLS", "SNP")
  
  #####################
  ### Main analysis ###
  #####################
  for (snp in 1:num_snps){
    name1 = names(training)[snp+15]
    
    # Quantile regression
    myform = as.formula(paste("Height ~ ", vars[snp_pos[snp]], " + Sex + poly(Age,2) + Array", sep = ""))
    mod_qr = tryCatch(summary(rq(formula = myform, data = training, tau = qs, method = "fn")), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    mod_ols = summary(lm(myform, training))
    
    results[snp,1] = mod_ols$coefficients[2,1]
    results[snp,2] = name1
    
    if(isTRUE(class(mod_qr)=="NULL")) { next } 
    else
      
      cat("snp",snp,":",name1, "\n")
    
    a = 1
    jez = as.data.frame(matrix(ncol = 2, nrow = length(qs)))
    for (i in 1:length(qs)) {
      mycoefs = mod_qr[[i]]
      jez[a,1] = mycoefs$coefficients[2,1]
      jez[a,2] = mycoefs$coefficients[2,2]
      a = a + 1
    }
    
    qr2 = rq(formula = myform, data = training, tau = qs, method = "fn")
    df2[c,2] = qr2$coef[1,]
    df2[c,3] = qr2$coef[2,]
    df2[c,5] = name1
    
    for (i in 1:length(qs)) {
      mycoefs = mod_qr[[i]]
      df2[j,4] = mycoefs$coefficients[2,2]
      j = j + 1
    }
    
    c = c + length(qs)
    
  }
  gc()
}

{hist(training$Height, breaks = 100)
hist(testing$Height, breaks = 100)

library(tidyverse)
testing = testing %>% mutate(quantile = ntile(Height, length(qs)))
testing$quantile = testing$quantile/10


score = as.data.frame(matrix(nrow = nrow(testing), ncol = 1))
names(score) = "PRS"

for(i in 1:nrow(testing)){
  
  PRS = 0
  
  for(j in 16:(ncol(testing)-1)){
    
    name = names(testing)[j]
    quant = testing[i,ncol(testing)]
    beta = df2$slope[which(df2$SNP == name & df2$x == quant)]
    if(length(beta > 1)){
    weighted = beta*testing[i,name]
    if(is.na(weighted) == TRUE){ next }
    else
      PRS = PRS + weighted
    }
  }
  
  score[i,1] = PRS
}

testing$PRS = score$PRS

hist(testing$PRS, breaks = 100)
cor.test(testing$PRS, testing$Height)
simple = lm(Height~PRS,testing)
#plot(simple)
summary(simple)}

mean(testing$PRS); sd(testing$PRS)

{score_ols = as.data.frame(matrix(nrow = nrow(testing), ncol = 1))
names(score_ols) = "PRS"

for(i in 1:nrow(testing)){
  
  PRS = 0
  
  for(j in 16:(ncol(testing)-1)){
    
    name = names(testing)[j]
    beta = results$Beta_OLS[which(results$SNP == name)]
    weighted = beta*testing[i,name]
    if(length(weighted) == 1){
      PRS = PRS + weighted
    }
  }
  
  score_ols[i,1] = PRS
}

testing$PRS_ols = score_ols$PRS
hist(testing$PRS_ols, breaks = 100)
simple = lm(Height~PRS,testing)
#plot(simple)
summary(simple)}



















