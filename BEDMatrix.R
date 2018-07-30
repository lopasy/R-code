### CQR GWAS
library(BEDMatrix); library(lawstat)
library(quantreg); library(metafor); library(ggplot2); library(ggpubr)
library(tidyr); library(foreach); library(doParallel); library(doSNOW)

missing = read.table("D:/CQR GWAS/missing.txt", quote="\"", comment.char="")
path = "F:/Full_ukb/Output/true"

ref = 3
num_snps = missing
qs = 1:19/20; qss = qs^2
num_boots = 1000
qr_df = as.data.frame(matrix(nrow = length(qs), ncol = 2))
names(qr_df) = c("tau","beta"); qr_df$tau = qs


m = BEDMatrix(path)
#num_snps = dim(m)[2]
master_true = read.csv("D:/master_true.txt", sep="")

{results        = as.data.frame(matrix(nrow = length(num_snps), ncol = 16))
  names(results) = c("SNP", "Beta_OLS", "LCI_OLS", "UCI_OLS", "P_OLS", "Linear_CQR", "Linear_P_CQR", 
                     "Quadratic_CQR", "Quadratic_P_CQR", "Linear_LCI_CQR", "Linear_UCI_CQR", 
                     "Quadratic_LCI_CQR", "Quadratic_UCI_CQR", "P_Levene_Observed", "P_Levene_Permuted", "MAF")
  
  df2 = as.data.frame(matrix(nrow = length(qs)*length(num_snps), ncol = 7))
  names(df2) = c("qs","cqr_int","cqr_slope","cqr_lci","cqr_uci", "N", "SNP")
  
  mr = as.data.frame(matrix(nrow = num_snps, ncol = 13))
  names(mr) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1", "BETA_CQR2", "P_CQR2", "LCI_CQR2", "UCI_CQR2")
  
  mr2 = as.data.frame(matrix(nrow = num_snps, ncol = 9))
  names(mr2) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1")
  
  cqr_res = 1:length(qs); cqr_place = 1; next_mr = 1}

e <- simpleError("test error")

for(i in num_snps){
  final = cbind(master_true, m[,i])
  name = paste(colnames(m)[i], sep ="")
  names(final)[50] = name
  myform = as.formula(paste("avMSE_jitter ~ ", name, " + UniEdu + Sex + poly(Age,2) + Array + 
                            PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", sep = ""))
  cat("snp",i,":",name, "\n")
  
  ##############################################################################################################
  
  # True CQR
  quantile_model = tryCatch(rq(formula = myform, data = final, tau = qs, method = "fn"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  if(length(quantile_model)>1){
    mod_sum = tryCatch(summary(quantile_model), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  if(isTRUE(class(mod_sum)=="NULL")) { next } 
  
  else
    # Get standard errors
    j = 1
  ses = as.data.frame(matrix(ncol = 2, nrow = length(qs)))
  rownames(ses) = qs
  for (i in 1:length(qs)) {
    mycoefs = mod_sum[[i]]
    ses[j,1] = mycoefs$coefficients[2,1]
    ses[j,2] = mycoefs$coefficients[2,2]
    j = j + 1
  }
  
  # CQR results
  df2[cqr_res,1] = qs*100
  df2[cqr_res,2] = quantile_model$coef[1,]
  df2[cqr_res,3] = quantile_model$coef[2,]
  df2[cqr_res,6] = length(quantile_model$residuals)/length(qs)
  df2[cqr_res,7] = name
  
  cqr_res = cqr_res + length(qs)
  
  for (snp in 1:length(qs)) {
    mycoefs = mod_sum[[snp]]
    df2[cqr_place,4] = mycoefs$coefficients[2,1] - (1.96*mycoefs$coefficients[2,2])
    df2[cqr_place,5] = mycoefs$coefficients[2,1] + (1.96*mycoefs$coefficients[2,2])
    cqr_place = cqr_place + 1
  }
  
  ##############################################################################################################
  
  # True MR
  res = try(rma(yi = ses[,1], sei = ses[,2], mods = ~ qs + qss, data = ses[1:length(qs),], control=list(stepadj=0.5)))
  
  if(isTRUE(class(res)=="try-error")) { next } 
  
  else
    
    # MR estimates
    mr[next_mr,2] = round(res$beta[1],3);mr[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr[next_mr,4] = round(res$ci.lb[1],3); mr[next_mr,5] = round(res$ci.ub[1],3)
  mr[next_mr,6] = round(res$beta[2],3);mr[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr[next_mr,8] = round(res$ci.lb[2],3); mr[next_mr,9] = round(res$ci.ub[2],3)
  mr[next_mr,10] = round(res$beta[3],3); mr[next_mr,11] = as.numeric(formatC(res$pval[3], digits = 2, format = "e")); mr[next_mr,12] = round(res$ci.lb[3],3); mr[next_mr,13] = round(res$ci.ub[3],3)
  mr[next_mr,1] = name
  
  # MR on 0.15, 0.3, 0.5 quantiles
  short = c("0.15", "0.3", "0.5"); reduced = c(0.15, 0.3, 0.5)
  res = rma(yi = ses[short,1], sei = ses[short,2], mods = ~ reduced, data = ses[short,])
  
  mr2[next_mr,2] = round(res$beta[1],3);mr2[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr2[next_mr,4] = round(res$ci.lb[1],3); mr2[next_mr,5] = round(res$ci.ub[1],3)
  mr2[next_mr,6] = round(res$beta[2],3);mr2[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr2[next_mr,8] = round(res$ci.lb[2],3); mr2[next_mr,9] = round(res$ci.ub[2],3)
  mr2[next_mr,1] = name
  
  next_mr = next_mr + 1
  
  ##############################################################################################################
  
  # True OLS
  results[i,1]  = name
  
  # MAF
  x             = sum(final[, name], na.rm = TRUE) / (2*sum(!is.na(final[,name])))
  results[i,16] = as.numeric(formatC(min(x, (1-x)), digits = 2, format = "f"))
  
  mod_sumM = summary(lm(formula = myform, data = final))
  beta_ols     = mod_sumM$coefficients[2,1]; se = mod_sumM$coefficients[2,2]
  lci      = beta_ols - (1.96*se); uci = beta_ols + (1.96*se)
  
  # Check if less than 50 copies
  copies        = min(table(final[, name]))
  if(copies > 50){
    results[i,2] = as.numeric(formatC(beta_ols, digits = 3, format = "f"))
    results[i,3] = as.numeric(formatC(lci, digits = 3, format = "f"))
    results[i,4] = as.numeric(formatC(uci, digits = 3, format = "f"))
    results[i,5] = as.numeric(formatC(mod_sumM$coefficients[2,4], digits = 2, format = "e"))
  } # End of if
  
  #qr_df$beta = beta_cqr
  #mod_lm     = lm(beta ~ poly(tau, 2), data = qr_df)
  #linear     = summary(mod_lm)$coefficients[2,1]
  #linear_p   = summary(mod_lm)$coefficients[2,4]
  #quadrat    = summary(mod_lm)$coefficients[3,1]
  #quadrat_p  = summary(mod_lm)$coefficients[3,4]
  
  #results[i,6] = as.numeric(formatC(linear,digits = 3, format = "f"))
  #results[i,7] = as.numeric(formatC(linear_p,digits = 3, format = "e"))
  #results[i,8] = as.numeric(formatC(quadrat,digits = 3, format = "f"))
  #results[i,9] = as.numeric(formatC(quadrat_p,digits = 3, format = "e"))
  #results[i,10] = as.numeric(formatC(summary(mod_lm)$coefficients[2,1] - 1.96*summary(mod_lm)$coefficients[2,2],digits = 3, format = "f"))
  #results[i,11] = as.numeric(formatC(summary(mod_lm)$coefficients[2,1] + 1.96*summary(mod_lm)$coefficients[2,2],digits = 3, format = "f"))
  #results[i,12] = as.numeric(formatC(summary(mod_lm)$coefficients[3,1] - 1.96*summary(mod_lm)$coefficients[3,2],digits = 3, format = "f"))
  #results[i,13] = as.numeric(formatC(summary(mod_lm)$coefficients[3,1] + 1.96*summary(mod_lm)$coefficients[3,2],digits = 3, format = "f"))
  
  # Observed Levene's test
  final4 = final[!is.na(final$avMSE),]
  mod_lev = levene.test(final4$avMSE, group = as.factor(final4[,name]), location = "mean", bootstrap = num_boots, correction.method = "none")
  results[i,14] = as.numeric(formatC(mod_lev$p.value, digits = 3, format = "e"))
}

write.csv(results, file = paste("D:/CQR GWAS/results_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr, file = paste("D:/CQR GWAS/mr_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr2, file = paste("D:/CQR GWAS/mr2_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(df2, file = paste("D:/CQR GWAS/df2_", ref, ".csv", sep = ""), row.names = F, quote = F)












