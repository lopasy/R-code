### CQR GWAS
library(BEDMatrix); library(lawstat); library(rms)
library(quantreg); library(metafor); library(ggplot2); library(ggpubr)
library(tidyr); library(foreach); library(doParallel); library(doSNOW)


ref = "AgeofOnset"
qs = 1:19/20; qss = qs^2; qsss = qs^3; qssss = qs^4
num_boots = 10000

master_true = read.csv("./cnv_all_details_included.txt", sep="")
Euro = read.table("./ukb_v2_10xSD_PC-derived_Europeans.keep", quote="\"", comment.char="")

master_true = master_true[which(master_true$FID %in% Euro$V1),]
master_true = master_true[which(master_true$het_miss_out == 0),]
master_true = master_true[which(master_true$in_kinship == 0),]
master_true = master_true[!is.na(master_true$AgeSpexWear),]

t1 = read.csv("./AoO_cluster1.txt", sep="")
master_true = master_true[which(master_true$FID %in% t1$FID),]

m = read.csv("./AoO.raw", sep="")
final = m[,c(1, 7:152,154:156)]; final = final[which(final$FID %in% master_true$FID),]; final = final[,-1]
num_snps = dim(final)[2]


# Define variables 
var_names   = names(final); vars = as.matrix(var_names)
num_vars    = nrow(vars); num_snps = 0
snp_pos     = matrix(nrow = num_vars, ncol = 1)

for (n in 1:num_vars){
  if (substr(vars[n],1,2) == "rs"){
    num_snps = num_snps + 1
    snp_pos[num_snps] = n
  }
}

# Recoding risk alleles in CREAM to match those in UKB
{data1        = read.csv(file = "./cream2017_ukbb_replicated.csv", header = TRUE)
  data2        = data1[, c("ukbSNP", "ukbCHR", "ukbPOS", "GENE.1", "ukbA1", "ukbA2", "ukbBETA", "ukbSE", "ukbP_bolt")]
  data2$EA     = ifelse(data2$ukbBETA < 0, as.character(data2$ukbA1), as.character(data2$ukbA2))
  
  data3        = as.data.frame(matrix(nrow = ncol(final), ncol = 1))
  data3$start  = var_names[snp_pos[1]:snp_pos[ncol(final)]]; data3$V1 = NULL
  data5        = separate(data = data3, col = start, into = c("ukbSNP", "plinkA1"), sep = "_")
  data6        = merge(data2, data5, "ukbSNP")
  data6$switch = ifelse(data6$EA == data6$plinkA1,0,1)}


# Ensure SNP is coded so that test allele produces a negative beta
for (n in 1:num_snps){
  if(data6$switch[n] == 1){
    old_snp  = paste(data6$ukbSNP[n],"_", data6$plinkA1[n], sep="")
    new_snp  = paste(data6$ukbSNP[n],"_", data6$EA[n], sep="")
    curr_col = which(colnames(final) == old_snp)
    final[, curr_col] = 2 - final[, curr_col]
    names(final)[curr_col] = new_snp
  }
}

# collect SNP names again after re-coding
var_names   = names(final); vars = as.matrix(var_names)
num_vars    = nrow(vars); num_snps = 0
snp_pos     = matrix(nrow = num_vars, ncol = 1)


for (n in 1:num_vars){
  if (substr(vars[n],1,2) == "rs"){
    num_snps = num_snps + 1
    snp_pos[num_snps] = n
  }}

{results = as.data.frame(matrix(nrow = length(num_snps), ncol = 7))
  names(results) = c("SNP", "Beta_OLS", "LCI_OLS", "UCI_OLS", "P_OLS", "P_Levene_Observed", "MAF")
  
  df2 = as.data.frame(matrix(nrow = length(qs)*length(num_snps), ncol = 22))
  names(df2) = c("qs","cqr_int","cqr_slope","cqr_lci","cqr_uci", "N", "SNP", "spline_pred", "spline_ci.lb", "spline_ci.ub",
                 "left_pred", "left_ci.lb", "left_ci.ub", "right_pred", "right_ci.lb", "right_ci.ub",
                 "left_quad_pred", "left_quad_ci.lb", "left_quad_ci.ub", "right_quad_pred", "right_quad_ci.lb", "right_quad_ci.ub")
  
  mr = as.data.frame(matrix(nrow = num_snps, ncol = 17))
  names(mr) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1",
                "BETA_CQR2", "P_CQR2", "LCI_CQR2", "UCI_CQR2", "BETA_CQR3", "P_CQR3", "LCI_CQR3", "UCI_CQR3")
  
  mr_left = as.data.frame(matrix(nrow = num_snps, ncol = 9))
  names(mr_left) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1")
  
  mr_right = as.data.frame(matrix(nrow = num_snps, ncol = 9))
  names(mr_right) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1")
  
  mr_left_quad = as.data.frame(matrix(nrow = num_snps, ncol = 13))
  names(mr_left_quad) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1", "BETA_CQR2", "P_CQR2", "LCI_CQR2", "UCI_CQR2")
  
  mr_right_quad = as.data.frame(matrix(nrow = num_snps, ncol = 13))
  names(mr_right_quad) = c("SNP", "Beta_MR", "P_MR", "LCI_MR", "UCI_MR", "BETA_CQR1", "P_CQR1", "LCI_CQR1", "UCI_CQR1", "BETA_CQR2", "P_CQR2", "LCI_CQR2", "UCI_CQR2")
  
  spline = list()
  
  cqr_res = 1:length(qs); cqr_place = 1; next_mr = 1; cqr_left = 1:8; cqr_right = 8:length(qs)}

e <- simpleError("test error")

m = final

for(i in 1:num_snps){
  final = cbind(master_true, m[,i])
  name = paste(colnames(m)[i], sep =""); names(final)[118] = name
  myform = as.formula(paste("AgeSpexWear ~ ", name, "* UniEdu +", name, " + UniEdu + Sex + Age + Geno_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", sep = ""))
  cat("snp",i,":",name, "\n")
  
  ##############################################################################################################
  
  # True CQR
  quantile_model = tryCatch(rq(formula = myform, data = final, tau = qs, method = "fn"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  if(isTRUE(class(quantile_model)=="NULL")) { next } 
  
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
    ses[j,1] = mycoefs$coefficients[2,1]; ses[j,2] = mycoefs$coefficients[2,2]
    j = j + 1
  }
  
  # CQR results
  df2[cqr_res,1] = qs*100
  df2[cqr_res,2] = quantile_model$coef[1,]; df2[cqr_res,3] = quantile_model$coef[2,]
  df2[cqr_res,6] = length(quantile_model$residuals)/length(qs)
  df2[cqr_res,7] = name
  
  
  for (snp in 1:length(qs)) {
    mycoefs = mod_sum[[snp]]
    df2[cqr_place,4] = mycoefs$coefficients[2,1] - (1.96*mycoefs$coefficients[2,2])
    df2[cqr_place,5] = mycoefs$coefficients[2,1] + (1.96*mycoefs$coefficients[2,2])
    cqr_place = cqr_place + 1
  }
  
  ##############################################################################################################
  
  # True MR
  knots = c(0.1,0.3,0.4,0.5,0.9)
  res = try(rma(ses[,1]~rcs(qs, knots), sei = ses[,2], data=ses, control=list(stepadj=0.5)))
  
  spline[[next_mr]] = smooth.spline(ses[,1]~qs, w = ses[,2], cv = T)
  
  predicted = predict(res)
  df2[cqr_res,8] = predicted$pred; df2[cqr_res,9] = predicted$ci.lb; df2[cqr_res,10] = predicted$ci.ub
  #df2[cqr_res,12] = predicted$cr.lb; df2[cqr_res,13] = predicted$cr.ub
  
  
  if(isTRUE(class(res)=="try-error")) { next } 
  
  else
    
    # MR estimates
    mr[next_mr,2] = round(res$beta[1],3);mr[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr[next_mr,4] = round(res$ci.lb[1],3); mr[next_mr,5] = round(res$ci.ub[1],3)
  mr[next_mr,6] = round(res$beta[2],3);mr[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr[next_mr,8] = round(res$ci.lb[2],3); mr[next_mr,9] = round(res$ci.ub[2],3)
  mr[next_mr,10] = round(res$beta[3],3); mr[next_mr,11] = as.numeric(formatC(res$pval[3], digits = 2, format = "e")); mr[next_mr,12] = round(res$ci.lb[3],3); mr[next_mr,13] = round(res$ci.ub[3],3)
  mr[next_mr,14] = round(res$beta[3],3); mr[next_mr,15] = as.numeric(formatC(res$pval[3], digits = 2, format = "e")); mr[next_mr,16] = round(res$ci.lb[3],3); mr[next_mr,17] = round(res$ci.ub[3],3)
  mr[next_mr,1] = name
  
  # MR for left AoO side
  short = c("0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4"); reduced = c(0.05, 0.1, 0.15,0.2,0.25,0.3,0.35,0.4)
  res = rma(yi = ses[short,1], sei = ses[short,2], mods = ~ reduced, data = ses[short,])
  
  predicted = predict(res)
  df2[cqr_left,11] = predicted$pred; df2[cqr_left,12] = predicted$ci.lb; df2[cqr_left,13] = predicted$ci.ub
  
  mr_left[next_mr,2] = round(res$beta[1],3);mr_left[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr_left[next_mr,4] = round(res$ci.lb[1],3); mr_left[next_mr,5] = round(res$ci.ub[1],3)
  mr_left[next_mr,6] = round(res$beta[2],3);mr_left[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr_left[next_mr,8] = round(res$ci.lb[2],3); mr_left[next_mr,9] = round(res$ci.ub[2],3)
  mr_left[next_mr,1] = name
  
  # MR for right AoO side
  short = c("0.4","0.45", "0.5", "0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95") 
  reduced = c(0.4,0.45, 0.5, 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  res = rma(yi = ses[short,1], sei = ses[short,2], mods = ~ reduced, data = ses[short,])
  
  predicted = predict(res)
  df2[cqr_right,14] = predicted$pred; df2[cqr_right,15] = predicted$ci.lb; df2[cqr_right,16] = predicted$ci.ub
  
  mr_right[next_mr,2] = round(res$beta[1],3);mr_right[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr_right[next_mr,4] = round(res$ci.lb[1],3); mr_right[next_mr,5] = round(res$ci.ub[1],3)
  mr_right[next_mr,6] = round(res$beta[2],3);mr_right[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr_right[next_mr,8] = round(res$ci.lb[2],3); mr_right[next_mr,9] = round(res$ci.ub[2],3)
  mr_right[next_mr,1] = name
  
  
  # Quadratic MR for left AoO side
  short = c("0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4"); reduced = c(0.05, 0.1, 0.15,0.2,0.25,0.3,0.35,0.4); reduceds = reduced^2
  res = rma(yi = ses[short,1], sei = ses[short,2], mods = ~ reduced + reduceds, data = ses[short,])
  
  predicted = predict(res)
  df2[cqr_left,17] = predicted$pred; df2[cqr_left,18] = predicted$ci.lb; df2[cqr_left,19] = predicted$ci.ub
  
  mr_left_quad[next_mr,2] = round(res$beta[1],3);mr_left_quad[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr_left_quad[next_mr,4] = round(res$ci.lb[1],3); mr_left_quad[next_mr,5] = round(res$ci.ub[1],3)
  mr_left_quad[next_mr,6] = round(res$beta[2],3);mr_left_quad[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr_left_quad[next_mr,8] = round(res$ci.lb[2],3); mr_left_quad[next_mr,9] = round(res$ci.ub[2],3)
  mr_left_quad[next_mr,10] = round(res$beta[3],3); mr_left_quad[next_mr,11] = as.numeric(formatC(res$pval[3], digits = 2, format = "e")); mr_left_quad[next_mr,12] = round(res$ci.lb[3],3); mr_left_quad[next_mr,13] = round(res$ci.ub[3],3)
  mr_left_quad[next_mr,1] = name
  
  # Quadratic MR for right AoO side
  short = c("0.4","0.45", "0.5", "0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95") 
  reduced = c(0.4,0.45, 0.5, 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95); reduceds = reduced^2
  res = rma(yi = ses[short,1], sei = ses[short,2], mods = ~ reduced + reduceds, data = ses[short,])
  
  predicted = predict(res)
  df2[cqr_right,20] = predicted$pred; df2[cqr_right,21] = predicted$ci.lb; df2[cqr_right,22] = predicted$ci.ub
  
  mr_right_quad[next_mr,2] = round(res$beta[1],3);mr_right_quad[next_mr,3] = as.numeric(formatC(res$pval[1], digits = 2, format = "e")); mr_right_quad[next_mr,4] = round(res$ci.lb[1],3); mr_right_quad[next_mr,5] = round(res$ci.ub[1],3)
  mr_right_quad[next_mr,6] = round(res$beta[2],3);mr_right_quad[next_mr,7] = as.numeric(formatC(res$pval[2], digits = 2, format = "e")); mr_right_quad[next_mr,8] = round(res$ci.lb[2],3); mr_right_quad[next_mr,9] = round(res$ci.ub[2],3)
  mr_right_quad[next_mr,10] = round(res$beta[3],3); mr_right_quad[next_mr,11] = as.numeric(formatC(res$pval[3], digits = 2, format = "e")); mr_right_quad[next_mr,12] = round(res$ci.lb[3],3); mr_right_quad[next_mr,13] = round(res$ci.ub[3],3)
  mr_right_quad[next_mr,1] = name
  
  cqr_res = cqr_res + length(qs)
  cqr_left = cqr_left + length(qs)
  cqr_right = cqr_right + length(qs)
  
  ##############################################################################################################
  
  # True OLS
  results[next_mr,1] = name
  
  # MAF
  x = sum(final[, name], na.rm = TRUE) / (2*sum(!is.na(final[,name])))
  results[next_mr,7] = as.numeric(formatC(min(x, (1-x)), digits = 4, format = "f"))
  
  mod_sumM = summary(lm(formula = myform, data = final))
  beta_ols = mod_sumM$coefficients[17,1]; se = mod_sumM$coefficients[17,2]
  lci      = beta_ols - (1.96*se); uci = beta_ols + (1.96*se)
  
  # Check if less than 50 copies
  copies        = min(table(final[, name]))
  if(copies > 50){
    results[next_mr,2] = as.numeric(formatC(beta_ols, digits = 3, format = "f"))
    results[next_mr,3] = as.numeric(formatC(lci, digits = 3, format = "f"))
    results[next_mr,4] = as.numeric(formatC(uci, digits = 3, format = "f"))
    results[next_mr,5] = as.numeric(formatC(mod_sumM$coefficients[17,4], digits = 2, format = "e"))
  } # End of if
  
  # Observed Levene's test
  final4 = final[!is.na(final$AgeSpexWear),]
  mod_lev = levene.test(final4$AgeSpexWear, group = as.factor(final4[,name]), location = "mean", bootstrap = num_boots, correction.method = "none")
  results[next_mr,6] = as.numeric(formatC(mod_lev$p.value, digits = 2, format = "e"))
  
  next_mr = next_mr + 1
}

write.csv(results, file = paste("D:/CQR GWAS/results_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr, file = paste("D:/CQR GWAS/mr_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr_left, file = paste("D:/CQR GWAS/mr_left_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr_right, file = paste("D:/CQR GWAS/mr_right_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr_left_quad, file = paste("D:/CQR GWAS/mr_left_quad_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(mr_right_quad, file = paste("D:/CQR GWAS/mr_right_quad_int_", ref, ".csv", sep = ""), row.names = F, quote = F)
write.csv(df2, file = paste("D:/CQR GWAS/df2_int_", ref, ".csv", sep = ""), row.names = F, quote = F)


single_cqr_spline_noedu = spline












