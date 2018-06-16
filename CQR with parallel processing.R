# Cross-validation for quantile regression
library(quantreg); library(metafor); library(ggplot2)
library(ggpubr); library(tidyr); library(foreach); library(doParallel); library(doSNOW)


master_true = read.csv("D:/master_true.txt", sep="")
CQR = read.csv("D:/CQR.raw", sep="")
final = merge(master_true, CQR, 1); final = final[,c(3:13,19:22,53:202)]; final = final[,-162]

ppi = 400
qs = 1:19/20

# Define variables 
var_names   = names(final)
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

# Recoding risk alleles in CREAM to match those in UKB
{data1        = read.csv(file = "D:/cream2017_ukbb_replicated.csv", header = TRUE)
  data2        = data1[, c("ukbSNP", "ukbCHR", "ukbPOS", "GENE.1", "ukbA1", "ukbA2", "ukbBETA", "ukbSE", "ukbP_bolt")]
  data2$EA     = ifelse(data2$ukbBETA < 0, as.character(data2$ukbA1), as.character(data2$ukbA2))
  
  data3        = as.data.frame(matrix(nrow = ncol(final)-15, ncol = 1))
  data3$start  = var_names[snp_pos[1]:snp_pos[ncol(final)-15]]; data3$V1 = NULL
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
var_names   = names(final)
vars        = as.matrix(var_names)
num_vars    = nrow(vars)
num_snps    = 0
snp_pos     = matrix(nrow = num_vars, ncol = 1)


for (n in 1:num_vars){
  if (substr(vars[n],1,2) == "rs"){
    num_snps = num_snps + 1
    snp_pos[num_snps] = n
  }}


stopwords = c("_A","_C","_T","_G")
final2 = final
names(final2) = gsub(paste0(stopwords, collapse = "|"), "", colnames(final))

rm(master_true, CQR, stopwords, var_names, num_vars, curr_col, old_snp, n, new_snp)
rm(data1, data2, data3, data5, data6)

# LM vs. CQR for SNP x UNI p-values
qr_GWAS = function(snp, num_boot){
  name = paste(names(final2[snp + 15]), sep ="")
  myform = as.formula(paste("avMSE ~ ", vars[snp_pos[snp]], "+ UniEdu + Sex + poly(Age,2) + Array + 
                            PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +", vars[snp_pos[snp]], ": UniEdu", sep = ""))
  cat("snp",snp,":",vars[snp_pos[snp]], "\n")
  
  mod_lm1            <- summary(lm(myform,final))
  
  cqr_est1           <- as.data.frame(matrix(ncol=19,nrow=18))
  cqr_se1            <- as.data.frame(matrix(ncol=19,nrow=18))
  cqr_p1             <- as.data.frame(matrix(ncol=19,nrow=18))
  
  row.names(cqr_est1)    <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")
  row.names(cqr_se1)     <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")
  row.names(cqr_p1)      <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")
  

    mod_qr                       <- try(rq(formula = myform, data = final, tau = qs))
    if(length(mod_qr)>1){
      mod_sum                      <- try(summary(mod_qr, se = "boot", R = num_boot, bsmethod = "mcmb")) # , se = "boot", R = num_boot, bsmethod = "mcmb"
      if(length(mod_sum)>1){
        
        for (i in 1:length(qs)) {
          mycoefs = mod_sum[[i]]
          cqr_est1[1:18,i] = mycoefs$coefficients[,1]
          cqr_se1[1:18,i] = mycoefs$coefficients[,2]
          cqr_p1[1:18,i] = mycoefs$coefficients[,4]
        }
      }
    }

  cqr_z1  <- as.numeric(cqr_est1[18,1:19]/cqr_se1[18,1:19])
  
  
  df = as.data.frame(matrix(nrow = length(qs),ncol = 0))
  df$qs = qs*100
  df$cqr = -log10(2*pnorm(-abs(cqr_z1)))
  df$ols = -log10(mod_lm1$coef[18,4])
  df$gwas = -log10(0.00000005)
  
  
  myplot = ggplot(data = df, aes(qs)) +
    geom_line(aes(y=cqr),size=1,colour="red") +
    geom_hline(colour="blue", aes(yintercept=ols)) +
    geom_hline(colour="black", aes(yintercept=gwas)) +
    scale_x_continuous(breaks = seq(0,100,by=20)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold", vjust = 1), 
          legend.title = element_text(), legend.position = "none",
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(face = "bold", color = "black", size = 5),
          axis.text.y = element_text(face = "bold", color = "black",  size = 5)) +
    ggtitle(name)
  
  list(cqr_est1[18,],cqr_se1[18,],df$cqr,myplot)
  
}


  
# Use parrallel processing
no_cores = detectCores()
cluster = makeCluster(no_cores-2, type = "SOCK")
registerDoSNOW(cluster)

test = foreach(snp = 100:40, .combine = list, .multicombine = TRUE, .packages = c("ggplot2", "quantreg")) %dopar% qr_GWAS(snp, 10)
stopImplicitCluster()

betas = as.data.frame(matrix(ncol=length(qs),nrow=num_snps)); names(betas) = qs
ses = as.data.frame(matrix(ncol=length(qs),nrow=num_snps)); names(ses) = qs
pvals = as.data.frame(matrix(ncol=length(qs),nrow=num_snps)); names(pvals) = qs
plots = list(1:num_snps)

for(i in 1:length(test)){
  betas[i,] = test[[i]][[1]]
}

for(i in 1:length(test)){
  ses[i,] = test[[i]][[2]]
}

for(i in 1:length(test)){
  pvals[i,] = test[[i]][[3]]
}

for(i in 1:length(test)){
  plots[[i]] = test[[i]][[4]]
}

  

plot_num = 1
for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                     plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], 
                     plots[[plot_num + 9]], plots[[plot_num + 10]], plots[[plot_num + 11]], plots[[plot_num + 12]], plots[[plot_num + 13]],
                     plots[[plot_num + 14]], plots[[plot_num + 15]], plots[[plot_num + 16]], plots[[plot_num + 17]], plots[[plot_num + 18]],
                     plots[[plot_num + 19]],
                     ncol = 4, nrow = 5)
  file_out = paste("D:/", "CQR_GWAS", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Refractive error percentile", color = "black", size = 10),
                             left = text_grob("Negative log10 (P-value))", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 20
}












