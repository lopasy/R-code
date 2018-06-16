# Cross-validation for quantile regression
library(quantreg); library(metafor); library(ggplot2); library(ggpubr)
master_true = read.csv("~/Documents/Useless/master_true.txt", sep="")
CQR = read.csv("~/Documents/CQR_17_04_2018/CQR.raw", sep="")
final = merge(master_true, CQR, 1); final = final[,c(3:13,19:22,53:202)]; final = final[,-162]


num_cv = 5
num_boot = 10
#left_out_sample1 = 1:14597; left_out_sample2 = 14598:29194
#left_out_sample3 = 29195:43791; left_out_sample4 = 43792:58388
#left_out_sample5 = 58389:72985


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

plots = list()

# LM vs. CQR for SNP x UNI p-values
for(snp in 1:num_snps){
  
myform = as.formula(paste("avMSE ~ ", vars[snp_pos[snp]], "+ UniEdu + Sex + poly(Age,2) + Array + 
                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +", vars[snp_pos[snp]], ": UniEdu", sep = ""))
  
mod_lm1            <- summary(lm(myform,final))

cqr_est1           <- as.data.frame(matrix(ncol=19,nrow=18))
cqr_se1            <- as.data.frame(matrix(ncol=19,nrow=18))
cqr_p1             <- as.data.frame(matrix(ncol=19,nrow=18))

row.names(cqr_est1)    <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")
row.names(cqr_se1)     <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")
row.names(cqr_p1)      <- c("Intercept","SNP","UniEdu","Sex","Age","Age2","Array","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SNPxEdu")

for (i in 1:19){
  q                            <- i/20
  mod_qr                       <- 0
  mod_qr                       <- try(rq(formula = myform, data = final, tau = q))
  if(length(mod_qr)>1){
    mod_sum                      <- 0
    mod_sum                      <- try(summary(mod_qr))
    if(length(mod_sum)>1){
      cqr_est1[1:18,i]               <- mod_sum$coefficients[,1]
      cqr_se1[1:18,i]                <- mod_sum$coefficients[,2]
      cqr_p1[1:18,i]                 <- mod_sum$coefficients[,4]
    } # end try if
  } # end try if
} # next quantile
cqr_p1

cqr_z1  <- as.numeric(cqr_est1[8,1:19]/cqr_se1[8,1:19])
cqr_zp1 <- 2*pnorm(-abs(cqr_z1))

myplot <- ggplot()+
  geom_line(aes(x=qs,y=-log10(as.numeric(cqr_zp1))),size=1,colour="red") +
  geom_point(aes(x=qs,y=-log10(as.numeric(cqr_zp1))),size=1,colour="red") +
  geom_hline(colour="blue", aes(yintercept=-log10(mod_lm1$coef[8,4]))) +
  scale_x_continuous(limits=c(0.05,0.95), breaks=seq(0.1,0.9,by=0.2), labels=seq(0.1,0.9,by=0.2)) +
  labs(x="Refractive error percentile",y="Negative log10 (P-value)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold", vjust = -6,
                                  margin = margin(t = -10, r = 20, b = 20, l = 0)), 
        legend.title = element_text(), legend.position = "none",
        axis.text.x = element_text(face = "bold", color = "black", size = 5),
        axis.text.y = element_text(face = "bold", color = "black",  size = 5)) +
  ggtitle(vars[snp_pos[snp]])

 plots[[plot_num]] = myplot
 #snp_title = snp_title + 1
 plot_num = plot_num + 1
 
}

plot_num = 1
for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                   plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], ncol = 3, nrow = 3)
  file_out = paste("/scratch/share_PR300/Alfred/Files/CQR/", "avMSE_plot", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 5*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Refractive error percentile", color = "black", size = 10),
                left = text_grob("Genetic effect size (diopters per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 9
}












