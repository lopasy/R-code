library(data.table)
{chr = list.files(pattern = "df2_A.", full.names = TRUE); chr = lapply(chr, fread)
df2 = do.call(rbind.data.frame, chr)

chr = list.files(pattern = "results_A.", full.names = TRUE); chr = lapply(chr, fread)
results = do.call(rbind.data.frame, chr)

chr = list.files(pattern = "mr_A.", full.names = TRUE); chr = lapply(chr, fread)
mr = do.call(rbind.data.frame, chr)


target = as.data.frame(results$SNP)
data1 = read.csv("D:/CQR GWAS/cream2017_ukbb_replicated.csv"); data2 = data1[,c(9,5)]

stopwords = c("_A","_C","_T","_G")
target[,1] = gsub(paste0(stopwords, collapse = "|"), "", target[,1])
df2$SNP = gsub(paste0(stopwords, collapse = "|"), "", df2$SNP)
results$SNP = gsub(paste0(stopwords, collapse = "|"), "", results$SNP)
mr$SNP = gsub(paste0(stopwords, collapse = "|"), "", mr$SNP)
snps = merge(target,data2,1)


{df2$beta_ols = 0; df2$lci_ols = 0; df2$uci_ols = 0
df2$beta_mr = 0; df2$lci_mr = 0; df2$uci_mr = 0
df2$cqr_beta_mr = 0; df2$cqr_lci_mr = 0; df2$cqr_uci_mr = 0
df2$cqr2_beta_mr = 0; df2$cqr2_lci_mr = 0; df2$cqr2_uci_mr = 0
df2$Gene = 0}

target = results$SNP

df2 = as.data.frame(df2); results = as.data.frame(results); mr = as.data.frame(mr)}

for(i in target){
  df2[which(df2$SNP == i), 8] = results[which(results$SNP == i), 2]
  df2[which(df2$SNP == i), 9] = results[which(results$SNP == i), 3]
  df2[which(df2$SNP == i), 10] = results[which(results$SNP == i), 4]
  df2[which(df2$SNP == i), 11] = mr[which(mr$SNP == i), 2]
  df2[which(df2$SNP == i), 12] = mr[which(mr$SNP == i), 4]
  df2[which(df2$SNP == i), 13] = mr[which(mr$SNP == i), 5]
  df2[which(df2$SNP == i), 14] = mr[which(mr$SNP == i), 6]
  df2[which(df2$SNP == i), 15] = mr[which(mr$SNP == i), 8]
  df2[which(df2$SNP == i), 16] = mr[which(mr$SNP == i), 9]
  df2[which(df2$SNP == i), 17] = mr[which(mr$SNP == i), 10]
  df2[which(df2$SNP == i), 18] = mr[which(mr$SNP == i), 12]
  df2[which(df2$SNP == i), 19] = mr[which(mr$SNP == i), 13]
}


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
df2$Gene = as.character(rep.row(snps$GENE, 19))


plots = list()
target = results$SNP
titles = 1
plot_num = 1
qs = 1:19

for(i in 1:length(target)){
  name1 = df2[qs, 20][1]
  name2 = paste("(", target[titles], ")", sep ="")
  
myplot = ggplot(data = df2[qs,],aes(qs)) + 
  geom_ribbon(aes(ymin = cqr_lci, ymax = cqr_uci), fill = "grey70", alpha = 0.8) +
  geom_line(aes(y = cqr_slope), colour = "black", size = 0.6) +
  geom_line(aes(y = beta_mr+cqr_beta_mr*(qs/100)+cqr2_beta_mr*(qs/100)^2), colour = "blue", size = 0.4) +
  geom_line(aes(y = lci_mr+cqr_lci_mr*(qs/100)+cqr2_lci_mr*(qs/100)^2), linetype = "longdash", size = 0.4, colour = "blue") + 
  geom_line(aes(y = uci_mr+cqr_uci_mr*(qs/100)+cqr2_uci_mr*(qs/100)^2), linetype = "longdash", size = 0.4, colour = "blue") +
  geom_hline(aes(yintercept = beta_ols), colour = "red", size = 0.4) +
  geom_hline(aes(yintercept = uci_ols), colour = "red", size = 0.4, linetype = 2) +
  geom_hline(aes(yintercept = lci_ols), colour = "red", size = 0.4, linetype = 2) +
  scale_x_continuous(breaks = seq(0,100, 20)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold", vjust = 1), 
        legend.title = element_text(), legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face = "bold", color = "black", size = 5),
        axis.text.y = element_text(face = "bold", color = "black",  size = 5)) +
  ggtitle(paste(name1,"\n", name2))

plots[[plot_num]] = myplot
titles = titles + 1
plot_num = plot_num + 1
qs = qs + 19
}


plot_num = 1
ppi = 400

for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                   plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], 
                   plots[[plot_num + 9]], plots[[plot_num + 10]], plots[[plot_num + 11]], plots[[plot_num + 12]], plots[[plot_num + 13]], 
                   plots[[plot_num + 14]], plots[[plot_num + 15]], ncol = 4, nrow = 4)
  file_out = paste("D:/", "AgeOfOnset_plot", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                left = text_grob("Genetic effect size (diopters per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 16
}



################################################################
################################################################
####################### Age Of Onset ###########################
################################################################
################################################################

results_edu = results_edu[,c(1:5, 14, 16)]
results_noedu = results_noedu[,c(1:5, 14, 16)]
results_int = results_int[,c(1:5, 14, 16)]

test = merge(results_edu, results_noedu, 1)
test = merge(test, results_int, 1)


  target = as.data.frame(test$SNP)
  data1 = read.csv("D:/CQR GWAS/cream2017_ukbb_replicated.csv"); data2 = data1[,c(9,5)]
  
  stopwords = c("_A","_C","_T","_G")
  target[,1] = gsub(paste0(stopwords, collapse = "|"), "", target[,1])
  #df2$SNP = gsub(paste0(stopwords, collapse = "|"), "", df2$SNP)
  #results$SNP = gsub(paste0(stopwords, collapse = "|"), "", results$SNP)
  #mr$SNP = gsub(paste0(stopwords, collapse = "|"), "", mr$SNP)
  snps = merge(target,data2,1)
  
  
#  {df2$beta_ols = 0; df2$lci_ols = 0; df2$uci_ols = 0
#    df2$beta_mr = 0; df2$lci_mr = 0; df2$uci_mr = 0
#    df2$cqr_beta_mr = 0; df2$cqr_lci_mr = 0; df2$cqr_uci_mr = 0
#    df2$cqr2_beta_mr = 0; df2$cqr2_lci_mr = 0; df2$cqr2_uci_mr = 0
#    df2$Gene = 0}
#  
#  target = results$SNP
#  
#  df2 = as.data.frame(df2); results = as.data.frame(results); mr = as.data.frame(mr)

#for(i in target){
#  df2[which(df2$SNP == i), 8] = results[which(results$SNP == i), 2]
#  df2[which(df2$SNP == i), 9] = results[which(results$SNP == i), 3]
#  df2[which(df2$SNP == i), 10] = results[which(results$SNP == i), 4]
#  df2[which(df2$SNP == i), 11] = mr[which(mr$SNP == i), 2]
#  df2[which(df2$SNP == i), 12] = mr[which(mr$SNP == i), 4]
#  df2[which(df2$SNP == i), 13] = mr[which(mr$SNP == i), 5]
#  df2[which(df2$SNP == i), 14] = mr[which(mr$SNP == i), 6]
#  df2[which(df2$SNP == i), 15] = mr[which(mr$SNP == i), 8]
#  df2[which(df2$SNP == i), 16] = mr[which(mr$SNP == i), 9]
#  df2[which(df2$SNP == i), 17] = mr[which(mr$SNP == i), 10]
#  df2[which(df2$SNP == i), 18] = mr[which(mr$SNP == i), 12]
#  df2[which(df2$SNP == i), 19] = mr[which(mr$SNP == i), 13]
#}


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

d1 = df2_edu[1:19,]; d2 = df2_noedu[1:19,]
df2 = cbind(test1, test2, 7)

df2$Gene = as.character(rep.row(snps$GENE, 19))
df2$Gene = "HIVEP3"

plots = list()
target = results_edu$SNP
titles = 1
plot_num = 1
qs = 1:19

for(i in 1:length(target)){
  name1 = df2[qs, 26][1]
  name2 = paste("(", target[titles], ")", sep ="")
  
  myplot = ggplot(data = df2[qs,],aes(qs)) + 
    geom_ribbon(aes(ymin = df2[,4], ymax = df2[,5]), fill = "deepskyblue1", alpha = 0.2) +
    geom_ribbon(aes(ymin = df2[,17], ymax = df2[,18]), fill = "red1", alpha = 0.2) +
    geom_line(aes(y = df2[,3]), colour = "deepskyblue1", size = 0.6) +
    geom_line(aes(y = df2[,16]), colour = "red1", size = 0.6) +
    scale_x_continuous(breaks = seq(0,100, 20)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold", vjust = 1), 
          legend.title = element_text(), legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(face = "bold", color = "black", size = 5),
          axis.text.y = element_text(face = "bold", color = "black",  size = 5)) +
    ggtitle(paste(name1,"\n", name2))
  
  plots[[plot_num]] = myplot
  titles = titles + 1
  plot_num = plot_num + 1
  qs = qs + 19
}


plot_num = 1
ppi = 400

for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                     plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], 
                     plots[[plot_num + 9]], plots[[plot_num + 10]], plots[[plot_num + 11]], plots[[plot_num + 12]], plots[[plot_num + 13]], 
                     plots[[plot_num + 14]], plots[[plot_num + 15]], ncol = 4, nrow = 4)
  file_out = paste("D:/", "AgeOfOnset_plot", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                             left = text_grob("Genetic effect size (diopters per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 16
}





