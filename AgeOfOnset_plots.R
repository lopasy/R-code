library(data.table)
################################################################
################################################################
####################### Age Of Onset ###########################
################################################################
################################################################
library(plyr)
test = merge(results_edu, results_noedu, 1)
test = merge(test, results_int, 1)


data1 = read.csv("~/Documents/CQR_17_04_2018/cream2017_ukbb_replicated.csv"); data2 = data1[,c(9,5)]
target = as.data.frame(data2$ukbSNP); snps = merge(target,data2,1,); names(snps)[1] = "SNP_edu"

stopwords = c("_A","_C","_T","_G")
AoO_cqr[,7] = gsub(paste0(stopwords, collapse = "|"), "", AoO_cqr[,7])

snps = join(AoO_cqr, snps, "SNP_edu")


plots = list()
target = unique(snps$SNP_edu)
titles = 1
plot_num = 1
qs = 1:19

for(i in 1:length(target)){
  name1 = snps[qs, 62][1]
  name2 = paste("(", target[titles], ")", sep ="")
  
  myplot = ggplot(data = snps[qs,],aes(qs)) + 
    geom_ribbon(aes(ymin = cqr_lci_edu, ymax = cqr_uci_edu), fill = "deepskyblue1", alpha = 0.2) +
    geom_ribbon(aes(ymin = cqr_lci_noedu, ymax = cqr_uci_noedu), fill = "red1", alpha = 0.2) +
    geom_line(aes(y = cqr_slope_edu), colour = "deepskyblue1", size = 0.6) +
    geom_line(aes(y = cqr_slope_noedu), colour = "red1", size = 0.6) +
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
                     plots[[plot_num + 14]], plots[[plot_num + 15]], plots[[plot_num + 16]], plots[[plot_num + 17]],
                     plots[[plot_num + 18]], plots[[plot_num + 19]]
                     , ncol = 4, nrow = 5)
  file_out = paste("D:/", "AgeOfOnset_plot_edu.v.noedu", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                             left = text_grob("Genetic effect size (years per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 20
}





{chr = list.files(pattern = "df2_A.", full.names = TRUE); chr = lapply(chr, fread)
  df2 = do.call(rbind.data.frame, chr)
  
  chr = list.files(pattern = "results_A.", full.names = TRUE); chr = lapply(chr, fread)
  results = do.call(rbind.data.frame, chr)
  
  chr = list.files(pattern = "mr_A.", full.names = TRUE); chr = lapply(chr, fread)
  mr = do.call(rbind.data.frame, chr)
  
  
  target = as.data.frame(mr$SNP)
  data1 = read.csv("D:/CQR GWAS/cream2017_ukbb_replicated.csv"); data2 = data1[,c(9,5)]
  
  stopwords = c("_A","_C","_T","_G")
  target[,1] = gsub(paste0(stopwords, collapse = "|"), "", target[,1])
  df2$SNP = gsub(paste0(stopwords, collapse = "|"), "", df2$SNP)
  results$SNP = gsub(paste0(stopwords, collapse = "|"), "", results$SNP)
  mr$SNP = gsub(paste0(stopwords, collapse = "|"), "", mr$SNP)
  snps = merge(target,data2,1)
  
  
  #df2 = snps
  
  {df2$beta_ols = 0; df2$lci_ols = 0; df2$uci_ols = 0}
  
  target = results_noedu$SNP
  
  df2 = as.data.frame(df2); results = as.data.frame(results); mr = as.data.frame(mr)}

for(i in target){
  df2[which(df2$SNP == i), 23] = results[which(results$SNP == i), 2]
  df2[which(df2$SNP == i), 24] = results[which(results$SNP == i), 3]
  df2[which(df2$SNP == i), 25] = results[which(results$SNP == i), 4]
}



plots = list()
target = results$SNP
titles = 1
plot_num = 1
qs = 1:9

for(i in 1:length(target)){
  #name1 = df2[qs, 62][1]
  #name2 = paste("(", target[titles], ")", sep ="")
  
  myplot = ggplot(data = df2[qs,],aes(qs)) + 
    geom_ribbon(aes(ymin = cqr_lci, ymax = cqr_uci), fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = cqr_slope), colour = "black", size = 0.6) +
    geom_line(aes(y = spline_pred), colour = "blue", size = 0.4) +
    geom_line(aes(y = spline_ci.lb), linetype = "longdash", size = 0.4, colour = "blue") + 
    geom_line(aes(y = spline_ci.ub), linetype = "longdash", size = 0.4, colour = "blue") +
    #geom_hline(aes(yintercept = beta_ols), colour = "red", size = 0.4) +
    #geom_hline(aes(yintercept = uci_ols), colour = "red", size = 0.4, linetype = 2) +
    #geom_hline(aes(yintercept = lci_ols), colour = "red", size = 0.4, linetype = 2) +
    scale_x_continuous(breaks = seq(0,100, 20)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 6, face = "bold", vjust = 1), 
          legend.title = element_text(), legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(face = "bold", color = "black", size = 5),
          axis.text.y = element_text(face = "bold", color = "black",  size = 5)) #+
    #ggtitle(paste(name1,"\n", name2))
  
  plots[[plot_num]] = myplot
  titles = titles + 1
  plot_num = plot_num + 1
  qs = qs + 9
}


plot_num = 1
ppi = 400

for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                     plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], 
                     plots[[plot_num + 9]], plots[[plot_num + 10]], plots[[plot_num + 11]], plots[[plot_num + 12]], plots[[plot_num + 13]], 
                     plots[[plot_num + 14]], plots[[plot_num + 15]], plots[[plot_num + 16]], plots[[plot_num + 17]],
                     plots[[plot_num + 18]], plots[[plot_num + 19]], ncol = 4, nrow = 5)
  file_out = paste("D:/", "AgeOfOnset_plot_cluster1_", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                             left = text_grob("Genetic effect size (years per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 20
}


#############################
####### Interaction #########
#############################

plots = list()
target = results$SNP
titles = 1
plot_num = 1
qs = 1:19

for(i in 1:length(target)){
  name1 = df2[qs, 62][1]
  name2 = paste("(", target[titles], ")", sep ="")
  
  myplot = ggplot(data = df2[qs,],aes(qs)) + 
    geom_ribbon(aes(ymin = cqr_lci_int, ymax = cqr_uci_int), fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = cqr_slope_int), colour = "black", size = 0.6) +
    geom_line(aes(y = spline_pred_int), colour = "blue", size = 0.4) +
    geom_line(aes(y = spline_ci.lb_int), linetype = "longdash", size = 0.4, colour = "blue") + 
    geom_line(aes(y = spline_ci.ub_int), linetype = "longdash", size = 0.4, colour = "blue") +
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
                     plots[[plot_num + 14]], plots[[plot_num + 15]], plots[[plot_num + 16]], plots[[plot_num + 17]],
                     plots[[plot_num + 18]], plots[[plot_num + 19]], ncol = 4, nrow = 5)
  file_out = paste("D:/", "AgeOfOnset_plot_int", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                             left = text_grob("Genetic effect size (years per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 20
}





