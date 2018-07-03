library(plyr); library(ggplot2); library(ggpubr)

target = as.data.frame(mr2$SNP)
data1 = read.csv("D:/CQR GWAS/cream2017_ukbb_replicated.csv"); data2 = data1[,c(9,5)]
names(data2)[1]="SNP"

stopwords = c("_A","_C","_T","_G")
target[,1] = gsub(paste0(stopwords, collapse = "|"), "", target[,1])
test$SNP = gsub(paste0(stopwords, collapse = "|"), "", test$SNP)
results$SNP = gsub(paste0(stopwords, collapse = "|"), "", results$SNP)
mr$SNP = gsub(paste0(stopwords, collapse = "|"), "", mr$SNP)
snps = merge(results,data2,"SNP")


#df2 = snps

{test$beta_ols = 0; test$lci_ols = 0; test$uci_ols = 0}

target = results$SNP
for(i in target){
  test[which(test$SNP == i), 13] = results[which(results$SNP == i), 2]
  test[which(test$SNP == i), 14] = results[which(results$SNP == i), 3]
  test[which(test$SNP == i), 15] = results[which(results$SNP == i), 4]
}


plots = list()
target = results$SNP
titles = 1
plot_num = 1
qs = 1:9
qs2 = 1342:1350

for(i in 1:length(target)){
  name1 = snps[which(snps$SNP == target[i]),8]
  name2 = paste("(", target[titles], ")", sep ="")
  
myplot = ggplot(data = test[c(qs,qs2),],aes(x = qs, group = age)) + 
  geom_ribbon(aes(ymin = cqr_lci, ymax = cqr_uci), fill = "grey70", alpha = 0.8) +
  geom_line(aes(y = cqr_slope), colour = "black", size = 0.6) +
  geom_line(aes(y = spline_pred), colour = "blue", size = 0.4) +
  geom_line(aes(y = spline_ci.lb), linetype = "longdash", size = 0.4, colour = "blue") + 
  geom_line(aes(y = spline_ci.ub), linetype = "longdash", size = 0.4, colour = "blue") +
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
  facet_wrap(~group, scales="free_x") +
  ggtitle(paste(name1,"\n", name2))

plots[[plot_num]] = myplot
titles = titles + 1
plot_num = plot_num + 1
qs = qs + 9
qs2 = qs2 + 9
}

plot_num = 1
ppi = 400

for(i in 1:30){
  figure = ggarrange(plots[[plot_num]], plots[[plot_num + 1]], plots[[plot_num + 2]], plots[[plot_num + 3]], 
                     plots[[plot_num + 4]], plots[[plot_num + 5]], plots[[plot_num + 6]], plots[[plot_num + 7]], plots[[plot_num + 8]], 
                     plots[[plot_num + 9]], plots[[plot_num + 10]], plots[[plot_num + 11]], ncol = 3, nrow = 4)
  file_out = paste("D:/", "AgeOfOnset_plot_grouped", i , ".png", sep = "")
  png(file_out, width = 5*ppi, height = 6*ppi, res = ppi)
  last_fig = annotate_figure(figure, bottom = text_grob("Age of Onset percentile", color = "black", size = 10),
                             left = text_grob("Genetic effect size (years per copy of the risk allele)", color = "black", rot = 90, size = 10))
  print(last_fig)
  dev.off()
  plot_num = plot_num + 12
}


