### GS CQR ###
# avMSE
{cream = read.csv("F:/Full_ukb/January/gg.raw", sep = "")
master = read.csv("D:/master_true.txt", sep = "")
final = merge(master, cream, "FID")
final = final[,c(3:13, 19:22, 53:202)]
final = final[,-162]
final3 = final[,16:164]

gene_scores = as.data.frame(rowSums(final3, na.rm = T)); names(gene_scores) = "GS"
final = cbind(final, gene_scores)

qs = seq(5,95,5)/100
n_boots = 2000
ols_avMSE = lm(avMSE ~ GS + Sex + poly(Age,2) + Array, data = final)
ols_mean  = summary(ols_avMSE)$coef[2,1]
ols_upper = summary(ols_avMSE)$coef[2,1] + (1.96*summary(ols)$coef[2,2])
ols_lower = summary(ols_avMSE)$coef[2,1] - (1.96*summary(ols)$coef[2,2])
p = summary(ols_avMSE)$coef[2,4]
cqr_avMSE = summary(rq(avMSE ~ GS + Sex + poly(Age,2) + Array, data = final, tau = qs), se = "boot", bsmethod = "mcmb", R = n_boots)
qr2 = rq(avMSE ~ GS + Sex + poly(Age,2) + Array, data = final, tau = qs)

meta = as.data.frame(matrix(nrow = 19, ncol = 4))
names(meta) = c("BETA", "SE", "P", "95% CI")
meta$quantile = rep(seq(0.05,0.95,0.05))

df2 = as.data.frame(matrix(nrow = 19, ncol = 5))
names(df2) = c("x","intercept","slope","lci","uci")
df2$x         = qs*100
df2$intercept = qr2$coef[1,]
df2$slope     = qr2$coef[2,]

g = 1
for (i in 1:19) {
  quantile = cqr_avMSE[[i]]
  meta[g, 1] = quantile$coefficients[2,1]
  meta[g, 2] = quantile$coefficients[2,2]
  meta[g, 3] = quantile$coefficients[2,4]
  meta[g, 4] = paste("[", quantile$coefficients[2,1] - (1.96*quantile$coefficients[2,2]), ";",
                     quantile$coefficients[2,1] + (1.96*quantile$coefficients[2,2]), "]")
  df2[g,4] = quantile$coefficients[2,1] - (1.96*quantile$coefficients[2,2])
  df2[g,5] = quantile$coefficients[2,1] + (1.96*quantile$coefficients[2,2])
  g = g + 1
}


res = rma(yi = meta[,1], sei = meta[,2], mods = ~ poly(quantile, 2), data = meta)
predicted = as.data.frame(predict(res, interval = 'confidence', level = 0.95))
df2 = cbind(df2, predicted)

file_out = paste("E:/2018/", "avMSE GS", ".png", sep = "")
png(file_out, width = 5*ppi, height = 5*ppi, res = ppi)
myplot = ggplot(data = df2,aes(x)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70", alpha = 0.8) +
  geom_line(aes(y = slope), colour = "black", size = 1.5) +
  geom_line(aes(y = pred), colour = "blue", size = 0.8) +
  geom_line(aes(y = ci.lb), linetype = "longdash", size = 0.8, colour = "blue") + 
  geom_line(aes(y = ci.ub), linetype = "longdash", size = 0.8, colour = "blue") +
  geom_hline(yintercept = ols_mean, colour = "red", size = 0.8) +
  geom_hline(yintercept = ols_upper, colour = "red", size = 0.8, linetype = 2) +
  geom_hline(yintercept = ols_lower, colour = "red", size = 0.8, linetype = 2) +
  scale_x_continuous(breaks = seq(0,100, 10)) +
  labs(title = "GS-avMSE", x = "\navMSE percentile", y = "Change in avMSE (D)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", vjust = -4,
                                  margin = margin(t = 20, r = 20, b = 0, l = 10)), 
        legend.title = element_text(), legend.position = "none",
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold",margin = margin(r = 10)), 
        axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black",  size = 10))

print(myplot)
dev.off()

miss = as.data.frame(matrix(nrow = nrow(final3), ncol = 1))
imp = ncol(final3)
for(i in 1:imp){
  miss[i,] = sum(is.na(final3[,i]))
}}


# Height
{cream = read.csv("D:/height.raw", sep = "")
master = read.csv("D:/master_true.txt", sep = "")
final = merge(master, cream, "FID")
final = final[,c(3:13, 19:22, 53:320)]
final3 = final[,16:283]

gene_scores = as.data.frame(rowSums(final3, na.rm = T)); names(gene_scores) = "GS"
final = cbind(final, gene_scores)

qs = seq(5,95,5)/100
n_boots = 1000
ols_avMSE = lm(avMSE ~ GS + Sex + poly(Age,1) + Array, data = final)
ols_mean  = summary(ols_avMSE)$coef[2,1]
ols_upper = summary(ols_avMSE)$coef[2,1] + (1.96*summary(ols)$coef[2,2])
ols_lower = summary(ols_avMSE)$coef[2,1] - (1.96*summary(ols)$coef[2,2])
p = summary(ols_avMSE)$coef[2,4]
summary(ols_avMSE)$coef[2,4]
cqr_avMSE = summary(rq(avMSE ~ GS + Sex + poly(Age, 1) + Array, data = final, tau = qs), se = "boot", bsmethod = "mcmb", R = n_boots)
qr2 = rq(avMSE ~ GS + Sex + poly(Age, 1) + Array, data = final, tau = qs)

meta = as.data.frame(matrix(nrow = 19, ncol = 4))
names(meta) = c("BETA", "SE", "P", "95% CI")
meta$quantile = rep(seq(0.05,0.95,0.05))

df2 = as.data.frame(matrix(nrow = 19, ncol = 5))
names(df2) = c("x","intercept","slope","lci","uci")
df2$x         = qs*100
df2$intercept = qr2$coef[1,]
df2$slope     = qr2$coef[2,]

g = 1
for (i in 1:19) {
  quantile = cqr_avMSE[[i]]
  meta[g, 1] = quantile$coefficients[2,1]
  meta[g, 2] = quantile$coefficients[2,2]
  meta[g, 3] = quantile$coefficients[2,4]
  meta[g, 4] = paste("[", quantile$coefficients[2,1] - (1.96*quantile$coefficients[2,2]), ";",
                     quantile$coefficients[2,1] + (1.96*quantile$coefficients[2,2]), "]")
  df2[g,4] = quantile$coefficients[2,1] - (1.96*quantile$coefficients[2,2])
  df2[g,5] = quantile$coefficients[2,1] + (1.96*quantile$coefficients[2,2])
  g = g + 1
}


res = rma(yi = meta[,1], sei = meta[,2], mods = ~ poly(quantile, 2), data = meta)
predicted = as.data.frame(predict(res, interval = 'confidence', level = 0.95))
df2 = cbind(df2, predicted)

file_out = paste("E:/2018/", "height GS", ".png", sep = "")
  png(file_out, width = 5*ppi, height = 5*ppi, res = ppi)
  myplot = ggplot(data = df2,aes(x)) + 
    geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = slope), colour = "black", size = 1.5) +
    geom_line(aes(y = pred), colour = "blue", size = 0.8) +
    geom_line(aes(y = ci.lb), linetype = "longdash", size = 0.8, colour = "blue") + 
    geom_line(aes(y = ci.ub), linetype = "longdash", size = 0.8, colour = "blue") +
    geom_hline(yintercept = ols_mean, colour = "red", size = 0.8) +
    geom_hline(yintercept = ols_upper, colour = "red", size = 0.8, linetype = 2) +
    geom_hline(yintercept = ols_lower, colour = "red", size = 0.8, linetype = 2) +
    scale_x_continuous(breaks = seq(0,100, 10)) +
    labs(title = "GS-Height", x = "\nHeight percentile", y = "Change in Height (cm)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", vjust = 0,
                                    margin = margin(t = 20, r = 20, b = 0, l = 10)), 
          legend.title = element_text(), legend.position = "none",
          axis.title.x = element_text(size=14,face="bold"),
          axis.title.y = element_text(size=14,face="bold",margin = margin(r = 10)), 
          axis.text.x = element_text(face = "bold", color = "black", size = 10),
          axis.text.y = element_text(face = "bold", color = "black",  size = 10))
  
  print(myplot)
  dev.off()


miss = as.data.frame(matrix(nrow = nrow(final3), ncol = 1))
imp = nrow(final3)
for(i in 1:imp){
  miss[i,] = sum(is.na(final3[i,]))/ncol(final3)
}}






