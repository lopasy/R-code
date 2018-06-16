library(car)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

{par(mfrow=c(2,4))
plot(final$BETA_INT, -log10(final$P), main='Volcano plot', ylab = substitute(-log[10](pval)),
     xlab = expression(hat(beta)[GxE]), cex.lab = 1.2)
hist(final$P, 200, col='black', xlab = "P-value", main = "Histogram of P values", cex.lab = 1.2)
qqPlot(final$P, main='QQ-plot', ylab = substitute(-log[10]("observed P")),
       xlab = substitute(-log[10]("expected P")))
plot(-log10(final$P), main='Manhattan plot', xlab = "GxG Comparison P-value",
     ylab = substitute(-log[10](pval)), cex.lab = 1.2)
abline(h = -log10(0.000000005), col="blue", lwd=4, lty=2)
text(20,8.4,expression("5x10"^"-9"))
abline(h = -log10(0.0000000005), col="red", lwd=3, lty=1)
text(24,9.39,expression("5x10"^"-10"))
arrows(720,9.19,800,9.33, length=0.1, lwd = 2)
text(680,9.10,"chr4 rs73244033\nchr6 rs17196517", cex = 0.9)
arrows(450,9.16,510,9.06, length=0.1, lwd = 2)
text(400,9.24,"Intron in DSE gene", cex = 0.9)
hist(table(final$SNP1), 20, main = "Frequency of SNP 1 ", cex.lab = 1.2, xlab = "Number of interactions")
hist(table(final$SNP2), 20, main = "Frequency of SNP 2 ", cex.lab = 1.2, xlab = "Number of interactions")}

par(mfrow=c(1,1))

final[which(final$CHR1 == 4 & final$CHR2 == 6),]

a = table(final$CHR1 == 4, final$CHR2) # test
sum(table(final$CHR1 == 18, final$CHR2))
plot(table(final$CHR1 == 4, final$CHR2))[2,]
# par(mfrow=c(1,1))

c = c(1,  3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22)
a = matrix(nrow = 21, ncol = 21)
for (i in 1:21){
  tryCatch({
  a[i,] = table(final$CHR1 == i, final$CHR2)[2,]
  if (i==14) stop("SKIP")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
a = cbind(a,rep(0))
a[,3:22] = a[,2:21]
a[,2] = 0
a[c(14,18,19,21),] = 0
a = rbind(a,rep(0))


palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)
heatmap.2(a, Rowv = NA, Colv = NA, col = heat.colors(256), scale = "column", margins = c(5,10), trace = "none",
          cexRow=1.4, cexCol=1.4)

b = melt(a)
colnames(b)[1] = "SNP1"
colnames(b)[2] = "SNP2"
b = ddply(b, .(SNP2), transform, rescale = rescale(value))
(p <- ggplot(b, aes(SNP2, SNP1)) + geom_tile(aes(fill = rescale),
                                                     colour = "white") + scale_fill_gradient(low = "white",
                                                     high = "steelblue"))
p



write.csv(a, file = "GxG_count.csv", col.names = F, row.names = F, quote = F)

