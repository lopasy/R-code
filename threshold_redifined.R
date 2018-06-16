library(data.table)
meta = fread("~/Desktop/meta.txt")
meta1 = fread("~/Desktop/meta.txt")
meta2 = fread("~/Desktop/meta.txt")
meta3 = fread("~/Desktop/meta.txt")


meta$mean=rowSums(meta[,c(2:6)], na.rm = T)
meta$mean=meta$mean/5
meta1$mean=rowSums(meta1[,c(2:6)], na.rm = T)
meta1$mean=meta1$mean/5
meta2$mean=rowSums(meta2[,c(2:6)], na.rm = T)
meta2$mean=meta2$mean/5
meta3$mean=rowSums(meta3[,c(2:6)], na.rm = T)
meta3$mean=meta3$mean/5

par(mfrow=c(1,2))

se <- function(x) sd(x)/sqrt(length(x))
meta$se = 0; meta1$se = 0; meta2$se = 0; meta3$se = 0
tt = 1:23
for (i in tt){
meta[i,8] = se(meta[i,c(2:6)])
}
for (i in tt){
  meta1[i,8] = se(meta1[i,c(2:6)])
}
for (i in tt){
  meta2[i,8] = se(meta2[i,c(2:6)])
}
for (i in tt){
  meta3[i,8] = se(meta3[i,c(2:6)])
}

ld$estimate = 0 
ld$lower = 0
ld$upper = 0
ld$se = 0
ld=as.matrix(ld)
for(i in tt){
ld[i,8] = ci( ld[i,c(2:7)] )[1]
ld[i,9] = ci( ld[i,c(2:7)] )[2]
ld[i,10] = ci( ld[i,c(2:7)] )[3]
ld[i,11] = ci( ld[i,c(2:7)] )[4]
}

meta1$Threshold = meta$Threshold + 0.05
meta2$Threshold = meta$Threshold + 0.1
meta3$Threshold = meta$Threshold + 0.15


{plot(meta$Threshold, meta$mean, type = "b", col = "firebrick", lwd = 2, main = "Mean FDR", 
      xlab = "-log(p)", ylab = "FDR", cex.lab = 1.5, ylim = c(0.5,1))
segments(meta$Threshold, meta$mean - meta$se,
         meta$Threshold, meta$mean + meta$se)
epsilon <- 0.025
segments(meta$Threshold-epsilon,meta$mean - meta$se,meta$Threshold+epsilon,meta$mean - meta$se)
segments(meta$Threshold-epsilon,meta$mean + meta$se,meta$Threshold+epsilon,meta$mean + meta$se)

lines(meta1$Threshold, meta1$mean, type = "b", col = "firebrick1", lwd = 2, main = "Mean FDR", 
     xlab = "-log(p)", ylab = "FDR", cex.lab = 1.5, ylim = c(0.5,1))
segments(meta1$Threshold, meta1$mean - meta1$se,
         meta1$Threshold, meta1$mean + meta1$se)
segments(meta1$Threshold-epsilon,meta1$mean - meta1$se,meta1$Threshold+epsilon,meta1$mean - meta1$se)
segments(meta1$Threshold-epsilon,meta1$mean + meta1$se,meta1$Threshold+epsilon,meta1$mean + meta1$se)

lines(meta2$Threshold, meta2$mean, type = "b", col = "deepskyblue", lwd = 2, main = "Mean FDR", 
      xlab = "-log(p)", ylab = "FDR", cex.lab = 1.5, ylim = c(0.5,1))
segments(meta2$Threshold, meta2$mean - meta2$se,
         meta2$Threshold, meta2$mean + meta2$se)
segments(meta2$Threshold-epsilon,meta2$mean - meta2$se,meta2$Threshold+epsilon,meta2$mean - meta2$se)
segments(meta2$Threshold-epsilon,meta2$mean + meta2$se,meta2$Threshold+epsilon,meta2$mean + meta2$se)

lines(meta3$Threshold, meta3$mean, type = "b", col = "darkblue", lwd = 2, main = "Mean FDR", 
      xlab = "-log(p)", ylab = "FDR", cex.lab = 1.5, ylim = c(0.5,1))
segments(meta3$Threshold, meta3$mean - meta3$se,
         meta3$Threshold, meta3$mean + meta3$se)
segments(meta3$Threshold-epsilon,meta3$mean - meta3$se,meta3$Threshold+epsilon,meta3$mean - meta3$se)
segments(meta3$Threshold-epsilon,meta3$mean + meta3$se,meta3$Threshold+epsilon,meta3$mean + meta3$se)

legend("bottomleft", c("Predicted continious", "Predicted case-control", "True continious", "True case-control"), bty = "n", col = c("firebrick","black","deepskyblue", "darkblue"),
       text.col = c("firebrick","black","deepskyblue", "darkblue"), lty = c(-1, -1, -1, -1),
       cex = 1.2,box.lwd = 0, horiz = F)
abline(h = 0.8, v = 5.45, lty = 2)}

