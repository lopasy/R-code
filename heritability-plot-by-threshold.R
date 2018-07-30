her = as.data.frame(her)
{par(mfrow=c(1,2))
  
  her$ho3 = her$Threshold - 0.05
  her$Total_ho3 = her$Threshold + 0.05
  
  plot(her$ho3, her$ho2, col = "deepskyblue", pch = 16, main = "Observed scale", 
       ylab = "Heritability", xlab = "avMSE in diopters", ylim = c(-0.05,0.2), xlim = c(-3,0.02), type = "b", lwd = 2,
       cex.axis = 1.35, cex.lab = 1.35, cex.main = 1.5)
  
  legend(x="topleft",col = c("deepskyblue","firebrick3", "gold3"), legend = c("G","GxE", "Combined"), 
         bty = "n", pch = c(16, 2, 8), lty = c(NA, NA, NA), lwd = 2, cex = 1.2)
  ################################################################################################## 
  segments(her$ho3, her$ho2 - her$ho2_se,
           her$ho3, her$ho2 + her$ho2_se)
  epsilon <- 0.025
  segments(her$ho3-epsilon,her$ho2 - her$ho2_se,her$ho3+epsilon,her$ho2 - her$ho2_se)
  segments(her$ho3-epsilon,her$ho2 + her$ho2_se,her$ho3+epsilon,her$ho2 + her$ho2_se)
  ##################################################################################################  
  points(her$Threshold, her$GxE_ho2, col="firebrick3", pch = 2, type = "b", lwd = 2)
  
  segments(her$Threshold, her$GxE_ho2 - her$GxE_ho2_se,
           her$Threshold, her$GxE_ho2 + her$GxE_ho2_se)
  epsilon <- 0.025
  segments(her$Threshold-epsilon,her$GxE_ho2 - her$GxE_ho2_se,her$Threshold+epsilon,her$GxE_ho2 - her$GxE_ho2_se)
  segments(her$Threshold-epsilon,her$GxE_ho2 + her$GxE_ho2_se,her$Threshold+epsilon,her$GxE_ho2 + her$GxE_ho2_se)
  ##################################################################################################  
  points(her$Total_ho3, her$Total_ho2, col="gold3", pch = 8, type = "b", lwd = 2)
  
  segments(her$Total_ho3, her$Total_ho2 - her$Totalh_o2_se,
           her$Total_ho3, her$Total_ho2 + her$Totalh_o2_se)
  epsilon <- 0.025
  segments(her$Total_ho3-epsilon,her$Total_ho2 - her$Totalh_o2_se,her$Total_ho3+epsilon,her$Total_ho2 - her$Totalh_o2_se)
  segments(her$Total_ho3-epsilon,her$Total_ho2 + her$Totalh_o2_se,her$Total_ho3+epsilon,her$Total_ho2 + her$Totalh_o2_se)
  ##################################################################################################
  
  
  plot(her$ho3, her$h_L2, col = "deepskyblue", pch = 16, main = "Liability scale", 
       ylab = "Heritability", xlab = "avMSE in diopters", ylim = c(-0.05,0.3), xlim = c(-3,0.02), type = "b", lwd = 2,
       cex.axis = 1.35, cex.lab = 1.35, cex.main = 1.5)
  
  legend(x="topright",col = c("deepskyblue","firebrick3", "gold3"), legend = c("G","GxE", "Combined"), 
         bty = "n", pch = c(16, 2, 8), lty = c(NA, NA, NA), lwd = 2, cex = 1.2)
  ##################################################################################################  
  segments(her$ho3, her$h_L2 - her$h_L2_se,
           her$ho3, her$h_L2 + her$h_L2_se)
  epsilon <- 0.025
  segments(her$ho3-epsilon,her$h_L2 - her$h_L2_se,her$ho3+epsilon,her$h_L2 - her$h_L2_se)
  segments(her$ho3-epsilon,her$h_L2 + her$h_L2_se,her$ho3+epsilon,her$h_L2 + her$h_L2_se)
  ##################################################################################################  
  points(her$Threshold, her$GxE_h_L2, col="firebrick3", pch = 2, type = "b", lwd = 2)
  
  segments(her$Threshold, her$GxE_h_L2 - her$GxE_h_L2se,
           her$Threshold, her$GxE_h_L2 + her$GxE_h_L2se)
  epsilon <- 0.025
  segments(her$Threshold-epsilon,her$GxE_h_L2 - her$GxE_h_L2se,her$Threshold+epsilon,her$GxE_h_L2 - her$GxE_h_L2se)
  segments(her$Threshold-epsilon,her$GxE_h_L2 + her$GxE_h_L2se,her$Threshold+epsilon,her$GxE_h_L2 + her$GxE_h_L2se)
  ################################################################################################## 
  points(her$Total_ho3, her$Total_h_L2, col="gold3", pch = 8, type = "b", lwd = 2)
  
  segments(her$Total_ho3, her$Total_h_L2 - her$Total_h_L2_se,
           her$Total_ho3, her$Total_h_L2 + her$Total_h_L2_se)
  epsilon <- 0.025
  segments(her$Total_ho3-epsilon,her$Total_h_L2 - her$Total_h_L2_se,her$Total_ho3+epsilon,her$Total_h_L2 - her$Total_h_L2_se)
  segments(her$Total_ho3-epsilon,her$Total_h_L2 + her$Total_h_L2_se,her$Total_ho3+epsilon,her$Total_h_L2 + her$Total_h_L2_se)
}

a = round(her$Threshold, 1)
{plot(her$Threshold, her$Prevalence, col = "darkblue", pch = 16, main = "Disease prevalence for different avMSE thresholds", 
     ylab = "Prevalence", xlab = "avMSE in diopters", ylim = c(0,0.5), xaxt = "n", type = "b",
     cex.main = 1.5, cex.lab = 1.35, cex.axis = 1.35)
axis(1, at = her$Threshold, cex.axis = 1.35)

pred_thr=c(0.43,0.4,0.38,0.35,0.32,0.3,0.26,0.23,0.2,0.17,0.13,0.1,0.071,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
pred_thr=as.data.frame(pred_thr)
her$Prev = pred_thr
points(her$Threshold, pred_thr$pred_thr, col = "red", pch = 16, main = "Disease prevalence under different avMSE thresholds", 
       ylab = "Prevalence", xlab = "Diopters", ylim = c(0,0.5),  type = "b")
legend(x="topleft",col = c("darkblue","red"), legend = c("True","Predicted"), 
       bty = "n", lwd = 2, cex = 1.35)
}







