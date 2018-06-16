library(ggplot2)

{EduYearsHigh$cc_0[EduYearsHigh$avMSE <= 0]= 1; EduYearsHigh$cc_0[EduYearsHigh$avMSE > 0]= 0
EduYearsHigh$cc_0.25[EduYearsHigh$avMSE <= -0.25]= 1; EduYearsHigh$cc_0.25[EduYearsHigh$avMSE > -0.25]= 0
EduYearsHigh$cc_0.5[EduYearsHigh$avMSE <= -0.5]= 1; EduYearsHigh$cc_0.5[EduYearsHigh$avMSE > -0.5]= 0
EduYearsHigh$cc_0.75[EduYearsHigh$avMSE <= -0.75]= 1; EduYearsHigh$cc_0.75[EduYearsHigh$avMSE > -0.75]= 0
EduYearsHigh$cc_1[EduYearsHigh$avMSE <= -1]= 1; EduYearsHigh$cc_1[EduYearsHigh$avMSE > -1]= 0
EduYearsHigh$cc_1.25[EduYearsHigh$avMSE <= -1.25]= 1; EduYearsHigh$cc_1.25[EduYearsHigh$avMSE > -1.25]= 0
EduYearsHigh$cc_1.5[EduYearsHigh$avMSE <= -1.5]= 1; EduYearsHigh$cc_1.5[EduYearsHigh$avMSE > -1.5]= 0
EduYearsHigh$cc_1.75[EduYearsHigh$avMSE <= -1.75]= 1; EduYearsHigh$cc_1.75[EduYearsHigh$avMSE > -1.75]= 0
EduYearsHigh$cc_2[EduYearsHigh$avMSE <= -2]= 1; EduYearsHigh$cc_2[EduYearsHigh$avMSE > -2]= 0
EduYearsHigh$cc_2.25[EduYearsHigh$avMSE <= -2.25]= 1; EduYearsHigh$cc_2.25[EduYearsHigh$avMSE > -2.25]= 0
EduYearsHigh$cc_2.5[EduYearsHigh$avMSE <= -2.5]= 1; EduYearsHigh$cc_2.5[EduYearsHigh$avMSE > -2.5]= 0
EduYearsHigh$cc_2.75[EduYearsHigh$avMSE <= -2.75]= 1; EduYearsHigh$cc_2.75[EduYearsHigh$avMSE > -2.75]= 0
EduYearsHigh$cc_3[EduYearsHigh$avMSE <= -3]= 1; EduYearsHigh$cc_3[EduYearsHigh$avMSE > -3]= 0
EduYearsHigh$cc_3.25[EduYearsHigh$avMSE <= -3.25]= 1; EduYearsHigh$cc_3.25[EduYearsHigh$avMSE > -3.25]= 0
EduYearsHigh$cc_3.5[EduYearsHigh$avMSE <= -3.5]= 1; EduYearsHigh$cc_3.5[EduYearsHigh$avMSE > -3.5]= 0
EduYearsHigh$cc_3.75[EduYearsHigh$avMSE <= -3.75]= 1; EduYearsHigh$cc_3.75[EduYearsHigh$avMSE > -3.75]= 0
EduYearsHigh$cc_4[EduYearsHigh$avMSE <= -4]= 1; EduYearsHigh$cc_4[EduYearsHigh$avMSE > -4]= 0
EduYearsHigh$cc_4.25[EduYearsHigh$avMSE <= -4.25]= 1; EduYearsHigh$cc_4.25[EduYearsHigh$avMSE > -4.25]= 0
EduYearsHigh$cc_4.5[EduYearsHigh$avMSE <= -4.5]= 1; EduYearsHigh$cc_4.5[EduYearsHigh$avMSE > -4.5]= 0
EduYearsHigh$cc_4.75[EduYearsHigh$avMSE <= -4.75]= 1; EduYearsHigh$cc_4.75[EduYearsHigh$avMSE > -4.75]= 0
EduYearsHigh$cc_5[EduYearsHigh$avMSE <= -5]= 1; EduYearsHigh$cc_5[EduYearsHigh$avMSE > -5]= 0
EduYearsHigh$cc_5.25[EduYearsHigh$avMSE <= -5.25]= 1; EduYearsHigh$cc_5.25[EduYearsHigh$avMSE > -5.25]= 0
EduYearsHigh$cc_5.5[EduYearsHigh$avMSE <= -5.5]= 1; EduYearsHigh$cc_5.5[EduYearsHigh$avMSE > -5.5]= 0
EduYearsHigh$cc_5.75[EduYearsHigh$avMSE <= -5.75]= 1; EduYearsHigh$cc_5.75[EduYearsHigh$avMSE > -5.75]= 0
EduYearsHigh$cc_6[EduYearsHigh$avMSE <= -6]= 1; EduYearsHigh$cc_6[EduYearsHigh$avMSE > -6]= 0}


ind = 25:49

prop = matrix(nrow = length(ind), ncol = 2)
for(i in ind){
  prop[i-(length(ind)+1),] = table(EduYearsHigh[,i])
}
prop = as.data.frame(prop)
colnames(prop) = c("Controls", "Cases")

Thresholds = c(0,-0.25,-0.5,-0.75,-1,-1.25,-1.5,-1.75,-2,-2.25,-2.5,-2.75,-3,-3.25,-3.5,-3.75,-4,-4.25,-4.5,-4.75,-5,-5.25,-5.5,-5.75,-6)
props = cbind(Thresholds, prop)


controls = true.heritability.prevalence[which(true.heritability.prevalence$Status == 0),]
cases = true.heritability.prevalence[which(true.heritability.prevalence$Status == 1),]
final = cbind(controls[,2], cases[,2])
colnames(final) <- c("Control", "Cases")
rownames(final) <- Thresholds


{predicted$cc_0[predicted$predicted_avMSE <= 0]= 1; predicted$cc_0[predicted$predicted_avMSE > 0]= 0
  predicted$cc_0.25[predicted$predicted_avMSE <= -0.25]= 1; predicted$cc_0.25[predicted$predicted_avMSE > -0.25]= 0
  predicted$cc_0.5[predicted$predicted_avMSE <= -0.5]= 1; predicted$cc_0.5[predicted$predicted_avMSE > -0.5]= 0
  predicted$cc_0.75[predicted$predicted_avMSE <= -0.75]= 1; predicted$cc_0.75[predicted$predicted_avMSE > -0.75]= 0
  predicted$cc_1[predicted$predicted_avMSE <= -1]= 1; predicted$cc_1[predicted$predicted_avMSE > -1]= 0
  predicted$cc_1.25[predicted$predicted_avMSE <= -1.25]= 1; predicted$cc_1.25[predicted$predicted_avMSE > -1.25]= 0
  predicted$cc_1.5[predicted$predicted_avMSE <= -1.5]= 1; predicted$cc_1.5[predicted$predicted_avMSE > -1.5]= 0
  predicted$cc_1.75[predicted$predicted_avMSE <= -1.75]= 1; predicted$cc_1.75[predicted$predicted_avMSE > -1.75]= 0
  predicted$cc_2[predicted$predicted_avMSE <= -2]= 1; predicted$cc_2[predicted$predicted_avMSE > -2]= 0
  predicted$cc_2.25[predicted$predicted_avMSE <= -2.25]= 1; predicted$cc_2.25[predicted$predicted_avMSE > -2.25]= 0
  predicted$cc_2.5[predicted$predicted_avMSE <= -2.5]= 1; predicted$cc_2.5[predicted$predicted_avMSE > -2.5]= 0
  predicted$cc_2.75[predicted$predicted_avMSE <= -2.75]= 1; predicted$cc_2.75[predicted$predicted_avMSE > -2.75]= 0
  predicted$cc_3[predicted$predicted_avMSE <= -3]= 1; predicted$cc_3[predicted$predicted_avMSE > -3]= 0}


ind = 24:36

predicteds = matrix(nrow = length(ind), ncol = 2)
for(i in ind){
  predicteds[i-(length(ind)+10),] = table(predicted[,i])
}
predicteds = as.data.frame(predicteds)
colnames(predicteds) = c("Controls", "Cases")

Thresholds = c(0,-0.25,-0.5,-0.75,-1,-1.25,-1.5,-1.75,-2,-2.25,-2.5,-2.75,-3)
propsed = cbind(Thresholds, predicteds)


controls2 = predicted.heritability.prevalence[which(predicted.heritability.prevalence$Status == 0),]
cases2 = predicted.heritability.prevalence[which(predicted.heritability.prevalence$Status == 1),]
final2 = cbind(controls2[,2], cases2[,2])
colnames(final2) <- c("Control", "Cases")
rownames(final2) <- Thresholds



{par(mfrow=c(1,2))
  par(mar=c(7, 7, 2, 0.5) + 0.2, mgp = c(5.8, 1, 0))
  barplot(t(final), beside = T, ylab = "Number of individuals", xlab = "True avMSE in diopters",
          cex.names = 1.35, las = 1, yaxt = "n", col = c("darkblue","red"),
          cex.lab = 1.55, cex.axis = 1.3)
  box(bty="l")
  axis(2, at = seq(0, 70000, 5000), las = 1, cex.axis = 1.55)
  legend("topleft", 
         legend = c("Controls", "Cases"), 
         fill = c("darkblue", "red"), bty = "n", cex = 1.35)
  
  barplot(t(final2), beside = T, xlab = "Predicted avMSE in diopters",
          cex.names = 1.35, las = 1, yaxt = "n", col = c("darkblue","red"),
          cex.lab = 1.55, cex.axis = 1.3)
  box(bty="l")
  axis(2, at = seq(0, 160000, 10000), las = 1, cex.axis = 1.55)
  legend("topleft", 
         legend = c("Controls", "Cases"), 
         fill = c("darkblue", "red"), bty = "n", cex = 1.35)}







