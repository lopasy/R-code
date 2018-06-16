# True vs Predicted coefficients EduYearsHigh
predicted = predicted[which(predicted$SNPID %in% true$SNPID),]
cor.test(predicted$BETA, true$BETA)

predicted_0.05 = predicted[which(predicted$P < 0.000005),] 
true_0.05 = true[which(true$P < 0.000005),] 

predicted_0.05 = predicted_0.05[which(predicted_0.05$SNPID %in% true_0.05$SNPID),]
#or
true_0.05 = true_0.05[which(true_0.05$SNPID %in% predicted_0.05$SNPID),]

cor.test(predicted_0.05$BETA, true_0.05$BETA)

plot(predicted$BETA, true$BETA)


