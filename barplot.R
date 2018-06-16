# Convert your output to chi-squared values
# For z-scores, just square the t statistic
chisq <- reports$STAT^2

# For p-values, calculate chi-squared statistic
chisq <- qchisq(1-chr6$pval,1)

# Calculate lambda gc (??gc)
median(chisq)/qchisq(0.5,1)








z = merge(EduYearsHigh,cv_0.00000001, by = "FID")
mean(z$SCORESUM)

cor.test(z$avMSE, z$SCORESUM)


plot(density(z$avMSE))
lines(density(z$SCORESUM))

summary(lm(avMSE~SCORESUM,z))
summary(lm(avMSE~SCORESUM+EduYearsHigh+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array,z))

sd(z$SCORESUM)
sd(z$SCORESUM)^2

summary(lm(avMSE~SCORESUM+EduYearsHigh+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array,z))

plot(density(predicted$predicted_avMSE))
lines(density(predicted$SCORESUM))



a = c(1e-01, 5e-02, 1e-02, 5e-03, 1e-03, 5e-04, 1e-04, 5e-05, 1e-05, 5e-06, 1e-06, 5e-07, 1e-07, 5e-08, 1e-08)
png("plot_w_rotated_axis_labels.png", height=3, width=6, units="in", res=400)
{par(mar = c(6, 4.1, 3, 2) + 0.2)
x = barplot(score$GRS, names.arg=score$Threshold, ylim=c(0,0.08), col = "dodgerblue4", space = 0.3, xaxt = "n",
            ylab = "% Variance Explained", main = "GRS", cex.axis = 1.2, cex.lab = 1.2)
options(scipen=999)
text(cex=1.2, x=x[,1], y=-0.003, labels = paste(a), xpd=TRUE, srt=45, adj = 1)}
dev.off()



