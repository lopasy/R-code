# Make scatterplot

par(mgp = c(3,1,0))
plot(a$new,a$Heritability, main = "MSE Heritability By Chromosome", 
		xlab = "Chromosome length (Mb)", ylab = '', 
		cex.lab = 1.6, col = "blue", cex = 3, bg = "blue", pch = 21, cex.main = 1.8, cex.axis = 1.3)
mtext(side = 2, text = "Additive variance explained", line = 3, cex = 1.5)
mtext(side = 2, text = "by each chromosome", line = 2, cex = 1.5)


# Fit regression model
z = lm(Heritability ~ new, a)

# Plot the line of best fit
abline(z, lwd=2)

# Add calibrate library (is not needed if used with text)
library(calibrate)

# Add text
textxy(X = a$new, Y = a$Heritability, labs = c(1:22), cex = 1.5, pos = 3)

# Calculate confidence intervals
x = predict(z, interval = "confidence")

# Plot intervals
lines(a$new, x[,3], lwd = 2, lty = "dashed", col = "red")
lines(a$new, x[,2], lwd = 2, lty = "dashed", col = "red")