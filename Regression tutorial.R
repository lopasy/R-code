summary(lm(avMSE~Age,master))$residuals
master = master[,c(3:10,19:22)]
master = data.frame(scale(master))
plot(avMSE~Age,master)
abline(0,1)
model = lm(avMSE~Age,master)
abline(coef(model), lty = 5)
cor.test(master$avMSE, master$Age)$estimate

cor(master)

#master$av = master$avMSE
#model2 = lm(avMSE~av,master)
#abline(coef(model2), lty = 5)
#cor.test(master$avMSE, master$av)

master$regressed = cor.test(master$avMSE, master$Age)$estimate*master$Age
master$model = coef(model)[1] + coef(model)[2]*master$Age

summary(model)$r.squared 
cor.test(master$avMSE, master$Age)$estimate^2
1-deviance (model) /sum ((master$avMSE-mean (master$avMSE))^2)
sd(unscaled$avMSE)^2

"Suppose you want to predict y. If you do not know x, then your best prediction is mean(y) 
but the variability in this prediction is high. If you do know x, then your prediction will 
be given by the regression fit. This prediction will be less variable provided there is some 
relationship between x and y. R2 is one minus the ratio of the sum of squares for these two 
predictions. Thus, for perfect predictions the ratio will be zero and R2 will be one."


# Residual standard error
sqrt(sum(summary(lm(avMSE~Age,master))$residuals^2)/summary(lm(avMSE~Age,master))$df[2])
# Residual sum of squares
sum(summary(lm(avMSE~Age,master))$residuals^2)
deviance(model)
# Mean squared error
sum(summary(lm(avMSE~Age,master))$residuals^2)/nrow(master)
mean(model$residuals^2)
anova(model)[[3]][2]
# Total sum of squares
var(master$avMSE) * (nrow(master)-1)
sum((master$avMSE-mean(master$avMSE))^2)
# Sum of squares regression
var(master$avMSE) * (nrow(master)-1) - sum(summary(lm(avMSE~Age,master))$residuals^2)
# R2
1 - (sum(summary(lm(avMSE~Age,master))$residuals^2)/(var(master$avMSE) * (nrow(master)-1)))
#df
nrow(master)-length(coef(model))
# corelation
cov(master$avMSE,master$Age)/sqrt(var(master$avMSE)*var(master$Age))
# if the data is scaled, correlaition = covariance for two variables
cov(master)
#standard error
se <- function(x) sd(x)/sqrt(length(x))
se(master$Age)
#standard error of the estimate given by Var^(??^)=??^2(X'X), where ??^2 is mean squared error
betaHat <- solve(t(master$Age) %*% master$Age) %*% t(master$Age) %*% master$avMSE # Estimate for age
var_betaHat <- anova(model)[[3]][2] * solve(t(master$Age) %*% master$Age) # Variance of age
sqrt(diag(var_betaHat)) # se for age

#########################################
### If needed to produce permutations ###
#########################################
summary(lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master))$fstat
fstats = numeric(4000)
for (i in 1:4000){
  ge = lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master)
  fstats [i] < - summary(ge)$fstat[1]
}
length(fstats[fstats > 306.9984])/4000



### Joint CI ###
library(ellipse)
plot(ellipse(lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master),c(2,4)),type="l",xlim=c(-0.5,0.5), ylim = c(-0.3,0.1))
points(0,0) 
points (coef (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master)) [2], 
        coef (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master)) [4], pch=18)
abline (v=confint (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master)) [2,], lty=2)
abline (h=confint (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master)) [4,], lty=2)

### Always check correlation of two or more variables
summary (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master), corr=TRUE) $corr
### Always check because the signs might differ


### Orthogonality
summary (lm(avMSE~Age+Sex+UniEdu+Array+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,master), corr=TRUE)
### Only controlled experiments can be truly orthogonal. Observational studies will not.



# Examples of OLS residual distribution violations
{par(mfrow=c (4, 3))
for (i in 1:3) plot (1:50, rnorm (50))
for (i in 1:3) plot (1:50, (1:50)*rnorm(50))
for (i in 1:3) plot (1:50, sqrt ((1:50))*rnorm(50))
for (i in 1:3) plot(1:50, cos ((1:50)*pi/25)+rnorm(50))}

# Check normality assumption
qqnorm(residuals(summary(lm(avMSE~poly(I(Age),1)+UniEdu+Sex+Array+PC1+PC2+PC3+PC4+PC5,master))), ylab="Residuals")
qqline(residuals(summary(lm(avMSE~Age,master))))

# Examples of OLS normality violations
{par(mfrow=c(4, 3) )
for (i in 1:3) qqnorm (rnorm (50) )
for (i in 1:3) qqnorm (exp (rnorm (50) ) )
for (i in 1:3) qqnorm (rcauchy (50) )
for (i in 1:3) qqnorm (runif (50) )}





###Unscaled data###
unscaled = unscaled[1:1000,c(3:10,19:22)]
plot(avMSE~Age,unscaled)
model = lm(avMSE~Age,unscaled)
abline(coef(model), lty = 5)
cor.test(unscaled$avMSE, unscaled$Age)

unscaled$model = coef(model)[1] + coef(model)[2]*unscaled$Age
