setwd("~/Documents")
# Load data and required packages.
library(nlme); library(foreign); library(gridExtra); library(ggplot2); library(reshape2); library(longpower); library(ellipse); library(readr); library(grid)
cream2017_longitudinal_2017_03_15_long = read.spss("~/cream2017_longitudinal_2017-03-15_long.sav")
cream2017_hits <- read_delim("~/cream2017_hits.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# Significance threshold adjusted for multiple comparisons. Different for loop combos allow for different number of SNPs.
significance = 0.05/146 
significance = 0.05/148 
significance = 0.05/151

# To calculate power or sample size (if needed).
lmmpower(your_model, 
         pct.change = -0.30, 
         t = seq(0,5,1), 
         power = 0.80)
# t is time points (I assume it's the number of visits in my case), pct.change is the effect size.

# Change cream data to data frame. Remove the original.
all = as.data.frame(cream2017_longitudinal_2017_03_15_long); rm(cream2017_longitudinal_2017_03_15_long)
##########################################################################################################################
# Linear model (for testing purposes).
#model1 = lm(mse_at_visit ~ sex + ReadingBinary + NumMyopicParents + rs11210537_A + ReadingBinary*rs11210537_A, data = all)
#summary(model1)
# Loop over all SNPs in 'snps'. Same fixed effects are used in each case. 
#my_model1 = lapply(1:151, function(x) lm(all[,10] ~ snps[,x] + all$sex 
#                                         + all$NumMyopicParents + all$ReadingBinary + snps[,x]*all$ReadingBinary))
# Obtain summary details of interest.
#summaries = lapply(my_model1, summary)
#details = lapply(summaries, function(x) x$coefficients[, c(1,2,3,4)])
#head(details) 
# If needed can add r-squared values.
#r_squared = sapply(summaries, function(x) c(r_sq = x$r.squared, 
#                                            adj_r_sq = x$adj.r.squared))
#head(r_squared)
##########################################################################################################################
# Descriptive statistics for different visits.
a = all[which(all$visit == 1),]; b = all[which(all$visit == 2),]; c = all[which(all$visit == 3),]
d = all[which(all$visit == 4),]; e = all[which(all$visit == 5),]
# Combine histograms for MSE at different visits.
par(mfrow = c(2,3))
hist(a$mse_at_visit, main = "MSE at visit 1", xlab = "MSE"); hist(b$mse_at_visit, main = "MSE at visit 2", xlab = "MSE")
hist(c$mse_at_visit, main = "MSE at visit 3", xlab = "MSE"); hist(d$mse_at_visit, main = "MSE at visit 4", xlab = "MSE")
hist(e$mse_at_visit, main = "MSE at visit 5", xlab = "MSE")

# To compute correlations for each visit, do the following:
x = a$mse_at_visit; z = b$mse_at_visit; y = c$mse_at_visit; w = d$mse_at_visit; q = e$mse_at_visit
x[is.na(x)] = 0; z[is.na(z)] = 0; y[is.na(y)] = 0; w[is.na(w)] = 0; q[is.na(q)] = 0; cor(x, z); cor(z, y); cor(y, w); cor(w, q)

# Number of non-missing individuals in each visit.
length(na.omit(a$mse_at_visit)); length(na.omit(b$mse_at_visit)); length(na.omit(c$mse_at_visit)); length(na.omit(d$mse_at_visit)); length(na.omit(e$mse_at_visit))

# Look at the gender proportions
#table(a$sex) # Not correct. Leads to overcounting.
table(a$total_visits)
# Combine histograms for age at different visits.
par(mfrow = c(2,3))
hist(a$age_at_visit, main = "Visit 1", xlab = "Age"); hist(b$age_at_visit, main = "Visit 2", xlab = "Age")
hist(c$age_at_visit, main = "Visit 3", xlab = "Age"); hist(d$age_at_visit, main = "Visit 4", xlab = "Age")
hist(e$age_at_visit, main = "Visit 5", xlab = "Age")

# Some additional basic plots.
par(mfrow = c(1,2))
plot (sort (all$mse_at_visit), pch=".", col = blues9, xlab = "Subject", ylab = "MSE");plot (mse_at_visit ~ visit,all, col = "red")
abline(a= 0.3412, b = -0.1371, col = "blue", lwd = 2)
#######################################################################################################################
# Fitting mixed effects models with different covariates.
# Change algorithm's convergence. Will take the model longer to run in order to reach local/global maximum.
ctrl = lmeControl(maxIter = 300, msMaxIter = 300, opt='optim')

# MY PREFERED MODELS! Maximum likelihood estimation methods.
result1 = lme(mse_at_visit ~ sex,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML" )
result2 = lme(mse_at_visit ~ sex + age_at_visit,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML" )
result3 = lme(mse_at_visit ~ sex + ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML" )
result4 = lme(mse_at_visit ~ age_at_visit + ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML" )
result5 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML" )
# Test if autoregressive covariance structure fits the data better. Spoiler alert (it does but very slightly).
# Changing correlation does not seem to affect estimates and predictions. Possibly because MSE correlation
# between different visits do not change much. Could maybe look at different structure types?
{result1 = lme(mse_at_visit ~ sex,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result2 = lme(mse_at_visit ~ sex + age_at_visit,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result3 = lme(mse_at_visit ~ sex + ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result4 = lme(mse_at_visit ~ age_at_visit + ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result5 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))}
# Example of model comparison. Model 4 fits the data the best. Sex seems to be not important.
anova(result4, result5)

# Check different interaction combinations using model 5.
{result1 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + sex:ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result2 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + age_at_visit:ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result3 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + sex:ReadingBinary +
              age_at_visit:ReadingBinary + sex:age_at_visit:ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
result6 = lme(mse_at_visit ~ ReadingBinary + age_at_visit + age_at_visit:ReadingBinary,
              random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
              na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))}
anova(result2, result6)
# Model 2 and model 6 are very similar but anova shows that model 6 is the best. I'll keep sex anyway because of previous studies.

# Do some diagnostics.
plot(result2)
plot(ranef(result2))
interaction.plot(all$ReadingBinary, all$age_at_visit, all$mse_at_visit) # Maybe in some universe it makes sense. 
# Could maybe look at one visit at a time. Grouping different ages somehow?
a$test = ifelse(a$age_at_visit < 7, 1, ifelse(a$age_at_visit >= 7 & a$age_at_visit < 8, 2, 3))
#Look at confidence intervals. Are they wide? Ideally should be as narrow as possible.
intervals(result2)
#plot(ACF(result2, form = ~ I(age_at_visit - 7.5) | alfred_ID1), alpha = 0.05)
#qqnorm(result2);qqline(result2) # Need to change

# Plot estimates for fixed effects.
{int = intervals(result2)
class(int$fixed)
kf <- dim(int$fixed)[1]
plot(int$fixed[,2], kf:1,
     xlab="Effect Size", ylab="", xlim=range(int$fixed), main = "Confidence Intervals for fixed effects",
     axes=FALSE, pch = 4, lwd = 2)
axis(1)
axis(2, kf:1, labels = c("Intercept", "sex", "Reading", "Age", "Interaction"))
box(which = "plot", bty = "l")
segments(int$fixed[,1], kf:1,
         int$fixed[,3], kf:1)}
#summary(result2)$tTable[,"p-value"] < 0.05
#summary(result2)$tTable[,"p-value"] < 0.01
#summary(result2)$tTable[,"p-value"] < 0.0001
#text(0.710379089,5.1,labels="*", family="mono", font=2, ps=8) # 0.0001
#text(-0.063740760,2.1,labels="*", family="mono", font=2, ps=8) # 0.0001
#text(-0.022966882,1.1,labels="*", family="mono", font=2, ps=8) # 0.0001
#text(0.133726559,3.1,labels="+", family="mono", font=2, ps=8) # 0.01
#legend("bottomright", legend=c("p-value < 0.01", "p-value < 0.0001"), 
#        bty = "n", pch = c("+", "*"))}








# NOT MY PREFERED MODELS (including polynomials)!

# In general, adding age polynomials leads to reduced residual variance. 
# Residual variance increases after adjusting for correlated data structure.
# But anova results show that models with autoregressive covariance structure are prefered.
{results1 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
results2 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
results3 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
results4 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3
               + age_at_visit:ReadingBinary,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
anova(results1, results2, results3, results4, results5, results6)}
# Results4 seems best.
# Estimates obtained using REML are identical.

#For predicting polynomials. 
{test0 = a[which(names(predict(results1)) %in% a$alfred_ID1),]
test1 = a[which(names(predict(results2)) %in% a$alfred_ID1),]
test2 = a[which(names(predict(results3)) %in% a$alfred_ID1),]
test3 = a[which(names(predict(results4)) %in% a$alfred_ID1),]
par(mfrow = c(3, 2))
plot(a$age_at_visit, a$mse_at_visit, main = "Visit 1", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-10, 8, by = 2))
box()
plot(test0$age_at_visit, predict(results1), col = "blueviolet", main = "Prediction with age^2", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-8, 8, by = 2))
box()
plot(test1$age_at_visit, predict(results2), col = "brown4", main = "Prediction with age^3", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-8, 8, by=2))
box()
plot(test2$age_at_visit, predict(results3), col = "mediumaquamarine", main = "Prediction with age^4", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-8, 8, by=2))
box()
plot(test3$age_at_visit, predict(results4), col = "orange1", main = "Prediction with age^4 + age:Reading", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-8, 8, by=2))
box()
plot(test0$age_at_visit, predict(results1), col = "blueviolet", main = "Comparison between 4 models", xlab = "Age", ylab = "MSE", axes = FALSE)
axis(side = 1, at = c(6:10))
axis(side = 2, at = seq(-8, 8, by=2))
box()
points(test1$age_at_visit, predict(results2), col = "brown4", xlab = "Age", ylab = "MSE")
points(test2$age_at_visit, predict(results3), col = "mediumaquamarine", xlab = "Age", ylab = "MSE")
points(test3$age_at_visit, predict(results4), col = "orange1", xlab = "Age", ylab = "MSE")}
dev.off()

# Plotting predicted polynomial lines.
{plot(1, type = "n", xlab = "Age", xlim = c(7, 19), ylim = c(-1, 1), ylab = "MSE", main = "Comparing polynomials", pch = 2)
test = predict(results3)
lines(test, col = "blueviolet")
test = predict(results4)
lines(test, col = "brown4")
test = predict(results5)
lines(test, col = "mediumaquamarine")
test = predict(results6)
lines(test, col = "orange1")
legend("bottomright", legend = c("age^2", "age^3", "age^4", "age:Read"), 
       col = c("blueviolet", "brown4", "mediumaquamarine", "orange"), bty = "n", lty = 1, lwd = 3)}

{plot(1, type = "n", xlab = "Age", xlim =  c(7, 60), ylim = c(-1, 4), ylab = "MSE", main = "Comparing polynomials", pch = 2)
test = predict(results3)
lines(test, col = "blueviolet")
test = predict(results4)
lines(test, col = "brown4")
test = predict(results5)
lines(test, col = "mediumaquamarine")
test = predict(results6)
lines(test, col = "orange1")
legend("bottomright", legend = c("age^2", "age^3", "age^4", "age:Read"), 
       col = c("blueviolet", "brown4", "mediumaquamarine", "orange"), bty = "n", lty = 1, lwd = 3)}

# Estimated effect sizes for each model.
{par(mfrow = c(1, 2), las = 1)
int = intervals(test7)
class(int$fixed)
kf = dim(int$fixed)[1]
plot(int$fixed[,2], kf:1,
     xlab = "Effect Size", ylab ="", xlim = range(int$fixed), main = "Confidence Intervals for fixed effects",
     axes = FALSE, pch = 4, lwd = 2)
axis(1)
axis(2, kf:1, labels = c("Intercept", "Read", "SNP", "Age", "Age^2", 
                         "Age^3", "Age^4", "int1", "int2", "int3"))
box(which = "plot", bty = "l")
segments(int$fixed[,1], kf:1,
         int$fixed[,3], kf:1)
int = intervals(test8)
class(int$fixed)
kf = dim(int$fixed)[1]
plot(int$fixed[,2], kf:1,
     xlab = "Effect Size", ylab = "", xlim=range(int$fixed), main = "Confidence Intervals for fixed effects",
     axes = FALSE, pch = 4, lwd = 2)
axis(1)
axis(2, kf:1, labels = c("Intercept", "Read", "SNP", "Age", "Age^2", 
                         "Age^3", "Age^4", "int1", "int2", "int3", "3way"))
box(which = "plot", bty = "l")
segments(int$fixed[,1], kf:1,
         int$fixed[,3], kf:1)
int = intervals(results3)
class(int$fixed)
kf = dim(int$fixed)[1]
plot(int$fixed[,2], kf:1,
     xlab = "Effect Size", ylab = "", xlim = range(int$fixed), main = "Confidence Intervals for fixed effects",
     axes = FALSE, pch = 4, lwd = 2)
axis(1)
axis(2, kf:1, labels = c("Intercept", "sex", "Reading", "Age", "Age^2", "Age^3", "Age^4"))
box(which = "plot", bty = "l")
segments(int$fixed[,1], kf:1,
         int$fixed[,3], kf:1)
int = intervals(results4)
class(int$fixed)
kf = dim(int$fixed)[1]
plot(int$fixed[,2], kf:1,
     xlab = "Effect Size", ylab = "", xlim = range(int$fixed), main = "Confidence Intervals for fixed effects",
     axes = FALSE, pch = 4, lwd = 2)
axis(1)
axis(2, kf:1, labels = c("Intercept", "sex", "Reading", "Age", "Age^2", "Age^3", "Age^4", "Age:Read"))
box(which = "plot", bty = "l")
segments(int$fixed[,1], kf:1,
         int$fixed[,3], kf:1)}

# Plotting polynomials.
{pol2 = function(Age) fixef(results1)[1] + fixef(results1)[2] + fixef(results1)[3] + fixef(results1)[4] + fixef(results1)[5]*Age^2
rol2 = function(Age) fixef(results2)[1] + fixef(results2)[2] + fixef(results2)[3] + fixef(results2)[4] + fixef(results2)[5]*Age^2 + fixef(results2)[6]*Age^3
tol2 = function(Age) fixef(results3)[1] + fixef(results3)[2] + fixef(results3)[3] + fixef(results3)[4] + fixef(results3)[5]*Age^2 + fixef(results3)[6]*Age^3 + fixef(results3)[7]*Age^4
sol2 = function(Age) fixef(results4)[1] + fixef(results4)[2] + fixef(results4)[3] + fixef(results4)[4] + fixef(results4)[5]*Age^2 + fixef(results4)[6]*Age^3 + fixef(results4)[7]*Age^4 + fixef(results4)[8]

plot(1, type="n", xlab="Age", xlim=c(7, 19), ylim=c(-10, 70), ylab = "MSE", main = "Comparing polynomials", pch = 2)
  curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE)
  curve(rol2, from = 7, to = 19, col="brown4", lwd=2, add = TRUE)
  curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE)
  curve(sol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE)
  legend("topleft", legend=c("Age^2", "Age^3", "Age^4", "Age:Read"), 
         col = c("blueviolet", "brown4", "mediumaquamarine", "orange1"), bty = "n", lty = 1, lwd = 3)}
dev.off()
{plot(1, type="n", xlab="Age", xlim=c(7, 19), ylim=c(-15, 1), ylab = "MSE", main = "Comparing polynomials", pch = 2)
  curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE)
  curve(rol2, from = 7, to = 19, col="brown4", lwd=2, add = TRUE)
  legend("topleft", legend=c("Age^2", "Age^3"), 
         col = c("blueviolet", "brown4"), bty = "n", lty = 1, lwd = 3)}
{plot(1, type="n", xlab="Age", xlim=c(7, 80), ylim=c(-25, 10), ylab = "MSE", main = "Comparing polynomials", pch = 2)
  curve(pol2, from = 7, to = 80, col="blueviolet", lwd=2, add = TRUE)
  curve(rol2, from = 7, to = 80, col="brown4", lwd=2, add = TRUE)
  legend("topleft", legend=c("Age^2", "Age^3"), 
         col = c("blueviolet", "brown4"), bty = "n", lty = 1, lwd = 3)}


# Different age and reading interaction combos.
{results7 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3
               + age_at_visit:ReadingBinary + placeholder:ReadingBinary,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
results8 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3
               + age_at_visit:ReadingBinary + placeholder:ReadingBinary + placeholder2:ReadingBinary,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
results9 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3
               + age_at_visit:ReadingBinary + placeholder:ReadingBinary + placeholder2:ReadingBinary + placeholder3:ReadingBinary,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
anova(results6, results7, results8, results9)} # age:reading seems to be the best
# Separated by reading exposure.
Reading_no = function(Age) fixef(results4)[1] + fixef(results4)[2] + fixef(results4)[4] + fixef(results4)[5]*Age^2 + fixef(results4)[6]*Age^3 + fixef(results4)[7]*Age^4
Reading_yes = function(Age) fixef(results4)[1] + fixef(results4)[2] + fixef(results4)[3] + fixef(results4)[4] + fixef(results4)[5]*Age^2 + fixef(results4)[6]*Age^3 + fixef(results4)[7]*Age^4 + fixef(results4)[8]
{plot(1, type="n", xlab="Age", xlim=c(7, 19), ylim=c(20, 70), ylab = "MSE", main = "Age:Read", pch = 2)
  curve(Reading_no, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE)
  curve(Reading_yes, from = 7, to = 19, col="mediumaquamarine", lwd=3, add = TRUE)
  legend("topleft", legend=c("Reading = Low Exposure", "Reading = High Exposure"), 
         col = c("blueviolet", "mediumaquamarine"), bty = "n", lty = 1, lwd = 3)}

results7 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + placeholder2 + placeholder3
               + age_at_visit:ReadingBinary + sex:ReadingBinary + age_at_visit:sex + sex:age_at_visit:ReadingBinary,
               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
anova(results6, results7)

# Look at correlation structure for best model.
{colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(length(cov2cor(vcov(results1)))) # or colors = colramp(100)
par(mfrow = c(2, 2))
my.plotcorr(cov2cor(vcov(results1)), col = colors[5*cov2cor(vcov(results1)) + 6], main='Predictor correlations',upper.panel="number")
colors = colramp(length(cov2cor(vcov(results2))))
my.plotcorr(cov2cor(vcov(results2)), col = colors[5*cov2cor(vcov(results2)) + 6], main='Predictor correlations',upper.panel="number")
colors = colramp(length(cov2cor(vcov(results3))))
my.plotcorr(cov2cor(vcov(results3)), col = colors[5*cov2cor(vcov(results3)) + 6], main='Predictor correlations',upper.panel="number")
colors = colramp(length(cov2cor(vcov(results4))))
my.plotcorr(cov2cor(vcov(results4)), col = colors[5*cov2cor(vcov(results4)) + 6], main='Predictor correlations',upper.panel="number")}
dev.off()









# Look at calculated values. Generic - needs adjustment for each model.
#test = summary(results4)
#test$...

# Some cool plots. Plots are made for different visits.
{par(mfrow = c(3, 2))
plot(a$age_at_visit[a$ReadingBinary == "1"], a$mse_at_visit[a$ReadingBinary == "1"], 
     col = "blue", xlab = "Age at visit 1", ylab = "MSE", main = "MSE vs Age, Reading", ylim = c(-6.5,8))
points(a$age_at_visit[a$ReadingBinary == "0"], a$mse_at_visit[a$ReadingBinary == "0"], col = "red", pch = 6)
legend("topright", legend=c("Reading High Exposure", "Reading Low Exposure"), col = c("blue", "red"), pch=c(1,6), bty = "n")

plot(b$age_at_visit[b$ReadingBinary == "1"], b$mse_at_visit[b$ReadingBinary == "1"], 
     col = "blue", xlab = "Age at visit 2", ylab = "MSE", main = "MSE vs Age, Reading")
points(b$age_at_visit[b$ReadingBinary == "0"], b$mse_at_visit[b$ReadingBinary == "0"], col = "red", pch = 6)
legend("topright", legend=c("Reading High Exposure", "Reading Low Exposure"), col = c("blue", "red"), pch=c(1,6), bty = "n")

plot(c$age_at_visit[c$ReadingBinary == "1"], c$mse_at_visit[c$ReadingBinary == "1"], 
     col = "blue", xlab = "Age at visit 3", ylab = "MSE", main = "MSE vs Age, Reading")
points(c$age_at_visit[c$ReadingBinary == "0"], c$mse_at_visit[c$ReadingBinary == "0"], col = "red", pch = 6)
legend("topright", legend=c("Reading High Exposure", "Reading Low Exposure"), col = c("blue", "red"), pch=c(1,6), bty = "n")

plot(d$age_at_visit[d$ReadingBinary == "1"], d$mse_at_visit[d$ReadingBinary == "1"], 
     col = "blue", xlab = "Age at visit 4", ylab = "MSE", main = "MSE vs Age, Reading")
points(d$age_at_visit[d$ReadingBinary == "0"], d$mse_at_visit[d$ReadingBinary == "0"], col = "red", pch = 6)
legend("topright", legend=c("Reading High Exposure", "Reading Low Exposure"), col = c("blue", "red"), pch=c(1,6), bty = "n")

plot(e$age_at_visit[e$ReadingBinary == "1"], e$mse_at_visit[e$ReadingBinary == "1"], 
     col = "blue", xlab = "Age at visit 5", ylab = "MSE", main = "MSE vs Age, Reading")
points(e$age_at_visit[e$ReadingBinary == "0"], e$mse_at_visit[e$ReadingBinary == "0"], col = "red", pch = 6)
legend("bottomright", legend=c("Reading High Exposure", "Reading Low Exposure"), col = c("blue", "red"), pch=c(1,6), bty = "n")}

# For all visits combined.
{plot(all$age_at_visit[all$ReadingBinary == "0"], all$mse_at_visit[all$ReadingBinary == "0"], 
     col = "blue", xlab = "Age for all visits", ylab = "MSE", main = "MSE vs Age, Reading", pch = 2)
points(all$age_at_visit[all$ReadingBinary == "1"], all$mse_at_visit[all$ReadingBinary == "1"], col = "red", pch = 6)
legend("bottomleft", legend=c("Reading Low Exposure", "Reading High Exposure"), col = c("red", "blue"), pch=c(2, 6), bty = "n")
abline(a = 0.7167944, b = -0.0637380, col = "black", lwd = 3)
abline(a = 0.8513595, b = -0.0867026, col = "blue", lwd = 3)}

# Separated by sex. should not be used because the model is not the best. Again points to the fact that sex is not significant.
#results10 = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit
#               + sex:age_at_visit:ReadingBinary + sex:age_at_visit + age_at_visit:ReadingBinary + sex:ReadingBinary,
#               random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
#               na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
#{plot(all$age_at_visit[all$ReadingBinary == "0" & all$sex == "1"], all$mse_at_visit[all$ReadingBinary == "0" & all$sex == "1"], 
#     col = "blue", xlab = "Age for all visits", ylab = "MSE", main = "MSE vs Age, Reading", pch = 2)
#points(all$age_at_visit[all$ReadingBinary == "1" & all$sex == "1"], all$mse_at_visit[all$ReadingBinary == "1" & all$sex == "1"], col = "red", pch = 6)
#points(all$age_at_visit[all$ReadingBinary == "0" & all$sex == "0"], all$mse_at_visit[all$ReadingBinary == "0" & all$sex == "0"], col = "green", pch = 4)
#points(all$age_at_visit[all$ReadingBinary == "1" & all$sex == "0"], all$mse_at_visit[all$ReadingBinary == "1" & all$sex == "0"], col = "black", pch = 8)
#par(mfrow = c(2, 2))
#plot(1, type="n", xlab="Age", xlim=c(7, 18), ylim=c(-1, 1), ylab = "MSE", main = "MSE vs Sex:Age:Reading", pch = 2)
#legend("bottomleft", legend=c("Reading = 0, sex = 1", "Reading = 1, sex = 1", "Reading = 0, sex = 0", "Reading = 1, sex = 0"), 
#       col = c("blue", "red", "green", "black"), fill = c("blue", "red", "green", "black"))
#abline(a = 0.7159207, b = -0.0656346, col = "blue", lwd = 1)
#abline(a = 0.806664, b = -0.0766685, col = "red", lwd = 1)
#abline(a = 0.6816053, b = -0.0656346, col = "green", lwd = 1)
#abline(a = 0.7723486, b = -0.0656346, col = "black", lwd = 1)
#plot(1, type="n", xlab="Age", xlim=c(7, 25), ylim=c(-2, 1), ylab = "MSE", main = "MSE vs Sex:Age:Reading", pch = 2)
#legend("bottomleft", legend=c("Reading = 0, sex = 1", "Reading = 1, sex = 1", "Reading = 0, sex = 0", "Reading = 1, sex = 0"), 
#       col = c("blue", "red", "green", "black"), fill = c("blue", "red", "green", "black"))
#abline(a = 0.7159207, b = -0.0656346, col = "blue", lwd = 1)
#abline(a = 0.806664, b = -0.0766685, col = "red", lwd = 1)
#abline(a = 0.6816053, b = -0.0656346, col = "green", lwd = 1)
#abline(a = 0.7723486, b = -0.0656346, col = "black", lwd = 1)
#plot(1, type="n", xlab="Age", xlim=c(7, 40), ylim=c(-3, 1), ylab = "MSE", main = "MSE vs Sex:Age:Reading", pch = 2)
#legend("bottomleft", legend=c("Reading = 0, sex = 1", "Reading = 1, sex = 1", "Reading = 0, sex = 0", "Reading = 1, sex = 0"), 
#       col = c("blue", "red", "green", "black"), fill = c("blue", "red", "green", "black"))
#abline(a = 0.7159207, b = -0.0656346, col = "blue", lwd = 1)
#abline(a = 0.806664, b = -0.0766685, col = "red", lwd = 1)
#abline(a = 0.6816053, b = -0.0656346, col = "green", lwd = 1)
#abline(a = 0.7723486, b = -0.0656346, col = "black", lwd = 1)
#plot(1, type="n", xlab="Age", xlim=c(7, 60), ylim=c(-5, 1), ylab = "MSE", main = "MSE vs Sex:Age:Reading", pch = 2)
#legend("bottomleft", legend=c("Reading = 0, sex = 1", "Reading = 1, sex = 1", "Reading = 0, sex = 0", "Reading = 1, sex = 0"), 
#       col = c("blue", "red", "green", "black"), fill = c("blue", "red", "green", "black"))
#abline(a = 0.7159207, b = -0.0656346, col = "blue", lwd = 1)
#abline(a = 0.806664, b = -0.0766685, col = "red", lwd = 1)
#abline(a = 0.6816053, b = -0.0656346, col = "green", lwd = 1)
#abline(a = 0.7723486, b = -0.0656346, col = "black", lwd = 1)}










" age, age^2, age^3, age^4, ignore NumMyopicParents?, add SNPs only after all jazz with interactions, Interaction with sex, 
Interaction with ReadingBinary, Interaction with Outdoors, Higher order interactions: age:reading:snp, sex:reading:snp etc.
Look at visit as possible random effect. Is corAR1 correlation structure used correct?"

### Adding SNP as fixed effect.
# MY WAY (NO POLYNOMIALS)!

snps = all[,c(11:26)] # That's for testing. SNP 159 doesn't work.
# A mixed effect models with SNPs and fixed effects. NO INTERACTIONS! 
{ctrl = lmeControl(maxIter = 300, msMaxIter = 300, opt='optim')
results2 = NULL
start.time = Sys.time()
for (i in 1:ncol(snps)) { 
  variable = snps[,i]
  results2[[i]] = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + variable, 
                                     random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                                     na.action = na.omit, control = ctrl, method = "ML")
}
end.time = Sys.time(); end.time - start.time}

# Use table function to check how many copies of each snp there are.
lapply(all[,11:161], table) # 1, 32, 38, 50, 62, 73, 142, 148, 150, 159 have only 2 possibilities

# Create a list with SNP names.
test = list("rs10910076", "rs11210537", "rs11589487", "rs1237670", "rs11802995", "rs1556867", "rs2225986", 
            "rs1858001", "rs2745953", "rs11118367", "rs12405776", "rs6753137", "rs28658452", "rs17032696", "rs41393947", "rs10187371")


# Plot some SNPs.
{par(mfrow = c(4,4), oma = c(3, 3, 3, 0), mar = c(1.8, 1.8, 0.5, 1)) # c(bottom, left, top, right) 
for (j in 1:length(test)){ 
  for (i in 1:length(results2)){
    plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1, 0.5), pch = 2, axes = FALSE,  las = 1)
    axis(side = 1, at = 1:19, labels = if (i %/% 13 == 1)levels(19) else FALSE)
    axis(side = 2, labels = (i %% 4 == 1))
    box(which = "plot", bty = "l")
    abline(a = fixef(results2[[i]])[1] + fixef(results2[[i]])[2] + fixef(results2[[i]])[4],
         b = fixef(results2[[i]])[3], col="blueviolet", lwd=2, lty = 2)
    abline(a = fixef(results2[[i]])[1] + fixef(results2[[i]])[2] + fixef(results2[[i]])[4] + fixef(results2[[i]])[5],
         b = fixef(results2[[i]])[3], col="orange1", lwd=2, lty = 3)
    abline(a = fixef(results2[[i]])[1] + fixef(results2[[i]])[2] + fixef(results2[[i]])[4] + 2*fixef(results2[[i]])[5],
         b = fixef(results2[[i]])[3], col="mediumaquamarine", lwd=2, lty = 4)  
    title(main = "Including SNPs as fixed/ No interactions/ No polynomials (example for first 16 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
        cex.main = 2, cex.lab = 2)
    legend("bottomleft", legend=c("0", "1", "2"), 
         col = c("blueviolet", "orange1", "mediumaquamarine"), bty = "n", lty = c(2, 3, 4),
         title = test[i], title.adj = 1.5, adj = 6)
    mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
  }
}}

# This is for the SNP(159) that can not be used in a loop.
test = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A, # This is SNP 159
           random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
           na.action = na.omit, control = ctrl, method = "ML")

plot(1, type="n", xlab="Age", xlim=c(7, 19), ylim=c(-10, 10), ylab = "MSE", main = "Including SNPs", pch = 2)
abline(a = fixef(test)[1] + fixef(test)[2] + fixef(test)[4], b = fixef(test)[3], col="blueviolet", lwd=2)
abline(a = fixef(test)[1] + fixef(test)[2] + fixef(test)[4] + fixef(test)[5], b = fixef(test)[3], col="brown4", lwd=2)
abline(a = fixef(test)[1] + fixef(test)[2] + fixef(test)[4] + 2*fixef(test)[5], b = fixef(test)[3], col="mediumaquamarine", lwd=2)

# Plotting CIs for different SNPs.
{par(mfrow = c(4, 4), oma = c(3, 3, 2.5, 0), mar = c(1.8, 1.8, 0.5, 1), las = 1)
  for (i in 1:length(results2)){
    int = intervals(results2[[i]])
    class(int$fixed)
    kf <- dim(int$fixed)[1]
    plot(int$fixed[,2], kf:1,
         xlab="Effect Size", ylab="", xlim=range(int$fixed),
         axes=FALSE, pch = 4, lwd = 2)
    axis(1, labels = if (i %/% 13 == 1)levels(19) else FALSE)
    axis(2, kf:1, labels = (i %% 4 == 1))
    box(which = "plot", bty = "l")
    segments(int$fixed[,1], kf:1,
             int$fixed[,3], kf:1)
    title(main = "CIs", outer = TRUE, ylab = "Fixed effect estimates", xlab = "Effect Size", line = 1,
          cex.main = 2, cex.lab = 2)
  }}





# Check which interaction models are the best.
{test1 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A +
           age_at_visit:rs2150458_A, # This is SNP 159
           random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
           na.action = na.omit, control = ctrl, method = "ML")
test2 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A +
            ReadingBinary:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")
test3 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A +
            sex:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")
test4 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A + age_at_visit:sex +
            age_at_visit:sex:rs2150458_A + age_at_visit:rs2150458_A + sex:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")
test5 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A + age_at_visit:ReadingBinary +
            age_at_visit:ReadingBinary:rs2150458_A + age_at_visit:rs2150458_A + ReadingBinary:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")
test6 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A + sex:ReadingBinary +
            sex:ReadingBinary:rs2150458_A + sex:rs2150458_A + ReadingBinary:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")
test7 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs2150458_A + age_at_visit:sex +
            age_at_visit:ReadingBinary + age_at_visit:rs2150458_A + sex:ReadingBinary + sex:rs2150458_A +
            ReadingBinary:rs2150458_A + age_at_visit:sex:ReadingBinary + age_at_visit:sex:rs2150458_A +
            sex:ReadingBinary:rs2150458_A + age_at_visit:ReadingBinary:rs2150458_A +
            age_at_visit:sex:ReadingBinary:rs2150458_A, # This is SNP 159
            random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
            na.action = na.omit, control = ctrl, method = "ML")}
# Test 4, 5 and 7 are better than the model with no interaction. 5 is arguably the best.
dev.off()
my.plotcorr(cov2cor(vcov(test5)), col = colors[5*cov2cor(vcov(test5)) + 6], main='Predictor correlations',upper.panel="number")



# Let's use test 5.
{snps = all[,c(12:20)] # SNP 11(1) can not be looped
tests3 = NULL
for (i in 1:ncol(snps)) { 
  variable = snps[,i]
  tests3[[i]] = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + variable + age_at_visit:ReadingBinary +
                    age_at_visit:ReadingBinary:variable + age_at_visit:variable + ReadingBinary:variable, # This is SNP 159
                    random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                    na.action = na.omit, control = ctrl, method = "ML")
}}

# Interaction plots.
test = list("rs10910076", "rs11210537", "rs11589487", "rs1237670", "rs11802995", "rs1556867", "rs2225986", 
            "rs1858001", "rs2745953")
{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
     for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1, 0.5), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      # Reading 0, SNP 0
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2],
             b = fixef(tests3[[i]])[3], col="blueviolet", lwd=2, lty = 2)
      # Reading 0, SNP 1
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[5],
             b = fixef(tests3[[i]])[3] + fixef(tests3[[i]])[7], col="orange1", lwd=2, lty = 3)
      # Reading 0, SNP 2
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + 2*fixef(tests3[[i]])[5],
             b = fixef(tests3[[i]])[3] + 2*fixef(tests3[[i]])[7], col="mediumaquamarine", lwd=2, lty = 4)
      # Reading 1, SNP 0
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4],
             b = fixef(tests3[[i]])[3] + fixef(tests3[[i]])[6], col="blueviolet", lwd=4, lty = 2)
      # Reading 1, SNP 1
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5] + fixef(tests3[[i]])[8],
             b = fixef(tests3[[i]])[3] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[9], col="orange1", lwd=4, lty = 3)
      # Reading 1, SNP 2
      abline(a = fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + 2*fixef(tests3[[i]])[5] + 2*fixef(tests3[[i]])[8],
             b = fixef(tests3[[i]])[3] + 2*fixef(tests3[[i]])[7] + 2*fixef(tests3[[i]])[9], col="mediumaquamarine", lwd=4, lty = 4)  
      title(main = "Interaction with Reading and Age (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "orange1", "mediumaquamarine", "blueviolet", "orange1", "mediumaquamarine"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "grey40")
    }
  }}







# Let's plot 2nd degree polynomials + Read:SNP.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:ReadingBinary,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ I(age_at_visit - 7.5) | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1.25, 0.25), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7]
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7]
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7]
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7]
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Second degree polynomial with SNP and reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
}}


# Let's plot 2nd degree polynomials + Read:Age:SNP.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:age_at_visit:ReadingBinary + age_at_visit:variable + age_at_visit:ReadingBinary + variable:ReadingBinary,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1.25, 0.25), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10]
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10]
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[8]
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10]
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10]
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Second degree polynomial with Age:SNP:Reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
  }}

# Let's plot 3rd degree polynomials + Read:Age:SNP.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:age_at_visit:ReadingBinary + age_at_visit:variable + age_at_visit:ReadingBinary + variable:ReadingBinary + placeholder2,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1.25, 0.25), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[7]*Age^3
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7]*Age^3 + fixef(tests3[[i]])[8]
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7]*Age^3 + 2*fixef(tests3[[i]])[8]
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[7]*Age^3
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7]*Age^3 + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7]*Age^3 + 2*fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + 2*fixef(tests3[[i]])[11]
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Third degree polynomial with Age:SNP:Reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
  }}

# Let's plot 2nd degree polynomials + Read:Age^2:SNP.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:placeholder:ReadingBinary + placeholder:variable + placeholder:ReadingBinary + variable:ReadingBinary,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-1.5, 0.5), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10]
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10]
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[8]
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10]
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10]
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Second degree polynomial with Age:SNP:Reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
  }}



# Let's include 3rd degree polynomials.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:age_at_visit:ReadingBinary + age_at_visit:variable + age_at_visit:ReadingBinary + variable:ReadingBinary + placeholder2,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-40, 20), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[11]*Age^3
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[11]*Age^3
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Third degree polynomial with Age:SNP:Reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
  }}


# What about 4th degree polynomials.
{snps = all[,c(12:20)] # 11(1) can not be looped
  tests3 = NULL
  for (i in 1:ncol(snps)) { 
    variable = snps[,i]
    tests3[[i]] = lme(mse_at_visit ~ sex + ReadingBinary + age_at_visit + placeholder + variable
                      + variable:age_at_visit:ReadingBinary + age_at_visit:variable + age_at_visit:ReadingBinary + variable:ReadingBinary + placeholder2 + placeholder3,
                      random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                      na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(form=~ 1 | alfred_ID1))
  }}

{par(mfrow = c(3,3), oma = c(3, 3, 3.2, 0), mar = c(1.8, 1.8, 0.6, 1),  las = 1) # c(bottom, left, top, right) 
  for (j in 1:length(test)){ 
    for (i in 1:length(tests3)){
      plot(1, type="n", xlab="", xlim=c(7, 19), ylim=c(-100, 100), pch = 4, axes = FALSE, main = test[i])
      axis(side = 1, at = 1:19, labels = if (i %/% 7 == 1)levels(19) else FALSE)
      axis(side = 2, labels = (i %% 3 == 1))
      box(which = "plot", bty = "l")
      pol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      rol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      qol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      sol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      tol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + fixef(tests3[[i]])[6] + fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + fixef(tests3[[i]])[9] + fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      uol2 <- function(Age) fixef(tests3[[i]])[1] + fixef(tests3[[i]])[2] + fixef(tests3[[i]])[3] + fixef(tests3[[i]])[4] + fixef(tests3[[i]])[5]*Age^2 + 2*fixef(tests3[[i]])[6] + 2*fixef(tests3[[i]])[7] + fixef(tests3[[i]])[8] + 2*fixef(tests3[[i]])[9] + 2*fixef(tests3[[i]])[10] + fixef(tests3[[i]])[11]*Age^3 + fixef(tests3[[i]])[12]*Age^4
      
      curve(pol2, from = 7, to = 19, col="blueviolet", lwd=2, add = TRUE, lty = 2)
      curve(rol2, from = 7, to = 19, col="mediumaquamarine", lwd=2, add = TRUE, lty = 3)
      curve(qol2, from = 7, to = 19, col="orange1", lwd=2, add = TRUE, lty = 4)
      curve(sol2, from = 7, to = 19, col="blueviolet", lwd=4, add = TRUE, lty = 2)
      curve(tol2, from = 7, to = 19, col="mediumaquamarine", lwd=4, add = TRUE, lty = 3)
      curve(uol2, from = 7, to = 19, col="orange1", lwd=4, add = TRUE, lty = 4)
      title(main = "Fourth degree polynomial with Age:SNP:Reading interaction (example for first 9 SNPs)", outer = TRUE, ylab = "MSE", xlab = "Age", line = 1,
            cex.main = 2, cex.lab = 2)
      legend("bottomleft", legend=c("0", "1", "2", "0", "1", "2"), 
             col = c("blueviolet", "mediumaquamarine", "orange1", "blueviolet", "mediumaquamarine", "orange1"), 
             bty = "n", lty = c(2, 3, 4, 2, 3, 4), lwd = c(2, 2, 2, 4, 4, 4),
             title = "Low High", title.adj = 0.34, adj = 5.8, ncol = 2)
      mtext(letters[i], side = 3, line = -1.5, adj = 1, cex = 1.5, col = "black")
    } 
  }}



















# Test for pairwise SNP:SNP effect.
snps = all[,c(11:161)]
pairwise_SNP = NULL
start.time = Sys.time()
plot(1, type="n", xlab="Age", xlim=c(7, 15), ylim=c(-1, 1), ylab = "MSE", main = "MSE vs Age:Reading:variable", pch = 2)
for (i in 1:ncol(snps)) { 
  for (j in 2:ncol(snps)) {
    snp1 = snps[,i]
    snp2 = snps[,j]
    pairwise_SNP[[i]] = lme(mse_at_visit ~ poly(I(age_at_visit - 7.5),4) + snp1 + snp2 + 
                      snp1:snp2 + I(age_at_visit - 7.5):snp1 + I(age_at_visit - 7.5):snp2 +
                      I(age_at_visit - 7.5):snp1:snp2, 
                      random=~I(age_at_visit - 7.5) | alfred_ID1, correlation = corCAR1(form = ~ visit| alfred_ID1), 
                      na.action = na.omit, method="ML", data=all)
    print(sapply(pairwise_SNP, function(x) summary(x)$tTable[11,5]))
  }
}
end.time = Sys.time(); end.time - start.time

test = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs11210537_A + rs11589487_A 
           , 
           random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
           na.action = na.omit, control = ctrl)

# What is 7.5 actually doing? Represents deviation from the mean at 7.5 years?
# Is it used to remove correlation?
# alfred_ID1 refers to grouping factor(i.e. subject). So age differs randomly between individuals. 
head(results2,2) # To look at some elements of the summary
lapply(results2, summary) # To look at the summary of all models
summary(results2[[2]]) # To extract specific model

# Plot residuals (example) (basic).
test = lapply(results2, plot.lme)  
do.call(grid.arrange, c(test, nrow=1))  #or# do.call(grid.arrange, c(test, nrow=2))

# Regression coefficients.
lapply(results2, fixef)
lapply(results2, ranef)
lapply(results2, coef) # Note in this case only intercept changes for each subject. Slopes remain the same.

# Extract p-values for interaction effects. 
test = sapply(pairwise_SNP, function(x) summary(x)$tTable[11,5]) # Row numbers depend on the number of variables fitted.
test[test < significance] = significance
length(which(test < significance))
length(which(test < 0.05))

# Confidence intervals. On average, deviations are not large except for some of the fixed effects (in my model).
lapply(results2, intervals)
#lapply(results2, plot(intervals))

# In the book, authors suggest using anova's statistics, instead of t-values and p-values from the mixed model,
# to assess fixed parameters, if "cell means" are not estimated. P.S. for interactions in my models there seem
# to be no big dfference. 
test = lapply(results2, anova.lme) 
test = sapply(test, '[[', 4) # To look only at p-values
#test[6,] # To look only at p-values for interaction. Row numbers depend on the number of variables fitted.
which(test[6,] < significance)
which(test < 0.05)

# Alternative. Not quite correct. Loops over.
#variable = snps[,i]
#my_model2 = lapply(1:151, function(x) lme(mse_at_visit ~ sex + NumMyopicParents + ReadingBinary, 
#                                          random = ~ 1|variable/ReadingBinary, data = all, na.action = na.omit))
#summaries = lapply(my_model2, summary)

# Another alternative. Not quite correct. Loops over.
#myfunc = function(x){
#  out = with(all, lme(mse_at_visit ~ sex + NumMyopicParents + ReadingBinary, 
#                       random = ~ 1|variable/ReadingBinary, na.action = na.exclude))
#}
#test = lapply(snps, myfunc)
#########################################################################################################################
# Next model is for fitting multiple SNPs as fixed in one model. Not complete. Could possibly add multiple interactions??
my_model2 = lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + rs11210537_A + rs11589487_A + 
                  rs1237670_G + rs11802995_C + rs1556867_T + rs2225986_A + rs1858001_G + rs2745953_T + rs11118367_T + rs12405776_T +
                  rs12451582_A + rs28488643_G + rs35879249_CA + rs4793501_C + rs6420484_A + rs10853531_A + rs12965607_G + rs7253703_C,
                  random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, na.action = na.omit, control = ctrl)
summary(my_model2)
#########################################################################################################################
# Scripts from the book: "Mixed models in S and S-plus"
system.file("scripts", package = "nlme")
list.files(system.file("scripts", package = "nlme"))
file.show(system.file("scripts", "ch03.R", package = "nlme"))

# Some codes that show residual expectations under different models.
# For normality.
par(mfrow=c(3, 3) )
for (i in 1:9) qqnorm (rnorm (50) ) # Normal
for (i in 1:9) qqnorm (exp (rnorm (50) ) ) # Lognormal
for (i in 1:9) qqnorm (rcauchy (50) ) # Cauchy
for (i in 1:9) qqnorm (runif (50) ) # Uniform
par (mfrow=c(1, 1) )

shapiro.test (residuals (result54))

# For variance.
par(mfrow=c (3, 3))
for (i in 1:9) plot (1:50, rnorm (50)) # Constant variance
for (i in 1:9) plot (1:50, (1:50)*rnorm(50)) # Strong nonconstant variance
for (i in 1:9) plot (1:50, sqrt ((1:50))*rnorm(50)) # Mild nonconstant variance
for (i in 1:9) plot(1:50, cos ((1:50)*pi/25)+rnorm(50)) # Nonlinearity
par(mfrow=c (1, 1))









# An attempt to plot standard errors. Doesn't quite work. Did not attempt to solve.
fit <- lme(mse_at_visit ~ sex + age_at_visit + ReadingBinary + placeholder +
                           sex:ReadingBinary:placeholder,
                     random = ~ I(age_at_visit - 7.5) | alfred_ID1, data = all, 
                     na.action = na.omit, control = ctrl, method = "ML", cor=corAR1(0.744447, form=~ 1 | alfred_ID1))
prd <- data.frame(age_at_visit = seq(from = range(all$age_at_visit)[1], to = range(all$age_at_visit)[2], length.out = 100, na.action(na.omit)))
err <- predict(fit, newdata = prd, se.fit = TRUE)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

ggplot(prd, aes(x = hp, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") +
  geom_point(data = mtcars, aes(x = hp, y = mpg))






# Function for creating correlation plots.
my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
  # this is a modified version of the plotcorr function from the ellipse package
  # this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
  # diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
  # digits specifies the number of digits after the . to round to
  # unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
  # modified by Esteban Buz
  if (!require('ellipse', quietly = TRUE, character = TRUE)) {
    stop("Need the ellipse library")
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}

