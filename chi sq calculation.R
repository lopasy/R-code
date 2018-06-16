# Convert your output to chi-squared values
# For z-scores, just square the t statistic
chisq <- reports$STAT^2

# For p-values, calculate chi-squared statistic
chisq <- qchisq(data$pval,1, lower.tail = F)

# Calculate lambda gc (??gc)
median(chisq)/qchisq(0.5,1)
