library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(data.table)

dataD            <- as.data.frame(fread(file="D:/SPSS/UKBB2017/2017-01-20/ukb_phens2017-02-12.csv"))

dataD$FID        <- dataD$IID
dataD$avMSE      <- ifelse(dataD$outlierMSE==1,NA,dataD$avMSE)

data2            <- dataD[ which( dataD$Ethnicity=="White" & dataD$AgeSpexWear>=0 & dataD$avMSE!="NA" ) , ]

# Loop to find optimum polynomial order for age
threshold                <- 0
poly_count               <- 0
while(threshold < 0.05){
  poly_count               <- poly_count + 1
  first_model              <- lm(avMSE ~ Sex + poly(Age, poly_count),     data=data2)
  second_model             <- lm(avMSE ~ Sex + poly(Age, poly_count + 1), data=data2)
  stat_summary             <- anova(first_model, second_model)
  first_model_summary      <- summary(first_model)
  threshold                <- unlist(stat_summary)[[12]] #this gets the p-value from the anova; I had to look-up the syntax on the web
}
optimum_age <- poly_count

# Loop to find optimum polynomial order for age-of-onset-spex
threshold                <- 0
poly_count               <- 0
graph_data               <- as.data.frame(matrix(nrow=2,ncol=30))
names(graph_data)        <- c("Polynomial_order", "Rsquared_of_model")
while(threshold < 0.05){
  poly_count               <- poly_count + 1
  first_model              <- lm(avMSE ~ Sex + poly(Age, optimum_age) + poly(AgeSpexWear, poly_count),     data=data2)
  second_model             <- lm(avMSE ~ Sex + poly(Age, optimum_age) + poly(AgeSpexWear, poly_count + 1), data=data2)
  stat_summary             <- anova(first_model, second_model)
  first_model_summary      <- summary(first_model)
  graph_data[poly_count,1] <- poly_count
  graph_data[poly_count,2] <- first_model_summary$adj.r.squared
  threshold                <- unlist(stat_summary)[[12]] #this gets the p-value from the anova; I had to look-up the syntax on the web
}
optimum_spex <- poly_count

plot(graph_data$Polynomial_order, graph_data$Rsquared_of_model, xlab="Polynomial order", ylab="R-squared of model")
lines(graph_data$Polynomial_order, graph_data$Rsquared_of_model)

best_model <- lm(avMSE ~ Sex + poly(Age, optimum_age) + poly(AgeSpexWear, optimum_spex), data=data2)
summary(best_model)

# Predict
dataD$predicted_avMSE <- predict(best_model,newdata=dataD)
dataD$predicted_avMSE <- ifelse((is.na(dataD$avMSE) & dataD$Ethnicity=="White"), dataD$predicted_avMSE, NA)
