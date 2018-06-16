quantile(seven2$predicted_avMSE, c(0.1,0.15,0.2,0.25,0.3,0.7,0.75,0.8,0.85,0.9))

seven_high = seven2[which(seven2$predicted_avMSE < -2.7732361),]
seven_low = seven2[which(seven2$predicted_avMSE > 1.2937264),]
seven_merged = rbind(seven_high,seven_low)
summary(lm(predicted_avMSE~UniEdu+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array, seven_merged))$r.squared

seven_high = seven2[which(seven2$predicted_avMSE < -2.3776335),]
seven_low = seven2[which(seven2$predicted_avMSE > 1.1136830),]
dim(seven_merged)
seven_merged = rbind(seven_high,seven_low)
dim(seven_merged)
summary(lm(predicted_avMSE~UniEdu+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array, seven_merged))$r.squared

seven_high = seven2[which(seven2$predicted_avMSE < -2.0054681),]
seven_low = seven2[which(seven2$predicted_avMSE > 0.9558076),]
seven_merged = rbind(seven_high,seven_low)
dim(seven_merged)
summary(lm(predicted_avMSE~UniEdu+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array, seven_merged))$r.squared

seven_high = seven2[which(seven2$predicted_avMSE < -1.6151796),]
seven_low = seven2[which(seven2$predicted_avMSE > 0.8116784),]
seven_merged = rbind(seven_high,seven_low)
dim(seven_merged)
summary(lm(predicted_avMSE~UniEdu+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array, seven_merged))$r.squared

seven_high = seven2[which(seven2$predicted_avMSE < -1.2012650),]
seven_low = seven2[which(seven2$predicted_avMSE > 0.6747276),]
seven_merged = rbind(seven_high,seven_low)
dim(seven_merged)
summary(lm(predicted_avMSE~UniEdu+Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Array, seven_merged))$r.squared









