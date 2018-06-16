seven$rounded = round(seven$Age, 0)
counts <- table(seven$Sex, seven$YearOfBirth)
barplot(counts, main="Distibution of sex as a function of year of birth",
        xlab="Year of Birth", ylab = "Count", ylim = c(0,5500), col=c("darkblue","red"),
        legend = c("Male", "Female"), beside=TRUE, bty = "n", cex.lab = 1.5, cex.main = 1.5)

table(seven$UniEdu, seven$Sex)
chisq.test(table(seven$UniEdu, seven$Sex))


s1 = seven[which(seven$EduYearsHigh == 0 & seven$Sex == 0),]
s2 = seven[which(seven$EduYearsHigh == 0 & seven$Sex == 1),]
s3 = seven[which(seven$EduYearsHigh == 1 & seven$Sex == 0),]
s4 = seven[which(seven$EduYearsHigh == 1 & seven$Sex == 1),]
t.test(s1$predicted_avMSE, s2$predicted_avMSE)
t.test(s3$predicted_avMSE, s4$predicted_avMSE)
t.test(s1$predicted_avMSE, s3$predicted_avMSE)
t.test(s2$predicted_avMSE, s4$predicted_avMSE)

ci(seven$predicted_avMSE)
ci(s1$predicted_avMSE)
ci(s2$predicted_avMSE)
ci(s3$predicted_avMSE)
ci(s4$predicted_avMSE)

labels = LETTERS[1:4]
par(mfrow=c(2,2), oma = c(2, 3, 2, 0), mar = c(2.6, 1.8, 0.7, 1)) # c(bottom, left, top, right) 
hist(s1$avMSE,breaks = 30);hist(s2$avMSE,breaks = 30);hist(s3$avMSE,breaks = 30);hist(s4$avMSE,breaks = 30)
par(mfrow=c(2,2), oma = c(1, 3, 1, 0), mar = c(3, 1.8, 1, 1))
qqnorm(s1$avMSE); qqline(s1$avMSE, col = 2);qqnorm(s2$avMSE); qqline(s2$avMSE, col = 2)
qqnorm(s3$avMSE); qqline(s3$avMSE, col = 2);qqnorm(s4$avMSE); qqline(s4$avMSE, col = 2)

vioplot(s1$avMSE,s2$avMSE,s3$avMSE,s4$avMSE, names=c("Male no UniEdu", "Female no UniEdu", "Male UniEdu", "Female UniEdu"), col="gold")
title("avMSE differences between sex and education")

eduageold=seven$EducationAgeOLD
eduageold=na.omit(eduageold)
length(eduageold)








