ID <- 1:6000
#disease <- rbinom(6000, 1, 0.2)

set.seed(1)
genetic_var1 <- rnorm(6000, 0, 0.02)
genetic_var2 <- rnorm(6000, 0, 0.02)
genetic_var3 <- rnorm(6000, 0, 0.02)
genetic_var4 <- rnorm(6000, 0, 0.02)
genetic_var5 <- rnorm(6000, 0, 0.1)
error <- rnorm(6000, 0, 1)
gender <- rbinom(6000, 1, 0.46)
disease <- genetic_var1 + genetic_var2 + genetic_var3 + genetic_var4 + genetic_var5 + gender + error

data <- as.data.frame(cbind(ID, disease, genetic_var1, error, gender))


a = lm(disease ~ genetic_var1 + genetic_var2 + genetic_var3 + genetic_var4 + genetic_var5, data = data)
summary(a)

par(mfrow=c(1,2))

plot(error)
plot(a$residuals)





x <- sample(0:2,6000,replace = T, prob = c(0.05,0.34,0.61))
