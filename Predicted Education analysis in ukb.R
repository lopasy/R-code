# Exploring UKB data. Modelling and finding new ways to visualize data.
library(ggplot2)

# Create dataframe for ggplots
z = as.data.frame(table(seven$YearOfBirth,seven$EduYearsHigh))
colnames(z)=c("Year","Years of Education","Individuals")

# One example of distribution of education with respect to Year of Birth.
ggplot(z, aes(x=Year, y=Individuals, fill=`Years of Education`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Distribution of Years of Education with respect to year of birth")



# Plot predicted_avMSE as a function of Year of Birth etc. for those who have Education degree and those who don't.
library(lattice)
xyplot(predicted_predicted_avMSE ~ YearOfBirth, group=UniEdu,data = seven)
xyplot(predicted_predicted_avMSE ~ BirthMonth, group=UniEdu,data = seven)
a = seven[which(seven$UniEdu == 0),]
b = seven[which(seven$UniEdu == 1),]
plot(density(a$predicted_predicted_avMSE), col="blue", lwd=5,main="Distribution of predicted_avMSE in Degree holders and non-holders");lines(density(b$predicted_predicted_avMSE),lwd=5,col="darkred")

# To recode MSE as categorical. One hyperopic, one mildly myopic, one pathogenical myopia.
seven$categorical[seven$predicted_avMSE >-6 & seven$predicted_avMSE <= -0.5] <- 1; seven$categorical[seven$predicted_avMSE <= -6] <- 2; seven$categorical[seven$predicted_avMSE > -0.5] <- 0 


z = as.data.frame(table(seven$YearOfBirth, seven$HighestQualification))
colnames(z)=c("Year","Highest Qualification","Individuals")
ggplot(z, aes(x=Year, y=Individuals, fill=`Highest Qualification`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Highest Qualification with respect to year of birth")

z = as.data.frame(table(seven$YearOfBirth, seven$EducationAgeOLD))
colnames(z)=c("Year","Years","Individuals")
ggplot(z, aes(x=Year, y=Individuals, fill=`Years`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black")) +
  ggtitle("EducationAgeOLD with respect to year of birth")


# Create YearOfBirth bins
{seven$bins[seven$YearOfBirth == 1936 | seven$YearOfBirth == 1937 | seven$YearOfBirth == 1938] <- 1
  seven$bins[seven$YearOfBirth == 1939 | seven$YearOfBirth == 1940 | seven$YearOfBirth == 1941] <- 2
  seven$bins[seven$YearOfBirth == 1942 | seven$YearOfBirth == 1943 | seven$YearOfBirth == 1944] <- 3
  seven$bins[seven$YearOfBirth == 1945 | seven$YearOfBirth == 1946 | seven$YearOfBirth == 1947] <- 4
  seven$bins[seven$YearOfBirth == 1948 | seven$YearOfBirth == 1949 | seven$YearOfBirth == 1950] <- 5
  seven$bins[seven$YearOfBirth == 1951 | seven$YearOfBirth == 1952 | seven$YearOfBirth == 1953] <- 6
  seven$bins[seven$YearOfBirth == 1954 | seven$YearOfBirth == 1955 | seven$YearOfBirth == 1956] <- 7
  seven$bins[seven$YearOfBirth == 1957 | seven$YearOfBirth == 1958 | seven$YearOfBirth == 1959] <- 8
  seven$bins[seven$YearOfBirth == 1960 | seven$YearOfBirth == 1961 | seven$YearOfBirth == 1962] <- 9
  seven$bins[seven$YearOfBirth == 1963 | seven$YearOfBirth == 1964 | seven$YearOfBirth == 1965] <- 10
  seven$bins[seven$YearOfBirth == 1966 | seven$YearOfBirth == 1967 | seven$YearOfBirth == 1968 | seven$YearOfBirth == 1969 | seven$YearOfBirth == 1970] <- 11}

{plot(1, type="n", xlab="", ylab="", xlim=c(-10, 7), ylim=c(0, 0.85), main = "Distribution of predicted_avMSE as a function of binned Years of Education")
  a = seven[which(seven$bins == 1),]; lines(density(a$predicted_avMSE), col = "deepskyblue", lwd = 3)
  a = seven[which(seven$bins == 2),]; lines(density(a$predicted_avMSE), col = "deepskyblue2", lwd = 3)
  a = seven[which(seven$bins == 3),]; lines(density(a$predicted_avMSE), col = "deepskyblue4", lwd = 3)
  a = seven[which(seven$bins == 4),]; lines(density(a$predicted_avMSE), col = "dodgerblue4", lwd = 3)
  a = seven[which(seven$bins == 5),]; lines(density(a$predicted_avMSE), col = "firebrick1", lwd = 3)
  a = seven[which(seven$bins == 6),]; lines(density(a$predicted_avMSE), col = "firebrick3", lwd = 3)
  a = seven[which(seven$bins == 7),]; lines(density(a$predicted_avMSE), col = "firebrick4", lwd = 3)
  a = seven[which(seven$bins == 8),]; lines(density(a$predicted_avMSE), col = "gold", lwd = 3)
  a = seven[which(seven$bins == 9),]; lines(density(a$predicted_avMSE), col = "gold2", lwd = 3)
  a = seven[which(seven$bins == 10),]; lines(density(a$predicted_avMSE), col = "goldenrod2", lwd = 3)
  a = seven[which(seven$bins == 11),]; lines(density(a$predicted_avMSE), col = "goldenrod4", lwd = 3)
  legend(x=-10, y=0.85, bty = "n", fill = c("deepskyblue","deepskyblue2","deepskyblue4","dodgerblue4","firebrick1","firebrick3",
                                            "firebrick4","gold","gold2","goldenrod2","goldenrod4"), 
         legend = c("1936-1938","1939-1941","1942-1944","1945-1947","1948-1950","1951-1953","1954-1956","1957-1959","1960-1962",
                    "1961-1965","1966-1970"))}

'z = as.data.frame(table(seven$bins, seven$predicted_avMSE))
colnames(z)=c("Year","Years","Individuals")
ggplot(z, aes(x=Years, colour=`Year`)) +
geom_density() +
theme_bw() +
theme(panel.border = element_rect(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
axis.line = element_line(colour = "black")) +
ggtitle("EducationAgeOLD with respect to year of birth")'

# Distribution of Sex as a f(YOB).
z = as.data.frame(table(seven$YearOfBirth,seven$Sex))
colnames(z)=c("Year","Sex","Individuals")
z$Sex = as.character(z$Sex)
z$Sex[z$Sex == 0] = "Male"; z$Sex[z$Sex == 1] = "Female"
ggplot(z, aes(x=Year, y = Individuals, fill=`Sex`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black"),legend.title=element_blank()) +
  ggtitle("Distribution of Sex as a function of year of birth")


# Split data by Year and calculate min, max, mean, median predicted_avMSE, SD etc.
years = 1937:1969
df2 = list()
fd2 = list()
for (i in years){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$YearOfBirth == i & seven$Sex == 0),]
  df2[[name]] = df
  fd = seven[which(seven$YearOfBirth == i & seven$Sex == 1),]
  fd2[[name]] = fd
}

# Function to calculate standard error
se <- function(x) sd(x)/sqrt(length(x))


predicted_avMSE = matrix(ncol = 5, nrow = 33)
colnames(predicted_avMSE) = c("Year", "Mean predicted_avMSE Male", "Mean predicted_avMSE Female", "SE Male", "SE Female")
i=1937
j=1
for (j in seq_along(1:33)){ 
  for (i in years){
    predicted_avMSE[j,1] = i
    name = paste('Year:',i,sep='')
    predicted_avMSE[j,2] = mean(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,3] = mean(fd2[[name]]$predicted_avMSE)
    predicted_avMSE[j,4] = se(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,5] = se(fd2[[name]]$predicted_avMSE)
    j=j+1
    i=i+1
  }
}
predicted_avMSE = as.data.frame(predicted_avMSE)
predicted_avMSE$YearMale = predicted_avMSE$Year-0.15
predicted_avMSE$YearFemale = predicted_avMSE$Year+0.15
par(mfrow=c(1,2))
{plot(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female`, col = "darkblue", pch = 16, main = "Decrease of predicted avMSE with respect to year of birth", 
      ylab = "avMSE in diopters", xlab = "Year", ylim = c(-2,1), type = "b",
      cex.main = 1.5, cex.axis = 1.35, cex.lab = 1.35)
  points(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male`, col="red", pch = 16, type = "b")
  legend(x="topright",fill = c("darkblue","red"), legend = c("Female","Male"), bty = "n", cex = 1.35)
  segments(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,
           predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`)
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  
  segments(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,
           predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`)
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
  counts <- table(seven$Sex, seven$YearOfBirth)
  barplot(counts, main="Distribution of sex with respect to year of birth",
          xlab="Year of Birth", ylab = "Number of participants", ylim = c(0,5500), col=c("red","darkblue"),
          beside=TRUE, bty = "n", cex.lab = 1.35, cex.main = 1.5, cex.axis = 1.35, cex.names=1.35)
  legend(x = "topright", legend = c("Female", "Male"), fill=c("darkblue","red"), cex = 1.35, bty = "n")}



# If use bins instead
years = 1:11
df2 = list()
fd2 = list()
for (i in years){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$bins == i & seven$Sex == 0),]
  df2[[name]] = df
  fd = seven[which(seven$bins == i & seven$Sex == 1),]
  fd2[[name]] = fd
}

predicted_avMSE = matrix(ncol = 5, nrow = 11)
colnames(predicted_avMSE) = c("Bin", "Mean predicted_avMSE Male", "Mean predicted_avMSE Female", "SE Male", "SE Female")
for (j in seq_along(1:33)){ 
  for (i in years){
    predicted_avMSE[j,1] = i
    name = paste('Year:',i,sep='')
    predicted_avMSE[j,2] = mean(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,3] = mean(fd2[[name]]$predicted_avMSE)
    predicted_avMSE[j,4] = se(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,5] = se(fd2[[name]]$predicted_avMSE)
    j=j+1
    i=i+1
  }
}
predicted_avMSE = as.data.frame(predicted_avMSE)
predicted_avMSE$YearMale = predicted_avMSE$Bin-0.15
predicted_avMSE$YearFemale = predicted_avMSE$Bin+0.15
{plot(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female`, col = "deepskyblue", pch = 16, main = "Decrase of predicted_avMSE as a function of Year of Birth binned", 
      ylab = "predicted_avMSE", xlab = "Bin", ylim = c(-1,1.5))
  points(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male`, col="firebrick3", pch = 16)
  legend(x=10,y=1.5,fill = c("deepskyblue","firebrick3"), legend = c("Female","Male"),bty = "n")
  segments(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,
           predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`)
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  
  segments(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,
           predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`)
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
}


# What about highest qualification
qual = 0:6
df2 = list()
fd2 = list()
for (i in qual){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$HighestQualification == i & seven$Sex == 0),]
  df2[[name]] = df
  fd = seven[which(seven$HighestQualification == i & seven$Sex == 1),]
  fd2[[name]] = fd
}

predicted_avMSE = matrix(ncol = 5, nrow = 7)
colnames(predicted_avMSE) = c("Year", "Mean predicted_avMSE Male", "Mean predicted_avMSE Female", "SE Male", "SE Female")
for (j in seq_along(0:6)){ 
  for (i in qual){
    predicted_avMSE[j,1] = i
    name = paste('Year:',i,sep='')
    predicted_avMSE[j,2] = mean(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,3] = mean(fd2[[name]]$predicted_avMSE)
    predicted_avMSE[j,4] = se(df2[[name]]$predicted_avMSE)
    predicted_avMSE[j,5] = se(fd2[[name]]$predicted_avMSE)
    j=j+1
    i=i+1
  }
}
predicted_avMSE = as.data.frame(predicted_avMSE)
predicted_avMSE$YearMale = predicted_avMSE$Year-0.15
predicted_avMSE$YearFemale = predicted_avMSE$Year+0.15
{plot(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female`, col = "deepskyblue", pch = 16, main = "Predicted_avMSE for individuals with different Education Qualifications", 
      ylab = "predicted avMSE", xlab = "Highest Qualification", ylim = c(-0.8,0.3), xlim = c(-0.3,6.3))
  points(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male`, col="firebrick3", pch = 16)
  legend(x=5.5,y=-0.7,fill = c("deepskyblue","firebrick3"), legend = c("Female","Male"), bty = "n")
  segments(predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,
           predicted_avMSE$YearFemale, predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` - predicted_avMSE$`SE Female`)
  segments(predicted_avMSE$YearFemale-epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`,predicted_avMSE$YearFemale+epsilon,predicted_avMSE$`Mean predicted_avMSE Female` + predicted_avMSE$`SE Female`)
  
  segments(predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,
           predicted_avMSE$YearMale, predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
  epsilon <- 0.1
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` - predicted_avMSE$`SE Male`)
  segments(predicted_avMSE$YearMale-epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`,predicted_avMSE$YearMale+epsilon,predicted_avMSE$`Mean predicted_avMSE Male` + predicted_avMSE$`SE Male`)
}



# Based on education status
years = 1937:1969
df2 = list();df3 = list()
fd2 = list();fd3 = list()
for (i in years){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$YearOfBirth == i & seven$Sex == 0),]
  noedu = df[which(df$EduYearsHigh == 0),]; edu = df[which(df$EduYearsHigh == 1),]
  df2[[name]] = noedu; df3[[name]] = edu
  fd = seven[which(seven$YearOfBirth == i & seven$Sex == 1),]
  noedus = fd[which(fd$EduYearsHigh == 0),]; edus = fd[which(fd$EduYearsHigh == 1),]
  fd2[[name]] = noedus; fd3[[name]] = edus 
}



avmse = matrix(ncol = 9, nrow = 33)
colnames(avmse) = c("Year", "Mean avMSE Male noedu", "Mean avMSE Male edu", "SE Male noedu", "SE Male edu",
                    "Mean avMSE Female noedu", "Mean avMSE Female edu", "SE Female noedu", "SE Female edu")
i=1937
j=1
for (j in seq_along(1:33)){ 
  for (i in years){
    avmse[j,1] = i
    name = paste('Year:',i,sep='')
    avmse[j,2] = mean(df2[[name]]$predicted_avMSE)
    avmse[j,3] = mean(df3[[name]]$predicted_avMSE)
    avmse[j,4] = se(df2[[name]]$predicted_avMSE)
    avmse[j,5] = se(df3[[name]]$predicted_avMSE)
    avmse[j,6] = mean(fd2[[name]]$predicted_avMSE)
    avmse[j,7] = mean(fd3[[name]]$predicted_avMSE)
    avmse[j,8] = se(fd2[[name]]$predicted_avMSE)
    avmse[j,9] = se(fd3[[name]]$predicted_avMSE)
    j=j+1
    i=i+1
  }
}
avmse = as.data.frame(avmse)
avmse$YearMale = avmse$Year-0.15
avmse$YearFemale = avmse$Year+0.15
par(mfrow=c(1,2))
{plot(avmse$YearFemale, avmse$`Mean avMSE Female noedu`, col = "darkblue", pch = 16, main = "Females", 
      ylab = "predicted avMSE in diopters", xlab = "Year of birth", ylim = c(-2.2,1.2), type = "b", cex.main = 1.5, cex.axis = 1.35, cex.lab = 1.35)
  points(avmse$YearMale, avmse$`Mean avMSE Female edu`, col="red", pch = 16, type = "b")
  legend(x="topright",fill = c("darkblue","red"), legend = c("No degree","Degree"),bty = "n", cex = 1.35)
  segments(avmse$YearFemale, avmse$`Mean avMSE Female noedu` - avmse$`SE Female noedu`,
           avmse$YearFemale, avmse$`Mean avMSE Female noedu` + avmse$`SE Female noedu`)
  epsilon <- 0.1
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female noedu` - avmse$`SE Female noedu`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female noedu` - avmse$`SE Female noedu`)
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female noedu` + avmse$`SE Female noedu`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female noedu` + avmse$`SE Female noedu`)
  
  segments(avmse$YearMale, avmse$`Mean avMSE Female edu` - avmse$`SE Female edu`,
           avmse$YearMale, avmse$`Mean avMSE Female edu` + avmse$`SE Female edu`)
  
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Female edu` - avmse$`SE Female edu`,avmse$YearMale+epsilon,avmse$`Mean avMSE Female edu` - avmse$`SE Female edu`)
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Female edu` + avmse$`SE Female edu`,avmse$YearMale+epsilon,avmse$`Mean avMSE Female edu` + avmse$`SE Female edu`)
  

plot(avmse$YearMale, avmse$`Mean avMSE Male noedu`, col = "darkblue", pch = 16, main = "Males", 
      ylab = "predicted avMSE in diopters", xlab = "Year of birth", ylim = c(-2.2,1.3), type = "b", cex.main = 1.5, cex.axis = 1.35, cex.lab = 1.35)
  points(avmse$YearFemale, avmse$`Mean avMSE Male edu`, col="red", pch = 16, type = "b")
  legend(x="topright",fill = c("darkblue","red"), legend = c("No degree","Degree"),bty = "n", cex = 1.35)
  segments(avmse$YearMale, avmse$`Mean avMSE Male noedu` - avmse$`SE Male noedu`,
           avmse$YearMale, avmse$`Mean avMSE Male noedu` + avmse$`SE Male noedu`)
  epsilon <- 0.15
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male noedu` - avmse$`SE Male noedu`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male noedu` - avmse$`SE Male noedu`)
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male noedu` + avmse$`SE Male noedu`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male noedu` + avmse$`SE Male noedu`)
  
  segments(avmse$YearFemale, avmse$`Mean avMSE Male edu` - avmse$`SE Male edu`,
           avmse$YearFemale, avmse$`Mean avMSE Male edu` + avmse$`SE Male edu`)
  
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Male edu` - avmse$`SE Male edu`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Male edu` - avmse$`SE Male edu`)
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Male edu` + avmse$`SE Male edu`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Male edu` + avmse$`SE Male edu`)
}

par(cex.lab=1.35, cex.axis=1.35) 
vioplot(s1$predicted_avMSE,s3$predicted_avMSE,s2$predicted_avMSE,s4$predicted_avMSE, 
        names=c("Male no degree", "Male with degree", "Female no degree", "Female with degree"), col="gold")
title(ylab = "predicted avMSE in diopters", cex.lab = 1.35, main = "Distribution of avMSE", cex.main = 1.5, cex.axis = 1.35, cex.lab = 1.35)

par(mfrow=c(2,2), oma = c(2, 3, 2, 0), mar = c(2.6, 1.8, 0.7, 1)) # c(bottom, left, top, right) 
hist(s1$predicted_avMSE,breaks = 30, main = "Male - Uni");hist(s2$predicted_avMSE,breaks = 30, main = "Female - Uni")
hist(s3$predicted_avMSE,breaks = 30, main = "Male + Uni");hist(s4$predicted_avMSE,breaks = 30, main = "Female + Uni")


s1 = seven[which(seven$Sex == 0),]
s2 = seven[which(seven$Sex == 1),]
ci(s1$predicted_avMSE)
ci(s2$predicted_avMSE)
ci(seven$predicted_avMSE)

s1 = seven[which(seven$EduYearsHigh == 0 & seven$Sex == 0),]
s2 = seven[which(seven$EduYearsHigh == 0 & seven$Sex == 1),]
s3 = seven[which(seven$EduYearsHigh == 1 & seven$Sex == 0),]
s4 = seven[which(seven$EduYearsHigh == 1 & seven$Sex == 1),]
