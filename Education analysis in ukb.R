# Exploring UKB data. Modelling and finding new ways to visualize data.
library(ggplot2)

# Create dataframe for ggplots
z = as.data.frame(table(seven$YearOfBirth,seven$EducationAgeOLD))
colnames(z)=c("Year of Birth","Years of Education","Individuals")

# One example of distribution of education with respect to Year of Birth.
ggplot(z, aes(x=Year, y=Individuals, fill=`Years of Education`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Distribution of Years of Education with respect to year of birth")



# Plot avMSE as a function of Year of Birth etc. for those who have Education degree and those who don't.
library(lattice)
xyplot(avMSE ~ YearOfBirth, group=UniEdu,data = seven)
xyplot(avMSE ~ BirthMonth, group=UniEdu,data = seven)
a = seven[which(seven$UniEdu == 0),]
b = seven[which(seven$UniEdu == 1),]
plot(density(a$avMSE), col="blue", lwd=5,main="Distribution of avMSE in Degree holders and non-holders");lines(density(b$avMSE),lwd=5,col="darkred")

# To recode MSE as categorical. One hyperopic, one mildly myopic, one pathogenical myopia.
seven$categorical[seven$avMSE >-6 & seven$avMSE <= -0.5] <- 1; seven$categorical[seven$avMSE <= -6] <- 2; seven$categorical[seven$avMSE > -0.5] <- 0 


z = as.data.frame(table(seven$YearOfBirth, seven$HighestQualification))
colnames(z)=c("Year of Birth","Highest Qualification","Individuals")
ggplot(z, aes(x=Year, y=Individuals, fill=`Highest Qualification`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Highest Qualification with respect to year of birth")

z = as.data.frame(table(seven$YearOfBirth, seven$EducationAgeOLD))
colnames(z)=c("Year of Birth","Years","Individuals")
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

{plot(1, type="n", xlab="", ylab="", xlim=c(-25, 20), ylim=c(0, 0.50), main = "Distribution of avMSE as a function of binned Years of Education")
  a = seven[which(seven$bins == 1),]; lines(density(a$avMSE), col = "darkblue", lwd = 3)
  a = seven[which(seven$bins == 2),]; lines(density(a$avMSE), col = "darkblue2", lwd = 3)
  a = seven[which(seven$bins == 3),]; lines(density(a$avMSE), col = "darkblue4", lwd = 3)
  a = seven[which(seven$bins == 4),]; lines(density(a$avMSE), col = "dodgerblue4", lwd = 3)
  a = seven[which(seven$bins == 5),]; lines(density(a$avMSE), col = "firebrick1", lwd = 3)
  a = seven[which(seven$bins == 6),]; lines(density(a$avMSE), col = "red", lwd = 3)
  a = seven[which(seven$bins == 7),]; lines(density(a$avMSE), col = "firebrick4", lwd = 3)
  a = seven[which(seven$bins == 8),]; lines(density(a$avMSE), col = "gold", lwd = 3)
  a = seven[which(seven$bins == 9),]; lines(density(a$avMSE), col = "gold2", lwd = 3)
  a = seven[which(seven$bins == 10),]; lines(density(a$avMSE), col = "goldenrod2", lwd = 3)
  a = seven[which(seven$bins == 11),]; lines(density(a$avMSE), col = "goldenrod4", lwd = 3)
  legend(x=-25,y=0.5,fill = c("darkblue","darkblue2","darkblue4","dodgerblue4","firebrick1","red",
                              "firebrick4","gold","gold2","goldenrod2","goldenrod4"), 
         legend = c("1936-1938","1939-1941","1942-1944","1945-1947","1948-1950","1951-1953","1954-1956","1957-1959","1960-1962",
                    "1961-1965","1966-1970"))}

'z = as.data.frame(table(seven$bins, seven$avMSE))
colnames(z)=c("Year of Birth","Years","Individuals")
ggplot(z, aes(x=Years, colour=`Year`)) +
geom_density() +
theme_bw() +
theme(panel.border = element_rect(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
axis.line = element_line(colour = "black")) +
ggtitle("EducationAgeOLD with respect to year of birth")'

# Distribution of Sex as a f(YOB).
z = as.data.frame(table(seven$YearOfBirth,seven$Sex))
colnames(z)=c("Year of Birth","Sex","Individuals")
z$Sex = as.character(z$Sex)
z$Sex[z$Sex == 0] = "Male"; z$Sex[z$Sex == 1] = "Female"
ggplot(z, aes(x=Year, y = Individuals, fill=`Sex`)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  theme(panel.border = element_rect(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5), 
        axis.line = element_line(colour = "black"),legend.title=element_blank()) +
  ggtitle("Distribution of Sex as a function of year of birth")

# POOLING
seven$YearOfBirth[seven$YearOfBirth == 1936] = 1940; seven$YearOfBirth[seven$YearOfBirth == 1937] = 1940
seven$YearOfBirth[seven$YearOfBirth == 1938] = 1940; seven$YearOfBirth[seven$YearOfBirth == 1939] = 1940
seven$YearOfBirth[seven$YearOfBirth == 1970] = 1969

# Split data by Year and calculate min, max, mean, median avMSE, SD etc.
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


avmse = matrix(ncol = 5, nrow = 33)
colnames(avmse) = c("Year of Birth", "Mean avMSE Male", "Mean avMSE Female", "SE Male", "SE Female")
i=1937
j=1
for (j in seq_along(1:33)){ 
  for (i in years){
    avmse[j,1] = i
    name = paste('Year:',i,sep='')
    avmse[j,2] = mean(df2[[name]]$avMSE)
    avmse[j,3] = mean(fd2[[name]]$avMSE)
    avmse[j,4] = se(df2[[name]]$avMSE)
    avmse[j,5] = se(fd2[[name]]$avMSE)
    j=j+1
    i=i+1
  }
}
avmse = as.data.frame(avmse)
avmse$YearMale = avmse$Year-0.15
avmse$YearFemale = avmse$Year+0.15
{plot(avmse$YearFemale, avmse$`Mean avMSE Female`, col = "darkblue", pch = 16, main = "Decrase of avMSE as a function of Year of Birth", 
      ylab = "avMSE in Dioptres", xlab = "Year of Birth", ylim = c(-1.2,3.3), type = "b", cex = 1.4, cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4)
  points(avmse$YearMale, avmse$`Mean avMSE Male`, col="red", pch = 16, type = "b", cex = 1.4)
  legend("topright",fill = c("darkblue","red"), legend = c("Female","Male"),bty = "n", cex = 1.4)
  segments(avmse$YearFemale, avmse$`Mean avMSE Female` - avmse$`SE Female`,
           avmse$YearFemale, avmse$`Mean avMSE Female` + avmse$`SE Female`)
  epsilon <- 0.1
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female` - avmse$`SE Female`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female` - avmse$`SE Female`)
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female` + avmse$`SE Female`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female` + avmse$`SE Female`)
  
  segments(avmse$YearMale, avmse$`Mean avMSE Male` - avmse$`SE Male`,
           avmse$YearMale, avmse$`Mean avMSE Male` + avmse$`SE Male`)
  epsilon <- 0.1
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male` - avmse$`SE Male`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male` - avmse$`SE Male`)
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male` + avmse$`SE Male`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male` + avmse$`SE Male`)
}

{x<-seq(1937,1969,1)
y = avmse$`Mean avMSE Female`

fit = lm(y~x+I(x^2))
prd <- data.frame(x = seq(1937, 2000, by = 1))
prd$pred = predict(fit, prd)
lines(sort(x), fitted(fit)[order(x)], col='blue', type='b') 
cor(y,fitted(fit)[order(x)])


y = avmse$`Mean avMSE Male`

fit = lm(y~x+I(x^2))
prd <- data.frame(x = seq(1937, 1969, by = 1))
prd$pred = predict(fit, prd)
lines(sort(x), fitted(fit)[order(x)], col='red', type='b') 
cor(y,fitted(fit)[order(x)])}



library(gridExtra)
pdata<-expand.grid(YearOfBirth=seq(1936,1967,by=1), Sex=seq(0,1,by=1), UniEdu=seq(0,1,by=1), avmse="1")
fit = lm(avMSE~poly(I(YearOfBirth),1)+Sex*UniEdu, seven)
pdata[,4]<-predict(fit,pdata,level=0)

uno = pdata[which(pdata$Sex == 0 & pdata$UniEdu == 0),]
dos = pdata[which(pdata$Sex == 0 & pdata$UniEdu == 1),]
tres = pdata[which(pdata$Sex == 1 & pdata$UniEdu == 0),]
quatro = pdata[which(pdata$Sex == 1 & pdata$UniEdu == 1),]
 
plot(1, type="n", xlab="", ylab="", xlim=c(1936, 1967), ylim=c(-2, 2))
lines(uno$YearOfBirth,uno$avmse, col = "red")
lines(dos$YearOfBirth,dos$avmse, col = "green")
lines(tres$YearOfBirth,tres$avmse, col = "blue")
lines(quatro$YearOfBirth,quatro$avmse, col = "black")

theme_fred <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      plot.margin = unit(c(0.25,0.25,0.5,0.75), "lines"),
      plot.background = element_rect(fill="white",colour="white"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1.3), vjust=-0.25),
      axis.title.y = element_text(colour = "black", size=rel(1.3), angle=90, vjust=1.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.key = element_rect(fill = "white",colour="white"),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size=rel(0.7), face="bold"),
      legend.text  = element_text(size=rel(0.7), face="plain"),
      legend.position=c(0.01,0), legend.justification=c(0.01,0),
      legend.box.just = "left"
    )   
}
theme_set(theme_fred())


ggplot(pdata, aes(YearOfBirth, avmse, colour=Sex, linetype = as.factor(Sex))) + 
  scale_linetype_manual(values=c(2,1)) +
  geom_line  (size=1.1, aes(colour=Sex, linetype = as.factor(Sex))) + 
  scale_shape_manual(values = c(21,21))  +
  scale_x_continuous(limits=c(1935,1968), breaks=seq(1936,1967, 1)) +
  scale_y_continuous(limits=c(-1.5, 1.25), breaks=seq(-1.5,1.25,0.25)) +
  theme(panel.background = element_rect(fill="white",colour="black") )



#from this graph find approximate starting values
'colnames(avmse)[2] = "avMSE_male"
{f = fitModel(avMSE_male ~ A*Year+B^2, data = avmse)
coef(f)
A=coef(f)[1];B=coef(f)[2]
m = nls(avMSE_male ~ A*Year+B^2, data = avmse, start=c(A=A,B=B))}
#get some estimation of goodness of fit
cor(avmse$avMSE_male, predict(m))
#plot the fit
lines(avmse$Year,predict(m),col="red",lty=2,lwd=3)'




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

avmse = matrix(ncol = 5, nrow = 11)
colnames(avmse) = c("Bin", "Mean avMSE Male", "Mean avMSE Female", "SE Male", "SE Female")
for (j in seq_along(1:33)){ 
  for (i in years){
    avmse[j,1] = i
    name = paste('Year:',i,sep='')
    avmse[j,2] = mean(df2[[name]]$avMSE)
    avmse[j,3] = mean(fd2[[name]]$avMSE)
    avmse[j,4] = se(df2[[name]]$avMSE)
    avmse[j,5] = se(fd2[[name]]$avMSE)
    j=j+1
    i=i+1
  }
}
avmse = as.data.frame(avmse)
avmse$YearMale = avmse$Bin-0.15
avmse$YearFemale = avmse$Bin+0.15
{plot(avmse$YearFemale, avmse$`Mean avMSE Female`, col = "darkblue", pch = 16, main = "Decrase of avMSE as a function of Year of Birth binned", 
      ylab = "avMSE in Dioptres", xlab = "Bin", ylim = c(-1,1.5), type = "b")
  points(avmse$YearMale, avmse$`Mean avMSE Male`, col="red", pch = 16, type = "b")
  legend(x=10,y=1.5,fill = c("darkblue","red"), legend = c("Female","Male"),bty = "n")
  segments(avmse$YearFemale, avmse$`Mean avMSE Female` - avmse$`SE Female`,
           avmse$YearFemale, avmse$`Mean avMSE Female` + avmse$`SE Female`)
  epsilon <- 0.1
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female` - avmse$`SE Female`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female` - avmse$`SE Female`)
  segments(avmse$YearFemale-epsilon,avmse$`Mean avMSE Female` + avmse$`SE Female`,avmse$YearFemale+epsilon,avmse$`Mean avMSE Female` + avmse$`SE Female`)
  
  segments(avmse$YearMale, avmse$`Mean avMSE Male` - avmse$`SE Male`,
           avmse$YearMale, avmse$`Mean avMSE Male` + avmse$`SE Male`)
  epsilon <- 0.1
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male` - avmse$`SE Male`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male` - avmse$`SE Male`)
  segments(avmse$YearMale-epsilon,avmse$`Mean avMSE Male` + avmse$`SE Male`,avmse$YearMale+epsilon,avmse$`Mean avMSE Male` + avmse$`SE Male`)
}




# What about highest qualification
qual = 0:6
df2 = list()
fd2 = list()
year = 1936:1967
for (i in qual){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$HighestQualification == i & seven$Sex == 0),]
  df2[[name]] = df
  fd = seven[which(seven$HighestQualification == i & seven$Sex == 1),]
  fd2[[name]] = fd
}

avMSE = matrix(ncol = 5, nrow = 7)
colnames(avMSE) = c("Year of Birth", "Mean avMSE Male", "Mean avMSE Female", "SE Male", "SE Female")
i = 0
j = 1
for (j in seq_along(0:6)){ 
  for (i in qual){
    avMSE[j,1] = i
    name = paste('Year:',i,sep='')
    avMSE[j,2] = mean(df2[[name]]$avMSE)
    avMSE[j,3] = mean(fd2[[name]]$avMSE)
    avMSE[j,4] = se(df2[[name]]$avMSE)
    avMSE[j,5] = se(fd2[[name]]$avMSE)
    j=j+1
    i=i+1
  }
}
avMSE = as.data.frame(avMSE)
avMSE$YearMale = avMSE$Year-0.15
avMSE$YearFemale = avMSE$Year+0.15
{plot(avMSE$YearFemale, avMSE$`Mean avMSE Female`, col = "darkblue", pch = 16, main = "avMSE for individuals with different Education Qualifications", 
      ylab = "predicted avMSE", xlab = "Highest Qualification", ylim = c(-1,1), xlim = c(-0.3,6.3), type = "b")
  points(avMSE$YearMale, avMSE$`Mean avMSE Male`, col="red", pch = 16, type = "b")
  legend(x=5.5,y=-0.7,fill = c("darkblue","red"), legend = c("Female","Male"), bty = "n")
  segments(avMSE$YearFemale, avMSE$`Mean avMSE Female` - avMSE$`SE Female`,
           avMSE$YearFemale, avMSE$`Mean avMSE Female` + avMSE$`SE Female`)
  epsilon <- 0.1
  segments(avMSE$YearFemale-epsilon,avMSE$`Mean avMSE Female` - avMSE$`SE Female`,avMSE$YearFemale+epsilon,avMSE$`Mean avMSE Female` - avMSE$`SE Female`)
  segments(avMSE$YearFemale-epsilon,avMSE$`Mean avMSE Female` + avMSE$`SE Female`,avMSE$YearFemale+epsilon,avMSE$`Mean avMSE Female` + avMSE$`SE Female`)
  
  segments(avMSE$YearMale, avMSE$`Mean avMSE Male` - avMSE$`SE Male`,
           avMSE$YearMale, avMSE$`Mean avMSE Male` + avMSE$`SE Male`)
  epsilon <- 0.1
  segments(avMSE$YearMale-epsilon,avMSE$`Mean avMSE Male` - avMSE$`SE Male`,avMSE$YearMale+epsilon,avMSE$`Mean avMSE Male` - avMSE$`SE Male`)
  segments(avMSE$YearMale-epsilon,avMSE$`Mean avMSE Male` + avMSE$`SE Male`,avMSE$YearMale+epsilon,avMSE$`Mean avMSE Male` + avMSE$`SE Male`)
}







# Based on education status
years = 1937:1969
df2 = list();df3 = list()
fd2 = list();fd3 = list()
for (i in years){
  name = paste('Year:',i,sep='')
  df = seven[which(seven$YearOfBirth == i & seven$Sex == 0),]
  noedu = df[which(df$UniEdu == 0),]; edu = df[which(df$UniEdu == 1),]
  df2[[name]] = noedu; df3[[name]] = edu
  fd = seven[which(seven$YearOfBirth == i & seven$Sex == 1),]
  noedus = fd[which(fd$UniEdu == 0),]; edus = fd[which(fd$UniEdu == 1),]
  fd2[[name]] = noedus; fd3[[name]] = edus 
}



avmse = matrix(ncol = 9, nrow = 33)
colnames(avmse) = c("Year of Birth", "Mean avMSE Male noedu", "Mean avMSE Male edu", "SE Male noedu", "SE Male edu",
                    "Mean avMSE Female noedu", "Mean avMSE Female edu", "SE Female noedu", "SE Female edu")
i=1937
j=1
for (j in seq_along(1:33)){ 
  for (i in years){
    avmse[j,1] = i
    name = paste('Year:',i,sep='')
    avmse[j,2] = mean(df2[[name]]$avMSE)
    avmse[j,3] = mean(df3[[name]]$avMSE)
    avmse[j,4] = se(df2[[name]]$avMSE)
    avmse[j,5] = se(df3[[name]]$avMSE)
    avmse[j,6] = mean(fd2[[name]]$avMSE)
    avmse[j,7] = mean(fd3[[name]]$avMSE)
    avmse[j,8] = se(fd2[[name]]$avMSE)
    avmse[j,9] = se(fd3[[name]]$avMSE)
    j=j+1
    i=i+1
  }
}
avmse = as.data.frame(avmse)
avmse$YearMale = avmse$`Year of Birth`-0.15
avmse$YearFemale = avmse$`Year of Birth`+0.15
par(mfrow=c(1,2))
{plot(avmse$YearFemale, avmse$`Mean avMSE Female noedu`, col = "darkblue", pch = 16, main = "Females", 
      ylab = "Dioptres", xlab = "Year of birth", ylim = c(-2,5), type = "b", cex.main = 1.5, cex.lab = 1.35, cex.axis = 1.35)
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
      ylab = "Dioptres", xlab = "Year of birth", ylim = c(-2,5), type = "b", cex.main = 1.5, cex.lab = 1.35, cex.axis = 1.35)
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


# For UniEdu
s1 = seven[which(seven$Sex == 0 & seven$UniEdu == 0),]
s2 = seven[which(seven$Sex == 1 & seven$UniEdu == 0),]
s3 = seven[which(seven$Sex == 0 & seven$UniEdu == 1),]
s4 = seven[which(seven$Sex == 1 & seven$UniEdu == 1),]

c1 = seven[which(seven$Sex == 0),]
c2 = seven[which(seven$Sex == 1),]

ci1 = seven[which(seven$UniEdu == 0),]
ci2 = seven[which(seven$UniEdu == 1),]

library(Rmisc)
CI(s1$avMSE)
CI(s2$avMSE)
CI(s3$avMSE)
CI(s4$avMSE)


# For EducationAgeOLD
s1 = seven[which(seven$Sex == 0),]
s2 = seven[which(seven$Sex == 1),]

ci(s1$avMSE)
ci(s2$avMSE)

library(vioplot)
par(cex.lab=1.35, cex.axis=1.35) 
vioplot(s1$avMSE,s3$avMSE,s2$avMSE,s4$avMSE, 
        names=c("Male no degree", "Male with degree", "Female no degree", "Female with degree"), col="gold", lty = )
title(ylab = "Dioptres", cex.lab = 1.35, main = "Distribution of avMSE", cex.main = 1.5)

library(ggplot2)
seven$strat[seven$Sex == 0 & seven$UniEdu == 0] = "Male no edu"
seven$strat[seven$Sex == 0 & seven$UniEdu == 1] = "Male edu"
seven$strat[seven$Sex == 1 & seven$UniEdu == 0] = "Female no edu"
seven$strat[seven$Sex == 1 & seven$UniEdu == 1] = "Female edu"


ggplot(seven, aes(x = factor(strat), y = avMSE, fill = strat)) + 
  geom_violin(trim=F)+
  scale_fill_brewer(palette="Blues")+
  labs(x="Group", y = "Dioptres")+
  geom_boxplot(width=0.18, fill="white")+
  theme_bw()+
  theme(axis.title.y = element_text(size=20),
        axis.text= element_text(size=20),
                axis.title.x = element_text(size=20))
    




