library(ggplot2); library(gridExtra); library(reshape2); library(grid)
t = as.data.frame(t)
colnames(t)[1] <- "Effect"; colnames(t)[2] <- "Sample size"; colnames(t)[3] <- "Power"
colnames(t)[4] <- "Effect"; colnames(t)[5] <- "Sample size"; colnames(t)[6] <- "Power"
colnames(t)[7] <- "Effect"; colnames(t)[8] <- "Sample size"; colnames(t)[9] <- "Power"
colnames(t)[10] <- "Effect"; colnames(t)[11] <- "Sample size"; colnames(t)[12] <- "Power"
colnames(t)[13] <- "Effect"; colnames(t)[14] <- "Sample size"; colnames(t)[15] <- "Power"
colnames(t)[16] <- "Effect"; colnames(t)[17] <- "Sample size"; colnames(t)[18] <- "Power"


Effect = c(0.1,0.12,0.14,0.16,0.18,0.2)
ggplot(data=t, aes_string(x = "`Sample size`", y = "Power")) +
  geom_line(aes(color = factor(Effect)), stat = "identity", size = 1.5)+
  geom_line(aes(color = factor(Effect)), stat = "identity", size = 1.5)+
  geom_line(aes(y=t[,9]),size=1.0)+
  geom_line(aes(y=t[,12]),size=1.0)+
  geom_line(aes(y=t[,15]),size=1.0)+
  geom_line(aes(y=t[,18]),size=1.0)+
  scale_x_continuous(breaks=seq(10000,100000,15000))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  geom_hline(yintercept = 0.80, linetype = 2)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold", margin = margin(t = 0, r = 20, b = 20, l = 0)), 
        legend.title = element_text(), 
        axis.text= element_text(size=12), axis.title=element_text(size=16,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))+
  labs(colour = "Effect size")+ 
  ggtitle("Power to detect different GxE effects, MAF 0.05")+
  guides(colour = guide_legend(reverse=T))
  
  
  
  
par(mfrow=c(2,2))  
options(scipen=5)  
ggplot(data=t, aes_string(x = "`Sample size`", y = "Power"))+
  geom_line(aes(color = factor(Effect)), stat = "identity", size = 1.5)+
  geom_text(x=10300, y=0.95, label="C",size=22)+
  scale_x_continuous(breaks=seq(10000,100000,10000))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  geom_hline(yintercept = 0.80, linetype = 2)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold", margin = margin(t = 0, r = 20, b = 20, l = 0)), 
        legend.title = element_text(size=20), legend.text=element_text(size=22),legend.key.size = unit(1.5, "cm"),
        axis.text= element_text(size=22), axis.title=element_text(size=26,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0)))+
  labs(colour = "Effect size")+ 
  guides(colour = guide_legend(reverse=T,keywidth=2))
  
  
  
options(scipen=5)  
{plot(t[,14], t[,15],type = "l",col = "deepskyblue4",lty = 1, lwd = 3, xlab = "Sample size", 
                  ylab = "Power", ylim = 0:1, xlim = c(10000,100000))
lines(t[,11], t[,12],type = "l",col = "seagreen",lty = 1, lwd = 3)
lines(t[,8], t[,9],type = "l",col = "mediumvioletred",lty = 1, lwd = 3)
lines(t[,5], t[,6],type = "l",col = "darkorange2",lty = 1, lwd = 3)
lines(t[,2], t[,3],type = "l",col = "darkred",lty = 1, lwd = 3)
abline(h = 0.8, lty = 2, lwd = 2.5)
title(main = "MAF 0.05")
legend(x = 80000, y = 0.2, legend = c("0.1 GxE","0.08 GxE","0.06 GxE","0.04 GxE","0.02 GxE"), 
       fill = c("deepskyblue4","seagreen","mediumvioletred","darkorange2","darkred"), pch = 1)}
legend(x = 80000, y = 0.2, legend = "0.08 GxE", col = "seagreen")
legend(x = 80000, y = 0.2, legend = "0.06 GxE", col = "mediumvioletred")
legend(x = 80000, y = 0.2, legend = "0.04 GxE", col = "darkorange2")
legend(x = 80000, y = 0.2, legend = "0.02 GxE", col = "darkred")



{data = table(edu$UniEdu, edu$edi)
names = c("No Degree", "Degree")
data = data.frame(cbind(data),names)
colnames(data)=c(13,14,15,16,17,18,19,20,21,22,23,24,25,26,"names")
data.m = melt(data, id.vars='names')
data = data[,-15]
data = t(data)
colnames(data)=c("No Degree", "Degree")

g = tableGrob(data)
lapply(1:length(g$grobs),function(i){g$grobs[[i]]$gp$fontsize <<- 15})

ggplot(data=data.m, aes(x = variable, value, fill = names)) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  scale_y_continuous(breaks=seq(0,25000,5000)) +
  scale_fill_discrete(name = "Education Status") +
  geom_text(x=1.5, y=15000, label="A",size=20) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold", margin = margin(t = 0, r = 20, b = 20, l = 0)), 
        legend.title = element_text(size=16), legend.text=element_text(size=16),legend.key.size = unit(0.8, "cm"),
        axis.text= element_text(size=20), axis.title=element_text(size=20,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0))) +
  labs(x = "Age at  education completion") + 
  labs(y = "Number of individuals") +
  annotation_custom(g, xmin=2, xmax=22.5, ymin=6000, ymax=15000)}




{data = table(avMSE2$UniEdu, avMSE2$EducationAgeOLD)
names = c("No Degree", "Degree")
data = data.frame(cbind(data),names)
colnames(data)=c(13,14,15,16,17,18,19,20,21,22,23,24,25,26,"names")
data.m = melt(data, id.vars='names')
data = data[,-15]
data = t(data)
colnames(data)=c("No Degree", "Degree")

g = tableGrob(data)
lapply(1:length(g$grobs),function(i){g$grobs[[i]]$gp$fontsize <<- 15})

ggplot(data=data.m, aes(x = variable, value, fill = names)) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  scale_y_continuous(breaks=seq(0,25000,5000)) +
  scale_fill_discrete(name = "Education Status") +
  geom_text(x=1.5, y=24500, label="B",size=20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold", margin = margin(t = 0, r = 20, b = 20, l = 0)), 
        legend.title = element_text(size=16), legend.text=element_text(size=16),legend.key.size = unit(0.8, "cm"),
        axis.text= element_text(size=20), axis.title=element_text(size=20,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 0, l = 0))) +
  labs(x = "Age at  education completion") + 
  labs(y = "Number of individuals") +
  annotation_custom(g, xmin=2, xmax=22.5, ymin=14000, ymax=20000)}
  


