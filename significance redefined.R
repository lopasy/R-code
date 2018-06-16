type1=.005
type1Power=0.05
type2=0.25
p=1-c(9000:9990)/10000
xbar = qnorm(1-p/2)
# alternative based on 80% POWER IN 5% TEST
muPower = qnorm(1-type2)+qnorm(1-type1Power/2)
bfPow = 0.5*(dnorm(xbar,muPower,1)+dnorm(xbar,-
                                           muPower,1))/dnorm(xbar,0,1)
muUMPBT = qnorm(0.9975)
bfUMPBT = 0.5*(dnorm(xbar,muUMPBT,1)+dnorm(xbar,-
                                             muUMPBT,1))/dnorm(xbar,0,1)
# two-sided "LR" bound
bfLR = 0.5/exp(-0.5*xbar^2)
bfLocal = -1/(2.71*p*log(p))
#coordinates for dashed lines
data = data.frame(p,bfLocal,bfLR,bfPow,bfUMPBT)
U_005 = max(data$bfLR[data$p=="0.005"])
L_005 = min(data$bfLocal[data$p=="0.005"])
U_05 = max(data$bfLR[data$p=="0.05"])
L_05 = min(data$bfUMPBT[data$p=="0.05"])
# Local bound; no need for two-sided adjustment
#plot margins
par(mai=c(0.8,0.8,.1,0.4))
par(mgp=c(2,1,0))
matplot(p,cbind(bfLR,-1/(2.71*p*log(p))),type='n',log='xy',
        xlab=expression(paste(italic(P) ,"-value")),
        ylab="Bayes Factor",
        ylim = c(0.3,100),
        bty="n",xaxt="n",yaxt="n")
lines(p,bfPow,col="red",lwd=2.5)
lines(p,bfLR,col="black",lwd=2.5)
lines(p,bfUMPBT,col="blue",lwd=2.5)
lines(p,bfLocal,col="green",lwd=2.5)
legend(0.015,100,c(expression(paste("Power")),"Likelihood Ratio
                   Bound","UMPBT",expression(paste("Local-",italic(H)[1],"
                                                   Bound"))),lty=c(1,1,1,1),
       lwd=c(2.5,2.5,2.5,2.5),col=c("red","black","blue","green"),
       cex = 0.8)
#text(0.062,65, "\u03B1", font =3, cex = 0.9)
#customizing axes
#x axis
16
axis(side=1,at=c(-
                   2,0.001,0.0025,0.005,0.010,0.025,0.050,0.100,0.14),
     labels =
       c("","0.0010","0.0025","0.0050","0.0100","0.0250","0.0500","0.1000",
         ""),lwd=1,
     tck = -0.01, padj = -1.1, cex.axis = .8)
#y axis on the left - main
axis(side=2,at=c(-0.2, 0.3,0.5,1,2,5,10,20,50,100),labels =
       c("","0.3","0.5","1.0","2.0","5.0","10.0","20.0","50.0","100.0"),lwd
     =1,las= 1,
     tck = -0.01, hadj = 0.6, cex.axis = .8)
#y axis on the left - secondary (red labels)
axis(side=2,at=c(L_005,U_005),labels = c(13.9,25.7),lwd=1,las= 1,
     tck = -0.01, hadj = 0.6, cex.axis = .6,col.axis="red")
#y axis on the right - main
axis(side=4,at=c(-0.2, 0.3,0.5,1,2,5,10,20,50,100),labels =
       c("","0.3","0.5","1.0","2.0","5.0","10.0","20.0","50.0","100.0"),lwd
     =1,las= 1,
     tck = -0.01, hadj = 0.4, cex.axis = .8)
#y axis on the right - secondary (red labels)
axis(side=4,at=c(L_05,U_05),labels = c(2.4,3.4),lwd=1,las= 1,
     tck = -0.01, hadj = 0.4, cex.axis = .6,col.axis="red")
###dashed lines
segments(x0 = 0.000011, y0= U_005, x1 = 0.005, y1 = U_005, col =
           "gray40", lty = 2)
segments(x0 = 0.000011, y0= L_005, x1 = 0.005, y1 = L_005, col =
           "gray40", lty = 2)
segments(x0 = 0.005, y0= 0.00000001, x1 = 0.005, y1 = U_005, col =
           "gray40", lty = 2)
segments(x0 = 0.05, y0= U_05, x1 = 0.14, y1 = U_05, col = "gray40",
         lty = 2)
segments(x0 = 0.05, y0= L_05, x1 = 0.14, y1 = L_05, col = "gray40",
         lty = 2)
segments(x0 = 0.05, y0= 0.00000001, x1 = 0.05, y1 = U_05, col =
           "gray40", lty = 2)

pow1=c(5:999)/1000 # power range for 0.005 tests
pow2=c(50:999)/1000 # power range for 0.05 tests
alpha=0.005 # test size
pi0=5/6 # prior probability
N=10^6 # doesn't matter
#graph margins
par(mai=c(0.8,0.8,0.1,0.1))
par(mgp=c(2,1,0))
plot(pow1,alpha*N*pi0/(alpha*N*pi0+pow1*(1-pi0)*N),type='n',ylim = c(0,1), xlim = c(0,1.5), xlab="Power", ylab="False positive rate", bty="n", xaxt="n", yaxt="n")
#grid lines
segments(x0 = -0.058, y0 = 0, x1 = 1, y1 = 0,lty=1,col = "gray92")
segments(x0 = -0.058, y0 = 0.2, x1 = 1, y1 = 0.2,lty=1,col =
           "gray92")
segments(x0 = -0.058, y0 = 0.4, x1 = 1, y1 = 0.4,lty=1,col =
           "gray92")
segments(x0 = -0.058, y0 = 0.6, x1 = 1, y1 = 0.6,lty=1,col =
           "gray92")
segments(x0 = -0.058, y0 = 0.8, x1 = 1, y1 = 0.8,lty=1,col =
           "gray92")
segments(x0 = -0.058, y0 = 1, x1 = 1, y1 = 1,lty=1,col = "gray92")
lines(pow1,alpha*N*pi0/(alpha*N*pi0+pow1*(1-
                                            pi0)*N),lty=1,col="blue",lwd=2)
odd_1_5_1 = alpha*N*pi0/(alpha*N*pi0+pow1[995]*(1-pi0)*N)
alpha=0.05
pi0=5/6
lines(pow2,alpha*N*pi0/(alpha*N*pi0+pow2*(1-
                                            pi0)*N),lty=2,col="blue",lwd=2)
odd_1_5_2 = alpha*N*pi0/(alpha*N*pi0+pow2[950]*(1-pi0)*N)
alpha=0.05
pi0=10/11
lines(pow2,alpha*N*pi0/(alpha*N*pi0+pow2*(1-
                                            pi0)*N),lty=2,col="red",lwd=2)
odd_1_10_2 = alpha*N*pi0/(alpha*N*pi0+pow2[950]*(1-pi0)*N)
alpha=0.005
pi0=10/11
lines(pow1,alpha*N*pi0/(alpha*N*pi0+pow1*(1-
                                            pi0)*N),lty=1,col="red",lwd=2)
odd_1_10_1 = alpha*N*pi0/(alpha*N*pi0+pow1[995]*(1-pi0)*N)
alpha=0.05
pi0=40/41

lines(pow2,alpha*N*pi0/(alpha*N*pi0+pow2*(1-
                                            pi0)*N),lty=2,col="green",lwd=2)
odd_1_40_2 = alpha*N*pi0/(alpha*N*pi0+pow2[950]*(1-pi0)*N)
alpha=0.005
pi0=40/41
lines(pow1,alpha*N*pi0/(alpha*N*pi0+pow1*(1-
                                            pi0)*N),lty=1,col="green",lwd=2)
odd_1_40_1 = alpha*N*pi0/(alpha*N*pi0+pow1[995]*(1-pi0)*N)
#customizing axes
axis(side=2,at=c(-0.5,0,0.2,0.4,0.6,0.8,1.0),labels =
       c("","0.0","0.2","0.4","0.6","0.8","1.0"),
     lwd=1,las= 1,tck = -0.01, hadj = 0.4, cex.axis = .8)
axis(side=1,at=c(-0.5,0,0.2,0.4,0.6,0.8,1.0),labels =
       c("","0.0","0.2","0.4","0.6","0.8","1.0"),
     lwd=1,las= 1, tck = -0.01, padj = -1.1, cex.axis = .8)
legend(1.05,1,c("Prior odds = 1:40","Prior odds = 1:10","Prior odds
= 1:5"),pch=c(15,15,15),
       col=c("green","red","blue"), cex = 1)
############### Use these commands to add brackets in Figure 2
library(pBrackets)
#add text and brackets
text(1.11,(odd_1_5_2+odd_1_40_2)/2, expression(paste(italic(P)," <
0.05 threshold")), cex = 0.9,adj=0)
text(1.11,(odd_1_5_1+odd_1_40_1)/2, expression(paste(italic(P)," <
0.005 threshold")), cex = 0.9,adj=0)
brackets(1.03, odd_1_40_1, 1.03, odd_1_5_1, h = NULL, ticks = 0.5,
         curvature = 0.7, type = 1,
         col = 1, lwd = 1, lty = 1, xpd = FALSE)
brackets(1.03, odd_1_40_2, 1.03, odd_1_5_2, h = NULL, ticks = 0.5,
         curvature = 0.7, type = 1,
         col = 1, lwd = 1, lty = 1, xpd = FALSE)
