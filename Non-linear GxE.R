#  Created by Yuehua Cui and Shujie Ma on 8/20/2011.
#  Copyright 2011 Michigan State University. All rights reserved.

#  This code R is coded for the paper in Ma et al. (2011) Bininformatics 27(15):2119-2126. 
#  Varying coefficient model for gene-environment interaction: a non-linear look. 


library(splines); library(BEDMatrix)
B=100 ### number of wild bootstrap samples. The larger the better, but depends on your computational power
path = "C:/Users/Alfred/Documents/CQR_17_04_2018/GWAS/directly_genotyped"
m = BEDMatrix(path)
num_snps = 50001:60000; m = m[,num_snps]

### read data file ###
master_true = read.csv("D:/master_true.txt", sep="")
edu = read.csv("~/Documents/ukb_epi_mse2017-02-16.csv"); edu = edu[,c(1,14)]
final = merge(master_true, edu, 1)
final = final[,c(1,20,50)]

final = cbind(final, m)
final = final[, -1]

k=ncol(final)-2 ### number of SNPs; the first two cols contain the y and x variables
SNPdata=final[,3:(k+2)]
data=final

### the following code is to remove SNPs with MAF<0.05 (you can also choose other threshold)
#index=matrix(0,k,1)
#for(i in 1:k)
#{
#	tab=table(SNPdata[,i])
#	if(length(table(tab))==2)
#	{
#	    PA=min(tab[1]/(tab[1]+tab[2]),tab[2]/(tab[1]+tab[2]))
#	}
#	if(length(table(tab))==3)
#	{
#	    PA=(tab[1]*2+tab[2])/((tab[1]+tab[2]+tab[3])*2)
#	}
#	index[i]=(PA<0.05)*i  
#}
#index=na.omit(index)
#index=index[index>0]
#SNPdata=SNPdata[,-index] ### delete SNPs with MAF<0.05 (you can also choose other cutoff)

### Power transformation function for response y ###
xlambda=function(x,lamb)  #box-cox transformation
{
  if (lamb != 0) {(x^(lamb)-1)/lamb}
  else {log(x)}
}

llambda=function(x,lamb)  #l-lambda
{
  n=length(x); m=length(lamb); ll=rep(0,m);
  for (i in 1:m)
  { 
    ll[i]=(lamb[i]-1)*n*mean(log(x)) - n/2*log( (n-1)/n*var(xlambda(x,lamb[i])) ) 
  }
  ll 
}
#########################

{pv0=matrix(1,k,1) ### matrix to save the p-values for the overall genetic test when fitting a VC model
pvc=matrix(1,k,1) ### matrix to save the p-values for the constant coefficient test
pvi=matrix(1,k,1) ### matrix to save the p-values for the linear coefficient test

pl0=matrix(1,k,1) ### matrix to save the p-values for the overall genetic test when fitting a regular linear model
plc=matrix(1,k,1)
pli=matrix(1,k,1)}

for(m in 1:k){
  d=cbind(data[,1],data[,2],SNPdata[,m])
  d=na.omit(d) ### eliminate obs. with missing values
  d=subset(d,subset=((d[,1]>=0.6)&(d[,2]<50))) ### eliminate extreme observations
  yy=d[,1]
  ### The following code is for the boxcox transformation to make the y variable more normal
  ### You don't have to do this if you beleive the data are pretty normal
  #lambda=seq(1.4,2.6,0.02)
  #lam=llambda(yy,lambda)
  #lambd=lambda[which(lam==max(lam))]
  #y=xlambda(yy,lambd)
  y=yy
  ########################
  
  xx=d[,2]; xx=as.matrix(xx)
  t1=d[,3]
  t11=t1-1 ### genetic variable for the additive effect
  t12=as.numeric(t1==1)-0.5  ### genetic variable for the dominance effect
  n=nrow(d)
  
  ### transform x to make it uniform ###
  x=as.matrix(pnorm((xx-mean(xx))/sd(xx)))
  ###########################
  ### reduced model without genetic effect
  xl0=cbind(matrix(1,c(n,1)),x) 
  betax0=solve(t(xl0)%*%xl0)
  bt0=betax0%*%t(xl0)
  betax0=bt0%*%y
  mhat0=xl0%*%betax0
  
  ### linear model without interaction
  xl1=cbind(matrix(1,c(n,1)),x,t11,t12)
  beta1=solve(t(xl1)%*%xl1)
  bt1=beta1%*%t(xl1)
  beta1=bt1%*%y
  mhat1=xl1%*%beta1
  
  #### linear model with interaction
  xli=cbind(matrix(1,c(n,1)),x,t11,t12,x*t11,x*t12)
  betaxi=solve(t(xli)%*%xli)
  bti=betaxi%*%t(xli)
  betaxi=bti%*%y
  mhati=xli%*%betaxi
  
  ###########################
  ### fit the varying coefficient model 
  ### BIC for knots and order selection ###
  BIC=matrix(0,4,4)
  for (o in 2:5){
    for (N in 1:4){
      xknots=bs(x,df=N+o+1,degree=o,intercept = TRUE)
      xknots=as.matrix(xknots)
      tt1=t11%*%array(1,c(1,ncol(xknots)))
      tt2=t12%*%array(1,c(1,ncol(xknots)))
      z1=xknots*tt1
      z2=xknots*tt2
      xt=cbind(matrix(1,c(n,1)),x,z1,z2)
      bt=solve(t(xt)%*%xt)%*%t(xt)
      beta=bt%*%y
      mhat=xt%*%beta
      s2hat=sum((y-mhat)^2)/n
      BIC[N,(o-1)]=log(s2hat)+ncol(xt)*log(n)/n
      #N=N+1
    }
    #o=o+1
  }
  N=which(BIC == min(BIC), arr.ind = TRUE)[1]
  o=which(BIC == min(BIC), arr.ind = TRUE)[2]+1
  
  #o=4
  #N=floor((n^(1/(2*o+1))))+2
  xknots=bs(x,df=(N+o+1),degree=o,intercept = TRUE)
  #xknots=ns(x,df=(N+2),intercept = TRUE)      
  xknots=as.matrix(xknots)
  tt1=t11%*%array(1,c(1,ncol(xknots)))
  tt2=t12%*%array(1,c(1,ncol(xknots)))
  z1=xknots*tt1
  z2=xknots*tt2
  xt=cbind(matrix(1,c(n,1)),x,z1,z2)
  beta=solve(t(xt)%*%xt)
  bt=beta%*%t(xt)
  beta=bt%*%y
  mhat=xt%*%beta
  beta12=beta[3:(ncol(xknots)+2)]
  beta22=beta[(ncol(xknots)+3):nrow(beta)]
  
  #Calculation of beta1(x) and beta2(x)
  beta12=xknots%*%beta12 ### additive effect over x
  beta22=xknots%*%beta22 ### dominance effect over x
  
  #plot(xx, beta12)
  #plot(xx, mhat)
  #########################
  #Test statistic
  Tnv=t(mhat-mhat0)%*%(mhat-mhat0)/n
  Tnc=t(mhat-mhat1)%*%(mhat-mhat1)/n
  Tni=t(mhat-mhati)%*%(mhat-mhati)/n
  
  #########################
  ### bootstrap pvalue; 
  ### If you believe the error follows a normal distribution, 
  ### then a likelihood ratio test can be applied to save computation time 
  
  NewB=matrix(1,B,3)
  for (j in 1:B)
  {     set.seed(j) 
    delta=runif(n)
    delta=matrix(as.numeric(delta<(1+sqrt(5))/(2*sqrt(5))),c(n,1))
    delta=sqrt(5)*(0.5-delta)+0.5
    
    ### For testing beta(x)=0
    ybv=mhat0+(y-mhat)*delta
    betab=bt%*%ybv
    mhatb=xt%*%betab
    
    betax0b=bt0%*%ybv
    mhat0b=xl0%*%betax0b
    
    Tnvb=t(mhatb-mhat0b)%*%(mhatb-mhat0b)/n
    pv=as.numeric(Tnvb>Tnv)
    
    ### For testing beta(x)=beta, constant coeffcient
    ybc=mhat1+(y-mhat)*delta
    betabc=bt%*%ybc
    mhatbc=xt%*%betabc
    
    betax1b=bt1%*%ybc
    mhat1b=xl1%*%betax1b
    
    Tnvc=t(mhatbc-mhat1b)%*%(mhatbc-mhat1b)/n
    pc=as.numeric(Tnvc>Tnc)
    
    ### For testing beta(x)=beta1+beta2*X, linear function
    ybi=mhati+(y-mhat)*delta
    betabi=bt%*%ybi
    mhatbi=xt%*%betabi
    
    betax1bi=bti%*%ybi
    mhat1bi=xli%*%betax1bi
    
    Tnvi=t(mhatbi-mhat1bi)%*%(mhatbi-mhat1bi)/n
    pi=as.numeric(Tnvi>Tni)
    
    NewB[j,]=c(pv, pc, pi)
  }
  
  ################
  pv0[m]=mean(NewB[,1])
  pvc[m]=mean(NewB[,2])
  pvi[m]=mean(NewB[,3])
  
  ################
  #pvalues for testing beta=0 and gamma=0 in linear regression
  l0=logLik(lm(y~x))  ### fit a linear model with only covariate X
  lm1=logLik(lm(y~x+t1))  ### fit a linear model with G and X
  l1=logLik(lm(y~x+t1+t1*x))  ### fit a linear model with G+X+G*X
  LR1=2*(lm1[1]-l0[1])  ### the likelihood ratio statistic for testing genetic effect
  LR2=2*(l1[1]-l0[1])  ### the likelihood ratio statistic for testing genetic effect including the interaction effect
  LRi=2*(l1[1]-lm1[1])  ### the likelihood ratio statistic for testing interaction effect effect
  
  
  pl0[m]=1-pchisq(LR1,1) ### the corresponding p-value for the above linear models
  plc[m]=1-pchisq(LR2,2)
  pli[m]=1-pchisq(LRi,1)
  
  cat(m,"\n")
}

p=data.frame(cbind(pv0,pvc,pvi,pl0,plc,pli))
colnames(p)<-c("pvc_0","pvc_constant","pvc_linear","pl_0","pl_noInter","pl_interaction")
write.table(p,"p-values.txt",sep="\t")

### The following code is for the calculation of the heritability H^2 ###
pA=table(t1)[3]/(table(t1)[1]+table(t1)[2]+table(t1)[3])
varT1=2*pA*(1-pA) ### variance for the additive variable
varT2=(1/4)*(1-(2*pA-1)^4) ### varaince for the dominance variable
covT12=2*pA*(pA-1)*(2*pA-1) ### covaraince between additive and dominance variables
res2=(y-mhat)*(y-mhat)
sigbeta=solve(t(xknots)%*%xknots)
sigbt=sigbeta%*%t(xknots)
sigbeta=sigbt%*%res2
sighat2=xknots%*%sigbeta
betavT=beta12^2*varT1+beta22^2*varT2+2*beta12*beta22*covT12
H2=betavT*((betavT+sighat2)^(-1))
plot(xx,H2)

yhat=(mhat*lambd+1)^(1/lambd) ### transform yhat to original scale
b=cbind(t1,xx,beta12,beta22,H2,yy,x)
write.table(b,"non_lin_GxE.txt",sep="\t")



