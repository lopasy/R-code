#some functions for converting between probabilities and odds
prob2odds <-function(p) p/(1-p)
odds2prob <-function(o) o/(1+o)

###FUNCTION TO GENERATE DATA
#arguments are;
#
#n        number of data points to generate
#snp.prob probability of the minor allele
#prob.AA  probability of disease in those with AA genotype
#beta.Aa  log odds ratio between Aa's and AA's
#beta.aa  log odds ratio between aa's and AA's
#
#(note that 2*beta.Aa = beta.aa, under the additive 'linear' model)

make.data<-function(n, snp.prob, prob.AA, beta.Aa, beta.aa){
  
  g <- rbinom(n, 2, snp.prob) # i.e. in Hardy-Weinburg
  
  odds.AA <- prob2odds(prob.AA) 
  prob.Aa <- odds2prob(odds.AA * exp(beta.Aa))
  prob.aa <- odds2prob(odds.AA * exp(beta.aa))
  
  yprobs <- (g==0)*prob.AA + (g==1)*prob.Aa + (g==2) * prob.aa
  
  y <- rbinom(n, 1, prob=yprobs)
  
  data.frame(snp=g, y=y)
}

###FUNCTIONS TO IMPLEMENT THE DIFFERENT REGRESSION MODELS
#additive, dominant and recessive models
#some special cases are taken care of in the first lines
#to prevent error messages later on

additive.p<-function(snp, y){
  model<-glm(y~as.numeric(snp), family=binomial, data = data)
  coef(summary(model))[2,4] # i.e. the p-value for association
}

dominant.p<-function(snp, y){
  if(all(snp>=1)) return(1.0) ## no common-allele homozygotes
  model<-glm(y~I(snp >= 1), family=binomial, data = data)
  coef(summary(model))[2,4]
}

recessive.p<-function(snp, y){
  if (all(snp<2)) return(1.0) ## no rare-allele homozygotes
  model<-glm(y~I(snp==2), family=binomial, data = data)
  coef(summary(model))[2,4]
}

one.p<-function(n, snp.prob, prob.AA, beta.Aa, beta.aa){
  data<-make.data(n, snp.prob, prob.AA, beta.Aa, beta.aa)
  c(add=additive.p(data$snp, data$y),
    dom=dominant.p(data$snp, data$y),
    rec=recessive.p(data$snp, data$y))
}

###FUNCTION TO DO MANY SIMULATED ANALYSES
#
#'LOTS' tells R how many simulations to do, the other parameters are described
#above

power<-function(LOTS, n, snp.prob, prob.AA, beta.Aa, beta.aa){
  
  pvalues<-replicate(LOTS, one.p(n, snp.prob, prob.AA, beta.Aa, beta.aa))
  
  apply(pvalues < 0.05, 1, mean)
  
}

###USING THE FUNCTIONS

# varying beta.aa, fixing the other parameters
# we write a 'wrapper' function for sapply;

example.power<-function(betas.aa, LOTS=500){
  
  sapply(betas.aa, function(beta.aa) power(LOTS, n=200, snp.prob=0.2, prob.AA=0.3, beta.Aa=0.5, beta.aa))
  
}

###A FAIRLY QUICK EXAMPLE (1-2 minutes)

set.seed(4)
pow1<-example.power(c(0,0.5,0.75, 1,1.25,1.5), LOTS=500)  
alarm() #beep!

beta.ratio <- c(0,0.5,0.75, 1,1.25,1.5)/0.5
plot(beta.ratio, pow1[1,], xlab=expression(beta[aa]/beta[Aa]),
     ylab="Power",type="l",col="forestgreen",lwd=2,
     ylim=c(0,max(pow1[1:3,]) ), main="n=200, p=0.2, Pr(Y|AA)=0.3, beta.Aa=0.5, LOTS=500" )

lines(beta.ratio, pow1[2,], col="purple", lwd=2)
lines(beta.ratio, pow1[3,], col="goldenrod",lwd=2)

legend("topleft", legend=c("additive","dominant","recessive"),
       col=c("forestgreen","purple","goldenrod"),lwd=2,lty=1, bty="n")

#it looks like the additive model is never the worst option
#- and it only loses out to the recessive model only when 
# there is a really big violation of the beta.aa = 2*beta.Aa

#let's repeat the analysis, but with more accuracy and a (few)
#more points on the x-axis

set.seed(5)
pow2<-example.power(seq(0,1.5, l=21), LOTS=2000)  
alarm()

beta.ratio <- seq(0,1.5, l=21)/0.5
plot(beta.ratio, pow2[1,], xlab=expression(beta[aa]/beta[Aa]),
     ylab="Power",type="l",col="forestgreen",lwd=2,
     ylim=c(0,max(pow2[1:3,]) ), main="n=200, p=0.2, Pr(Y|AA)=0.3, beta.Aa=0.5, LOTS=2000" )
lines(beta.ratio, pow2[2,], col="purple", lwd=2)
lines(beta.ratio, pow2[3,], col="goldenrod",lwd=2)
legend("topleft", legend=c("additive","dominant","recessive"),
       col=c("forestgreen","purple","goldenrod"),lwd=2,lty=1)

pdf("SpecExlarge.pdf") #see file
dev.off()