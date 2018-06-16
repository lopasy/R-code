# Mixture distribution toy examples
N <- 100000

components <- sample(1:3,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
mus <- c(0,10,3)
sds <- sqrt(c(1,1,0.1))
samples <- rnorm(n=N,mean=mus[components],sd=sds[components])



x = seq(-20,20,.1)
truth = .3*dnorm(x,0,1) + .5*dnorm(x,10,1) + .2*dnorm(x,3,.1)
plot(density(samples),main="Density Estimate of the Mixture Model",ylim=c(0,.2),lwd=2)
lines(x,truth,col="red",lwd=2)
legend("topleft",c("True Density","Estimated Density"),col=c("red","black"),lwd=2)



# R code for examples in Lecture 20

# Data preparation
snoqualmie <- read.csv("http://www.stat.cmu.edu/~cshalizi/402/lectures/16-glm-practicals/snoqualmie.csv",header=FALSE)
snoqualmie.vector <- na.omit(unlist(snoqualmie))
snoq <- snoqualmie.vector[snoqualmie.vector > 0]


### Figure 1
plot(hist(snoq,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="Precipitation (1/100 inch)",main="Precipitation in Snoqualmie Falls")
lines(density(snoq),lty=2)


# Two-component Gaussian mixture
library(mixtools)
snoq.k2 <- normalmixEM(snoq,k=2,maxit=100,epsilon=0.01)

# Function to add Gaussian mixture components, vertically scaled, to the
# current plot
# Presumes the mixture object has the structure used by mixtools
plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

### Figure 2
plot(hist(snoq,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="Precipitation (1/100 inch)",main="Precipitation in Snoqualmie Falls")
lines(density(snoq),lty=2)
sapply(1:2,plot.normal.components,mixture=snoq.k2)


# Function to calculate the cumulative distribution function of a Gaussian
# mixture model
# Presumes the mixture object has the structure used by mixtools
# Doesn't implement some of the usual options for CDF functions in R, like
# returning the log probability, or the upper tail probability
pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}


#### Figure 3
# Distinct values in the data
distinct.snoq <- sort(unique(snoq))
# Theoretical CDF evaluated at each distinct value
tcdfs <- pnormmix(distinct.snoq,mixture=snoq.k2)
# Empirical CDF evaluated at each distinct value
# ecdf(snoq) returns an object which is a _function_, suitable for application
# to new vectors
ecdfs <- ecdf(snoq)(distinct.snoq)
# Plot them against each other
plot(tcdfs,ecdfs,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
# Main diagonal for visual reference
abline(0,1)


# Probability density function for a Gaussian mixture
# Presumes the mixture object has the structure used by mixtools
dnormalmix <- function(x,mixture,log=FALSE) {
  lambda <- mixture$lambda
  k <- length(lambda)
  # Calculate share of likelihood for all data for one component
  like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  # Create array with likelihood shares from all components over all data
  likes <- sapply(1:k,like.component,x=x)
  # Add up contributions from components
  d <- rowSums(likes)
  if (log) {
    d <- log(d)
  }
  return(d)
}

# Log likelihood function for a Gaussian mixture, potentially on new data
loglike.normalmix <- function(x,mixture) {
  loglike <- dnormalmix(x,mixture,log=TRUE)
  return(sum(loglike))
}

# Evaluate various numbers of Gaussian components by data-set splitting
# (i.e., very crude cross-validation)
n <- length(snoq)
data.points <- 1:n
data.points <- sample(data.points) # Permute randomly
train <- data.points[1:floor(n/2)] # First random half is training
test <- data.points[-(1:floor(n/2))] # 2nd random half is testing
candidate.component.numbers <- 2:10
loglikes <- vector(length=1+length(candidate.component.numbers))
# k=1 needs special handling
mu<-mean(snoq[train]) # MLE of mean
sigma <- sd(snoq[train])*sqrt((n-1)/n) # MLE of standard deviation
loglikes[1] <- sum(dnorm(snoq[test],mu,sigma,log=TRUE))
for (k in candidate.component.numbers) {
  mixture <- normalmixEM(snoq[train],k=k,maxit=400,epsilon=1e-2)
  loglikes[k] <- loglike.normalmix(snoq[test],mixture=mixture)
}

### Figure 4
plot(x=1:10, y=loglikes,xlab="Number of mixture components",
     ylab="Log-likelihood on testing data")

### Figure 5
snoq.k9 <- normalmixEM(snoq,k=9,maxit=400,epsilon=1e-2)
plot(hist(snoq,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="Precipitation (1/100 inch)",main="Precipitation in Snoqualmie Falls")
lines(density(snoq),lty=2)
sapply(1:9,plot.normal.components,mixture=snoq.k9)

### Figure 6
# Assigments for distinct.snoq and ecdfs are redundant if you've already
# made Figure 3
distinct.snoq <- sort(unique(snoq))
tcdfs <- pnormmix(distinct.snoq,mixture=snoq.k9)
ecdfs <- ecdf(snoq)(distinct.snoq)
plot(tcdfs,ecdfs,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

### Figure 7
plot(0,xlim=range(snoq.k9$mu),ylim=range(snoq.k9$sigma),type="n",
     xlab="Component mean", ylab="Component standard deviation")
points(x=snoq.k9$mu,y=snoq.k9$sigma,pch=as.character(1:9),
       cex=sqrt(0.5+5*snoq.k9$lambda))

### Figure 8
plot(density(snoq),lty=2,ylim=c(0,0.04),
     main=paste("Comparison of density estimates\n",
                "Kernel vs. Gaussian mixture"),
     xlab="Precipitation (1/100 inch)")
curve(dnormalmix(x,snoq.k9),add=TRUE)


# Do the classes of the Gaussian mixture make sense as annual weather patterns?
# Most probable class for each day:
day.classes <- apply(snoq.k9$posterior,1,which.max)
# Make a copy of the original, structured data set to edit
snoqualmie.classes <- snoqualmie
# Figure out which days had precipitation
wet.days <- (snoqualmie > 0) & !(is.na(snoqualmie))
# Replace actual precipitation amounts with classes
snoqualmie.classes[wet.days] <- day.classes
# Problem: the number of the classes doesn't correspond to e.g. amount of
# precipitation expected.  Solution: label by expected precipitation, not by
# class number.
snoqualmie.classes[wet.days] <- snoq.k9$mu[day.classes]

### Figure 9
plot(0,xlim=c(1,366),ylim=range(snoq.k9$mu),type="n",xaxt="n",
     xlab="Day of year",ylab="Expected precipiation (1/100 inch)")
axis(1,at=1+(0:11)*30)
for (year in 1:nrow(snoqualmie.classes)) {
  points(1:366,snoqualmie.classes[year,],pch=16,cex=0.2)
}

# Next line is currently (5 April 2011) used to invoke a bug-patch kindly
# provided by Dr. Derek Young; the patch will be incorporated in the next
# update to mixtools, so should not be needed after April 2011
source("http://www.stat.cmu.edu/~cshalizi/402/lectures/20-mixture-examples/bootcomp.R")
snoq.boot <- boot.comp(snoq,max.comp=10,mix.type="normalmix",
                       maxit=400,epsilon=1e-2)
# Running this takes about 5 minutes

### Figure 10
# automatically produced as a side-effect of running boot.comp()

### Figure 11
library(mvtnorm)
x.points <- seq(-3,3,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- c(1,1)
sigma <- matrix(c(2,1,1,1),nrow=2)
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),mean=mu,sigma=sigma)
  }
}
contour(x.points,y.points,z)
# Using expand.grid, as in Lecture 6, would be more elegant than this double
# for loop


# Intro to mixtures of normal distributions


# Why mixtures of normal distributions?

# -- typically for linguists because the subjects and / or items come from 2 or more different groups and contribute different random effects depending on the group they come from
# -- these groups could have different sizes in the population, e.g., 20% of the population strongly prefers to give indefinites narrow scope, while the other 80% of the population has a mild preference for wide scope; their contribution to the experimental results should therefore be weighted by the relative sizes of the 2 groups



# Example 1: a mixture of 2 normals

# -- we choose one or the other normal; their weights are 0.3 and 0.7, i.e., when we generate the mixture distribution, we choose the first normal 30% of the time and the second normal 70% of the time

w <- c(0.3, 0.7)
t(rmultinom(1, size=1, prob=w))

# -- the first normal distribution: mu[1]=-2, sigma[1]=1
# -- the second normal distribution: mu[2]=4, sigma[2]=1.5

mu <- c(-2, 4)
sigma <- c(1, 1.5)

(sample.draw <- t(rmultinom(1, size=1, prob=w)))
as.numeric(sample.draw %*% mu)
as.numeric(sample.draw %*% sigma)

# -- we generate 10000 observations from the mixture distribution:

n.obs <- 10000
y <- vector(length=n.obs)
str(y)

for (i in 1:n.obs) {
  sample.draw <- t(rmultinom(1, size=1, prob=w))
  y[i] <- rnorm(1, as.numeric(sample.draw %*% mu), as.numeric(sample.draw %*% sigma))
}

hist(y, col="lightblue", freq=FALSE, breaks=30, ylim=range(0, 0.45))
lines(density(y), col="blue", lwd=2)
temp <- seq(min(y), max(y), length.out=1000)
lines(temp, w[1]*dnorm(temp, mu[1], sigma[1]), col="red", lty=2, lwd=3) # the first component
lines(temp, w[2]*dnorm(temp, mu[2], sigma[2]), col="green", lty=2, lwd=3) # the second component
lines(temp, (w[1]*dnorm(temp, mu[1], sigma[1]))+(w[2]*dnorm(temp, mu[2], sigma[2])), col="black", lty=3, lwd=3) # the theoretical density



# Example 2: a mixture of 3 normals

# -- their weights are 0.2, 0.3 and 0.5, i.e., when we generate the mixture distribution, we choose the first normal 20% of the time, the second normal 30% of the time and the third normal 50% of the time

w <- c(0.2, 0.3, 0.5)
t(rmultinom(1, size=1, prob=w))

# -- the first normal distribution: mu[1]=-2, sigma[1]=1
# -- the second normal distribution: mu[2]=4, sigma[2]=1.5
# -- the third normal distribution: mu[3]=1, sigma[3]=0.5

mu <- c(-2, 4, 1)
sigma <- c(1, 1.5, 0.5)

(sample.draw <- t(rmultinom(1, size=1, prob=w)))
as.numeric(sample.draw %*% mu)
as.numeric(sample.draw %*% sigma)

# -- we generate 10000 observations from the mixture distribution:

n.obs <- 10000
y <- vector(length=n.obs)
str(y)

for (i in 1:n.obs) {
  sample.draw <- t(rmultinom(1, size=1, prob=w))
  y[i] <- rnorm(1, as.numeric(sample.draw %*% mu), as.numeric(sample.draw %*% sigma))
}

hist(y, col="lightblue", freq=FALSE, breaks=30, ylim=range(0, 0.85))
lines(density(y), col="blue", lwd=2)
temp <- seq(min(y), max(y), length.out=1000)
lines(temp, w[1]*dnorm(temp, mu[1], sigma[1]), col="red", lty=2, lwd=3) # the first component
lines(temp, w[2]*dnorm(temp, mu[2], sigma[2]), col="green", lty=2, lwd=3) # the second component
lines(temp, w[3]*dnorm(temp, mu[3], sigma[3]), col="darkred", lty=2, lwd=3) # the third component
lines(temp, (w[1]*dnorm(temp, mu[1], sigma[1]))+(w[2]*dnorm(temp, mu[2], sigma[2]))+(w[3]*dnorm(temp, mu[3], sigma[3])), col="black", lty=3, lwd=3) # the theoretical density



# Example 3 (really just a preliminary to example 4): a mixture of 4 normals with randomly assigned weights

# -- the weights need to sum to 1
# -- we can think of the probability from 0 to 1 as a stick and of the weights as a way of breaking the stick into 4 parts

# We break the stick as follows:
# 1. we randomly choose a proportion, i.e., a number between 0 and 1; we do this with a random draw from a Beta(1, alpha) distribution, where alpha is a parameter that we are free to set

alpha <- 5
rbeta(1, 1, alpha)
temp <- seq(0, 1, length.out=100)
plot(temp, dbeta(temp, 1, alpha), type="l", col="blue", lwd=2)

# -- the bigger the alpha, the smaller the proportion

alpha <- c(0.2, 0.7, 1, 2, 5, 10, 20, 40)
par(mfrow=c(2, 4))
temp <- seq(0, 1, length.out=100)
for (i in 1:length(alpha)) {
  plot(temp, dbeta(temp, 1, alpha[i]), type="l", col="blue", lwd=2, main=paste("alpha =", alpha[i]))
}
par(mfrow=c(1, 1))

# 2. we break a part of the stick proportional to our randomly drawn proportion and take that to be the first weight

alpha <- 5
(q1 <- rbeta(1, 1, alpha))

w <- vector(length=4)

remaining.stick <- 1
(w[1] <- q1*remaining.stick)
(remaining.stick <- (1-q1)*remaining.stick)

# 3. iterate this procedure 2 more times on the remaining stick

(q2 <- rbeta(1, 1, alpha))
(w[2] <- q2*remaining.stick)
(remaining.stick <- (1-q2)*remaining.stick)

(q3 <- rbeta(1, 1, alpha))
(w[3] <- q3*remaining.stick)
(remaining.stick <- (1-q3)*remaining.stick)

# 4. we set the final weight to the remaining stick or equivalently to 1 minus the sum of the previous weights

remaining.stick
(w[4] <- 1-sum(w[1:3]))

# -- the resulting randomly drawn weights are guaranteed to sum to 1 by construction
w
sum(w)

# -- just as before, we can randomly draw one of the 4 normal distributions based on these weights
t(rmultinom(1, size=1, prob=w))

# -- the first normal distribution: mu[1]=-2, sigma[1]=1
# -- the second normal distribution: mu[2]=4, sigma[2]=1.5
# -- the third normal distribution: mu[3]=1, sigma[3]=0.5
# -- the fourth normal distribution: mu[4]=2.5, sigma[4]=0.75

mu <- c(-2, 4, 1, 2.5)
sigma <- c(1, 1.5, 0.5, 0.75)

(sample.draw <- t(rmultinom(1, size=1, prob=w)))
as.numeric(sample.draw %*% mu)
as.numeric(sample.draw %*% sigma)

# -- we generate 10000 observations from the mixture distribution:

n.obs <- 10000
y <- vector(length=n.obs)
str(y)

for (i in 1:n.obs) {
  sample.draw <- t(rmultinom(1, size=1, prob=w))
  y[i] <- rnorm(1, as.numeric(sample.draw %*% mu), as.numeric(sample.draw %*% sigma))
}

hist(y, col="lightblue", freq=FALSE, breaks=30, ylim=range(0, 0.85))
lines(density(y), col="blue", lwd=2)
temp <- seq(min(y), max(y), length.out=1000)
lines(temp, w[1]*dnorm(temp, mu[1], sigma[1]), col="red", lty=2, lwd=3) # the first component
lines(temp, w[2]*dnorm(temp, mu[2], sigma[2]), col="green", lty=2, lwd=3) # the second component
lines(temp, w[3]*dnorm(temp, mu[3], sigma[3]), col="darkred", lty=2, lwd=3) # the third component
lines(temp, w[4]*dnorm(temp, mu[4], sigma[4]), col="darkgreen", lty=2, lwd=3) # the fourth component
lines(temp, (w[1]*dnorm(temp, mu[1], sigma[1]))+(w[2]*dnorm(temp, mu[2], sigma[2]))+(w[3]*dnorm(temp, mu[3], sigma[3]))+(w[4]*dnorm(temp, mu[4], sigma[4])), col="black", lty=3, lwd=3) # the theoretical density


# Example 4: a mixture of normals with a random number of components and with randomly assigned weights to these components

# The parameter alpha simultaneously controls both how many groups there are and their relative weights:
# -- small alpha: few groups, one-two of which have large weights
# -- large alpha: many groups, most of them have small weights

# We can keep breaking the stick until the remaining stick is below 0.0001, let's say, at which point we can stop because the remaining components of the mixture are going to have neglijible contributions to the mixture distribution

alpha <- 5
w <- vector()
remaining.stick <- 1

while (remaining.stick >= 0.0001) {
  current.q <- rbeta(1, 1, alpha)
  current.w <- current.q*remaining.stick
  w <- c(w, current.w)
  remaining.stick <- (1-current.q)*remaining.stick
}

w
sum(w)

w <- c(w, 1-sum(w))
w
sum(w)
remaining.stick

plot(w, pch=20, col="blue")


# -- just as before, we can randomly draw one of the normal distributions based on these weights

t(rmultinom(1, size=1, prob=w))


# We can draw the weights many times to see how many weights we obtain on average and the corresponding dispersion:

n.draws <- 10000
length.weights <- vector(length=n.draws)

alpha <- 5

for (i in 1:n.draws) {
  w <- vector()
  remaining.stick <- 1
  while (remaining.stick >= 0.0001) {
    current.q <- rbeta(1, 1, alpha)
    current.w <- current.q*remaining.stick
    w <- c(w, current.w)
    remaining.stick <- (1-current.q)*remaining.stick
  }
  w <- c(w, 1-sum(w))
  length.weights[i] <- length(w)
}

summary(length.weights)
plot(prop.table(table(length.weights)), col="blue", lwd=2)


# We can also plot the weights themselves to see the range of variability

n.draws <- 10000

alpha <- 5

plot(1, 1, type="n", xlim=c(1, 80), ylim=c(0, 1), pch=20, col="blue", xlab="", ylab="")
for (i in 1:n.draws) {
  w <- vector()
  remaining.stick <- 1
  while (remaining.stick >= 0.0001) {
    current.q <- rbeta(1, 1, alpha)
    current.w <- current.q*remaining.stick
    w <- c(w, current.w)
    remaining.stick <- (1-current.q)*remaining.stick
  }
  w <- c(w, 1-sum(w))
  points(1:length(w), w, col="blue", pch=19, cex=.75)
}



# Let's see the range of variability for multiple values of alpha

alpha <- c(0.2, 0.7, 1, 2, 5, 10, 20, 40)
n.draws <- 5000

par(mfrow=c(4, 4))

for (i in 1:length(alpha)) {
  plot(1, 1, type="n", xlim=c(1, 150), ylim=c(0, 1), pch=20, col="blue", xlab="", ylab="", main=paste("weights for alpha =", alpha[i]))
  length.weights <- vector(length=n.draws)
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha[i])
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    points(1:length(w), w, col="blue", pch=".") 
    length.weights[j] <- length(w)
  }
  plot(prop.table(table(length.weights)), col="blue", lwd=1.5, xlab="", ylab="", main=paste("no of weights for alpha =", alpha[i]))
}

par(mfrow=c(1, 1))



# Let us know generate mixtures of normals for alpha = 0.5, which will induce relatively few components most of the time

alpha <- 0.5

w <- vector()
remaining.stick <- 1
while (remaining.stick >= 0.0001) {
  current.q <- rbeta(1, 1, alpha)
  current.w <- current.q*remaining.stick
  w <- c(w, current.w)
  remaining.stick <- (1-current.q)*remaining.stick
}
w <- c(w, 1-sum(w))
w
sum(w)

length(w)

plot(w, type="h", col="blue", lwd=3)

# For simplicity, we will assume all components of the mixture have the same sd

sigma <- 2

# For every component, we need a mean; we generate these means as random draws from an underlying normal distribution centered at 0 and with a large sd as the components

(mu <- rnorm(length(w), 0, 5))


# Just as before, we can generate 10000 observations from the mixture distribution:

n.obs <- 10000
y <- vector(length=n.obs)
str(y)

for (i in 1:n.obs) {
  sample.draw <- t(rmultinom(1, size=1, prob=w))
  y[i] <- rnorm(1, as.numeric(sample.draw %*% mu), sigma)
}

hist(y, col="lightblue", freq=FALSE, breaks=30, ylim=range(0, 0.85))
lines(density(y), col="blue", lwd=2)
temp <- seq(min(y), max(y), length.out=1000)
mixture.density <- 0
for (i in 1:length(w)) {
  lines(temp, w[i]*dnorm(temp, mu[i], sigma), col=paste("brown", i, sep=""), lty=2, lwd=3)
  mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
}
lines(temp, mixture.density, col="black", lty=3, lwd=3) # the theoretical density


# Let's explore the range of variability that we can induce given alpha = 0.5, sigma = 2 and the underlying distribution for the means of the components being N(0, 5)

alpha <- 0.5
sigma <- 2
temp <- seq(-20, 20, length.out=1000)

par(mfrow=c(2, 2))
for (k in 1:4) {
  
  plot(temp, dnorm(temp, 0, 5), ylim=c(0, 0.5), type="l", col="grey", lwd=3, xlab="", ylab="", main=paste("alpha=", alpha, ", sigma=", sigma, sep=""))
  
  n.draws <- 5
  
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha)
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    mu <- rnorm(length(w), 0, 5)
    mixture.density <- 0
    for (i in 1:length(w)) {
      mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
    }
    lines(temp, mixture.density, col=j, lty=2, lwd=2)    
  }
  
}
par(mfrow=c(1, 1))


# We can change alpha to see what happens -- the bigger the alpha, the closer the mixtures get to the underlying distribution:

# -- alpha = 1

alpha <- 1
sigma <- 2
temp <- seq(-20, 20, length.out=1000)

par(mfrow=c(2, 2))
for (k in 1:4) {
  
  plot(temp, dnorm(temp, 0, 5), ylim=c(0, 0.5), type="l", col="grey", lwd=3, xlab="", ylab="", main=paste("alpha=", alpha, ", sigma=", sigma, sep=""))
  
  n.draws <- 5
  
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha)
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    mu <- rnorm(length(w), 0, 5)
    mixture.density <- 0
    for (i in 1:length(w)) {
      mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
    }
    lines(temp, mixture.density, col=j, lty=2, lwd=2)    
  }
  
}
par(mfrow=c(1, 1))


# -- alpha = 5

alpha <- 5
sigma <- 2
temp <- seq(-20, 20, length.out=1000)

par(mfrow=c(2, 2))
for (k in 1:4) {
  
  plot(temp, dnorm(temp, 0, 5), ylim=c(0, 0.5), type="l", col="grey", lwd=3, xlab="", ylab="", main=paste("alpha=", alpha, ", sigma=", sigma, sep=""))
  
  n.draws <- 5
  
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha)
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    mu <- rnorm(length(w), 0, 5)
    mixture.density <- 0
    for (i in 1:length(w)) {
      mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
    }
    lines(temp, mixture.density, col=j, lty=2, lwd=2)    
  }
  
}
par(mfrow=c(1, 1))

# -- alpha = 15

alpha <- 15
sigma <- 2
temp <- seq(-20, 20, length.out=1000)

par(mfrow=c(2, 2))
for (k in 1:4) {
  
  plot(temp, dnorm(temp, 0, 5), ylim=c(0, 0.5), type="l", col="grey", lwd=3, xlab="", ylab="", main=paste("alpha=", alpha, ", sigma=", sigma, sep=""))
  
  n.draws <- 5
  
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha)
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    mu <- rnorm(length(w), 0, 5)
    mixture.density <- 0
    for (i in 1:length(w)) {
      mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
    }
    lines(temp, mixture.density, col=j, lty=2, lwd=2)    
  }
  
}
par(mfrow=c(1, 1))

# -- alpha = 40

alpha <- 40
sigma <- 2
temp <- seq(-20, 20, length.out=1000)

par(mfrow=c(2, 2))
for (k in 1:4) {
  
  plot(temp, dnorm(temp, 0, 5), ylim=c(0, 0.5), type="l", col="grey", lwd=3, xlab="", ylab="", main=paste("alpha=", alpha, ", sigma=", sigma, sep=""))
  
  n.draws <- 5
  
  for (j in 1:n.draws) {
    w <- vector()
    remaining.stick <- 1
    while (remaining.stick >= 0.0001) {
      current.q <- rbeta(1, 1, alpha)
      current.w <- current.q*remaining.stick
      w <- c(w, current.w)
      remaining.stick <- (1-current.q)*remaining.stick
    }
    w <- c(w, 1-sum(w))
    mu <- rnorm(length(w), 0, 5)
    mixture.density <- 0
    for (i in 1:length(w)) {
      mixture.density <- mixture.density+(w[i]*dnorm(temp, mu[i], sigma))
    }
    lines(temp, mixture.density, col=j, lty=2, lwd=2)    
  }
  
}
par(mfrow=c(1, 1))