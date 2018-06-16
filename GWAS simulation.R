# Simulating GWAS number 1
library(rmutil)

spike <- function(x, A=0.02, mu=0, theta=.015, phi=2000)
  A*log((1 - theta)^(abs(x-mu)/theta) * phi + 1)/log(phi)
chr <- rep( 1:22, round(seq(87, 17, len = 22)))
chr_start <- c(which(!duplicated(chr)), length(chr))
chr_start <- (chr_start[-length(chr_start)] + chr_start[-1])/2
P <- length(chr) ## nunmber of features
loc <- 1:P
centers <- c(0.05,0.25,0.55,.8,0.94)*P
As <- c(0.1,0.1,0.125,0.1,0.05)
tmp <- sapply(seq_along(centers), function(i) spike(loc, mu = centers[i], A=As[i]))
effect <- 1 + rowSums(tmp)
set.seed(1)
baseline_p <- rmutil::rbetabinom(P, 100, 0.25, 100) / 100
N <- 5000
disease <- sapply(effect*baseline_p, function(p) rbinom(N, 1, p = p))
control <- sapply(baseline_p, function(p) rbinom(N, 1, p = p))
X <- rbind(control, disease)
if(any(colMeans(X) ==0 | colMeans(X)==1)) stop("One column with MAF=0")
y <- rep(c(0,1), each=nrow(X)/2)

# Running GWAS
res <- apply(X, 2, function(x){
  tab<-table(x,y)
  c(tab[1,1]*tab[2,2]/(tab[1,2]*tab[2,1]),
    chisq.test(tab)$p.value)
})
plot(loc, -log10(res[2,]))
plot(loc, -log10(res[2,]), pch=16, xaxt="n", xlab="Chromosome", ylab="-log (base 10) p-value", ylim=c(0, -log10(0.05/P)))
axis(1, chr_start, seq_along(chr_start), tick=FALSE)
abline(h=-log10(0.05/P), lty=2)

# Smoothing odds ratio
logodds <- log(res[1,])
fit <- predict(loess(logodds~loc, span = 0.1), se=TRUE)
mat <- fit$fit + cbind(-2*fit$se.fit,0,2*fit$se.fit)
matplot(loc, mat, type="l", col=c("grey","black","grey"), lty=1, ylim=max(mat)*c(-1,1), ylab="Log odds", xlab="Chromosome", xaxt="n")
axis(1, chr_start, seq_along(chr_start), tick=FALSE)
abline(h=0, lty=2)




#GWAS simulation number 2
disease2 <- sapply(effect*baseline_p, function(p) rbinom(N, 1, p = p))
control2 <- sapply(baseline_p, function(p) rbinom(N, 1, p = p))
X <- rbind(control, control2, disease, disease2)
if(any(colMeans(X) ==0 | colMeans(X)==1)) stop("One column with MAF=0")
y <- rep(c(0,1), each=nrow(X)/2)

res2 <- apply(X, 2, function(x){
  tab<-table(x,y)
  c(tab[1,1]*tab[2,2] / (tab[1,2]*tab[2,1]),chisq.test(tab)$p.value)
})
plot(loc, -log10(res2[2,]), pch=16, xaxt="n", xlab="Chromosome", ylab="-log (base 10) p-value")
axis(1, chr_start, seq_along(chr_start), tick=FALSE)
abline(h=-log10(0.05/P), lty=2)

plot(loc, effect, col=chr, pch=16, xaxt="n", xlab="Chromosome", ylab="Effect")
axis(1, chr_start, seq_along(chr_start), tick=FALSE)

par(mfrow=c(2,1))
plot(loc, -log10(res2[2,]), col=ifelse(effect>1.01, "black", "grey"), pch=16, xaxt="n",
     xlab="Chromosome", ylab="-log (base 10) p-value", main="Double the sample size")
axis(1, chr_start, seq_along(chr_start), tick=FALSE)
abline(h=-log10(0.05/P), lty=2)
o <- order(res[2,])[1:25] ##top 25 variants
plot(loc, -log10(res[2,]), col=ifelse(effect>1.01, "black", "grey"), pch=16, ylim=c(0, -log10(0.05/P)), xaxt="n", xlab="Chromosome", ylab="-log (base 10) p-value", main="Lower the threshold")
axis(1, chr_start, seq_along(chr_start), tick=FALSE)
abline(h = -log10(max(res[2,o])), lty=2)
