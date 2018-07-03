library(rms)
library(SemiPar)


### plot data
with(ses, plot(qs, ses[,1], pch=19))
xs <- seq(0,25,by=1)

### model with restricted cubic spline
test = spm(ses[,1]~qs,spar.method = "REML")

knots <- c(0.1,0.3,0.4,0.5,0.9)
res <- rma(ses[,1]~rcs(qs, knots), sei = ses[,2], data=ses)
lines(qs, predict(res)$pred, col="red", lwd=2)
lines(qs, predict(res)$ci.lb, col="blue", lwd=2)
lines(qs, predict(res)$ci.ub, col="blue", lwd=2)


test = smooth.spline(ses[,1]~qs, w = ses[,2], cv = T)
lines(predict(xz),lwd=2,col="purple")

