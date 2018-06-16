library(locfdr)
library(sfsmisc)

EBay.PRS.linear <- function(zval, samp, nulltype=0, df=7, bre=120) {
  
  ztoVg.func <-function(Z,samp){
    Z^2 /(samp-2 + Z^2)
  }
  
  bias.corr.kernel<- function (zz,xout,bw="nrd0")
  {
    density.obj=density(zz,bw=bw)
    fz.func = splinefun(density.obj$x, density.obj$y)
    Psi.z = log(  fz.func(zz) /dnorm(zz)  )
    truez = D1ss(x= zz, y=Psi.z ,xout=xout)
    return(truez)
  }
  
  # Compute Tweedie formula adjusted z-values and beta
  zall.corr.ker =  bias.corr.kernel(zval,xout=zval)
  Vg.corr.ker = ztoVg.func(zall.corr.ker, samp  )
  beta.corr= sqrt( Vg.corr.ker )*sign(zval) 
  
  #Compute local fdr
  fdrobj = locfdr(zval, nulltype=nulltype, df=df, bre=bre, plot=0)
  fdr = fdrobj$fdr
  
  beta.standardized = ztoVg.func(zval, samp)
  beta.corr.fdr = (1-fdr)*beta.standardized
  beta.corr.fdrXTweedie = (1-fdr)*beta.corr
  return(list(beta.corr.fdrXTweedie = beta.corr.fdrXTweedie, 
              beta.corr.Tweedie = beta.corr,
              beta.corr.fdr = beta.corr.fdr))
  
}


zval =rnorm(100)
EBay.PRS.linear(zval, samp=200)
