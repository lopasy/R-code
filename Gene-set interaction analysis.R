set.seed(1)
n <- 100000
p <- 200
Y <- matrix(rnorm(n), ncol=1)
Z <- matrix(rbinom(n*p, 2, 0.3), nrow=n)
E <- matrix(rnorm(n))
X <- matrix(rnorm(n*2), nrow=n)
set.seed(2)
Ybinary <- matrix(rbinom(n, 1,0.5), ncol=1)
GESAT(Z, Y, E, X)
GESAT(Z, Ybinary, E, X, out_type="D")


# Jan 29, 2014
# v22: change from "Davies" to "davies" and "Liu" to "liu"
# Jan 27, 2014
# v21 is modified from v20
# but v21 only keeps the functions required for logistic regression
# since most of the functions are shared between linear and logistic
# removes old functions, burden tests

#library(penalized)

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------
GxEscore.logistic.GCV <- function(Y, Xtilde, Z, V, type="davies",
                                  lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL){
  
  # Y (n x 1 matrix):  binary outcome variable
  # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
  # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
  # V (n x p matrix): GxE terms that we are testing
  # n = no. of samples
  # lower cannot be zero
  # to use a fixed lambda, set nintervals=1
  # NB: if nintervals=1, upper is used and lower is ignored
  
  # v20: add three new arguments scale.Z=T, weights.Z=NULL, weights.V=NULL
  # v20: if nintervals=1, upper is used and lower is ignored
  # v19: gives Liu p-value if Davies doesn't converge
  # modified in v16 to introduce davies method
  # v16 added 3 new arguments  lower, upper and type, otherwise unchanged
  # v16 added some initial checks
  # old version (< v16) is similar to GxEscore.logistic(),
  # except the ridge parameter is selected by GCV
  # v16 also removed p-df LRT
  # v16 has davies method, so isn't similar to  GxEscore.logistic()
  
  if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- drop(nrow(Y))
  
  #---------------------------------------------------------------------------
  # fit ridge regression model under the null
  # Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
  #---------------------------------------------------------------------------
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  Z <- scale(Z, center=T, scale=scale.Z)
  W <- cbind(rep(1,n), scale(Xtilde), Z)         # all the variables in the null model
  
  if(nintervals>1){
    lambdahat <- chooseridge.logistic(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                                      intervals=nintervals, plot=plotGCV, file=plotfile, center.Z=F, scale.Z=F, weights.Z=NULL)
  }else{
    lambdahat <- upper
  }
  
  model <- ridge.logistic(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
  Yhat <- model$Yhat
  sqrtvarhat_vec <- sqrt(Yhat*(1-Yhat))
  lambda <- model$lambda
  #varhat_vec <- Yhat*(1-Yhat)
  #varhat <- diag(varhat_vec)
  #sigmahat <- solve(varhat)
  
  
  #---------------------------------------------------------------------------
  # Score statistic
  #---------------------------------------------------------------------------
  if(is.null(weights.V)==F){
    V <- t(t(V) * (weights.V))
  }
  
  #Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
  Q1 <- t(Y-Yhat) %*% V  #v8c
  Q <- Q1 %*% t(Q1)      #v8c
  
  
  p <- ncol(Z)
  qtilde <- ncol(Xtilde) + 1
  # +1 is for intercept,
  # no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
  }
  
  #Wvarhat <- W * varhat_vec                         #v8c,varhat %*% W = W * varhat_vec
  #inverse <- solve(t(W) %*% Wvarhat + temp)         #v8c
  #inverse <- solve(t(W) %*% varhat %*% W + temp)
  # NB: sqrtvarhat_vec = sqrt(varhat_vec)
  M1 <- sqrtvarhat_vec * W                      # M1 = sqrt(varhat) %*% W
  M0 <- t(sqrtvarhat_vec * V)                   # M0 = V' %*% sqrt(varhat)
  #inverse <- solve(t(M1) %*% M1 +temp)         # t(M1) %*% M1 = t(W) %*% varhat %*% W
  #M2 <- M0 - M0 %*% M1 %*% inverse %*% t(M1)
  transM1 <- t(M1)
  invM1 <- solve(transM1 %*% M1 +temp, transM1) # invM1 = inverse %*% t(M1) = solve(t(M1) %*% M1 +temp) %*% t(M1)
  M2 <- M0 - M0 %*% M1 %*% invM1
  M3 <- M2 %*% t(M2)
  
  
  #---------------------------------------------------------------------------
  # p-value from non-central chi-square approximation
  #---------------------------------------------------------------------------
  if(type=="liu"){
    #H <- W  %*% inverse %*% t(Wvarhat)       #v8c
    #H <- W  %*% inverse %*% t(W) %*% varhat
    
    #R1 <- t(diag(n)-H) %*% (V * varhat_vec)  #v8c
    #R <- R1 %*% t(R1)
    #v8c, this is correct, old version shdnt have t(diag(n)-H) but shd be (diag(n)-H)
    #R <- t(diag(n)-H) %*% varhat %*% V %*% t(V) %*% varhat %*% t(diag(n)-H)
    
    #A <- R %*% sigmahat            #v8c
    #A2 <- A %*% A                  #v8c
    #kappa1 <-  mtrace(A)           #v8c
    #kappa2 <-  2*mtrace(A2)        #v8c
    #kappa3 <-  8*sum(A * t(A2))    #v8c
    #kappa4 <-  48*sum(A2 * t(A2))  #v8c
    
    #kappa1 <-  mtrace(R %*% sigmahat)
    #kappa2 <-  2*mtrace(R %*% sigmahat %*% R %*% sigmahat)
    #kappa3 <-  8*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)
    #kappa4 <-  48*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)
    
    M4 <- M3 %*% M3
    kappa1 <-  mtrace(M3)             # = tr(M3)
    kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
    kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
    kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
    
    approx <- noncentralapproxdirect(kappa2, kappa3, kappa4)
    Q.Norm <-((Q - kappa1)/sqrt(kappa2))*approx$sigmaX +  approx$muX
    pvalue <- pchisq(Q.Norm, df=approx$df, ncp = approx$ncp,  lower.tail=F)
    Is_converge <- 1
    
  }else{
    
    
    #---------------------------------------------------------------------------
    # p-value from davies
    #---------------------------------------------------------------------------
    # sqrt(varhat) %*% W = W * sqrt(varhat_vec)
    # W' %*% sqrt(varhat) = t(W*sqrt(varhat_vec))
    # M1 <- t(sigmahat * sqrt(varhat_vec)) - W  %*% inverse %*% t(W * sqrt(varhat_vec))
    # M2 <- t(V *varhat_vec) %*% M1
    # M3 <- M2 %*% t(M2)
    
    daviesout <- Get_PValue_GESAT(M3, Q)
    pvalue <- daviesout$p.value
    Is_converge <- daviesout$is_converge
    
    # new in v19: do nchisq/Liu if Davies doesn't converge
    if(Is_converge<=0){
      M4 <- M3 %*% M3
      kappa1 <-  mtrace(M3)             # = tr(M3)
      kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
      kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
      kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
      approx <- noncentralapproxdirect(kappa2, kappa3, kappa4)
      Q.Norm <-((Q - kappa1)/sqrt(kappa2))*approx$sigmaX +  approx$muX
      pvalue <- pchisq(Q.Norm, df=approx$df, ncp = approx$ncp,  lower.tail=F)
    }
    
  }
  
  return(list(pvalue=pvalue, Is_converge=Is_converge, lambda=drop(model$lambda)))
  
}



ridge.logistic <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
  # in v20, functions were renamed
  # ridge.logistic() and ridge.select.logistic() were modified from ridge.logistic.lambda() in v19
  # modified in v20 to add three more arguments center.Z=T, scale.Z=T, weights.Z=NULL
  # ridge.logistic() and ridge.select.logistic() are similar but return different things
  # ridge.select.logistic() returns only lambda, GCV, effective.df
  # ridge.logistic() returns only lambda, Yhat, thetahat
  # note that Z should always be centered as scale() behaves weirdly for scaling if not centered
  
  #============================================================
  # Method 2: Use penalized package - do not center Y
  #============================================================
  # when unpenalized is specified as a matrix, intercept has to be specified explicitly
  # ridge.logistic.lambda() includes intercept, so do not add intercept as part of Xtilde
  # uses penalized package (load it externally)
  # Xtilde is a n*qtilde matrix,
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of binrary outcomes
  # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for Yhat
  
  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  
  n <- nrow(Y)
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  
  Z <- scale(Z, center=center.Z, scale=scale.Z)
  test2 <- penalized(response=Y, penalized=Z,  unpenalized=cbind(rep(1,n), scale(Xtilde)), lambda1=0, lambda2=lambda, model="logistic")
  if (test2@converged==F){ # new in v8b
    thetahat <- NA
    Yhat <- NA
  }else{
    thetahat <- coefficients(test2, "all")
    Yhat <- fitted(test2)
  }
  return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat))
}



ridge.select.logistic <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
  # in v20, functions were renamed
  # ridge.logistic() and ridge.select.logistic() were modified from ridge.logistic.lambda() in v19	
  # modified in v20 to add three more arguments center.Z=T, scale.Z=T, weights.Z=NULL
  # ridge.logistic() and ridge.select.logistic() are similar but return different things
  # ridge.select.logistic() returns only lambda, GCV, effective.df
  # ridge.logistic() returns only lambda, Yhat, thetahat
  # note that Z should always be centered as scale() behaves weirdly for scaling if not centered
  
  #============================================================
  # Method 2: Use penalized package - do not center Y
  #============================================================
  # when unpenalized is specified as a matrix, intercept has to be specified explicitly
  # ridge.logistic.lambda() includes intercept, so do not add intercept as part of Xtilde
  # uses penalized package (load it externally)
  # Xtilde is a n*qtilde matrix,
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of binrary outcomes
  # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for Yhat
  
  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  
  n <- nrow(Y) 
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  
  Z <- scale(Z, center=center.Z, scale=scale.Z)
  test2 <- penalized(response=Y, penalized=Z,  unpenalized=cbind(rep(1,n), scale(Xtilde)), lambda1=0, lambda2=lambda, model="logistic")
  if (test2@converged==F){ # new in v8b
    thetahat <- NA
    Yhat <- NA
    GCV <- NA
    effective.df <- NA
  }else{
    thetahat <- coefficients(test2, "all")
    Yhat <- fitted(test2)
    sqrtvarhat_vec <- sqrt(Yhat*(1-Yhat))
    
    if(sum(sqrtvarhat_vec<=0)>0){
      GCV <- NA
      effective.df <- NA
    }else{
      W <- cbind(rep(1,n), scale(Xtilde), Z)         # all the variables in the null model
      #varhat_vec <- Yhat*(1-Yhat)                           #v8c
      #varhat <- diag(varhat_vec)                            #v8c
      p <- ncol(Z)
      qtilde <- ncol(Xtilde) + 1
      # +1 is for intercept, no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
      if(p==1){
        temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
      }else{
        temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
      }
      
      #Wvarhat <- W * varhat_vec
      #v8c, varhat %*% W =  W * varhat_vec since varhat is diagonal
      #inverse <- try(solve(t(W)  %*% Wvarhat + temp))          #v8c
      #inverse <- try(solve(t(W) %*% varhat %*% W + temp))     # try fn new in v8B
      
      M1 <- sqrtvarhat_vec * W             	        # M1 = sqrt(varhat) %*% W
      MM1<- t(M1) %*% M1			# MM1 = M1' M1 = W' sqrt(varhat) sqrt(varhat) W = W' varhat W
      inverse <- try(solve(MM1+temp))
      
      error.inverse <- 0                                       # new in v8b
      if (is(inverse, "try-error")) error.inverse <- 1         # new in v8b
      if(error.inverse==1){
        GCV <- NA
        effective.df <- NA
      }else{
        
        #H <- W  %*% inverse %*% t(Wvarhat)              	                #v8c
        #H <- W  %*% inverse %*% t(W) %*% varhat
        #effective.df <- mtrace(H)-1                      	                # subtract 1 for the intercept
        equivH <- MM1 %*% inverse		 		        # tr(H) = tr(equivH)
        effective.df <- mtrace(equivH)-1               		            # subtract 1 for the intercept
        #GCV <- (t(Y-Yhat) %*% sigmahat %*% (Y-Yhat))/(n*(1-effective.df/n)^2)  #new GCV in v8
        GCV <- (sum(((Y-Yhat)/sqrtvarhat_vec)^2))/(n*(1-effective.df/n)^2)			
      }
      
    }	
  }
  return(list(lambda=lambda, GCV = GCV, effective.df=effective.df))
  #return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat, GCV = GCV, effective.df=effective.df))
}


chooseridge.logistic <- function(Y, Xtilde, Z, lambdastart=1e-20, lambdaend=sqrt(nrow(Y))/log(nrow(Y)),
                                 intervals=5, plot=F, file=NA, center.Z=T, scale.Z=T, weights.Z=NULL)
{
  # renamed chooseridge.logistic in v20 (used to be chooseridge)
  # modified in v20 to add center.Z=T, scale.Z=T, weights.Z=NULL as arguments
  # modified in v16 to add lambdastart and lambdaend as arguments, otherwise unchanged
  # need lambdastart, lambdaend >=0, lambdastart < lambdaend
  
  lambda <- c(exp(seq(log(lambdastart),log(lambdaend),length=intervals)))
  output <- c()
  for(ii in 1:length(lambda)){
    temp <- ridge.select.logistic(Y, Xtilde, Z, lambda[ii], center.Z, scale.Z, weights.Z)$GCV
    output <- c(output, temp)
    rm(temp)
  }
  
  lambdafinal <- lambda[which.min(output)]
  if(plot==T&is.na(file)==T){
    plot(lambda, output, xlab="lambda", ylab="GCV")
    abline(v=lambdafinal, col="red")
  }
  
  if(plot==T&is.na(file)==F){
    pdf(file)
    plot(lambda, output, xlab="lambda", ylab="GCV")
    abline(v=lambdafinal, col="red")
    dev.off()
  }
  
  return(lambdafinal)
}

#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to logistic
#-------------------------------------------------------------------------------------------------------------


# Jan 29, 2014
# v12: change from "Davies" to "davies" and "Liu" to "liu"
# Version 11: Jan 27, 2014
# identical to V10, except BT functions removed
# Verstion 10: Oct 1, 2013
# identical to Version 9, only SKAT_davies() function modified so that .qfc() calls from iSKAT
# Version 9: Sept 25, 2013
# identical to Version 8, only SKAT_davies() function modified so that .qfc() calls from SKAT
# this is to resolve namespace problems.
# Version 8: Jan 28, 2013
# change to GxEscore.linear.GCV() such that varhat can accomodate larger p
# Version 7: Dec 20, 2012
# change to ridge.select.linear() and ridge.linear()
# such that an intercept is now included but Y is not centered
# Version 6: Dec 17, 2012
# the only change is to add a try() function in ridge.select.linear()
# such that the GCV will only select a model that converges/ matrix is invertible
# Version 5: March 15, 2012
# Changes in V5:
# GxEscore() renamed to GxEscore.linear.GCV() for consistency compared to logistic model
# chooseridge() renamed to chooseridge.linear()
# ridge() renamed to ridge.linear()
# ridge.select() renamed to ridge.select.linear()
# GxEscore.linear.GCV() implements Davies method from v5 onwards
# GxEscore.linear.GCV have changes made to speed up computation, adds 3 arugments related to scale and weights
# chooseridge.linear() modified slightly, add 3 other arguments center.Z, scale.Z, and weights.Z
# ridge.linear and ridge.select.linear() modified to add 3 arguments center.Z, scale.Z, and weights.Z, matrix computation made faster 
# v5 adds 3 helper functions to get davies p-value which are called by GxEscore.linear.GCV()
# SKAT_davies()
# Get_Lambda() # modified to include only.values=T in eigen()
# Get_PValue_GESAT() 
# All 3 helper functions were copied unmodified from GxE-scoretest-logistic-snpset-v19.R except Get_Lambda() #modified to include only.values=T in eigen()

# Version 4 sets the upper limit of lambda to be 9*floor(sqrt(n)), where n=sample size, modification made in GxEscore()
# Version 3 sets the upper limit of lambda to be 3*floor(sqrt(n)), where n=sample size, modification made in GxEscore() and chooseridge()
# Version 2 was modified so that GxEscore() calculates min p-value as well, changed the input to function
# Version 1 adapted from GxE-scoretest-v9.R
# converts it into a function
# chooseridge() was modified such that lambda=0 is not an option
# p-value calculated without perturbations

#----------------------------------------------------------------------------------------------------------
# Note: these are the functions for linear (v5) and logistic regression (v20) respectively:
# GxEscore.logistic.GCV()                            : main function     
# ridge.logistic()                                   : fit final null model
# chooseridge.logistic() and ridge.select.logistic() : select ridge parameter
# Burdentests.GE.logistic                            : Burden tests

# GxEscore.linear.GCV()                              : main function
# ridge.linear()                                     : fit final null model     
# chooseridge.linear() and ridge.select.linear()     : select ridge parameter
# Burdentests.GE.linear                              : Burden tests

# The analogous functions in each take in exactly the same arguments in linear and logistic codes
# ridge.logistic() returns slightly different values from ridge.linear()
# other analogous functions return the same values in both linear and logistic codes
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------
GxEscore.linear.GCV <- function(Y, Xtilde, Z, V, type="davies",
                                lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL){
  
  # Y (n x 1 matrix):  continuous outcome variable
  # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
  # Do NOT include intercept as part of Xtilde (Y is centered in ridge, so no intercept needed)
  # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
  # V (n x p matrix): GxE terms that we are testing
  # n = no. of samples
  # lower cannot be zero
  # to use a fixed lambda, set nintervals=1
  # NB: if nintervals=1, upper is used and lower is ignored
  
  if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- drop(nrow(Y))
  
  #---------------------------------------------------------------------------
  # fit ridge regression model under the null
  # Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
  #---------------------------------------------------------------------------
  if(nintervals>1){
    if(is.null(weights.Z)==F & scale.Z==F){
      Z <- t(t(Z) * (weights.Z))
    }
    Z <-  scale(Z, center=T, scale=scale.Z) 
    lambdahat <- chooseridge.linear(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                                    intervals=nintervals, plot=plotGCV, file=plotfile, 
                                    center.Z=F, scale.Z=F, weights.Z=NULL)
    ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
    Yhat <- ridgemodel$Yhat
  }else{
    lambdahat <- upper
    ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=T, scale.Z=scale.Z, weights.Z=weights.Z)
    Yhat <- ridgemodel$Yhat
  }
  #---------------------------------------------------------------------------
  # Score statistic
  #---------------------------------------------------------------------------
  
  if(is.null(weights.V)==F){
    V <- t(t(V) * (weights.V))
  }
  
  # Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
  Q1 <- t(Y-Yhat) %*% V
  Q <- Q1 %*% t(Q1)
  
  
  #varhat <- var(Y-Yhat)                             # Change in GxE-scoretest-v8.R
  df1 <- sum(ridgemodel$W * t(ridgemodel$invW))      # Change in GxE-scoretest-v8.R
  varhat <- var(Y-Yhat) * (n-1) / (n - df1)          # Change in GxE-scoretest-v8.R
  
  # R <- t(diag(n)-H) %*% V %*% t(V) %*% (diag(n)-H)
  # R1 <- t(diag(n)-H) %*% V
  # R <- R1 %*% t(R1)
  #M1 <- t(V) - t(V) %*% W %*% inverse %*% t(W)
  M1 <- t(V) - t(V) %*% ridgemodel$W %*% ridgemodel$invW     # W = ridgemodel$W = all variables under the null with appropriate scaling, centering, invW = inverse %*% t(W)
  M2 <- M1 %*% t(M1)
  M3 <- drop(varhat)*M2
  
  #---------------------------------------------------------------------------
  # p-value from non-central chi-square approximation
  #---------------------------------------------------------------------------
  if(type=="liu"){
    
    #sigmahat <- drop(varhat)*diag(n)
    #A <- R %*% sigmahat
    #A2 <- A %*% A
    #kappa1 <-  mtrace(A)
    #kappa2 <-  2*mtrace(A2)
    #kappa3 <-  8*sum(A * t(A2))
    #kappa4 <-  48*sum(A2 * t(A2))
    
    # kappa1 <-  mtrace(R %*% sigmahat)
    # kappa2 <-  2*mtrace(R %*% sigmahat %*% R %*% sigmahat)
    # kappa3 <-  8*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)
    # kappa4 <-  48*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)
    
    M4 <- M3 %*% M3
    kappa1 <-  mtrace(M3)	          # = tr(M3)
    kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
    kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
    kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
    
    approx2 <- noncentralapproxdirect(kappa2, kappa3, kappa4)
    Q.Norm2 <-((Q - kappa1)/sqrt(kappa2))*approx2$sigmaX +  approx2$muX
    pvalue <- pchisq(Q.Norm2, df=approx2$df, ncp = approx2$ncp,  lower.tail=F)
    Is_converge <- 1
    
  }else{
    
    #---------------------------------------------------------------------------
    # p-value from davies
    #---------------------------------------------------------------------------
    #M3 <- drop(varhat)*R
    daviesout <- Get_PValue_GESAT(M3, Q)
    pvalue <- daviesout$p.value
    Is_converge <- daviesout$is_converge
    
    if(Is_converge<=0){
      #sigmahat <- drop(varhat)*diag(n)
      #A <- R %*% sigmahat
      #A2 <- A %*% A
      #kappa1 <-  mtrace(A)
      #kappa2 <-  2*mtrace(A2)
      #kappa3 <-  8*sum(A * t(A2))
      #kappa4 <-  48*sum(A2 * t(A2))
      
      M4 <- M3 %*% M3
      kappa1 <-  mtrace(M3)             # = tr(M3)
      kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
      kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
      kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
      
      approx2 <- noncentralapproxdirect(kappa2, kappa3, kappa4)
      Q.Norm2 <-((Q - kappa1)/sqrt(kappa2))*approx2$sigmaX +  approx2$muX
      pvalue <- pchisq(Q.Norm2, df=approx2$df, ncp = approx2$ncp,  lower.tail=F)
    }
    
  }
  
  return(list(pvalue=pvalue, Is_converge=Is_converge, lambda=drop(lambdahat)))
}


ridge.linear <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
  # modified in v7, intercept is included, Y is not centered
  # modified in v5 to return W, invW instead of H and GCV, matrix computations made faster, 
  # add 3 arguments center.Z=T, scale.Z=T, weights.Z=NULL
  # computes the ridge estimator for each value of the ridge parameter, lambda
  # Xtilde is a n*qtilde matrix, 
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of outcomes
  # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for yhat, (qtilde+p)*(qtilde+p) matrix for Hat matrix
  # when lambda=0, Yhat is the same as that from predict(lm(Y~Xtilde+Z))
  # NB: if scale.Z=T, regardless of what weights.Z is, this corresponds to beta(MAF; 0.5, 0.5) weight
  # to use weights specified as weights.Z, set scale.Z=F
  
  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  #if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- nrow(Y)
  qtilde <- ncol(Xtilde) + 1 		 # +1 is for intercept, new in v7
  p <- ncol(Z)	
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
  }
  
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  W <- cbind(rep(1,n), scale(Xtilde), scale(Z, center=center.Z, scale=scale.Z)) # doesn't matter if Xtilde is scaled or not
  #Ys <- Y - mean(Y)
  
  #inverse <- solve(t(W) %*% W + temp)	
  #invW <- inverse %*% t(W)
  transW <- t(W)
  invW <- solve(transW %*% W + temp, transW) #invW <- solve(t(W) %*% W + temp, t(W))
  #thetahat <- inverse %*% t(W) %*% Ys
  #thetahat <- invW %*% Ys
  thetahat <- invW %*% Y
  
  Yhat <- W %*% thetahat
  #Yshat <- W %*% thetahat
  #Yhat <- Yshat + mean(Y)
  #H <-  W %*% inverse %*% t(W)
  #GCV <- sum((Ys-Yshat)^2)/(n*(1-sum(diag(H))/n)^2)
  #equivH <- invW %*% W          # tr(H) = tr(equivH)
  #GCV <- sum((Ys-Yshat)^2)/(n*(1-sum(diag(equivH))/n)^2)
  
  return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat, W=W, invW=invW))
  #return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat, inverse=inverse, W=W, invW=invW))
}


ridge.select.linear <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
  # modified in v7 to include intercept, but not center Y
  # modified in v6 to add a try() function to matrix inversion
  # modified in v5 to add 3 arguments center.Z=T scale.Z=T, weights.Z=NULL and matrix computations made faster
  # function is identical to ridge.linear(), except it only returns GCV and lambda
  # computes the ridge estimator for each value of the ridge parameter, lambda
  # Xtilde is a n*qtilde matrix, 
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of outcomes
  # when lambda=0, Yhat is the same as that from predict(lm(Y~Xtilde+Z))
  # NB: if scale.Z=T, regardless of what weights.Z is, this corresponds to beta(MAF; 0.5, 0.5) weight
  # to use weights specified as weights.Z, set scale.Z=F
  
  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  #if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- nrow(Y)
  qtilde <- ncol(Xtilde) + 1			 # +1 is for intercept, new in v7
  p <- ncol(Z)	
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
  }
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  W <- cbind(rep(1,n), scale(Xtilde), scale(Z, center=center.Z, scale=scale.Z)) # doesn't matter if Xtilde is scaled or not
  #Ys <- Y - mean(Y)
  
  #inverse <- solve(t(W) %*% W + temp)
  #invW <- inverse %*% t(W)
  transW <- t(W)
  
  # change in v6 to add a try function	 
  invW <- try(solve(transW %*% W + temp, transW)) #invW <- solve(t(W) %*% W + temp, t(W))
  if(is(invW, "try-error")){
    GCV <- NA
    effective.df <- NA
  }else{
    
    #thetahat <- inverse %*% t(W) %*% Ys
    #thetahat <- invW %*% Ys
    thetahat <- invW %*% Y
    Yhat <- W %*% thetahat
    #Yshat <- W %*% thetahat
    #Yhat <- Yshat + mean(Y)
    #H <-  W %*% inverse %*% t(W)
    #GCV <- sum((Ys-Yshat)^2)/(n*(1-sum(diag(H))/n)^2)
    equivH <- invW %*% W         					# tr(H) = tr(equivH)
    effective.df <- mtrace(equivH)-1  			 	# -1 for intercept
    #GCV <- sum((Ys-Yshat)^2)/(n*(1-effective.df/n)^2)
    GCV <- sum((Y-Yhat)^2)/(n*(1-effective.df/n)^2)
  }
  
  return(list(lambda=lambda, GCV = GCV, effective.df=effective.df))  
}


chooseridge.linear <- function(Y, Xtilde, Z, lambdastart=1e-20, lambdaend=sqrt(nrow(Y))/log(nrow(Y)), intervals=5, plot=F, file=NA, center.Z=T, scale.Z=T, weights.Z=NULL){
  # modified in v5
  # need lambdastart, lambdaend >=0, lambdastart <= lambdaend
  lambda <- c(exp(seq(log(lambdastart),log(lambdaend),length=intervals)))
  
  output <- c()
  for(ii in 1:length(lambda)){
    temp <- ridge.select.linear(Y, Xtilde, Z, lambda[ii], center.Z, scale.Z, weights.Z)$GCV
    output <- c(output, temp)
    rm(temp)
  }
  lambdafinal <- lambda[which.min(output)]
  
  if(plot==T&is.na(file)==T){
    plot(lambda, output, xlab="lambda", ylab="GCV")
    abline(v=lambdafinal, col="red")
  }
  
  if(plot==T&is.na(file)==F){
    pdf(file)
    plot(lambda, output, xlab="lambda", ylab="GCV")
    abline(v=lambdafinal, col="red")
    dev.off()
  }
  
  return(lambdafinal)
}

#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# START: common functions shared by both linear and logistic models
# NB: Functions common to both have same names in both scripts
# NB: Common functions are identical in both GxE-scoretest-logistic-snpset-v19.R and GxE-scoretest-snpset-v5.R
#-------------------------------------------------------------------------------------------------------
Beta.Weights<-function(MAF,weights.beta){
  # copied without modification from Beta.Weights() SKAT 0.71
  
  n<-length(MAF)
  weights<-rep(0,n)	
  IDX_0<-which(MAF == 0)
  if(length(IDX_0) == n){
    stop("No polymorphic SNPs")
  } else if( length(IDX_0) == 0){
    weights<-dbeta(MAF,weights.beta[1],weights.beta[2])
  } else {
    weights[-IDX_0]<-dbeta(MAF[-IDX_0],weights.beta[1],weights.beta[2])
  }
  
  #print(length(IDX_0))
  #print(weights[-IDX_0])
  return(weights)
  
}


Get_PValue_GESAT <- function(K,Q){
  # copied without modification from GxE-scoretest-logistic-snpset-v20.R
  # new in v5, helper function to get davies p-value
  # modified from Get_PValue() from SKAT 0.71
  
  lambda <- Get_Lambda(K)
  out <- SKAT_davies(Q, lambda, acc=10^(-6))
  
  p.val <- out$Qq
  is_converge <- 1
  
  # check p-value
  if(p.val > 1 || p.val< 0 ){
    p.val <- NA
    is_converge <- -2
  }
  
  
  # check convergence
  if(length(lambda) == 1){
    p.val <-  NA
    is_converge <- -1
  } else if(out$ifault != 0){
    p.val <- NA
    is_converge <- 0
  }
  
  
  return(list(p.value=p.val, is_converge=is_converge))
  
}


Get_Lambda <- function(K){
  # copied without modification from GxE-scoretest-logistic-snpset-v19.R except to include only.values=T in eigen()
  # identical to Get_Lambda() in GxE-scoretest-logistic-snpset-v20.R
  # new in v5, helper function to get davies p-value
  # copied without modification from SKAT 0.71
  
  out.s <- eigen(K,symmetric=TRUE, only.values=T)
  #print(out.s$values)
  
  lambda1 <- out.s$values
  IDX1<-which(lambda1 >= 0)
  
  # eigenvalue bigger than sum(eigenvalues)/1000        
  IDX2 <- which(lambda1 > mean(lambda1[IDX1])/100000)
  
  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda <- lambda1[IDX2]
  lambda
  
}


SKAT_davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  # copied without modification from GxE-scoretest-logistic-snpset-v19.R
  # new in v5, helper function to get davies p-value
  # copied without modification from SKAT 0.71
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0))
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}


qqplots <- function(pvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.png"){
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  png(file = filename, width = 800, height = 800)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19, 
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,	
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


qqplotspdf <- function(pvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.pdf"){
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  pdf(file = filename, width = 6, height = 6)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19, 
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,	
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


qqplotspdfNA <- function(inputpvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.pdf"){
  # similar to qqplotspdf except it first removes NA values
  pvalues <- inputpvalues[is.na(inputpvalues)==F]
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  pdf(file = filename, width = 6, height = 6)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19,
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


skewness <-  function(x){
  return((mean((x-mean(x))^3))/(sd(x)^3))
}


kurtosis <- function(x){  
  return((mean((x-mean(x))^4))/(sd(x)^4)) 
}


noncentralapproxdirect <- function(k2, k3, k4){
  # k2, k3, k4 are the 2nd, 3rd and 4th cumulants
  
  s1 <- k3/((k2^1.5)*sqrt(8))
  s2 <- k4/(k2*k2*12)
  
  if(s1^2>s2){
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  }else{
    a <- 1/s1
    delta <- 0
    l <- 1/s1^2
  }
  return(list(df=l, ncp=delta, muX=l+delta, sigmaX=sqrt(2)*a))
}


mtrace <- function(X){
  return(sum(diag(X)))
}


mafcall <- function(x){
  return(min(sum(x, na.rm=T), 2*sum(is.na(x)==F)-sum(x, na.rm=T))/(2*sum(is.na(x)==F)))
}


checkpolymorphic <- function(x){
  return(length(unique(x))>1)
}
#-------------------------------------------------------------------------------------------------------
# END: common functions shared by both linear and logistic models
#
#-------------------------------------------------------------------------------------------------------


# Version 4: June 18, 2015
# change one line of code so that it can handle small p-values better
# changed for Get_Liu_PVal.MOD.Lambda_iSKAT()
# so that instead of using 1- pchisq(), use pchisq(lower.tail=F)

# Version 3: Jan 29, 2014
# Version 3 removed all the functions repeated in GxE-scoretest-snpset-v11.R
# Remaining functions are identical to those in Version 2.
# Version 3 adds a function: checkvariation2()


# Version 2: Feb 2, 2013
# (1) Added three new functions (shared between iSKAT.logistic() and iSKAT.linear()):
# Get_PValue.Lambda_iSKAT()
# Get_Liu_PVal.MOD.Lambda_iSKAT()
# Get_Liu_Params_Mod_Lambda_iSKAT()
# These three functions are only used for the optimal test iSKAT.logistic() and iSKAT.linear()
# Added one line to the function iSKAT_Optiaml_Each_Q()
# Details of the four changes in v2 (all changes apply to both logistic and linear):
# (i)   one-liner change in iSKAT_Optiaml_Each_Q()
# calls
# (ii)  Get_PValue.Lambda_iSKAT()
# which calls
# (iii) Get_Liu_PVal.MOD.Lambda_iSKAT()
# which calls
# (iv)  Get_Liu_Params_Mod_Lambda_iSKAT()
#
# (2) Changes to four iSKAT.linear() and GxEscore.linear.GCV() specific functions
# Changes are made so that functions for linear regression are consistent between GxE-scoretest-snpset-v8.R, iSKAT-Linear-v2.R
# (i) Change to iSKAT.linear() so that varhat can accomodate larger p
# The following three functions are also changed so that they are identical to those in GxE-scoretest-snpset-v8.R
# Version 8: Jan 28, 2013
# (ii) change to GxEscore.linear.GCV() such that varhat can accomodate larger p
# Version 7: Dec 20, 2012
# (iii) change to ridge.select.linear() and ridge.linear()
# such that an intercept is now included but Y is not centered
# Version 6: Dec 17, 2012
# (iv) the only change is to add a try() function in ridge.select.linear()
# such that the GCV will only select a model that converges/ matrix is invertible


# Version 1: April 16, 2012
# iSKAT-Linear
# modified from GxE-scoretest-snpset-v5.R
# Functions that are retained from GxE-scoretest-snpset-v5.R are identical
# Thus any GESAT related functions are the same in both GxE-scoretest-snpset-v5.R and iSKAT-Linear-v1.R
# iSKAT-Linear adds the function iSKAT.linear() and all other functions it needs
# All functions added to iSKAT-Linear-v1.R from GxE-scoretest-snpset-v5.R have names starting with iSKAT

#----------------------------------------------------------------------------------------------------------
# Note: these are the functions for linear (GxE-scoretest-snpset-v5.R, iSKAT-Linear-v1.R) 
# and logistic regression (GxE-scoretest-logistic-snpset-v20.R, iSKAT-Logistic-v1.R) respectively:
# GxEscore.logistic.GCV(), iSKAT.logistic()          : main function     
# ridge.logistic()                                   : fit final null model
# chooseridge.logistic() and ridge.select.logistic() : select ridge parameter
# Burdentests.GE.logistic                            : Burden tests

# GxEscore.linear.GCV(), iSKAT.linear()              : main function
# ridge.linear()                                     : fit final null model     
# chooseridge.linear() and ridge.select.linear()     : select ridge parameter
# Burdentests.GE.linear                              : Burden tests

# The analogous functions in each take in exactly the same arguments in linear and logistic codes
# ridge.logistic() returns slightly different values from ridge.linear()
# other analogous functions return the same values in both linear and logistic codes
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------
iSKAT.linear <- function(Y, Xtilde, Z, V, type="davies",
                         lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL, r.corr=(0:10)/10){
  
  # Y (n x 1 matrix):  continuous outcome variable
  # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
  # Do NOT include intercept as part of Xtilde (Y is centered in ridge, so no intercept needed)
  # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
  # V (n x p matrix): GxE terms that we are testing
  # n = no. of samples
  # lower cannot be zero
  # to use a fixed lambda, set nintervals=1
  # NB: if nintervals=1, upper is used and lower is ignored
  
  # iSKAT.linear() Modified from GxEscore.linear.GCV() v5 of code
  # added 1 argument: r.corr(), type is changed to lower-case
  # parts on fitting the null model are identical and unchanged
  
  # iSKAT.linear() v2 changed varhat estimate to make it consistent with GxEscore.linear.GCV() v8
  
  if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- drop(nrow(Y))
  
  #---------------------------------------------------------------------------
  # fit ridge regression model under the null
  # Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
  #---------------------------------------------------------------------------
  if(nintervals>1){
    if(is.null(weights.Z)==F & scale.Z==F){
      Z <- t(t(Z) * (weights.Z))
    }
    Z <-  scale(Z, center=T, scale=scale.Z) 
    lambdahat <- chooseridge.linear(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                                    intervals=nintervals, plot=plotGCV, file=plotfile, center.Z=F, scale.Z=F, weights.Z=NULL)
    ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
    Yhat <- ridgemodel$Yhat
  }else{
    lambdahat <- upper
    ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=T, scale.Z=scale.Z, weights.Z=weights.Z)
    Yhat <- ridgemodel$Yhat
  }
  #---------------------------------------------------------------------------
  # Score statistic
  #---------------------------------------------------------------------------
  
  if(is.null(weights.V)==F){
    V <- t(t(V) * (weights.V))
  }
  
  
  #---------------------------------------------------------------------------
  # new in iSKAT (everything above this line is unchanged from GESAT v5, except for changing "Davies" to "davies" and "Liu" to "liu")
  #---------------------------------------------------------------------------
  res <- Y - Yhat
  # Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
  # Q1 <- t(res) %*% V
  # Q <- Q1 %*% t(Q1)
  
  #varhat <- drop(var(res))                                # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
  df1 <- drop(sum(ridgemodel$W * t(ridgemodel$invW)))      # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
  varhat <- drop(var(res)) * (n-1) / (n - df1)             # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
  V1 <- V - ridgemodel$W %*% ridgemodel$invW %*% V         # V1=t(diag(n)-H) %*% V, H = t(H) for linear regression only
  results <- iSKAT_Optimal_Linear(res, V, X1=NULL, kernel=NULL, weights = NULL, s2=varhat, method=type, res.out=NULL,     
                                  n.Resampling =0, r.all=r.corr, V1)
  
  
  return(list(pvalue=results$p.value, param=results$param, lambda=drop(lambdahat)))
}



iSKAT_Optimal_Linear = function(res, V, X1=NULL, kernel=NULL, weights = NULL, s2, method=NULL
                                , res.out=NULL, n.Resampling =0, r.all, V1){
  
  # Modified from SKAT v0.73 SKAT_Optimal_Linear()
  # Note that arguments like X1, kernel, weights, res.out, n.Resampling are never used
  # they are kept for consistency with the corresponding SKAT function
  # all the Z's in original function are renamed V
  # added 1 argument V1
  
  # if r.all >=0.999 ,then r.all = 0.999. It is just for computation.
  IDX<-which(r.all >= 0.999)
  if(length(IDX) > 0){
    r.all[IDX]<-0.999	
  }
  
  n<-dim(V)[1]
  p.m<-dim(V)[2]	
  n.r<-length(r.all)
  
  
  ###########################################
  # Compute Q.r and Q.r.res
  ##########################################
  out.Q<-iSKAT_Optimal_Get_Q(V, res, r.all, n.Resampling, res.out)
  Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2
  
  ##################################################
  # Compute P-values 
  #################################################
  
  out<-iSKAT_Optimal_Get_Pvalue(Q.all, V1 / sqrt(2), r.all, method)
  
  param<-list(p.val.each=NULL,q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)
  
  id_temp<-which(param$p.val.each == min(param$p.val.each))
  id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
  if(length(id_temp1) > 0){
    param$rho[id_temp1] = 1
  }
  param$rho_est<-param$rho[id_temp]
  
  
  p.value<-out$p.value[1]
  
  p.value.resampling<-NULL
  if(n.Resampling > 0){
    p.value.resampling<-out$p.value[-1]
    #param$pval.each.resample<-out$p.val.each[-1]
  }
  
  re<-list(p.value = p.value, p.value.resampling = p.value.resampling
           , Test.Type = method, Q = NA, param=param )  
  
  return(re)	
  
}
#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# START: common functions shared by both linear and logistic models
# NB: Functions common to both have same names in both scripts
# NB: Common functions are identical in both GxE-scoretest-logistic-snpset-v19.R and GxE-scoretest-snpset-v5.R
#-------------------------------------------------------------------------------------------------------
iSKAT_Optimal_Param<-function(Z1,r.all){
  
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  # Function get parameters of optimal test
  n<-dim(Z1)[1]
  p.m<-dim(Z1)[2]
  r.n<-length(r.all)
  
  z_mean<-rowMeans(Z1)
  Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
  cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)
  
  Z.item1<-Z_mean %*% diag(cof1)	
  Z.item2<-Z1 - Z.item1
  
  # W3.2 Term : mixture chisq
  W3.2.t<-t(Z.item2) %*% Z.item2
  lambda<-Get_Lambda(W3.2.t)
  
  # W3.3 Term : variance of remaining ...
  W3.3.item<-sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
  
  # Mixture Parameters
  MuQ<-sum(lambda)
  VarQ<-sum(lambda^2) *2 + W3.3.item
  KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
  Df<-12/KerQ
  
  # W3.1 Term : tau1 * chisq_1
  tau<-rep(0,r.n)
  for(i in 1:r.n){
    r.corr<-r.all[i]
    #term1<-p.m*r.corr + cof1^2 * (1-r.corr)
    term1<-p.m^2*r.corr + sum(cof1^2) * (1-r.corr)
    tau[i]<-sum(term1) *  sum(z_mean^2)
  }
  
  out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau)
  return(out)
}



iSKAT_Optiaml_Each_Q<-function(param.m, Q.all, r.all, lambda.all){
  
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  #       Function get SKAT statistics with given rho
  #               Q.all is a matrix with n.q x n.r
  
  n.r<-length(r.all)
  c1<-rep(0,4)
  n.q<-dim(Q.all)[1]
  
  pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
  pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
  param.mat<-NULL
  
  for(i in 1:n.r){
    Q<-Q.all[,i]
    r.corr<-r.all[i]
    lambda.temp<-lambda.all[[i]] 
    c1[1]<-sum(lambda.temp)
    c1[2]<-sum(lambda.temp^2)
    c1[3]<-sum(lambda.temp^3)
    c1[4]<-sum(lambda.temp^4)
    param.temp<-iSKAT_Get_Liu_Params_Mod(c1)
    
    muQ<-param.temp$muQ
    varQ<-param.temp$sigmaQ^2
    df<-param.temp$l
    
    # get pvalue
    Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
    pval[,i]<- 1-pchisq(Q.Norm,  df = df)
    
    #----------------------------
    # one-line addition in iSKAT v2
    pval[,i]<-Get_PValue.Lambda_iSKAT(lambda.temp,Q)$p.value
    #---------------------------- 
    param.mat<-rbind(param.mat,c(muQ,varQ,df))
  }
  
  pmin<-apply(pval,1,min)
  for(i in 1:n.r){
    
    muQ<-param.mat[i,1]
    varQ<-param.mat[i,2]
    df<-param.mat[i,3]
    
    q.org<-qchisq(1-pmin,df=df)
    q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
    pmin.q[,i]<-q.q
    
  }
  
  out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
  return(out)
  
}


iSKAT_Optimal_Integrate_Func_Davies<-function(x,pmin.q,param.m,r.all){
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  
  n.r<-length(r.all)
  n.x<-length(x)
  
  temp1<-param.m$tau %x% t(x)
  
  temp<-(pmin.q - temp1)/(1-r.all)
  temp.min<-apply(temp,2,min)
  
  re<-rep(0,length(x))
  for(i in 1:length(x)){
    #a1<<-temp.min[i]
    min1<-temp.min[i]
    if(min1 > sum(param.m$lambda) * 10^4){
      temp<-0
    } else {
      min1.temp<- min1 - param.m$MuQ			
      sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
      min1.st<-min1.temp *sd1 + param.m$MuQ
      
      dav.re<-SKAT_davies(min1.st,param.m$lambda,acc=10^(-6))
      temp<-dav.re$Qq
      if(dav.re$ifault != 0){
        stop("dav.re$ifault is not 0")
      }
    }
    if(temp > 1){
      temp=1
    }
    #lambda.record<<-param.m$lambda
    #print(c(min1,temp,dav.re$ifault,sum(param.m$lambda)))
    re[i]<-(1-temp) * dchisq(x[i],df=1)
  }
  return(re)
  
}


iSKAT_Optimal_PValue_Davies<-function(pmin.q,param.m,r.all){
  
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  re<-try(integrate(iSKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15), silent = TRUE)
  if(class(re) == "try-error"){
    re<-iSKAT_Optimal_PValue_Liu(pmin.q,param.m,r.all)
    return(re)
  } 
  
  pvalue<-1-re[[1]]
  if(pvalue < 0){
    pvalue=0
  }
  return(pvalue)
  
}


iSKAT_Optimal_Integrate_Func_Liu<-function(x,pmin.q,param.m,r.all){
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT	
  #x<-1
  #print(length(x))
  #print(x)
  #X1<<-x
  #x<-X1
  
  n.r<-length(r.all)
  n.x<-length(x)
  
  temp1<-param.m$tau %x% t(x)
  
  temp<-(pmin.q - temp1)/(1-r.all)
  temp.min<-apply(temp,2,min)
  
  temp.q<-(temp.min - param.m$MuQ)/sqrt(param.m$VarQ)*sqrt(2*param.m$Df) + param.m$Df
  re<-pchisq(temp.q ,df=param.m$Df) * dchisq(x,df=1)
  return(re)
  
}


iSKAT_Optimal_PValue_Liu<-function(pmin.q,param.m,r.all){
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  
  re<-integrate(iSKAT_Optimal_Integrate_Func_Liu, lower=0, upper=30, subdivisions=500
                ,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15)
  
  pvalue<-1-re[[1]]
  return(pvalue)
  
}


iSKAT_Optimal_Get_Q<-function(Z1, res, r.all, n.Resampling, res.out, res.moments=NULL){
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  
  n.r<-length(r.all)
  p.m<-dim(Z1)[2]
  
  Q.r<-rep(0,n.r)
  Q.r.res<-NULL
  Q.sim<-NULL	
  
  temp<-t(res) %*% Z1
  for(i in 1:n.r){
    r.corr<-r.all[i]
    Q1<-(1-r.corr) * rowSums(temp^2)
    Q2<-r.corr * p.m^2 * rowMeans(temp)^2
    Q.r[i]<-Q1 + Q2
  }
  Q.r = Q.r /2
  if(n.Resampling > 0){
    
    temp<-t(res.out) %*% Z1
    Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
    for(i in 1:n.r){
      r.corr<-r.all[i]
      Q1<-(1-r.corr) * rowSums(temp^2)
      Q2<-r.corr * p.m^2 * rowMeans(temp)^2
      Q.r.res[,i]<-Q1 + Q2
    }
    Q.r.res = Q.r.res/2
  }
  
  if(!is.null(res.moments)){
    
    temp<-t(res.moments) %*% Z1
    n.moments<-dim(res.moments)[2]
    Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
    for(i in 1:n.r){
      r.corr<-r.all[i]
      Q1<-(1-r.corr) * rowSums(temp^2)
      Q2<-r.corr * p.m^2 * rowMeans(temp)^2
      Q.sim[,i]<-Q1 + Q2
    }
    Q.sim = Q.sim/2
    
  }
  
  re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.sim=Q.sim)
  return(re)
  
  
}


iSKAT_Optimal_Get_Pvalue<-function(Q.all, Z1, r.all, method){
  # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
  
  n.r<-length(r.all)
  n.q<-dim(Q.all)[1]
  p.m<-dim(Z1)[2]
  
  lambda.all<-list()
  for(i in 1:n.r){
    r.corr<-r.all[i]
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Z2<- Z1 %*% t(L)
    K1<-t(Z2) %*% Z2
    
    lambda.all[[i]]<-Get_Lambda(K1)
    
  }
  
  # Get Mixture param 
  param.m<-iSKAT_Optimal_Param(Z1,r.all)
  Each_Info<-iSKAT_Optiaml_Each_Q(param.m, Q.all, r.all, lambda.all)
  pmin.q<-Each_Info$pmin.q
  pval<-rep(0,n.q)
  
  if(method == "davies" || method=="optimal"){
    
    for(i in 1:n.q){
      pval[i]<-iSKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all)
    }
    
    
  } else if(method =="liu" || method =="liu.mod"){
    
    for(i in 1:n.q){
      pval[i]<-iSKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all)
    }
    
  } else {
    stop("Invalid Method!")
  }
  return(list(p.value=pval,p.val.each=Each_Info$pval))
  
}


iSKAT_Get_Liu_Params_Mod<-function(c1){
  # copied without modification from SKAT 0.73, except to prefix of name to iSKAT
  
  ## Helper function for getting the parameters for the null approximation
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}



Get_Liu_PVal.MOD.Lambda_iSKAT <- function(Q.all, lambda){
  # new in iSKAT v2, copied without modification from SKAT v0.81, renamed to append iSKAT
  param<-Get_Liu_Params_Mod_Lambda_iSKAT(lambda)
  
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  # change in version 4 the line below
  p.value<- pchisq(Q.Norm1,  df = param$l, ncp=param$d, lower.tail=FALSE)
  
  return(p.value)
  
}

Get_Liu_Params_Mod_Lambda_iSKAT <- function(lambda){
  ## new in iSKAT v2, copied without modification from SKAT v0.81, renamed to append iSKAT
  ## Helper function for getting the parameters for the null approximation
  
  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }
  
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}


Get_PValue.Lambda_iSKAT<-function(lambda,Q){
  # new in iSKAT v2, copied with a one-line modification from SKAT v0.81, renamed to append iSKAT 
  #print(lambda)
  n1<-length(Q)
  
  p.val<-rep(0,n1)
  p.val.liu<-rep(0,n1)
  is_converge<-rep(0,n1)
  p.val.liu<-Get_Liu_PVal.MOD.Lambda_iSKAT(Q, lambda)
  
  for(i in 1:n1){
    out<-SKAT_davies(Q[i],lambda,acc=10^(-6))
    
    p.val[i]<-out$Qq
    #p.val.liu[i]<-SKAT_liu(Q[i],lambda)
    
    is_converge[i]<-1
    
    # check convergence
    if(length(lambda) == 1){
      p.val[i]<-p.val.liu[i]
    } else if(out$ifault != 0){
      is_converge[i]<-0
      p.val[i]<-p.val.liu[i]  # this one-line modification is different from SKAT v0.81
    }
    
    # check p-value
    if(p.val[i] > 1 || p.val[i] <= 0 ){
      is_converge[i]<-0
      p.val[i]<-p.val.liu[i]
    }
  }
  
  return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge))
  
}



checkvariation <- function(x){
  # added in iSKAT v1
  # on top of checking that x is polymorphic, checks that there are at least 2 indiv in each level of x
  x <- x[is.na(x)==F]
  return(min(summary(as.factor(x)))>1&length(unique(x))>1)
}

checkvariation2 <- function(x){
  # added in iSKAT v3: to allow for dosage data, function returns TRUE only if there are more than 2 levels 
  # or if only 2 levels, there are no singletons
  # on top of checking that x is polymorphic, checks that there are at least 2 indiv in each level of x
  x <- x[is.na(x)==F]
  return((min(summary(as.factor(x)))>1&length(unique(x))==2)|length(unique(x))>2)
}

#-------------------------------------------------------------------------------------------------------
# END: common functions shared by both linear and logistic models
#
#-------------------------------------------------------------------------------------------------------


# Jan 29, 2014
# Version 3 keeps only functions not repeated in iSKAT-linear-v3.R and GxE-scoretest-snpset-v12.R 
# remaining functions are identical to those in v2
# Version 2: Feb 1, 2013
# Added three new functions (shared between iSKAT.logistic() and iSKAT.linear()):
# Get_PValue.Lambda_iSKAT()
# Get_Liu_PVal.MOD.Lambda_iSKAT()
# Get_Liu_Params_Mod_Lambda_iSKAT()
# These three functions are only used for the optimal test iSKAT.logistic() and iSKAT.linear()
# Added one line to the function iSKAT_Optiaml_Each_Q()
# Details of the four changes in v2 (all changes apply to both logistic and linear):
# (i)   one-liner change in iSKAT_Optiaml_Each_Q()
# calls
# (ii)  Get_PValue.Lambda_iSKAT()
# which calls
# (iii) Get_Liu_PVal.MOD.Lambda_iSKAT()
# which calls
# (iv)  Get_Liu_Params_Mod_Lambda_iSKAT()
#
# Version 1: April 23, 2012
# iSKAT-Logistic
# modified from GxE-scoretest-logistic-snpset-v20.R and iSKAT-Linear-v1.R
# Functions that are retained from GxE-scoretest-logistic-snpset-v20.R and iSKAT-Linear-v1.R are identical
# Thus any GESAT related functions are the same in  GxE-scoretest-logistic-snpset-v20.R/iSKAT-Linear-v1.R and iSKAT-Logistic-v1.R
# iSKAT-Linear adds the function iSKAT.linear() and all other functions it needs
# All functions added to iSKAT-Linear-v1.R from GxE-scoretest-snpset-v5.R have names starting with iSKAT

#----------------------------------------------------------------------------------------------------------
# Note: these are the functions for linear (GxE-scoretest-snpset-v5.R, iSKAT-Linear-v1.R) 
# and logistic regression (GxE-scoretest-logistic-snpset-v20.R, iSKAT-Logistic-v1.R) respectively:
# GxEscore.logistic.GCV(), iSKAT.logistic()          : main function     
# ridge.logistic()                                   : fit final null model
# chooseridge.logistic() and ridge.select.logistic() : select ridge parameter
# Burdentests.GE.logistic                            : Burden tests

# GxEscore.linear.GCV(), iSKAT.linear()              : main function
# ridge.linear()                                     : fit final null model     
# chooseridge.linear() and ridge.select.linear()     : select ridge parameter
# Burdentests.GE.linear                              : Burden tests

# The analogous functions in each take in exactly the same arguments in linear and logistic codes
# ridge.logistic() returns slightly different values from ridge.linear()
# other analogous functions return the same values in both linear and logistic codes
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------
iSKAT.logistic <- function(Y, Xtilde, Z, V, type="davies",
                           lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL, r.corr=(0:10)/10){
  
  # Y (n x 1 matrix):  binary outcome variable
  # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
  # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
  # V (n x p matrix): GxE terms that we are testing
  # n = no. of samples
  # lower cannot be zero
  # to use a fixed lambda, set nintervals=1
  # NB: if nintervals=1, upper is used and lower is ignored
  
  # iSKAT.logistic() Modified from GxEscore.logistic.GCV() v20 of code
  # added 1 argument: r.corr(), type is changed to lower-case
  # parts on fitting the null model are identical and unchanged
  
  if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  n <- drop(nrow(Y))
  
  #---------------------------------------------------------------------------
  # fit ridge regression model under the null
  # Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
  #---------------------------------------------------------------------------
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  Z <- scale(Z, center=T, scale=scale.Z)
  W <- cbind(rep(1,n), scale(Xtilde), Z)         # all the variables in the null model
  
  if(nintervals>1){
    lambdahat <- chooseridge.logistic(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                                      intervals=nintervals, plot=plotGCV, file=plotfile, center.Z=F, scale.Z=F, weights.Z=NULL)
  }else{
    lambdahat <- upper
  }
  
  model <- ridge.logistic(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
  Yhat <- model$Yhat
  sqrtvarhat_vec <- sqrt(Yhat*(1-Yhat))
  lambda <- model$lambda
  #varhat_vec <- Yhat*(1-Yhat)
  #varhat <- diag(varhat_vec)
  #sigmahat <- solve(varhat)
  
  
  #---------------------------------------------------------------------------
  # Score statistic
  #---------------------------------------------------------------------------
  if(is.null(weights.V)==F){
    V <- t(t(V) * (weights.V))
  }
  
  
  #---------------------------------------------------------------------------
  # new in iSKAT (everything above this line is unchanged from GESAT v20, except for changing "Davies" to "davies" and "Liu" to "liu")
  #---------------------------------------------------------------------------
  res <- Y - Yhat
  
  #Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
  #Q1 <- t(Y-Yhat) %*% V  #v8c
  #Q <- Q1 %*% t(Q1)      #v8c
  
  
  p <- ncol(Z)
  qtilde <- ncol(Xtilde) + 1
  # +1 is for intercept,
  # no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
  }
  
  M1 <- sqrtvarhat_vec * W                      # M1 = sqrt(varhat) %*% W
  M00 <-sqrtvarhat_vec * V		      # M00 = sqrt(varhat) %*% V
  transM1 <- t(M1)
  invM1 <- solve(transM1 %*% M1 +temp, transM1) # invM1 = inverse %*% t(M1) = solve(t(M1) %*% M1 +temp) %*% t(M1)
  V1 <- M00 -  M1 %*% invM1 %*% M00
  
  results <- iSKAT_Optimal_Logistic(res, V, X1=NULL, kernel=NULL, weights = NULL, pi_1=NULL, method=type, 
                                    res.out=NULL, n.Resampling =0, r.all=r.corr, V1)
  
  
  return(list(pvalue=results$p.value, param=results$param, lambda=drop(lambdahat)))
  
  
  
}







iSKAT_Optimal_Logistic  = function(res, V, X1=NULL, kernel=NULL, weights = NULL, pi_1=NULL , method = NULL
                                   , res.out=NULL, n.Resampling =0, r.all, V1){
  
  
  # Modified from SKAT v0.73 SKAT_Optimal_Logistic()
  # Note that arguments like X1, kernel, weights, pi_1, res.out, n.Resampling are never used
  # they are kept for consistency with the corresponding SKAT function
  # all the Z's in original function are renamed V
  # added 1 argument V1
  
  # if r.all >=0.999 ,then r.all = 0.999
  IDX<-which(r.all >= 0.999)
  if(length(IDX) > 0){
    r.all[IDX]<-0.999	
  }
  
  n<-dim(V)[1]
  p.m<-dim(V)[2]
  n.r<-length(r.all)
  
  
  ###########################################
  # Compute Q.r and Q.r.res
  ##########################################
  out.Q<-iSKAT_Optimal_Get_Q(V, res, r.all, n.Resampling, res.out)
  Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) 
  
  ##################################################
  # Compute P-values 
  #################################################
  
  out<-iSKAT_Optimal_Get_Pvalue(Q.all, V1 / sqrt(2), r.all, method)
  
  param<-list(p.val.each=NULL,q.val.each=NULL)
  param$p.val.each<-out$p.val.each[1,]
  param$q.val.each<-Q.all[1,]
  param$rho<-r.all
  param$minp<-min(param$p.val.each)
  
  id_temp<-which(param$p.val.each == min(param$p.val.each))
  id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
  if(length(id_temp1) > 0){
    param$rho[id_temp1] = 1
  }
  
  param$rho_est<-param$rho[id_temp]
  
  
  p.value<-out$p.value[1]
  p.value.resampling<-NULL
  if(n.Resampling > 0){
    p.value.resampling<-out$p.value[-1]
    #param$pval.each.resample<-out$p.val.each[-1]
  }
  
  re<-list(p.value = p.value, p.value.resampling = p.value.resampling
           , Test.Type = method, Q = NA, param=param )  
  
  return(re)	
  
}
#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------


# Last Edited: October 14, 2015
# Add a check for missingness in Z
# if is_check_genotype==FALSE and is_dosage== FALSE,
# then Z cannot have missigneness
# iSKAT/GESAT core functions cannot have missingness in all variables (the functions work assuming no missingness)
# and the checks for missingness occur in main.R
# in main.R, imputation for missing values takes place only IF is_check_genotype = TRUE OR is_dosage==TRUE
# thus if both of these are false, Z must have no missingness
# updated main.R such that if Z has missingness and is_check_genotype==FALSE and is_dosage== FALSE, it throws an error
#
#
# Last Edited: January 29, 2014
# v1.* of package onwards names the common variant method as GESAT and the rare variant method as iSKAT
# Four changes to GESAT()
# also, fixed two bugs related to weights.V in GESAT() main function call:
# (1)a one line bug to if(ncol(Z)!= length(weights.V)) 
# (2) weights.V are handled by directly multiplying it to Z, before GE matrix=V is formed
# previously if weights.V were specified, there wd be problems due to bugs
# (3)also now all colliner columns with Z are removed added "&(apply(abs(cor(V,Z))>0.9999999999,1,sum,na.rm=T)==0)"
# (4)also, changed to "davies" and "liu"
# iSKAT_Main_Check_Z() is also modified to add a check: if(is.na(sd1)==FALSE)
# Last Edited: January 27, 2014 
# v0.6 onwards allow for binary outcome

GESAT <- function(Z, Y, E, X=NULL, type="davies",
                  lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=FALSE, plotfile=NA, scale.Z=TRUE, weights.Z=NULL, weights.V=NULL, 
                  out_type="C", impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, SetID=NULL){
  
  #-------------------------------------
  # check outcome type
  if(out_type != "C" && out_type != "D"){
    stop("Invalid out_type!. Please use only \"C\" for continous outcome and \"D\" for dichotomous outcome.")
  }
  
  
  #-------------------------------------
  # check dimensions
  if(nrow(E)!= nrow(Y)) stop("Dimensions of E and Y don't match.")
  if(nrow(Z)!= nrow(Y)) stop("Dimensions of Z and Y don't match.")
  if(is.null(X)==FALSE){
    if(nrow(X)!= nrow(Y)) stop("Dimensions of X and Y don't match.")
    if(class(X)!= "matrix") stop("X is not a matrix.")
  }
  if(class(Z)!= "matrix") stop("Z is not a matrix.")
  if(class(E)!= "matrix") stop("E is not a matrix.")
  if(class(Y)!= "matrix") stop("Y is not a matrix.")
  
  #----------------------------------------- added on Oct 14, 2015 (start)
  if(is_check_genotype==FALSE & is_dosage== FALSE){               
    if(sum(is.na(Z))!= 0) stop("Z cannot have any missing values if is_check_genotype = FALSE and is_dosage=FALSE.")
  }
  #----------------------------------------- added on Oct 14, 2015 (end)
  
  if(is.null(weights.Z)==FALSE){
    if(ncol(Z)!= length(weights.Z)) stop("Dimensions of Z and weights.Z don't match.")
    if(sum(weights.Z<0)!=0) stop("weights.Z have to be non-negative and cannot be missing.")
  }
  
  if(is.null(weights.V)==FALSE){
    if(ncol(Z)!= length(weights.V)) stop("dimensions of V and weights.V don't match")  #change!!!
    if(sum(weights.V<0)!=0) stop("weights.V have to be non-negative and cannot be missing")
  }
  
  
  #-------------------------------------
  # check other arguments
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  
  
  #-------------------------------------
  # check for missing values in Y, E, X
  if(is.null(X)==FALSE){
    if(sum(is.na(X))!= 0) stop("X cannot have any missing values.")
  }
  
  if(sum(is.na(E))!= 0) stop("E cannot have any missing values.")
  if(sum(is.na(Y))!= 0) stop("Y cannot have any missing values.")
  
  
  #-------------------------------------
  # check that X doesn't have intercept
  if(is.null(X)==FALSE){
    
    if(ncol(X)==1){
      if(checkpolymorphic(X)==FALSE) stop("X should not include intercept and must have more than one level.")	
    }else{
      if(sum(apply(X, 2, checkpolymorphic))!= ncol(X)) stop("X should not include intercept and must have more than one level.")
      
    }
  }
  
  
  #-------------------------------------
  # check that E has more than one levels
  if(ncol(E)==1){
    if(checkpolymorphic(E)==FALSE) stop("E must have more than one level.")	
  }else{
    if(sum(apply(E, 2, checkpolymorphic))!= ncol(E)) stop("E must have more than one level.")
    
  }
  
  
  #-------------------------------------
  # check Z and impute
  if(is_dosage ==TRUE){
    impute.method="fixed"
  }
  
  if(is_check_genotype==TRUE | is_dosage==TRUE){
    Z.out <- iSKAT_MAIN_Check_Z(Z=Z, SetID=SetID, weights.Z=weights.Z, weights.V=weights.V, impute.method=impute.method,  									missing_cutoff=missing_cutoff)
    Z <- as.matrix(Z.out$Z.test)
    weights.Z <- as.vector(Z.out$weights.Z.test)
    weights.V <- as.vector(Z.out$weights.V.test)
  }
  
  if(is.null(Z)==TRUE){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has no SNPs." )
    } else {
      msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=0, n.GE.test=NA))
  }
  
  if(ncol(Z)<2){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has fewer than 2 SNPs. Do not need penalized procedures." )
    } else {
      msg <- sprintf("In %s, the Z matrix has fewer than 2 SNPs. Do not need penalized procedures.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=NA))
  }
  
  
  
  #-------------------------------------
  # check V
  Ztemp <- Z
  if(is.null(weights.V)==FALSE){
    Ztemp <- t(t(Z) * (weights.V))
  } 
  
  if(ncol(E)==1){
    V <- as.matrix(drop(E)*Ztemp)
  }else{
    V <- NULL
    for(hhh in 1:ncol(E)){
      V <- as.matrix(cbind(V, drop(E[,hhh])*Ztemp))
    }
  }
  
  Idx.V.check <- apply(V, 2, checkpolymorphic)&(apply(abs(cor(V,Z))>0.9999999999,1,sum,na.rm=T)==0)
  
  if(sum(Idx.V.check)==0){
    
    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix is empty. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix is empty.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))	}
  
  
  if(sum(Idx.V.check)==1){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix has fewer than 2 SNPs. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix has fewer than 2 SNPs.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=1))
  }
  
  V <- as.matrix(V[,Idx.V.check])
  
  
  
  #-------------------------------------
  # Run GESAT
  Xtilde <- as.matrix(cbind(X,E))
  if(out_type == "C"){
    iSKAT.out <- GxEscore.linear.GCV(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                     lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL)
  } 
  
  if(out_type == "D"){	   		
    iSKAT.out <- GxEscore.logistic.GCV(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                       lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL)
  } 	   		   
  
  return(list(pvalue=iSKAT.out$pvalue, Is_converge=iSKAT.out$Is_converge, 
              lambda=iSKAT.out$lambda, n.G.test=drop(ncol(Z)), n.GE.test=drop(ncol(V))))
}














#
#	Check the Z, and do imputation
#     modified significantly from SKAT_MAIN_Check_Z from V0.78
#
iSKAT_MAIN_Check_Z <- function(Z, SetID, weights.Z=NULL, weights.V=NULL, impute.method,  missing_cutoff){
  
  # check.Z.error = 0 : no snps removed, but some snps possibly imputed
  # check.Z.error = 1 : all snps removed, returns NULL matrix for Z
  # check.Z.error = 2 : some snps removed, remainder snps may have been imputed
  
  check.Z.error <- 0
  n <- nrow(Z)
  ##############################################
  # Recode Missing to NA
  
  IDX_MISS <- union(which(is.na(Z)), which(Z == 9))
  if(length(IDX_MISS) > 0){
    Z[IDX_MISS] <- NA
  } 
  
  ###################################################
  # Check missing rates and exclude any SNPs with missing rate > missing_cutoff
  # Also exclude non-polymorphic SNPs
  m <- ncol(Z)
  ID_INCLUDE_SNP <- NULL
  for(i in 1:m){
    missing.ratio <- length(which(is.na(Z[,i])))/n
    sd1 <- sd(Z[,i], na.rm=TRUE)
    if(is.na(sd1)==FALSE){  # change in v1.0 of iSKAT package, Jan 29, 2014
      if(missing.ratio < missing_cutoff && sd1 > 0){
        ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
      }
    }
  }
  
  if(length(ID_INCLUDE_SNP) == 0){
    
    if(is.null(SetID)){
      msg <- sprintf("ALL SNPs have either high missing rates or no-variation. ")
    } else {
      msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation. ",SetID )
    }
    warning(msg, call.=FALSE)
    
    re <- list(Z.test=NULL, weights.Z.test=NULL, weights.V.test=NULL, check.Z.error=1) 
    
    
  } else{
    
    if(m - length(ID_INCLUDE_SNP) > 0 ){
      
      if(is.null(SetID)){
        msg <- sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
      } else {
        msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
      }
      
      warning(msg, call.=FALSE)	
      Z <- as.matrix(Z[,ID_INCLUDE_SNP])
      if(is.null(weights.Z)==FALSE) weights.Z <- weights.Z[ID_INCLUDE_SNP]
      if(is.null(weights.V)==FALSE) weights.V <- weights.V[ID_INCLUDE_SNP]
      check.Z.error <- 2
      IDX_MISS <- which(is.na(Z))
    }
    
    
    
    ##########################################
    # Missing Imputation
    if(length(IDX_MISS) > 0){
      if(is.null(SetID)){
        msg <- sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
      } else {
        msg <- sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
      }
      
      warning(msg,call.=FALSE)
      Z <- Impute_iSKAT(Z,impute.method)
    } 
    re <- list(Z.test=Z, weights.Z.test=weights.Z, weights.V.test=weights.V, check.Z.error=check.Z.error)
    
  }
  
  return(re)
}




# copied without modifcation from SKAT V0.78, renamed Impute_iSKAT()
# Simple Imputation
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing has to be NA: a missing genotype value.

Impute_iSKAT <-function(Z, impute.method){
  
  p <- dim(Z)[2]
  
  if(impute.method =="random"){
    for(i in 1:p){
      IDX <- which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1 <- mean(Z[-IDX,i])/2
        Z[IDX,i] <- rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1 <- mean(Z[-IDX,i])/2
        Z[IDX,i] <- 2 * maf1
      }
    }
  } else {
    stop("Error: Imputation method shoud be either \"fixed\" or \"random\"! ")
  }
  
  return(Z)
}


# Last Edited: October 14, 2015
# Add a check for missingness in Z
# if is_check_genotype==FALSE and is_dosage== FALSE,
# then Z cannot have missigneness
# iSKAT/GESAT core functions cannot have missingness in all variables (the functions work assuming no missingness)
# and the checks for missingness occur in main.R
# in main.R, imputation for missing values takes place only IF is_check_genotype = TRUE OR is_dosage==TRUE
# thus if both of these are false, Z must have no missingness
# updated main.R such that if Z has missingness and is_check_genotype==FALSE and is_dosage== FALSE, it throws an error
#
#
# Last Edited: January 29, 2014 
# v1.* of package onwards names the common variant method as GESAT and the rare variant method as iSKAT
# iSKAT_MAIN_Check_Z_and_V() is new in v1.* of iSKAT package
# also removes common variants

iSKAT <- function(Z, Y, E, X=NULL, type="davies",
                  lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=FALSE, plotfile=NA, scale.Z=FALSE, weights.Z=NULL, weights.V=NULL, 
                  out_type="C", impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, 
                  r.corr=(0:10)/10, weights.beta=c(1,25), MAF_cutoff=0.05, SetID=NULL){
  
  # if weights.beta is specified and weights.Z is not, weights.Z will be set to beta weights
  # if weights.beta is specified and weights.V is not, weights.V will be set to beta weights
  # if weights.Z is specified, weights.beta is ignored for Z
  # if weights.V is specified, weights.beta is ignored for V
  # if scale.Z is T, weights.beta is ignored for Z and weights.Z is also ignored
  # NB: in v1.* package onwards, iSKAT main function has three more arguments r.corr,  weights.beta, MAF_cutoff than GESAT main function
  # the default on scale.Z is also different.
  # in iSKAT the default is to NOT scale Z and use weights.beta=c(1,25) weights
  
  #-------------------------------------
  # check outcome type
  if(out_type != "C" && out_type != "D"){
    stop("Invalid out_type!. Please use only \"C\" for continous outcome and \"D\" for dichotomous outcome.")
  }
  
  
  #-------------------------------------
  # check dimensions
  if(nrow(E)!= nrow(Y)) stop("Dimensions of E and Y don't match.")
  if(nrow(Z)!= nrow(Y)) stop("Dimensions of Z and Y don't match.")
  if(is.null(X)==FALSE){
    if(nrow(X)!= nrow(Y)) stop("Dimensions of X and Y don't match.")
    if(class(X)!= "matrix") stop("X is not a matrix.")
  }
  if(class(Z)!= "matrix") stop("Z is not a matrix.")
  if(class(E)!= "matrix") stop("E is not a matrix.")
  if(class(Y)!= "matrix") stop("Y is not a matrix.")
  
  
  #----------------------------------------- added on Oct 14, 2015 (start)
  if(is_check_genotype==FALSE & is_dosage== FALSE){
    if(sum(is.na(Z))!= 0) stop("Z cannot have any missing values if is_check_genotype = FALSE and is_dosage=FALSE.")
  }
  #----------------------------------------- added on Oct 14, 2015 (end)
  
  if(is.null(weights.Z)==FALSE){
    if(ncol(Z)!= length(weights.Z)) stop("Dimensions of Z and weights.Z don't match.")
    if(sum(weights.Z<0)!=0) stop("weights.Z have to be non-negative and cannot be missing.")
  }
  
  if(is.null(weights.V)==FALSE){
    if(ncol(Z)!= length(weights.V)) stop("dimensions of Z and weights.V don't match")
    if(sum(weights.V<0)!=0) stop("weights.V have to be non-negative and cannot be missing")
  }
  
  
  #-------------------------------------
  # check other arguments
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored, beta.weights are also ignored for Z! To use weights as weights.Z or beta.weights for Z, set scale.Z=F")
  
  
  #-------------------------------------
  # check for missing values in Y, E, X
  if(is.null(X)==FALSE){
    if(sum(is.na(X))!= 0) stop("X cannot have any missing values.")
  }
  
  if(sum(is.na(E))!= 0) stop("E cannot have any missing values.")
  if(sum(is.na(Y))!= 0) stop("Y cannot have any missing values.")
  
  
  #-------------------------------------
  # check that X doesn't have intercept
  
  if(is.null(X)==FALSE){
    
    if(ncol(X)==1){
      if(checkpolymorphic(X)==FALSE) stop("X should not include intercept and must have more than one level.")	
    }else{
      if(sum(apply(X, 2, checkpolymorphic))!= ncol(X)) stop("X should not include intercept and must have more than one level.")
      
    }
  }
  
  
  #-------------------------------------
  # check that E has more than one levels
  if(ncol(E)==1){
    if(checkpolymorphic(E)==FALSE) stop("E must have more than one level.")	
  }else{
    if(sum(apply(E, 2, checkpolymorphic))!= ncol(E)) stop("E must have more than one level.")
    
  }
  
  
  #-------------------------------------
  # check Z and impute
  if(is_dosage ==TRUE){
    impute.method="fixed"
  }
  
  if(is_check_genotype==TRUE | is_dosage==TRUE){
    Z.out <- iSKAT_MAIN_Check_Z_and_V(Z=Z, SetID=SetID, weights.Z=weights.Z, weights.V=weights.V, impute.method=impute.method,  										missing_cutoff=missing_cutoff, MAF_cutoff=MAF_cutoff)
    Z <- as.matrix(Z.out$Z.test)
    weights.Z <- as.vector(Z.out$weights.Z.test)
    weights.V <- as.vector(Z.out$weights.V.test)
    include.nonsingleton.Z <- as.vector(Z.out$include.nonsingleton.Z)
  }
  
  if(is.null(Z)==TRUE){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has no SNPs." )
    } else {
      msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=0, n.GE.test=NA))
  }
  
  if(ncol(Z)<2){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has fewer than 2 SNPs. Do not need penalized procedures." )
    } else {
      msg <- sprintf("In %s, the Z matrix has fewer than 2 SNPs. Do not need penalized procedures.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=NA))
  }
  
  
  #-------------------------------------
  # check and assign weights
  MAF_for_betaweights <- colMeans(as.matrix(Z), na.rm=T)/2
  MAF_for_betaweights <- pmin(MAF_for_betaweights, 1-MAF_for_betaweights)
  
  if(is.null(weights.beta)==FALSE&&is.null(weights.Z)==TRUE) weights.Z <- Beta.Weights(MAF_for_betaweights, weights.beta)
  if(is.null(weights.beta)==FALSE&&is.null(weights.V)==TRUE) weights.V <- Beta.Weights(MAF_for_betaweights, weights.beta)
  
  
  #-------------------------------------
  # check V 
  if(is.null(include.nonsingleton.Z)==TRUE){
    
    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix is empty.")
    } else {
      msg <- sprintf("In %s, the GxE matrix is empty.", SetID)
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))
    
  }
  
  Ztemp <- NULL
  Ztemp <- as.matrix(Z[,include.nonsingleton.Z])
  if(is.null(weights.V)==FALSE){
    weights.V.temp <- weights.V[include.nonsingleton.Z]
    Ztemp <- t(t(Ztemp) * (weights.V.temp))
  } 
  
  
  if(ncol(E)==1){
    V <- as.matrix(drop(E)*Ztemp)
  }else{
    V <- NULL
    for(hhh in 1:ncol(E)){
      V <- as.matrix(cbind(V, drop(E[,hhh])*Ztemp))
    }
  }
  
  Idx.V.check <- apply(V, 2, checkpolymorphic)&(apply(abs(cor(V,Z))>0.9999999999,1,sum,na.rm=T)==0)
  
  if(sum(Idx.V.check)==0){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix is empty. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix is empty.", SetID )
    }
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))
  }
  
  
  if(sum(Idx.V.check)==1){ 
    
    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix has fewer than 2 SNPs. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix has fewer than 2 SNPs.", SetID )
    }
    
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=1))
  }
  
  V <- as.matrix(V[,Idx.V.check])
  
  
  #-------------------------------------
  # Run iSKAT
  Xtilde <- as.matrix(cbind(X,E))
  if(out_type == "C"){
    iSKAT.out <- iSKAT.linear(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                              lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL, r.corr=r.corr)
  } 
  
  if(out_type == "D"){	   		
    iSKAT.out <- iSKAT.logistic(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL, r.corr=r.corr)
  } 	   		   
  
  return(list(pvalue=iSKAT.out$pvalue, param=iSKAT.out$param, 
              lambda=iSKAT.out$lambda, n.G.test=drop(ncol(Z)), n.GE.test=drop(ncol(V))))
}









#
#	Check the Z, and do imputation

#
iSKAT_MAIN_Check_Z_and_V <- function(Z, SetID, weights.Z=NULL, weights.V=NULL, impute.method,  missing_cutoff, MAF_cutoff){
  
  # check.Z.error = 0 : no snps removed, but some snps possibly imputed
  # check.Z.error = 1 : all snps removed, returns NULL matrix for Z
  # check.Z.error = 2 : some snps removed, remainder snps may have been imputed
  
  check.Z.error <- 0
  n <- nrow(Z)
  ##############################################
  # Recode Missing to NA
  
  IDX_MISS <- union(which(is.na(Z)), which(Z == 9))
  if(length(IDX_MISS) > 0){
    Z[IDX_MISS] <- NA
  } 
  
  ###################################################
  # Check missing rates and exclude any SNPs with missing rate > missing_cutoff
  # Also exclude non-polymorphic SNPs
  ########################################## modification #1
  m <- ncol(Z)
  ID_INCLUDE_SNP <- NULL
  for(i in 1:m){
    missing.ratio <- length(which(is.na(Z[,i])))/n
    sd1 <- sd(Z[,i], na.rm=TRUE)
    maf1 <- min(mean(Z[,i], na.rm=TRUE)/2, 1- mean(Z[,i], na.rm=TRUE)/2)
    if(is.na(sd1)==FALSE){
      if(missing.ratio < missing_cutoff && sd1 > 0 && maf1<MAF_cutoff){
        ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
      }
    }
  }
  
  
  if(length(ID_INCLUDE_SNP) == 0){
    
    if(is.null(SetID)){
      msg <- sprintf("ALL SNPs have either high missing rates or no-variation or exceed MAF_cutoff. ")
    } else {
      msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation or exceed MAF_cutoff. ",SetID )
    }
    warning(msg, call.=FALSE)
    
    re <- list(Z.test=NULL, weights.Z.test=NULL, weights.V.test=NULL, check.Z.error=1, include.nonsingleton.Z=NULL) 
    
    
  } else{
    
    if(m - length(ID_INCLUDE_SNP) > 0 ){
      
      if(is.null(SetID)){
        msg <- sprintf("%d SNPs with either high missing rates or no-variation or exceed MAF_cutoff are excluded!", m - length(ID_INCLUDE_SNP))
      } else {
        msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation or exceed MAF_cutoff are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
      }
      
      ########################################## end of modification #1
      
      warning(msg, call.=FALSE)	
      Z <- as.matrix(Z[,ID_INCLUDE_SNP])
      if(is.null(weights.Z)==FALSE) weights.Z <- weights.Z[ID_INCLUDE_SNP]
      if(is.null(weights.V)==FALSE) weights.V <- weights.V[ID_INCLUDE_SNP]
      check.Z.error <- 2
      IDX_MISS <- which(is.na(Z))
    }
    
    ########################################## modification #2
    # Check for non-singleton Z's before imputation
    
    ID_INCLUDE_SNP_GE <- NULL
    
    for(iii in 1:ncol(Z)){
      nonsingleton <- checkvariation2(Z[,iii])
      
      if(nonsingleton==TRUE){
        ID_INCLUDE_SNP_GE <- c(ID_INCLUDE_SNP_GE,iii)
      }
    }
    
    
    ########################################## end of modification #2
    
    ##########################################
    # Missing Imputation
    if(length(IDX_MISS) > 0){
      if(is.null(SetID)){
        msg <- sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
      } else {
        msg <- sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
      }
      
      warning(msg,call.=FALSE)
      Z <- Impute_iSKAT(Z,impute.method)
    } 
    re <- list(Z.test=Z, weights.Z.test=weights.Z, weights.V.test=weights.V, check.Z.error=check.Z.error, 
               include.nonsingleton.Z=ID_INCLUDE_SNP_GE)
    
  }
  
  return(re)
}

# Jan 29, 2014
# in v1.0 of iSKAT package there are two sets of SSD functions, one for GESAT(), one for iSKAT()

# modified from SKAT V0.78 SKAT.SSD.OneSet()
GESAT.SSD.OneSet <- function(SSD.INFO, SetID, ...){
  
  id1 <- which(SSD.INFO$SetInfo$SetID == SetID)
  if(length(id1) == 0){
    MSG <- sprintf("Error: cannot find set id [%s] from SSD!", SetID)
    stop(MSG)
  }	
  Set_Index <- SSD.INFO$SetInfo$SetIndex[id1]
  
  Z <- Get_Genotypes_SSD(SSD.INFO, Set_Index)
  re <- GESAT(Z, ...)
  
  return(re)
}


# modified from SKAT V0.78 SKAT.SSD.OneSet_SetIndex()
GESAT.SSD.OneSet_SetIndex <- function(SSD.INFO, SetIndex, ...){
  
  id1 <- which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG <- sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }	
  SetID <- SSD.INFO$SetInfo$SetID[id1]
  
  
  Z <- Get_Genotypes_SSD(SSD.INFO, SetIndex)
  re <- GESAT(Z, ...)
  return(re)
}


# modified from SKAT V0.78 SKAT.SSD.All()
# how GESAT.SSD.ALL handles errors differently from SKAT.SSD.All()
GESAT.SSD.All <- function(SSD.INFO, ...){
  
  N.Set <- SSD.INFO$nSets
  OUT.Pvalue <- rep(NA, N.Set)
  OUT.n.G.test <- rep(NA, N.Set)
  OUT.n.GE.test <- rep(NA, N.Set)
  OUT.Error <- rep(-1, N.Set)
  OUT.lambda <- rep(NA, N.Set)
  OUT.converge <- rep(NA, N.Set)
  
  
  for(i in 1:N.Set){
    
    Is.Error <- TRUE
    try1 <- try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
    
    if(class(try1) != "try-error"){
      Z <- try1
      Is.Error <- FALSE
    } else {
      err.msg <- geterrmessage()
      msg <- sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
    }
    
    
    if(!Is.Error){
      Is.Error <- TRUE
      try2 <- try(GESAT(Z, ...), silent = TRUE)
      
      if(class(try2) != "try-error"){
        re <- try2
        Is.Error <- FALSE
      } else {
        
        err.msg <- geterrmessage()
        msg <- sprintf("Error to run GESAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
        warning(msg,call.=FALSE)
      }
    }
    
    
    if(!Is.Error){
      
      OUT.Pvalue[i] <- re$pvalue
      OUT.converge[i] <- re$Is_converge
      OUT.lambda[i] <- re$lambda
      OUT.n.G.test[i] <- re$n.G.test
      OUT.n.GE.test[i] <- re$n.GE.test
      OUT.Error[i] <- 0
      
    }else{
      
      OUT.Pvalue[i] <- NA
      OUT.converge[i] <- NA
      OUT.lambda[i] <- NA
      OUT.n.G.test[i] <- NA
      OUT.n.GE.test[i] <- NA
      OUT.Error[i] <- 1
      
    }
  }
  
  
  out.tbl <- data.frame(SetID=SSD.INFO$SetInfo$SetID, pvalue=OUT.Pvalue, Is_converge=OUT.converge, lambda=OUT.lambda, 
                        n.G.test=OUT.n.G.test, n.GE.test=OUT.n.GE.test, Error=OUT.Error)
  
  
  
  return(out.tbl)	
}

# Jan 29, 2014
# in v1.0 of iSKAT package there are two sets of SSD functions, one for GESAT(), one for iSKAT()

# modified from SKAT V0.78 SKAT.SSD.OneSet()
iSKAT.SSD.OneSet <- function(SSD.INFO, SetID, ...){
  
  id1 <- which(SSD.INFO$SetInfo$SetID == SetID)
  if(length(id1) == 0){
    MSG <- sprintf("Error: cannot find set id [%s] from SSD!", SetID)
    stop(MSG)
  }	
  Set_Index <- SSD.INFO$SetInfo$SetIndex[id1]
  
  Z <- Get_Genotypes_SSD(SSD.INFO, Set_Index)
  re <- iSKAT(Z, ...)
  
  return(re)
}


# modified from SKAT V0.78 iSKAT.SSD.OneSet_SetIndex()
iSKAT.SSD.OneSet_SetIndex <- function(SSD.INFO, SetIndex, ...){
  
  id1 <- which(SSD.INFO$SetInfo$SetIndex == SetIndex)
  if(length(id1) == 0){
    MSG <- sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
    stop(MSG)
  }	
  SetID <- SSD.INFO$SetInfo$SetID[id1]
  
  
  Z <- Get_Genotypes_SSD(SSD.INFO, SetIndex)
  re <- iSKAT(Z, ...)
  return(re)
}


# modified from SKAT V0.78 SKAT.SSD.All()
# how iSKAT.SSD.ALL handles errors differently from SKAT.SSD.All()
iSKAT.SSD.All <- function(SSD.INFO, ...){
  
  N.Set <- SSD.INFO$nSets
  OUT.Pvalue <- rep(NA, N.Set)
  OUT.n.G.test <- rep(NA, N.Set)
  OUT.n.GE.test <- rep(NA, N.Set)
  OUT.Error <- rep(-1, N.Set)
  OUT.lambda <- rep(NA, N.Set)
  OUT.converge <- rep(NA, N.Set)
  
  
  for(i in 1:N.Set){
    
    Is.Error <- TRUE
    try1 <- try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
    
    if(class(try1) != "try-error"){
      Z <- try1
      Is.Error <- FALSE
    } else {
      err.msg <- geterrmessage()
      msg <- sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
    }
    
    
    if(!Is.Error){
      Is.Error <- TRUE
      try2 <- try(iSKAT(Z, ...), silent = TRUE)
      
      if(class(try2) != "try-error"){
        re <- try2
        Is.Error <- FALSE
      } else {
        
        err.msg <- geterrmessage()
        msg <- sprintf("Error to run iSKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
        warning(msg,call.=FALSE)
      }
    }
    
    
    if(!Is.Error){
      
      OUT.Pvalue[i] <- re$pvalue
      OUT.converge[i] <- re$param$rho_est[1]
      OUT.lambda[i] <- re$lambda
      OUT.n.G.test[i] <- re$n.G.test
      OUT.n.GE.test[i] <- re$n.GE.test
      OUT.Error[i] <- 0
      
    }else{
      
      OUT.Pvalue[i] <- NA
      OUT.converge[i] <- NA
      OUT.lambda[i] <- NA
      OUT.n.G.test[i] <- NA
      OUT.n.GE.test[i] <- NA
      OUT.Error[i] <- 1
      
    }
  }
  
  
  out.tbl <- data.frame(SetID=SSD.INFO$SetInfo$SetID, pvalue=OUT.Pvalue, rho_est=OUT.converge, lambda=OUT.lambda, 
                        n.G.test=OUT.n.G.test, n.GE.test=OUT.n.GE.test, Error=OUT.Error)
  
  
  
  return(out.tbl)	
}






