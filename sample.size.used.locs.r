sample.size.used.locs = function(alpha = alpha,
                                 p_val = p_val,
                                 HSF.True = HSF.True,
                                 S = S,
                                 beta = beta
                                 ){
  
# Derive statistical 'significance' inputs
#  zalpha <- qnorm(1-alpha)
#  zp <- qnorm(1-p_val/2)
  zalpha <- qnorm(1-alpha/2)
  zp <- qnorm(1-p_val)
  
  #Next sets of code taken from Street et al. 2021
  # Calculate exp(beta*R(x)) for each x, and the integral 
  expbetaRx <- HSF.True
  intexpbetaRx<-sum(values(expbetaRx))/(ncol(S)*nrow(S))
  
  # Calculate R(x)*exp(beta*R(x)) for each x, and the integral 
  RexpbetaRx<-values(S)*expbetaRx
  intRexpbetaRx<-sum(values(RexpbetaRx))/(ncol(S)*nrow(S))
  
  # Calculate R^2(x)*exp(beta*R(x)) for each x, and the integral
  RsqexpbetaRx<-values(S)*RexpbetaRx
  intRsqexpbetaRx<-sum(values(RsqexpbetaRx))/(ncol(S)*nrow(S))
  
  # Calculate the variance of R(X_beta).  This is Eqn 4 in the manuscript
  variance<-intRsqexpbetaRx/intexpbetaRx - (intRexpbetaRx/intexpbetaRx)^2
  variance
  
  # Calculate N_{alpha,p}(beta).  This is Equation (3) in the Main Text
  Nalphapbetas = ((zalpha+zp)^2/variance)*(beta^(-2))
  
  #This is the sample size of locations required
  list(Nalphapbetas=Nalphapbetas, variance=variance,
       zalpha=zalpha,zp=zp)
  
}
