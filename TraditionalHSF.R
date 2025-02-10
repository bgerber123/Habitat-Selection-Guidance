##   #TOC {
##     max-width: fit-content;
##     white-space: nowrap;
##   }
## 
##   div:has(> #TOC) {
##     display: flex;
##     flex-direction: row-reverse;
## }

## ----setup, include=FALSE-----------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----packages, results='hide',message=FALSE, class.source = "fold-hide"-------------------
#Load packages
  library(ResourceSelection)
  library(glmmTMB)
  library(plotrix)
  library(broom.mixed)
  library(knitr)
  library(amt)
  library(remotes)
  library(ggplot2)
  library(raster)

# A github repository to install
#  remotes::install_github("Pakillo/grateful")
  library(grateful)

# Source power analysis function
  source("./functions/sample.size.used.locs.r")

# Source bootstrapping function
  source("./functions/mean_ci_boot.r")

# Load spatial covariate data (used for power analysis)
  load("./data/Covs")


## ----setup.parameters, class.source = "fold-hide"-----------------------------------------
# Number of Sampled individuals (e.g., tracked via GPS telemetry)
  n.indiv = 20

# Define the true population-level coefficients
# and use these to simulate individual-level coefficients

  #Intercept - Make them all the same, except the last individual
  # This difference is used below.
    beta0 = c(rep(1.5,(n.indiv-1)),0.2)
    
  # Selection against at population-level 
  # Low variation among individuals 
    beta1.mu = -1
    beta1.sd = 0.2
    
    set.seed(543543)
    beta1 = rnorm(n.indiv, 
                  mean = beta1.mu,
                  sd = beta1.sd
                  )
  # No selection at population-level
  # Wide variation among individuals  
    beta2.mu = 0
    beta2.sd = 1
    
    set.seed(5435431)
    beta2 = rnorm(n.indiv, 
                  mean = beta2.mu,
                  sd = beta2.sd
                  )
  
  # Selection for this feature at population-level 
  # Low variation among individuals
    beta3.mu = 1
    beta3.sd = 0.2
    
    set.seed(1543543)
    beta3 = rnorm(n.indiv, 
                  mean = beta3.mu,
                  sd = beta3.sd
                  )

# Combine coefficients and plot these values
  betas = data.frame(b0 = beta0, 
                     b1 = beta1,
                     b2 = beta2,
                     b3 = beta3
                     )


## ----visualize.true.coefficients, fig.height=8,fig.width=6, class.source = "fold-hide"----
  par(mfrow=c(3,1))
  hist(betas[,2],xlab=bquote(beta[1]),xlim=c(-3,3),main="True Individual Coefficient Values",breaks=10,freq=FALSE)
  abline(v=beta1.mu,lwd=2,col=2)
  curve(dnorm(x,beta1.mu,beta1.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))
  
  hist(betas[,3],xlab=bquote(beta[2]),xlim=c(-3,3),main="",breaks=20,freq=FALSE)
  abline(v=beta2.mu,lwd=2,col=2)
  curve(dnorm(x,beta2.mu,beta2.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))

  hist(betas[,4],xlab=bquote(beta[3]),xlim=c(-3,3),main="",breaks=20,freq=FALSE)
  abline(v=beta3.mu,lwd=2,col=2)
  curve(dnorm(x,beta3.mu,beta3.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))



## ----simualte.data, cache=TRUE, class.source = "fold-hide"--------------------------------
# How many used samples per individual   
  n.used = 500

# Ratio of used:available
  m = 10 
  
#Create available samples  
  n=n.used*m 

# Create spatial covariate variables
# dist.dev = Euclidean distance to nearest development
# forest = percent forest cover at a relevant spatial scale
# shrub = shrub cover at a relevant spatial scale
  
# These are continuous variables that are scaled to a mean of 0 and stdev of 1
  set.seed(51234) 
  x = data.frame(dist.dev=rnorm(n),forest=runif(n),shrub=runif(n)) 

# Create one dataset per each individual with equal used samples per individual
# Note that in the simulation function (simulateUsedAvail) we are specifyng 
# a link of 'log' indicating we will be fitting a model to estimate 
# RELATIVE habitat selection rather than absolute selection.


  sims = apply(betas,1,FUN=function(b){
               ResourceSelection::simulateUsedAvail(data=x,
                                                    parms=b,
                                                    n.used=n.used,
                                                    m=m,
                                                    link="log"
                                                    ) 
    })



## ----data.examine-------------------------------------------------------------------------
# Check the dimensions of the data  
# This should have the same length as n.indiv  
  length(sims)

# Look at one individual's dataset
  dim(sims[[1]])
  head(sims[[1]])
  table(sims[[1]]$status)


## ----organize. simualted.data-------------------------------------------------------------
# Combine data into a single data.frame
  sims2 = do.call("rbind", sims)
  head(sims2)
  dim(sims2)

#Create ID vector for individual's data
  ID = rep(1:n.indiv, each = nrow(sims[[1]]))
  sims2$ID = ID
  head(sims2)
  dim(sims2)
  
# Create a vector to weight the 0's by 1000 and the 1's by 1
# This is discussed below  (see manuscript section 10. The model being approximated)
  sims2$weight = 1000^(1-sims2$status)
  head(sims2)


## ----first.fit.inference------------------------------------------------------------------
# We will use data from individual #20
  indiv.data20 = sims[[20]]

# Let's look at the data to remind us what it is
  head(indiv.data20)


## ----first.fit.inference2-----------------------------------------------------------------
# We have related our response variable of the used and available sample (status) to our covariates. This is an additive model, as opposed to having interactions.  
  model1 = glmmTMB::glmmTMB(status ~ dist.dev + forest + shrub, 
                            family = binomial(), 
                            data = indiv.data20 
                            )

# Look at estimated coefficients
  summary(model1)


## ----first.fit.inference.glm--------------------------------------------------------------
  model1.glm = stats::glm(status ~ dist.dev + forest + shrub, 
                          family = binomial(), 
                          data = indiv.data20 
                          )
  summary(model1.glm)


## ----amt----------------------------------------------------------------------------------
amt::fit_rsf


## ----RS-----------------------------------------------------------------------------------
# Get estimates without the intercept
  coef = fixef(model1)[[1]][-1]

  RS = exp(-2*coef[1] + 2*coef[2] + 2*coef[3]) / 
       exp(2*coef[1] + 2*coef[2] + -2*coef[3])
  RS


## ----RSS----------------------------------------------------------------------------------
# Coefficient for dist.dev
  exp(coef[2])


## ----predict------------------------------------------------------------------------------
# Force the intercept to be zero without any error
  model1$fit$par[1] = 0
  model1$sdr$par.fixed[1] = 0
  model1$sdr$cov.fixed[1,] = 0
  model1$sdr$cov.fixed[,1] = 0

# Create the dataframe for the combinations of variables for prediction  
  newdata=data.frame(dist.dev = seq(-2,2, by=0.1),
                     forest = 0,
                     shrub = 0
                     )
# Next, predict 
  preds = predict(model1,
                  type = "link", # NOTE: this is where we tell the 'predict' function to use the same link function as the habitat selection function we fit above
                  newdata = newdata, 
                  se.fit = TRUE
                  )
# Use the Normal Approximation to get confidence intervals (other approachces are possible too!)  
  preds$LCL = preds$fit-1.96*preds$se.fit
  preds$UCL = preds$fit+1.96*preds$se.fit

  preds = data.frame(preds)  
  
# Exponentiate predicted values
  preds.exp = exp(preds[,-2])  


## ----predict.plot, class.source = "fold-hide"---------------------------------------------
# Plot
  plot(newdata$dist.dev,preds.exp$fit,type="l",lwd=3,xlab="Distance to Development (Scaled)",
       ylab="Relative Intensity of Selection")
  lines(newdata$dist.dev,preds.exp$LCL,lwd=3,col=2,lty=4)
  lines(newdata$dist.dev,preds.exp$UCL,lwd=3,col=2,lty=4)
  abline(h=1, lwd=4, col=4)
  text(1.5,1.4, "Selection")
  text(1.5,0.6, "Avoidance")


## ----assumption1--------------------------------------------------------------------------
# Create a new (small) set of available locations 
  n.avail.new = 50
  set.seed(215464)
  avail.new = data.frame(status = rep(0,n.avail.new),
                         dist.dev = rnorm(n.avail.new,-0.5,3),
                         forest = rnorm(n.avail.new,0.5,3),
                         shrub = rnorm(n.avail.new,1.5,3)
                        )

# Append the new data to the original data
  indiv.data20.appended = rbind(indiv.data20,avail.new)

# Fit the same model with appended data
  model1.appended = glmmTMB(status ~ dist.dev + forest + shrub, 
                            family = binomial(), 
                            data = indiv.data20.appended 
                            )
# Re-fit the original data
    model1 = glmmTMB::glmmTMB(status ~ dist.dev + forest + shrub, 
                              family = binomial(), 
                              data = indiv.data20 
                            )
  
# Compare coefficients using the appended data and the original data
  coef.df = rbind(fixef(model1)[[1]][-1],
                  fixef(model1.appended)[[1]][-1]
                  )
  colnames(coef.df) = c("dist.dev","forest","shrub")
  rownames(coef.df) = c("Original Data","Data with Additional Locations")
  knitr::kable(coef.df)  


## ----assumption2--------------------------------------------------------------------------

# Lets use indidual 7's data
  indiv.data7 = sims[[7]]

# Find the used locations
  index1 = which(indiv.data7$status==1)

# Replicate the used locations 
  indiv.data7.replicated = rbind(indiv.data7,
                                 indiv.data7[index1,]
                                )

  fit.replicated = glmmTMB(status ~ dist.dev + forest + shrub  , 
                           family = binomial(), 
                           data = indiv.data7.replicated, 
                           )

# Fit the original data
    model1 = glmmTMB::glmmTMB(status ~ dist.dev+forest+shrub, 
                              family = binomial(), 
                              data = indiv.data7 
                            )
  
# Compare coefficients when using the replicated data and the original data
  coef.df[1,] = fixef(model1)[[1]][-1]
  coef.df[2,] = fixef(fit.replicated)[[1]][-1]
  rownames(coef.df)[2] = "Non-Independent Data"
  knitr::kable(coef.df)


## ----assumption2.2------------------------------------------------------------------------
# Compare measures of uncertainty
  summary(model1) # independent data
  summary(fit.replicated) # replicated (dependent) data


## ----sensitiivity, cache=TRUE-------------------------------------------------------------
# Grab individual one's data
  indiv.dat1 = sims[[1]]

# index the used (1) and available (0) data
  index.0 = which(indiv.dat1==0)
  index.1 = which(indiv.dat1==1)

# Create a sequence of the number of available locations to consider
  n.available = seq(10,length(index.0),by=100)

# Create storage objects  
  coef.save = NULL

# Loop through and modify the data and fit each model with increasing 
# number of available samples
  for(i in 1:length(n.available)){
    dat.temp = rbind(indiv.dat1[index.1,],
                     indiv.dat1[index.0[1:n.available[i]],]
                     )
  
  #Create the weighting vector (a 1 for each 1 and a 1000 for each 0)
  dat.temp$weight = 1000^(1-dat.temp$status)
  
  model.fit =  glmmTMB(status ~ dist.dev + forest + shrub, 
                        family = binomial(), 
                        data =  dat.temp,
                        weight = weight # add the weighting vector
                        )
  

  coef.save = rbind(coef.save,
                    fixef(model.fit)[[1]]
                  )
}

coef.save = data.frame(n.available,coef.save)
colnames(coef.save) = c("N.Avail","Intercept","beta1","beta2","beta3")


## ----plotting.sensitivity, fig.height=8, fig.width=6, class.source = "fold-hide"----------
par(mfrow=c(2,1))
plot(coef.save$N.Avail, coef.save$beta1,lwd=3,type="l",col=2,ylim=c(-1.5,1.5),
     main="Slope Coefficients",
     xlab="Number of Available Samples",ylab="Coefficient Estimate")
lines(coef.save$N.Avail, coef.save$beta2,lwd=3,col=2)
lines(coef.save$N.Avail, coef.save$beta3,lwd=3,col=2)
plot(coef.save$N.Avail, coef.save$Intercept,lwd=3,type="l",col=1,main="Intercept",
     xlab="Number of Available Samples",ylab="Coefficient Estimate")


## ----sensitiivity2, cache=TRUE,warnings=FALSE---------------------------------------------
# Size of available samples per used locaiton
  n.avail2 = c(10,50,100,1000,2000,4000)
  
# Number of sample replications
 n.sample = 20

#Save coefficients
  coef.save=NULL
  
#Loop across available sample sizes  
  for(i in 1:length(n.avail2)){
  
  #re-sample each available sample size 20 times  
  for (z in 1:n.sample){
  #Loop through each used location and reduce 
  #the number of available samples and do this n.sample times
    ind.dat = NULL

    dat.temp = rbind(indiv.dat1[index.1,],
                     indiv.dat1[sample(index.0,n.avail2[i]),]
                     )
    
  #Create the weighting vector (a 1 for each 1 and a 1000 for each 0)
  dat.temp$weight = 1000^(1-dat.temp$status)
  
  model.fit =  glmmTMB(status ~ dist.dev + forest + shrub, 
                        family = binomial(), 
                        data =  dat.temp,
                        weight = weight # add the weighting vector
                        )
  
  coef.save= rbind(coef.save, 
                   model.fit$fit$par
                   )
  }#end z sample loop
}#end i loop
  


## ----plot.sensitive.org2, class.source = "fold-hide"--------------------------------------
one = data.frame(N.Avail = rep(n.avail2,each = n.sample),
                 Sample = rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,1]
                )
two = data.frame(N.Avail = rep(n.avail2,each = n.sample),
                 Sample = rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,2]
                )
three = data.frame(N.Avail = rep(n.avail2,each = n.sample),
                 Sample = rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,3]
                )
four = data.frame(N.Avail = rep(n.avail2,each = n.sample),
                 Sample = rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,4]
                )

plot.data = rbind(one,two,three,four)
plot.data$Name = c(rep("Intercept",nrow(one)),
                   rep("dist.dev",nrow(one)),
                   rep("forest",nrow(one)),
                   rep("shrub",nrow(one))
                 )

plot.data$N.Avail = as.factor(plot.data$N.Avail)


colnames(plot.data) = c("N.Available.Sample","Sample","Coefficient.Estimate","Name")

ggplot(plot.data, aes(N.Available.Sample, Coefficient.Estimate, fill=factor(Name))) +
  theme_bw()+
  geom_boxplot()


## ----pooling,cache=TRUE-------------------------------------------------------------------
# Pooled model - no consideration of individual variability
  model.pool = glmmTMB(status ~ dist.dev + forest + shrub , 
                       family = binomial(), 
                       data = sims2, 
                    )

# Separate models to each individual
  indiv.fit = lapply(sims, FUN = function(x){
                     glmmTMB(status ~ dist.dev + forest + shrub ,
                             family = binomial(), 
                             data = x, 
                            )
                      })
  
# Look at estimates when pooling the data
  summary(model.pool)


## ----separate-----------------------------------------------------------------------------
summary(indiv.fit[[2]])


## ----separate2----------------------------------------------------------------------------
summary(indiv.fit[[20]])


## ----estimate.pool.separate---------------------------------------------------------------
# Extract pooled estimates with Confidence intervals
pool.effect = broom.mixed::tidy(model.pool, effects = "fixed", conf.int = TRUE)

# Extract separate individual estimates
indiv.coef = lapply(indiv.fit,FUN=function(x){
                      temp=broom.mixed::tidy(x, effects = "fixed", 
                                             conf.int = TRUE
                                             )
                      data.frame(temp[-1,c(4,8,9)])
})

# These are all coefs (1:3) for each individual 1:n.indiv
  estimates = sapply(indiv.coef,"[[",1)
  LCI = sapply(indiv.coef,"[[",2)
  UCI = sapply(indiv.coef,"[[",3)


## ----plot.pool.separate, class.source = "fold-hide"---------------------------------------
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE),
       width=c(2,1)
       )
plotCI(1:n.indiv,y=estimates[1,],
       ui=UCI[1,],
       li=LCI[1,],
       ylab="beta",
       xlab="Individual",
       lwd=2,
       ylim=c(-2.5,4)
)
 legend("topright",lwd=3,col=c(1,2,3),legend=c("dist.dev","forest","shrub"),box.lty=0,
        y.intersp=0.8)
 
plotCI(1:n.indiv,y=estimates[2,],
       ui=UCI[2,],
       li=LCI[2,],
       ylab="beta_1",
       xlab="Individual",
       lwd=2,
       col=2,
       add=TRUE
       )
plotCI(1:n.indiv,y=estimates[3,],
       ui=UCI[3,],
       li=LCI[3,],
       ylab="beta_1",
       xlab="Individual",
       lwd=2,
       col=3,
       add=TRUE
       )
plotCI(c(1,1,1),y=pool.effect$estimate[2:4],
       xlim=c(0.95,1.05),
       ui=pool.effect$conf.high[2:4],
       li=pool.effect$conf.low[2:4],
       ylab="beta",
       xlab="Population",
       lwd=2,
       ylim=c(-2.5,4),
       col=c(1,2,3),
       xaxt='n'
       )


## ----power`-------------------------------------------------------------------------------
S = covs[[1]]
values(S) = scale(values(S))  


## ----power2, fig.height=4, fig.width=6----------------------------------------------------
# Create habitat selection function values
# We will use a slope of -0.2 for our spatial covariate
  beta = matrix(c(1,-0.2),
                 nrow = 2,
                 ncol = 1
                ) 
# Design matrix
  X = cbind(1,values(S))  
  
#Intensity of selection  
  lambda = exp(X%*%beta) 

# Create spatial layer with HSF values  
  S.HSF = S
  values(S.HSF) = lambda

# Plots
  par(mfrow=c(1,2),mar=c(4,4,4,4))
  plot(S, main = "Spatial Covariate")
  plot(S.HSF, main = "Relative Intensity of Selection")


## ----power3-------------------------------------------------------------------------------
# Define Inputs
  alpha = 0.05
  p_val = 0.05


## ----power4-------------------------------------------------------------------------------
N.used = sample.size.used.locs(alpha = alpha, 
                               p_val = p_val,
                               HSF.True = S.HSF,
                               S = S,
                               beta = beta[2]
                               )
N.used$Nalphapbetas


## ----power5-------------------------------------------------------------------------------
  n.combs =  20
  alpha = seq(0.0001,0.3,length.out = n.combs)
  N.used = rep(0,n.combs)
  
  #Loop through   
  for(i in 1:n.combs){
  temp = sample.size.used.locs(alpha=alpha[i],
                                    p_val = p_val,
                                    HSF.True = S.HSF,
                                    S = S,
                                    beta = beta[2])
  N.used[i] = temp$Nalphapbetas
  }


## ----power.plot1, class.source = "fold-hide"----------------------------------------------
  plot(alpha,N.used,col=2,type="l",lwd=3,xlab="Type I Error Rate",ylab="Sample Size of Used Locations Needed")


## ----power7-------------------------------------------------------------------------------
  n.combs =  40
  N.used = rep(0,n.combs)
  beta.combs=seq(-2,2,length.out=n.combs)
  
  for(i in 1:n.combs){
  
    beta = matrix(c(-1,beta.combs[i]),
                   nrow=,ncol=1)  
    X = cbind(1,values(S))  
    lambda = exp(X%*%beta) 
  
    HSF.True = S
    values(HSF.True)= lambda
    temp = sample.size.used.locs(alpha = 0.05,
                                      p_val = 0.05,
                                      HSF.True = HSF.True,
                                      S = S,
                                      beta = beta.combs[i]
                                      )
    N.used[i] = temp$Nalphapbetas
  }


## ----power.plot2, class.source = "fold-hide"----------------------------------------------
  plot(beta.combs,N.used,col=2,type="b",lwd=3,xlab="Coefficient",ylab="Sample Size of Used Locations Needed")


## ----power9, cache=TRUE-------------------------------------------------------------------
  n.combs = 20
  n.reps = 30
  sd.combs = seq(0.1,0.5,length.out = n.combs)
  N.used = N.variance = NULL
  save.hsf = vector("list",20)

# Loop over sd.combs
  for(i in 1:n.combs){
  
    # Replicates
    for(j in 1:n.reps){  
      # Simulate Spatial support of point process
      S = raster(xmn=0,xmx=1.5,ymn=0,ymx=1,res=0.05)
    
      # Random spatial covariate
      phi = 0.7 # spatial decay parameter
      
      n = ncell(S)
      xy = xyFromCell(S, 1:n)
      d = as.matrix(dist(xy))		
      Sigma = exp(-d/phi)
      Sigma = t(chol(Sigma))
      values(S) = Sigma %*% rnorm(n,0,sd=sd.combs[i])
  
      values(S) = scale(values(S))
      # Smooth resource covariate
      S = disaggregate(S,fact=2,method="bilinear")
      
      beta = matrix(c(-1,-0.2),nrow=2,ncol=1)  #coefs
      X = cbind(1,values(S))  # design matrix
      lambda = exp(X%*%beta) #intensity of selection
  
      HSF.True=S
      values(HSF.True)= lambda
      N.temp = sample.size.used.locs(alpha = 0.05,
                                     p_val = 0.05,
                                     HSF.True = HSF.True,
                                     S = S,
                                     beta = beta[2]
                                     )
      
      N.used = c(N.used,N.temp$Nalphapbetas)
      N.variance = c(N.variance, N.temp$variance)
    }
  }


## ----power.plot3, class.source = "fold-hide"----------------------------------------------
  N.used.plot = data.frame(N.used = N.used,
                           N.variance = N.variance,
                           sd = rep(sd.combs,each=n.reps)
                          )
  N.used.plot = N.used.plot[order(N.used.plot$N.variance),]
plot(N.used.plot$N.variance,N.used.plot$N.used,lwd=3,col=2,type="l",xlab="Landscape Complexity",ylab="Sample Size of Used Locations Needed")


## ----boot.setup---------------------------------------------------------------------------
  coef.list = lapply(indiv.fit,FUN=function(x){fixef(x)[[1]]})
  
  coef.df=do.call(rbind.data.frame, coef.list)
  colnames(coef.df)=names(fixef(indiv.fit[[1]])[[1]])

# Remove intercept
  coef.df = coef.df[,-1]

# Add name (could be modified to keep track if animals have an ID)
  coef.df$name = 1:nrow(coef.df)

# Frequency could be greater than one if there were multiple estimates of a single individual, perhaps across years
  coef.df$freq = 1


## ----boot---------------------------------------------------------------------------------
#How many bootstraps to do? More will lead to results with less error
  nboot = 1000

#Which columns have coefficients in coef.df
  col.locs=1:3

  boot=list()
  for(i in 1:nboot){
    boot[[i]]=apply(
                     coef.df[sample(nrow(coef.df), nrow(coef.df), 
                                    replace=TRUE, 
                                    prob=coef.df$freq),
                             col.locs 
                             ],
                     2, 
                     median, 
                     na.rm=TRUE
                     ) 
  }


## ----boot.summary-------------------------------------------------------------------------
  boot.pop.esimates=mean_ci_boot(boot)
  knitr::kable(boot.pop.esimates)


## ----random.intercept, cache=TRUE---------------------------------------------------------
# Setup HSF with logistic regression approximation - random intercept model
# Do not fit the model
   re_int = glmmTMB(status ~ dist.dev + forest + shrub  + (1|ID), 
                    family = binomial(), 
                    data = sims2, 
                    doFit = FALSE
                    )

# Make the intercept variance large and fixed (i.e., do not estimate).
  re_int$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the intercept  
  re_int$mapArg = list(theta=factor(c(NA)))
  
# Now ready to fit the model  
  re_int = glmmTMB:::fitTMB(re_int)


## ----examine.random.intercept.model-------------------------------------------------------
# The fixed components - these are the population level (across individual) effects. 
  fixef(re_int)
  
# The random components- these are the effect differences for each individual from the population mean estimate (referred to as the fixed effect in the 'conditional model' statement)
  ranef(re_int)


## ----examine.random.intercept.model2------------------------------------------------------
# Summarize results  
  summary(re_int)


## ----alt.RE.intercept---------------------------------------------------------------------
# Make ID a factor
  sims2$ID=as.factor(sims2$ID)

# Setup a design matrix with individual ID's
  X = model.matrix(~ID,
                   data = sims2,
                   contrasts.arg =list(ID="contr.sum")
                   )

# Fixed effect intercepts that vary by individual
  fix_int = glmmTMB(status ~ 0 + X + dist.dev + forest + shrub, 
                   family = binomial(), 
                   data = sims2
  )


## ----fix.re.compare1----------------------------------------------------------------------
kable(summary(re_int)[[6]][[1]][2:4,])


## ----fix.re.compare2----------------------------------------------------------------------
kable(summary(fix_int)[[6]][[1]][(n.indiv+1):(n.indiv+1+2),])


## ----fix.re.compare3----------------------------------------------------------------------
int.re = ranef(re_int)[[1]][[1]]

int.fe = c(summary(fix_int)[[6]][[1]][1,1]+
           summary(fix_int)[[6]][[1]][2:20],
           summary(fix_int)[[6]][[1]][1,1]
           )
kable(round(data.frame(individual = 1:20,int.fe=int.fe,int.re=int.re),digits=2))


## ----model2, cache=TRUE-------------------------------------------------------------------
#Fit RSF with intercept and slopes with random effect
  re_int_slopes =  glmmTMB(status ~ dist.dev + forest + shrub  + 
                                     (1|ID) + (0+dist.dev|ID) +
                                     (0+forest|ID) + (0+shrub|ID), 
                            family = binomial(), 
                            data = sims2, 
                            doFit = FALSE, 
                            weights = weight
                            )
# make the intercept large and fixed
  re_int_slopes$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the intercept, but 
# to do so for the other three variables
  re_int_slopes$mapArg = list(theta=factor(c(NA,1:3)))
  
# Now ready to fit the model    
  re_int_slopes = glmmTMB:::fitTMB(re_int_slopes)


## ----RE.results.fixed---------------------------------------------------------------------
  fixef(re_int_slopes)
  broom.mixed::tidy(re_int_slopes, effects = "fixed", conf.int = TRUE)[-1,c(3,4,8,9)]


## ----RE.results.random--------------------------------------------------------------------
  ranef(re_int_slopes)
  broom.mixed::tidy(re_int_slopes, effects = "ran_vals", conf.int = TRUE)[-c(1:20),c(5,6,8,9)]



## ----RE.results---------------------------------------------------------------------------
summary(re_int_slopes)


## ----RE.est.forest.indiv.pop--------------------------------------------------------------
  indiv.coef.popModel = broom.mixed::tidy(re_int_slopes, 
                                          effects = "ran_vals", 
                                          conf.int = TRUE)
  indiv.forest.popModel = data.frame(indiv.coef.popModel[41:60,c(5,6,8,9)])

# Add back mean to get individual estimates  
  indiv.forest.popModel[,2:4] = indiv.forest.popModel[,2:4]+
                                rep(fixef(re_int_slopes)[[1]][3],each=20)
    
# Extract population mean and uncertainty for forest effect
  pop.coef.popModel = data.frame(broom.mixed::tidy(re_int_slopes, 
                                                   effects = "fixed",
                                                   conf.int = TRUE)
                                 )
  pop.coef.popModel=pop.coef.popModel[3,c(4,8,9)]


## ----RE.plot.forest.indiv.pop, class.source = "fold-hide"---------------------------------
#Plot
  plotCI(x = 1:20,
         y = indiv.forest.popModel$estimate,
         ui = indiv.forest.popModel$conf.high,
         li = indiv.forest.popModel$conf.low,
         lwd = 3,col=1,
         xlab="Individual",ylab="Forest Coefficient")
  abline(h=pop.coef.popModel$estimate,lwd=3,col=1,lty=4)
  abline(h=pop.coef.popModel$conf.low,lwd=3,col=2,lty=4)
  abline(h=pop.coef.popModel$conf.high,lwd=3,col=2,lty=4)


## ----power.population1,cache=TRUE---------------------------------------------------------
#Number of individuals
  M = 30
  
# Population-level coefs
  beta_mu = 0.2
  beta_s = 0.5 # note the amount of variation
  beta_s2 = beta_s^2
  
# Individual level coefs
  set.seed(4254)
  beta.ind <- rnorm(M, beta_mu,beta_s)
  
# Need to calculate landscape complexity or Var[R(X_beta)] for each individual, given the covariate and the individual coefficient

N.variance = NULL    
for(i in 1:M){
  
      beta = matrix(c(-1,beta.ind[i]),nrow=2,ncol=1)  #coefs
      X = cbind(1,values(covs[[1]]))  # design matrix
      lambda = exp(X%*%beta) #intensity of selection
  
      HSF.True=covs[[1]]
      values(HSF.True)= lambda
      N.temp = sample.size.used.locs(alpha = 0.05,
                                     p_val = 0.05,
                                     HSF.True = HSF.True,
                                     S = covs[[1]],
                                     beta = beta[2]
                                     )
      
      N.variance = c(N.variance, N.temp$variance)
}  
  
# Assume the same number of locations per individual
  N = rep(1000,M)

# Equation 5 of Street et al. 2021
  sigma <- 1/sqrt(N*N.variance)

# Get from function  
  zalpha = N.temp$zalpha
  zp = N.temp$zp

# We need M to be greater than or equal to this value; (Equation 8; Street et al. 2021):
  
  Eq8 = (beta_s2*(zalpha+zp)^2 + sqrt((beta_s^4)*(zalpha+zp)^4 + 
              4*(beta_mu^2)*(zalpha+zp)^2*sum(sigma))) / (2*beta_mu^2)

paste("M (tracked individuals) should be greater than or equal to ", round(Eq8,digits=0))


## ----power.population2,cache=TRUE, class.source = "fold-hide"-----------------------------
  M = 30
  
  beta_mu = seq(0.1,0.75,length.out=10)
  beta_s = 0.5 
  beta_s2 = beta_s^2

  min.indiv = rep(NA, length(beta_mu))

for(z in 1:length(beta_mu)){    
# Individual level coefs
  set.seed(4254)
  beta.ind <- rnorm(M, beta_mu[z],beta_s)
  
  N.variance = NULL    
  for(i in 1:M){
    
        beta = matrix(c(-1,beta.ind[i]),nrow=2,ncol=1)  
        X = cbind(1,values(covs[[1]]))  
        lambda = exp(X%*%beta) 
    
        HSF.True=covs[[1]]
        values(HSF.True)= lambda
        N.temp = sample.size.used.locs(alpha = 0.05,
                                       p_val = 0.05,
                                       HSF.True = HSF.True,
                                       S = covs[[1]],
                                       beta = beta[2]
                                       )
        
        N.variance = c(N.variance, N.temp$variance)
  }  #end i loop
    
    N = rep(1000,M)
    sigma <- 1/sqrt(N*N.variance)
    zalpha = N.temp$zalpha
    zp = N.temp$zp
  
    
    Eq8 = (beta_s2*(zalpha+zp)^2 + sqrt((beta_s^4)*(zalpha+zp)^4 + 
                4*(beta_mu[z]^2)*(zalpha+zp)^2*sum(sigma))) / (2*beta_mu[z]^2)
    
    min.indiv[z] = Eq8
  
}#end z loop


## ----power.pop.plot-----------------------------------------------------------------------
plot(beta_mu,min.indiv,ylab="Minimum Number of Individuals Tracked",lwd=3,col=4,type="b")


## ----power.population3, class.source = "fold-hide"----------------------------------------
  M = 30
  
  beta_mu = 0.2
  beta_s = seq(0.01,0.5,length.out=10) 
  beta_s2 = beta_s^2

  min.indiv = rep(NA, length(beta_mu))

for(z in 1:length(beta_s)){    
  set.seed(4254)
  beta.ind <- rnorm(M, beta_mu,beta_s)
  
  N.variance = NULL    
  for(i in 1:M){
    
        beta = matrix(c(-1,beta.ind[i]),nrow=2,ncol=1)  
        X = cbind(1,values(covs[[1]]))  
        lambda = exp(X%*%beta) 
    
        HSF.True=covs[[1]]
        values(HSF.True)= lambda
        N.temp = sample.size.used.locs(alpha = 0.05,
                                       p_val = 0.05,
                                       HSF.True = HSF.True,
                                       S = covs[[1]],
                                       beta = beta[2]
                                       )
        
        N.variance = c(N.variance, N.temp$variance)
  }  #end i loop
    
    N = rep(1000,M)
    sigma <- 1/sqrt(N*N.variance)
    zalpha = N.temp$zalpha
    zp = N.temp$zp
  
    
    Eq8 = (beta_s2[z]*(zalpha+zp)^2 + sqrt((beta_s[z]^4)*(zalpha+zp)^4 + 
                4*(beta_mu^2)*(zalpha+zp)^2*sum(sigma))) / (2*beta_mu^2)
    
    min.indiv[z] = Eq8
  
}#end z loop


## ----power.pop.plot2----------------------------------------------------------------------
plot(beta_s2,min.indiv,ylab="Minimum Number of Individuals Tracked",lwd=3,col=4,type="b")


## ----behav.sim, class.source = "fold-hide"------------------------------------------------

# Setup simulation input where there is more data from the second behavior
 n.used.behav1 = 200
 n.used.behav2 = 200

# Ratio of used:available
  m = 10 
  
#Create available samples  
  n1 = n.used.behav1*m
  n2 = n.used.behav2*m 

# Consider a single categorical variable of forest and not-forest
  set.seed(51234) 
  x1 = data.frame(forest = rbinom(n1,1,0.4))
  x2 = data.frame(forest = rbinom(n2,1,0.4)) 

# Intercept and effect of forest
  parms.behav1 = c(0.5,2) # the 2 indicates selection
  parms.behav2 = c(0.5,-2) # the -0 indicates avoidance
  
# simulate data  
  sim.behav1 = ResourceSelection::simulateUsedAvail(data = x1,
                                                   parms = parms.behav1,
                                                   n.used = n.used.behav1,
                                                   m = m,
                                                   link = "log"
                                                   ) 
  sim.behav2 = ResourceSelection::simulateUsedAvail(data = x2,
                                                   parms = parms.behav2,
                                                   n.used = n.used.behav2,
                                                   m = m,
                                                   link = "log"
                                                   ) 
# Combine the data  
  data.ignore = rbind(sim.behav1,
                      sim.behav2
                      )  


## ----behav.fit----------------------------------------------------------------------------
  model.ignore.behav = glmmTMB::glmmTMB(status ~ forest, 
                                        family = binomial(), 
                                        data = data.ignore 
                                        )
  summary(model.ignore.behav)


## ----packages.cite, eval=TRUE, echo=FALSE-------------------------------------------------
  pkgs = cite_packages(output = "table", out.dir = ".",out.format ="pdf")
  knitr::kable(pkgs)

