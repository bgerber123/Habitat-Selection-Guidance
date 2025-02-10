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
  library(geoR)
  library(circular)
  library(glmmTMB)
  library(raster)
  library(sp)
  library(survival)
  library(Rfast)
  library(remotes)
  library(plotrix)
  library(ggplot2)
  library(knitr)

# Install github repository  
#    remotes::install_github("Pakillo/grateful")

  library(grateful)

# Source simulation function
  source("./functions/sim.ind.movement.hsf.r")

# Source bootstrapping function
  source("./functions/mean_ci_boot.r")

# Load spatial covariates stored in a save R object
  load("./data/Covs")


## ----setup.parameters, class.source = "fold-hide"-----------------------------------------
# Number of Sampled individuals (e.g., tracked via GPS telemetry)
  n.indiv = 30

# Define the true population-level coefficients
# and use these to simulate individual-level coefficients

# Selection against at population-level 
# Low variation among individuals 
  beta1.mu = -1
  beta1.sd = 0.2
  
  set.seed(543543)
  beta1=rnorm(n.indiv, 
              mean = beta1.mu,
              sd = beta1.sd
  )
# No selection at population-level
# Wide variation among individuals  
    beta2.mu = 0
    beta2.sd = 1
  
    set.seed(5435431)
    beta2=rnorm(n.indiv, 
                mean = beta2.mu,
                sd = beta2.sd
    )
  
  # Selection for this feature at population-level 
  # Low variation among individuals
    beta3.mu = 1
    beta3.sd = 0.2
  
    set.seed(1543543)
    beta3=rnorm(n.indiv, 
                mean = beta3.mu,
                sd = beta3.sd
    )
    
# Combine coefficients and plot these values
  betas=data.frame(b1 = beta1,
                   b2 = beta2,
                   b3 = beta3)    


## ----visualize.true.coefficients, fig.height=6,fig.width=8, class.source = "fold-hide"----
  par(mfrow=c(3,1))
  hist(betas[,1],xlab=bquote(beta[1]),xlim=c(-3,3),main="True Individual Coefficient Values",breaks=10,freq=FALSE)
  abline(v=beta1.mu,lwd=2,col=2)
  curve(dnorm(x,beta1.mu,beta1.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))
  
  hist(betas[,2],xlab=bquote(beta[2]),xlim=c(-3,3),main="",breaks=20,freq=FALSE)
  abline(v=beta2.mu,lwd=2,col=2)
  curve(dnorm(x,beta2.mu,beta2.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))

  hist(betas[,3],xlab=bquote(beta[3]),xlim=c(-3,3),main="",breaks=20,freq=FALSE)
  abline(v=beta3.mu,lwd=2,col=2)
  curve(dnorm(x,beta3.mu,beta3.sd),lwd=3,col=3,add=TRUE)
  legend("topright",lwd=3,col=c("gray","red","green"),legend=c("Indiv. Coefs", "Pop. Mean","True Distribution"))



## ----covs.plot, class.source = "fold-hide",fig.height=6,fig.width=8-----------------------
# Combine into 1 raster stack
  covs = stack(covs[[1]], 
               covs[[2]],
               covs[[3]]
               ) 

# Change the extent to be larger to accommodate more realistic movements 
#  (not required but makes me feel better)
  extent(covs) = c(0,4000,0,4000)
  
# Names of variables
  names(covs[[1]]) = 'dist.dev'
  names(covs[[2]]) = 'forest'
  names(covs[[3]]) = 'shrub'
  
par(mfrow=c(2,3))
  plot(covs[[1]], main='dist.dev')
  plot(covs[[2]], main='forest')
  plot(covs[[3]], main='shrub')
  hist(values(covs[[1]]), main='dist.dev')
  hist(values(covs[[2]]), main='forest')
  hist(values(covs[[3]]), main='shrub')


## ----hsf.true-----------------------------------------------------------------------------
  hsf.true = vector("list",n.indiv)
  for (i in 1:n.indiv){
    hsf.true[[i]] = (covs[[1]]*betas[i,1]+
                     covs[[2]]*betas[i,2]+
                     covs[[3]]*betas[i,3]
                     )
  }


## ----simulate.data, cache=TRUE------------------------------------------------------------
# Set number of available samples per used
  n.avail = 100

# Number of movements
  n.time = 400	

# Number of possible steps to choose from at each iteration- these are not available locations, this is for the simulation of the movement path  
  n.step = 400	
  
# Population (across individual) turning angle parameters for von mises distribution
  angle.mean = 0
  angle.kappa = 0.0001  

# Step length rate of exponential distribution
  step.rate = 0.01
  
# Simulate individual-level movement-based habitat selection data
  sims =  sim.ind.movement.hsf(hsf.true = hsf.true,
                               n.time = n.time,
                               n.step = n.step,
                               n.avail = n.avail,
                               angle.mean = angle.mean,
                               angle.kappa = angle.kappa,
                               step.rate = step.rate
                               )


## ----simulate.data2-----------------------------------------------------------------------
  names(sims)

# Individual 1 data    
  head(sims$locs[[1]])
  head(sims$indiv[[1]])
  table(sims$indiv[[1]]$status)


## ----plot.tracks, class.source = "fold-hide"----------------------------------------------
#Plot the true HSF with locations for individual k
  k=10
  par(mfrow=c(1,1))
  plot(hsf.true[[k]])
  points(sims$locs[[k]]$use.x,
         sims$locs[[k]]$use.y,
         pch=16
         )


## ----organize. simualted.data-------------------------------------------------------------
# Combine data into a single data.frame
  sims2 = do.call("rbind", sims$indiv)
  head(sims2)
  dim(sims2)


## ----organize2. simualted.data------------------------------------------------------------
# Create ID vector for individual's data
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

  ID = rep(LETTERS702[1:n.indiv],each=nrow(sims$indiv[[1]]))
  sims2$ID=ID

# Create new strata
  sims2$indiv.id.strata = paste0(sims2$ID,
                                 sims2$strata
                                 )
# The number of unique identifiers needed
# n.step*n.indiv
  length(unique(sims2$indiv.id.strata)) == n.step*n.indiv
  
  dim(sims2)  
  head(sims2)


## ----first.fit.inference------------------------------------------------------------------
# We will use data from individual 1
  indiv.data1 = sims$indiv[[1]]

# Let's look at the data to remind us what it is
  head(indiv.data1)


## ----first.fit.inference2-----------------------------------------------------------------
  model1 = clogit(status ~ dist.dev + forest + shrub + strata(strata), 
                  data = indiv.data1
                  )	

# Look at estimated coefficients
  summary(model1)


## ----amt function-------------------------------------------------------------------------
amt::fit_ssf


## ----amt----------------------------------------------------------------------------------
  model1.amt1 = amt::fit_ssf(data = indiv.data1,
                             formula = status ~ dist.dev + forest + 
                                                shrub + strata(strata),
                             model=TRUE
                            )

  model1.amt2 = amt::fit_clogit(data = indiv.data1,
                                formula = status ~ dist.dev+forest+ 
                                                   shrub + strata(strata),
                                model=TRUE
                                )
  
  coef(model1)
  coef(model1.amt1)
  coef(model1.amt2)


## ----first.fit.inference.poisson.tmb,warnings=FALSE---------------------------------------
# Fit the model with fixed effect stratum-specific intercepts
	model1.tmb = glmmTMB(status ~ dist.dev + forest + 
	                              shrub + strata(strata), 
	                     data = indiv.data1,
	                     family = poisson
	                     )

# Or using a random effect with fixed variance
	model1.tmb2 = glmmTMB(status ~ dist.dev + forest + 
	                               shrub + (1| strata), 
	                      data = indiv.data1,
	                      family = poisson,
	                      doFit = FALSE
	                      )

# Make the intercept variance large and fixed
  model1.tmb2$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the variance of the intercept  
  model1.tmb2$mapArg = list(theta=factor(c(NA)))
  
# Now ready to fit the model  
  model1.tmb2 = glmmTMB:::fitTMB(model1.tmb2)


## ----poisson.clogit.comparison1-----------------------------------------------------------
knitr::kable(summary(model1.tmb2)[[6]]$cond[2:4,],digits=3)


## ----poisson.clogit.comparison2-----------------------------------------------------------
knitr::kable(summary(model1.tmb)[[6]]$cond[2:4,],digits=3)


## ----poisson.clogit.comparison3-----------------------------------------------------------
knitr::kable(summary(model1)[[7]][,-2],digits=3)


## ----RS-----------------------------------------------------------------------------------
# Get estimates
  coef = coef(model1)

# Relative Selection
  RS = exp(-2*coef[1] + 2*coef[2] + 2*coef[3]) / exp(2*coef[1] + 2*coef[2] + -2*coef[3])
  RS


## ----RSS----------------------------------------------------------------------------------
# Coefficient for forest
  exp(coef[2])


## ----preds.amt----------------------------------------------------------------------------
# Make a new data.frame for s1
  s1 = data.frame(
          dist.dev = seq(from = -2, to = 2, length.out = 100),
          forest = 0,
          shrub = 0,
          strata = 1 # needed, but the number is irrelevant
          )

# data.frame for s2
  s2 = data.frame(
         dist.dev = 0, # mean of dist.dev
         forest = 0,
         shrub = 0,
         strata = 1 # needed, but the number is irrelevant
    )

# Calculate log-RSS
  lr2 = amt::log_rss(model1.amt1, 
                     s1, 
                     s2, 
                     ci = "se", 
                     ci_level = 0.95
                     )


## ----plot.rss, class.source = "fold-hide"-------------------------------------------------
  plot(lr2,lwd=3)
  lines(lr2$df$dist.dev_x1,lr2$df$lwr,lty=2)
  lines(lr2$df$dist.dev_x1,lr2$df$upr,lty=2)
  abline(h=0,lwd=2,col=2)


## ----sensitiivity, cache=TRUE,warnings=FALSE----------------------------------------------
# Grab individual one's data
  indiv.dat1 = sims$indiv[[1]]

#Size of available samples per used locaiton
  n.avail2 = seq(2,100,by=2)
  
#Save coefficients
  coef.save = NULL
  
#Loop across available sample sizes  
  for(i in 1:length(n.avail2)){
  #Loop through each used location and reduce the number of available samples
    ind.dat=NULL
  for(j in 1:n.step){
    index = which(indiv.dat1$strata==j)
    #Get the used locaiton and the sub-sampled availables
    ind.dat.temp = indiv.dat1[index,][1:(n.avail2[i]+1),]
    ind.dat = rbind(ind.dat,ind.dat.temp)
  }
    
# fit model with data
	model.temp = clogit(status ~ dist.dev + forest + shrub + strata(strata), 
                      data = ind.dat
                      )	
  coef.save = rbind(coef.save, 
                    coef(model.temp)
                    )
}

coef.save = data.frame(n.avail2,coef.save)
colnames(coef.save) = c("N.Avail","beta1","beta2","beta3")


## ----plotting.sensitiivty, class.source = "fold-hide"-------------------------------------
par(mfrow=c(1,1))
plot(coef.save$N.Avail, coef.save$beta1,lwd=3,type="l",col=2,ylim=c(-01,1),
     main="Slope Coefficients",
     xlab="Number of Available Samples",ylab="Coefficient Estimate")
lines(coef.save$N.Avail, coef.save$beta2,lwd=3,col=2)
lines(coef.save$N.Avail, coef.save$beta3,lwd=3,col=2)

#knitr::kable(coef.save,digits=3)


## ----sensitiivity2, cache=TRUE,warnings=FALSE---------------------------------------------
# Grab individual one's data
  indiv.dat1 = sims$indiv[[1]]

# Size of available samples per used locaiton
  n.avail2 = c(2,10,20,40,60,80,100)
  
# Number of sample replications
 n.sample = 20
  
# Save coefficients
  coef.save = NULL
  
#Loop across available sample sizes  
  for(i in 1:length(n.avail2)){
  
  #re-sample each available sample size 20 times  
  for (z in 1:n.sample){  
  #Loop through each used location and reduce the number of available samples
    ind.dat = NULL
  for(j in 1:n.step){
    index = which(indiv.dat1$strata==j)
    #Get the used locaiton and the sub-sampled availables
    rand.sample = sample(index[-1],n.avail2[i])
    ind.dat.temp = indiv.dat1[c(index[1],rand.sample),]
    ind.dat = rbind(ind.dat,ind.dat.temp)
  } # end j step loop
    
# fit model with data
	model.temp = clogit(status ~ dist.dev + forest + shrub + strata(strata), 
                      data = ind.dat
                      )	
  coef.save = rbind(coef.save, 
                   coef(model.temp)
                  )
  }#end z sample loop
}#end i loop


## ----plot.sensitive, class.source = "fold-hide"-------------------------------------------
one = data.frame(N.Avail=rep(n.avail2,each=n.sample),
                 Sample=rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,1]
                )
two = data.frame(N.Avail=rep(n.avail2,each=n.sample),
                 Sample=rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,2]
                )
three = data.frame(N.Avail=rep(n.avail2,each=n.sample),
                 Sample=rep(1:n.sample,length(n.avail2)),
                 beta = coef.save[,3]
                )
plot.data = rbind(one,two,three)
plot.data$Name = c(rep("dist.dev",nrow(one)),
                   rep("forest",nrow(one)),
                   rep("shrub",nrow(one))
                 )

plot.data$N.Avail = as.factor(plot.data$N.Avail)


colnames(plot.data) = c("N.Available.Sample","Sample","Coefficient.Estimate","Name")

ggplot2::ggplot(plot.data, aes(N.Available.Sample, Coefficient.Estimate, fill=factor(Name))) +
  theme_bw()+
  geom_boxplot()


## ----pooling, cache=TRUE, warnings=FALSE, error=FALSE, message=FALSE----------------------
# Pooled model - no consideration of individual variability
  model.pool = glmmTMB(status ~ dist.dev + forest + shrub +
                               (1|indiv.id.strata), 
	                     data=sims2,
	                     family=poisson,
	                     doFit = FALSE
	                     )
  model.pool$parameters$theta[1] = log(1e3)
  model.pool$mapArg = list(theta=factor(c(NA)))
  options(warn=-1)
  model.pool = glmmTMB:::fitTMB(model.pool)
  options(warn=0)
# Look at estimates when pooling the data
  summary(model.pool)


## ----indiv.separate.models,cache=TRUE, warnings=FALSE, error=FALSE, message=FALSE---------
# Separate models to each individual
  indiv.fit = lapply(sims[[1]],FUN=function(x){
                  temp=glmmTMB(status ~ dist.dev+forest+shrub + (1|strata), 
	                             data=x,
	                             family=poisson,
                               doFit = FALSE
	                             )

              temp$parameters$theta[1] = log(1e3)
              temp$mapArg = list(theta=factor(c(NA)))
              options(warn=-1)
              temp = glmmTMB:::fitTMB(temp)
              options(warn=0)
              temp
              })


## ----separate-----------------------------------------------------------------------------
summary(indiv.fit[[2]])


## ----separate2----------------------------------------------------------------------------
summary(indiv.fit[[10]])


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
       width=c(2,1))
plotCI(1:n.indiv,
       y=estimates[1,],
       ui=UCI[1,],
       li=LCI[1,],
       ylab="beta",
       xlab="Individual",
       lwd=2,
       ylim=c(-2.5,4)
)
 legend("topright",lwd=3,col=c(1,2,3),
        legend=c("dist.dev","forest","shrub"),box.lty=0,y.intersp=0.8)
 
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
# Source summary function
  boot.pop.esimates=mean_ci_boot(boot)
  rownames(boot.pop.esimates)=c("dist.dev","forest","shrub")
  knitr::kable(boot.pop.esimates,digits=3)


## ----random.intercept, cache=TRUE---------------------------------------------------------
# Setup HSF with Poisson regression approximation - random intercept model
# indiv.id.strata indicates individuals and strata and there the variance is fixed so there is no shrinkage
  re_int = glmmTMB(status ~ dist.dev + forest + shrub  +
                             (1|indiv.id.strata), 
                    family=poisson ,
                    data = sims2, 
                    doFit=FALSE
                    )

# Make the intercept variance large and fixed (i.e., do not estimate).
  re_int$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the variance of the intercept  
  re_int$mapArg = list(theta=factor(c(NA)))
  
# Now ready to fit the model  
  options(warn=-1)
  re_int = glmmTMB:::fitTMB(re_int)
  options(warn=0)


## ----examine.random.intercept.model-------------------------------------------------------
  fixef(re_int)


## ----examine2-----------------------------------------------------------------------------
  head(ranef(re_int)[[1]][[1]])
  summary(re_int)


## ----random.int.slope, cache=TRUE---------------------------------------------------------
#Fit RSF with intercept and slopes with random effect
  re_int_slopes =  glmmTMB(status ~ -1 + dist.dev + forest + shrub  +
                                    (1|indiv.id.strata) +
                                    (0+dist.dev|ID) + (0+forest|ID) + 
                                    (0+shrub|ID), 
                            family = poisson, 
                            data = sims2, 
                            doFit = FALSE
                            )
# make the intercept variance large and fixed (i.e., do not estimate).
  re_int_slopes$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the variance of intercept, but 
# to do so for the other three variables
  re_int_slopes$mapArg = list(theta=factor(c(NA,1:3)))
  
# Now ready to fit the model 
  options(warn=-1)
  re_int_slopes = glmmTMB:::fitTMB(re_int_slopes)
  options(warn=0)


## ----RE.results.fixed1--------------------------------------------------------------------
  fixef(re_int_slopes)


## ----RE.results.fixed2--------------------------------------------------------------------
  broom.mixed::tidy(re_int_slopes, 
                    effects = "fixed", 
                    conf.int = TRUE)[,c(3,4,8,9)]


## ----RE.results.random--------------------------------------------------------------------
  ranef.out = ranef(re_int_slopes)
  ranef.out$cond$ID
  ranef.out.ci = broom.mixed::tidy(re_int_slopes, 
                                   effects = "ran_vals", 
                                   conf.int = TRUE)

#Non-intercept estimates per indivual  
  ranef.out.ci[-which(ranef.out.ci$term=="(Intercept)"),-c(1,2,3,4)]


## ----RE.results---------------------------------------------------------------------------
  summary(re_int_slopes)


## ----RE.est.forest.indiv.pop--------------------------------------------------------------
  indiv.coef.popModel = broom.mixed::tidy(re_int_slopes, effects = "ran_vals", conf.int = TRUE)
  index=which(indiv.coef.popModel$term=="forest")
  indiv.forest.popModel = data.frame(indiv.coef.popModel[index,c(5,6,8,9)])

# Add back mean to get individual estimates  
  indiv.forest.popModel[,2:4] = indiv.forest.popModel[,2:4]+rep(fixef(re_int_slopes)[[1]][2],each=20)
    
# Extract population mean and uncertainty for forest effect
  pop.coef.popModel = data.frame(broom.mixed::tidy(re_int_slopes, effects = "fixed", conf.int = TRUE))
  pop.coef.popModel=pop.coef.popModel[,c(4,8,9)]


## ----RE.plot.forest.indiv.pop, class.source = "fold-hide"---------------------------------
#Plot
  plotCI(x = 1:n.indiv,
         y = indiv.forest.popModel$estimate,
         ui = indiv.forest.popModel$conf.high,
         li = indiv.forest.popModel$conf.low,
         lwd = 3,col=1,
         xlab="Individual",ylab="Forest Coefficient")
  abline(h=pop.coef.popModel$estimate[2],lwd=3,col=1,lty=4)
  abline(h=pop.coef.popModel$conf.low[2],lwd=3,col=2,lty=4)
  abline(h=pop.coef.popModel$conf.high[2],lwd=3,col=2,lty=4)


## ----behav.sim----------------------------------------------------------------------------
  hsf.true.behav = vector("list",2)
  hsf.true.behav[[1]] = covs[[2]]*2
  hsf.true.behav[[2]] = covs[[2]]*-2

  sim.behav1 =  sim.ind.movement.hsf(hsf.true = hsf.true.behav,
                                     n.time = n.time,
                                     n.step = n.step,
                                     n.avail = n.avail,
                                     angle.mean = angle.mean,
                                     angle.kappa = angle.kappa,
                                     step.rate = step.rate
                                     )

# Combine the data  
  data.ignore= rbind(sim.behav1$indiv[[1]],
                     sim.behav1$indiv[[2]]
                     )  
  
# Create ID vector for the different behavioral data
  ID=rep(1:2,each = nrow(sim.behav1$indiv[[1]]))
  data.ignore$ID = ID
# Create ID vector for unique strata within individual
  data.ignore$indiv.id.strata=unique(data.ignore$ID+data.ignore$strata*100)


## ----behav.fit----------------------------------------------------------------------------
  model.ignore.behav =  clogit(status ~ dist.dev + strata(indiv.id.strata), 
                               data = data.ignore
                               )	
  summary(model.ignore.behav)


## ----packages.cite, eval=TRUE, echo=FALSE-------------------------------------------------
  pkgs = cite_packages(output = "table", out.dir = ".",out.format ="pdf")
  knitr::kable(pkgs)

