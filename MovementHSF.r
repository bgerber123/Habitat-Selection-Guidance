# Setup Environment
  rm(list=ls())
  library(geoR)
  library(circular)
  library(glmmTMB)
  library(raster)
  library(sp)
  library(survival)
  library(Rfast)
  source("sim.ind.movement.hsf.r")
  load("Covs")
##############################
#Number of individuals to simulate data
  n.indiv = 20

# Population-level coefficients  
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
  
###############################
# Look at covariates

par(mfrow=c(2,3))
  plot(covs[[1]])
  plot(covs[[2]])
  plot(covs[[3]])
  hist(values(covs[[1]]))
  hist(values(covs[[2]]))
  hist(values(covs[[3]]))

# Combine into 1 raster stack
  covs=stack(covs[[1]], covs[[2]],covs[[3]]) 

# Change the extent to be larger to accommodate more realistic movements 
#  (not required but makes me feel better)
  extent(covs)=c(0,1000,0,1000)

# create the linear combination of the true HSF as a raster
# for each individual
# Pull out the max hsf for use in scaling
  
  hsf.true=vector("list",n.indiv)
  maxhsf=rep(NA,n.indiv)
  for (i in 1:n.indiv){
    hsf.true[[i]]=(covs[[1]]*betas[i,1]+covs[[2]]*betas[i,1]+covs[[3]]*betas[i,3])
    maxhsf[i]=cellStats(hsf.true[[i]], max)
  }
  
  
  
# Simulate individual-level movement-based habitat selection data
  sim.data =  sim.ind.movement.hsf(hsf.true)

  
#Plot the true HSF with locations for individual k
  k=10
  par(mfrow=c(1,1))
  plot(hsf.true[[k]])
  points(sim.data$locs[[k]]$use.x,
         sim.data$locs[[k]]$use.y,
         pch=16
         )
  
  
	
	
	# Fit with conditional logistic
model.fit=clogit(use~x1+x2+x3+x4+x5+x6+x7+x8+x2:x5+x6:x8+x1:x5+strata(strat), data=data)	
	beta.save=rbind(beta.save,summary(model.fit)$coefficients[,1])
	

	# Fit with TM model builder using Poisson 
	
	# TMBstruc=glmmTMB(use~x5+(1|strat),family=poisson, data=data, doFit=FALSE)
	# TMBstruc=glmmTMB(use~x1+x2+x3+x4+x5+x6+x7+x8+x2:x5+x6:x8+x1:x5+(1|strat),family=poisson, data=data, doFit=FALSE)
	# TMBstruc$parameters$theta[1]=log(1e6)
	# TMBstruc$mapArg = list(theta=factor(c(NA)))
	# mod=glmmTMB:::fitTMB(TMBstruc)
	# beta.save=rbind(beta.save,summary(mod)$coefficients$cond[,1])
	
	
	