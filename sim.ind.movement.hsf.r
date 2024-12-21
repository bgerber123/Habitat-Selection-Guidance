sim.ind.movement.hsf = function(hsf.true=hsf.true,
                                n.time=n.time,
                                n.step=n.step,
                                n.avail=n.avail,
                                angle.mean=angle.mean,
                                angle.kappa=angle.kappa,
                                step.rate=step.rate
                                ){
  
  n.indiv=length(hsf.true)  
  
#Constrained random starting of each kth simuation (contrained within the sample space)
  set.seed(343554)
  start.x=runif(n.indiv,400,700)
  start.y=runif(n.indiv,400,700)
  
  # Create empty vectors and the use stratification variable
  avail.x=NULL
  avail.y=NULL
  strat.u=seq(1, n.time, 1)
  strat.a=NULL


# Object to save individual data  
  locs.data=indiv.data = vector("list",n.indiv)
  
  for(k in 1:n.indiv){
    step.save=NULL
    #Random starting point b/w 400 and 700
    x=start.x[k]
    y=start.y[k]
    
    # simulate movement path
    for(j in 1:n.time){
      
      # Draw n.step turn angles and step lengths
      step=rexp(n.step, rate=step.rate)
      options(warn=-1)
      turn=rvonmises(n = n.step,
                     m = angle.mean,
                     k = angle.kappa
                     )
      options(warn=0)
      # Calculate the difference in x and y direction
      delta.y=sin(turn)*step
      delta.x=cos(turn)*step
      
      # From that calculate new x and y
      new.x=as.numeric(x[j]+delta.x)
      new.y=as.numeric(y[j]+delta.y)
      
      # Extract the true underlying hsf and exponentiate
      hsf.ext=exp(extract(hsf.true[[k]],
                          cbind(new.x,new.y)
                          )
                  )
      
      # Scale so that all values are <=1 but still same relative magnitude
      maxhsf=raster::cellStats(hsf.true[[k]], max)
      hsf.10=hsf.ext/maxhsf
      
      # Select a single point from these
      indic=rmultinom(1, size=1, prob=hsf.10)
      
      # keep the x and y coordinates at that location
      x[j+1]=new.x[indic==1]
      y[j+1]=new.y[indic==1]
      
      # Save the step length
      step.save[j]=step[indic==1]
    }
    
    # # Plot	
    # par(mfrow=c(1,1))
    # plot(hsf.true[[i]])
    # points(x,y,pch=16)
    
    # Calculate the mean step length for use in drawing availables
      mean.step=mean(step.save)
    
    # Drop first location as it is not endpoint of step
      use.x=x[2:length(x)]
      use.y=y[2:length(y)]
    
    # Empty vectors to save things
    
      step.cov=NULL
      avail.x=NULL
      avail.y=NULL
      strat.a=NULL
    
    # Draw availables
    
    for(i in 1:length(use.x)){
      
      # For each used location draw n.avail availables
      # For step draw from exponential with rate=1/(2SL) where SL = mean step length 
      options(warn=-1)
      step=rexp(n.avail, rate=1/(2*mean.step))
      turn=rvonmises(n.avail,
                     angle.mean,
                     angle.kappa
                     )
      options(warn=0)
      # Calculate change in x and y
      delta.y=as.numeric(sin(turn)*step)
      delta.x=as.numeric(cos(turn)*step)
      
      # Calculate new x and y- because i dropped the first x and y above, can just use i to indicate which x and y the steps depart from
      new.x=as.numeric(x[i]+delta.x)
      new.y=as.numeric(y[i]+delta.y)
      
      # Append onto available and create strata
      avail.x=c(avail.x, new.x)
      avail.y=c(avail.y, new.y)
      strat.a=c(strat.a, rep(i, n.avail))
      step.cov=c(step.cov, step)
    }
    
    # Now extract covariate the covariates
    covs.u=extract(covs, cbind(use.x, use.y))
    covs.a=extract(covs, cbind(avail.x, avail.y))
    
    # Set up the rest of the data file
    use=c(rep(1, length(use.x)), rep(0, length(avail.x)))
    strat=c(strat.u, strat.a)
    covs.x=rbind(covs.u, covs.a)
    step.x=c(step.save, step.cov)
    
    data=data.frame(use, strat, covs.x, step.x)
    locs=data.frame(use.x,use.y)
    names(data)=c("status","strata","dist.dev","forest","shrub","step")
    
    data$dist.dev = scale(data$dist.dev)
    data$forest = scale(data$forest)
    data$shrub = scale(data$shrub)
    
    indiv.data[[k]] = data[order(data$strata),]
    locs.data[[k]] = locs
  } #loop across individuals  
  
  list(indiv=indiv.data,locs=locs.data)
  
} #end function
