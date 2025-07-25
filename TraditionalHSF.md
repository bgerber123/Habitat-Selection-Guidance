---
title: "Traditional Habitat Selection <br> Guidance"
author: "Brian D. Gerber, Casey Setash, Jacob S. Ivan. and Joseph M. Northrup"
date: "2025-06-17"
output: 
  html_document:
    toc: true
    toc_float: true
    number_sections: true
    code_folding: show
    keep_md: true
---

```{=html}
<style type="text/css">
body, td {
   font-size: 16px;
}
code.r{
  font-size: 16px;
}
pre {
  font-size: 16px
}
</style>
```

<style type="text/css">
  #TOC {
    max-width: fit-content;
    white-space: nowrap;
  }
  
  div:has(> #TOC) {
    display: flex;
    flex-direction: row-reverse;
}
</style>



## Introduction

This vignette is associated with the manuscript titled, 'A plain language review and guidance for modeling animal habitat-selection'. We will be demonstrating some of the points from the manuscript on fitting models to estimate traditional habitat selection functions (HSF) at the third-order of selection (e.g. selection within the home range).

*Note:* some of the code-chunks are not automatically displayed. To show the code, select 'show'.

<br>

## Setup

### Environment


```{.r .fold-hide}
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
```

<br>

### Simulate data

We will consider a habitat selection analysis of individuals within a broad geographic area (e.g., home range; third-order of selection). Within each individual's home range the effects of predictors (i.e., $\beta_1$ and $\beta_2$) on habitat selection are assumed to come from a distribution of effects that can be characterized by a mean and standard deviation (i.e., random effect).


```{.r .fold-hide}
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
```


```{.r .fold-hide}
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
```

![](TraditionalHSF_files/figure-html/visualize.true.coefficients-1.png)<!-- -->

We are now ready to simulate individual-level data. We will do this with the `simulateUsedAvail` function from the R package `ResourceSelection`. Because of this, we are not simulating spatial (x-y) data from explicit spatial layers, but we will connect how these data represent typical GPS data and setup when working with spatial variables.


```{.r .fold-hide}
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
```

For more information on why a log link indicates relative selection, see Fieberg et al. (2021):

"Instead of focusing on $p_i$, as is typical in applications to presence–absence data, logistic regression applied to use-availability data should simply be viewed as a
convenient tool for estimating coefficients in a habitat-selection function, $w(X(s); \beta) = \text{exp}(X_{1}(s)\beta_1 + ... X_{k}(s)\beta_k)$  (Boyce & McDonald, 1999; Boyce et al., 2002), where we have X(s) to highlight that the predictors correspond to measurements at specific point locations in geographical space, s. As we will see in the next section, this expression is equivalent to the intensity function of an IPP model but with the intercept (the log of the baseline intensity) removed; the baseline intensity gives the expected density of points when all covariates are 0. Because habitat-selection functions do not include this baseline intensity, they are said to measure ‘relative probabilities of use’, or alternatively, said to be ‘proportional to the probability of use’ (Manly et al., 2002). "


Before moving forward lets understand the data for one individual.


``` r
# Check the dimensions of the data  
# This should have the same length as n.indiv  
  length(sims)
```

```
## [1] 20
```

``` r
# Look at one individual's dataset
  dim(sims[[1]])
```

```
## [1] 5500    4
```

``` r
  head(sims[[1]])
```

```
##      status    dist.dev      forest     shrub
## 1969      1  0.14004930 0.162122994 0.2868476
## 3827      1 -0.43224008 0.006257207 0.8402856
## 3266      1 -3.18972176 0.905023379 0.2473926
## 2556      1 -1.73908829 0.714554879 0.6718038
## 2241      1 -0.73270288 0.738071773 0.5666879
## 498       1  0.09237958 0.246186710 0.8294609
```

``` r
  table(sims[[1]]$status)
```

```
## 
##    0    1 
## 5000  500
```

We see from this individuals' data frame that we have the columns status, dist.dev, forest, shrub. The column `status` is our response variable. A `1` indicates an animal location (used point) and a `0` indicates a potentially used location (i.e., available location) within the individual's home range. Each 1 and 0 have corresponding spatial variables, which we have called `dist.dev` (Euclidean distance to nearest development), `forest` (percent of landcover that is forested), and `shrub` (percent of landcover that is covered in shrub). Each of these is a continuous spatial variable and has been standardized to have a mean of 0 and standard deviation of 1 (also called Normalizing). This is a common setup in statistical models to put all the variables on an equivalent scale, i.e., the scale of one standard deviation. Note that we have a large available sample, much more than the used sample. This is correct. Remember that the available sample (the 0's) is what is going to allow us to use logistic regression to approximate the IPP model.

We now want to package our simulated data into a single data frame. We will use the objects `sims` and `sims2` below when we want to fit models to each individual separately and together, respectively.


``` r
# Combine data into a single data.frame
  sims2 = do.call("rbind", sims)
  head(sims2)
```

```
##      status    dist.dev      forest     shrub
## 1969      1  0.14004930 0.162122994 0.2868476
## 3827      1 -0.43224008 0.006257207 0.8402856
## 3266      1 -3.18972176 0.905023379 0.2473926
## 2556      1 -1.73908829 0.714554879 0.6718038
## 2241      1 -0.73270288 0.738071773 0.5666879
## 498       1  0.09237958 0.246186710 0.8294609
```

``` r
  dim(sims2)
```

```
## [1] 110000      4
```

``` r
#Create ID vector for individual's data
  ID = rep(1:n.indiv, each = nrow(sims[[1]]))
  sims2$ID = ID
  head(sims2)
```

```
##      status    dist.dev      forest     shrub ID
## 1969      1  0.14004930 0.162122994 0.2868476  1
## 3827      1 -0.43224008 0.006257207 0.8402856  1
## 3266      1 -3.18972176 0.905023379 0.2473926  1
## 2556      1 -1.73908829 0.714554879 0.6718038  1
## 2241      1 -0.73270288 0.738071773 0.5666879  1
## 498       1  0.09237958 0.246186710 0.8294609  1
```

``` r
  dim(sims2)
```

```
## [1] 110000      5
```

``` r
# Create a vector to weight the 0's by 1000 and the 1's by 1
# This is discussed below  (see manuscript section 10. The model being approximated)
  sims2$weight = 1000^(1-sims2$status)
  head(sims2)
```

```
##      status    dist.dev      forest     shrub ID weight
## 1969      1  0.14004930 0.162122994 0.2868476  1      1
## 3827      1 -0.43224008 0.006257207 0.8402856  1      1
## 3266      1 -3.18972176 0.905023379 0.2473926  1      1
## 2556      1 -1.73908829 0.714554879 0.6718038  1      1
## 2241      1 -0.73270288 0.738071773 0.5666879  1      1
## 498       1  0.09237958 0.246186710 0.8294609  1      1
```

<br>

## Manuscript sections

We are now ready to fit models to estimate our habitat selection parameters and demonstrate points made in the manuscript. We have organized the sections below to generally match with the sections of the manuscript.

<br>

### What is habitat?

Definitions only.

<br>

### What is habitat selection?

Definitions only.

<br>

### What is a habitat-selection function?

Definitions only.

We want to remind the reader of the nuance between a habitat selection function and the statistical model we will be fitting. The model we are fitting (as coded) and the true underlying model (Inhomogeneous Poisson Point process; IPP) is not the habitat selection function. The habitat selection function is a component of the IPP model. We will be using a logistic regression modeling (or Generalized Binomial Regression Model) to approximate the IPP model, in which our estimated parameters will be the parameters within the habitat selection function. To understand the model, see equation 1 of [Gerber and Northrup, 2020](https://doi.org/10.1002/ecy.2953).

It is also important to differentiate a habitat selection function from a habitat selection analysis. As mentioned above, a habit selection function is a component of the model we want to fit, while a habitat analysis refers to all the steps that are being considered; this includes, for example, how the data are setup, the model fit, the algorithms used, the inference or predictions from the model, and how the model is evaluated. 

<br>

### ‘Habitat selection function’ or ‘resource selection function’? 

Definitions only.

<br>

### What about scale?

Definitions only.

In our data setup and analysis we will be estimating habitat selection at the third-order of the scales defined in Johnson (1980). Most commonly, we can think of this as selection of habitat within a defined home-range. Our data generating process above was not specific about the spatial bounds of each individual's home-range, but this is how you can think of the data. The animal telemetry locations were used to create the 1's in our data and the corresponding spatial covariates (e.g., dist.dev) as well as the home range boundary. We then make a grid of points within the home range to extract the spatial covariate values and we assign these to the response of 0 and are the available sample. These two sets of values (1's and covariate values and 0's and covariate values) are appended together in a single data frame.

There are many home-range estimators in the literature. Choose one that meets the requirements of your data, for example, whether you have assumed independence between consecutive spatial locations (e.g., kernel density estimation) or whether spatial autocorrelation is due to a high-frequency of sampling (e.g., short fix interval between relocations relative to the animals speed; movement-based kernel estimator).

<br>


### Considering objectives and data collection

This is the hardest part of the whole habitat selection study. How do you decide on the what, where, and how many individuals and relocations to answer the question at hand. Talk to many people. Ask questions to many people.

The modeling mentioned here is about differences between study goals, such as inference and prediction. Modeling for inference versus prediction is not a straightforward distinction. There are lots of opinions from the statistical folks (which are great and everyone should read them, e.g., [Shmueli, 2010](https://doi.org/10.1214/10-STS330) and [Scholz and Burner, 2022](https://arxiv.org/abs/2210.06927)) and the science philosophy folks.

We reference the manuscript [Gerber and Northrup, 2020](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2953) in regards to when the study goal is prediction. Associated with this manuscript is code in the Supporting Information file (ecy2953-sup-0004-DataS1.Zip) that pertains to optimizing for predicting spatial distributions of habitat selection (i.e., making a map). This process can jeopardize inference, e.g., make the interpretation of estimated effects unreliable. In contrast, if inference is sought you should think hard about a single model that includes the most important variables for the species and for your question that you want to consider so that estimated effects and measures of uncertainty are reliable ([Bolker 2023](https://doi.org/10.32942/X2Z01P), [Tredennick et al. 2021](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3336)).

<br>

###  What is the scientific goal of implementing a habitat-selection function?

In this section, we discuss some mechanics of inference and prediction, as in interpreting coefficients and predicting relative habitat selection. We will walk through the basics here, but refer the reader to full in depth coverage in [Fieberg et al. 2021](https://doi.org/10.1111/1365-2656.13441) and [Avgar et al. 2017](https://doi.org/10.1002/ece3.3122). We also mention the goal of using habitat selection to predict animal abundance. Demonstrating this is beyond our work. We suggest starting with the reading of section 1.6, "When is species density a reliable reflection of habitat suitability?" of [Matthiopoulos, Fieberg, and Aarts, 2023](http://hdl.handle.net/11299/217469).

<br>

#### Inference

Let's fit a model to one individual's used and available data and discuss the coefficients.


``` r
# We will use data from individual #20
  indiv.data20 = sims[[20]]

# Let's look at the data to remind us what it is
  head(indiv.data20)
```

```
##      status   dist.dev    forest     shrub
## 4750      1 -1.7175434 0.8433439 0.5885524
## 2483      1 -1.4853814 0.2001264 0.7956416
## 4802      1  0.6492396 0.2110496 0.4721530
## 2363      1 -1.3261416 0.9622746 0.9805034
## 1674      1 -0.7316487 0.6683905 0.1342194
## 4756      1 -2.3376280 0.9587345 0.1744920
```

We are now ready to approximate our true model using the generalized linear modeling machinery of logistic regression. We will do this using two different packages and functions. First, we will use the `glmmTMB` function in package `glmmTMB`.


``` r
# We have related our response variable of the used and available sample (status) to our covariates. This is an additive model, as opposed to having interactions.  
  model1 = glmmTMB::glmmTMB(status ~ dist.dev + forest + shrub, 
                            family = binomial(), 
                            data = indiv.data20 
                            )

# Look at estimated coefficients
  summary(model1)
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: indiv.data20
## 
##      AIC      BIC   logLik deviance df.resid 
##   2817.6   2844.1  -1404.8   2809.6     5496 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -3.4354     0.1518 -22.634  < 2e-16 ***
## dist.dev     -1.1067     0.0542 -20.418  < 2e-16 ***
## forest        0.1244     0.1709   0.728    0.467    
## shrub         0.8778     0.1781   4.929 8.26e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Let's first consider the `Intercept`. As mentioned in the section '12. Population Inference', this parameter is largely the ratio of used to available sample. If we increase our available sample, this estimate will get smaller. Its value has a lot to do with how we setup the data to approximate the true underlying model. It has no biological interpretation, so we can skip it.

Next, we have estimated the effect of `dist.dev` as -1.11. This estimate is negative, indicating that as the value of `dist.dev` increases from its mean (defined at 0) while the values of `forest` and `shrub` (the other variables in the model) are held at their mean that habitat selection decreases. In other words, habitat use relative to what is available to this individual (as we defined it!) decreases the further from development. That's a lot of qualifiers to understand this estimate. These are important though. In terms of evidence of an effect, we can say that at a defined probability of Type I error of 0.05, we have a [statistically clear](https://doi.org/10.1111/2041-210X.13159) effect and that that it is not likely zero. The evidence of this is the very small p-value of 1.1558137\times 10^{-92}.

Next, we have estimated the effect of `forest` as 0.12. This estimate is positive, indicating that as the value of `forest` increases from its mean (defined at 0) while the values of `dist.dev` and `shrub` (the other variables in the model) are held at their mean that habitat selection increases relative to what is available to this individual. At our defined Type I error we do not have statistical clarity and there is a reasonable chance that it is zero. (p-value = 0.4665834).

Next, we have estimated the effect of `shrub` as 0.88. This estimate is positive, indicating that as the value of `shrub` increases from its mean (defined at 0) while the values of `dist.dev` and `forest` (the other variables in the model) are held at their mean that habitat selection increases relative to what is available to this individual. At our defined Type I error we do have statistical clarity of an effect that is not zero (p-value = 8.2627753\times 10^{-7}).

<br>

##### **Other functions**

We can estimate the parameters with other functions that implement the logistic regression model. For example, the `glm` function in the `stats` package as:


``` r
  model1.glm = stats::glm(status ~ dist.dev + forest + shrub, 
                          family = binomial(), 
                          data = indiv.data20 
                          )
  summary(model1.glm)
```

```
## 
## Call:
## stats::glm(formula = status ~ dist.dev + forest + shrub, family = binomial(), 
##     data = indiv.data20)
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -3.4354     0.1518 -22.634  < 2e-16 ***
## dist.dev     -1.1067     0.0542 -20.418  < 2e-16 ***
## forest        0.1244     0.1709   0.728    0.467    
## shrub         0.8778     0.1781   4.929 8.26e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 3351.0  on 5499  degrees of freedom
## Residual deviance: 2809.6  on 5496  degrees of freedom
## AIC: 2817.6
## 
## Number of Fisher Scoring iterations: 6
```

Notice that we estimate the exact same quantities.

The package `amt` can do the same thing using the function `fit_rsf` which is a wrapper function for using the `stats glm` function.


``` r
amt::fit_rsf
```

```
## function (data, formula, ...) 
## {
##     m <- stats::glm(formula, data = data, family = stats::binomial(link = "logit"), 
##         ...)
##     m <- list(model = m)
##     class(m) <- c("fit_logit", "glm", class(m))
##     m
## }
## <bytecode: 0x00000237297814d8>
## <environment: namespace:amt>
```

<br>

##### **Relative Selection**

How do we quantitatively evaluate two locations in terms of habitat selection using our model results? We can do so using Equation 8 of [Fieberg et al. 2021](https://doi.org/10.1111/1365-2656.13441).

Perhaps we want to compare a location that is very near development (dist.dev = -2) at high forest and shrub cover (forest = 2, shrub = 2) with that of a location far from development (dist.dev = 2) and also at high forest but low shrub cover (forest =2, shrub = -2).


``` r
# Get estimates without the intercept
  coef = fixef(model1)[[1]][-1]

  RS = exp(-2*coef[1] + 2*coef[2] + 2*coef[3]) / 
       exp(2*coef[1] + 2*coef[2] + -2*coef[3])
  RS
```

```
## dist.dev 
## 2801.048
```

If given equivalent availability between these two types of sites ( 1) near development at high forest and shrub cover, 2) far from development at high forest and shrub cover ), this individual would relatively select the first location by a factor of 2801.05.


Why are we exponentiating our linear combinations of terms (additions of covariates times estimated coefficients)? This is because in this setup, we have assumed the habitat selection process follows an exponential form (see section 10. The model being approximated) and equation 1 of [Gerber and Northrup, 2020](https://doi.org/10.1002/ecy.2953).

<br>

##### **Relative Selection Strength**

[Avgar et al. 2017](https://doi.org/10.1002/ece3.3122) refers to the relative selection strength (RSS) as the exponentiation of each coefficient.

For example,


``` r
# Coefficient for dist.dev
  exp(coef[2])
```

```
##   forest 
## 1.132473
```

Given two locations that differ by 1 standard deviation of dist.dev, but are otherwise equal, this individual would be 1.1324727 as likely to choose the location with higher dist.dev (or, equivalently, 0.8830235 times more likely to choose the location with the lower dist.dev.

<br>

#### Prediction

Lets combine our covariates now to predict relative habitat selection over more scenarios. To do this, we need to remember we used generalized linear modeling functions to approximate our model. We need to take care when predicting. First, we need to remember to drop the intercept. Second, we need to remember our true link function for relative habitat selection is the log-link (inverse-link is the exponential function) and not the logit-link, which is the default when fitting logistic regression models.


``` r
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
```


```{.r .fold-hide}
# Plot
  plot(newdata$dist.dev,preds.exp$fit,type="l",lwd=3,xlab="Distance to Development (Scaled)",
       ylab="Relative Intensity of Selection")
  lines(newdata$dist.dev,preds.exp$LCL,lwd=3,col=2,lty=4)
  lines(newdata$dist.dev,preds.exp$UCL,lwd=3,col=2,lty=4)
  abline(h=1, lwd=4, col=4)
  text(1.5,1.4, "Selection")
  text(1.5,0.6, "Avoidance")
```

![](TraditionalHSF_files/figure-html/predict.plot-1.png)<!-- -->

We see that for this individual, given equal availability and `forest` and `shrub` at their mean values, avoid areas far from development (high values on x-axis) and select for areas close to development (low values on x-axis).

We can use this same process to extract values across spatial layers to predict relative habitat selection that is mapped.

<br>

###  Traditional habitat selection function (HSF) or step-selection function (SSF)?

In fitting these data using a traditional HSF, we are making two important assumptions. First, we are assuming that the available locations for an individual are accessible and could be used if the individual chooses. Note that in our data setup, we have not assumed the set of available locations are the same for each individual. Each individual has a different set, based on their home-range (remember this is not explicitly how we sampled, but it is implied). Second, we are assuming that each used location is independent from each other. Lets consider the implications of violating these assumptions.

<br>

#### Critical Assumption 1: Accessibility of habitat

We can demonstrate this effect by comparing results from individual 20 (above) with that of the same model but additional available locations are added to the available sample. Essentially, we are going to consider the effect of assuming these new locations are available to this individual.


``` r
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
```



|                               |   dist.dev|    forest|     shrub|
|:------------------------------|----------:|---------:|---------:|
|Original Data                  | -1.1066673| 0.1244035| 0.8777699|
|Data with Additional Locations | -0.9722748| 0.0645601| 0.2185216|

**What did we find?** Simply, that our estimates are quite different when considering these additional available locations that were inaccessible to the individual. That is a challenge for us when setting up the data. Our estimates can change quite a bit. This is an important reminder that our estimated effects highly depend on our assumptions of what is available to each animal. We want to avoid including available locations that are in fact inaccessible to the animal because it will change our inference in unknown ways.

<br>

#### Critical Assumption 2: Independence of location data

We can evaluate this assumption by including non-independent data. We can easily do this by replicating our used animal location data. It's as if we were tracking an individual and took two locations every 4 hours, where the 2 location fixes were only separated by 5 minutes. While 4 hours may be enough for the individual to travel throughout its home range (if it chooses), 2 minutes is certainly not enough time. Thus, we have two data points really close in space and time for every 4 hour interval.


``` r
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
```



|                     |   dist.dev|     forest|    shrub|
|:--------------------|----------:|----------:|--------:|
|Original Data        | -0.6883261| -0.2891435| 1.135597|
|Non-Independent Data | -0.6905494| -0.2946979| 1.147856|

**What did we find?** Well, the point estimates look similar. That's good. Now let's look at measures of uncertainty.


``` r
# Compare measures of uncertainty
  summary(model1) # independent data
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: indiv.data7
## 
##      AIC      BIC   logLik deviance df.resid 
##   3110.5   3136.9  -1551.2   3102.5     5496 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -3.01865    0.14274 -21.147  < 2e-16 ***
## dist.dev    -0.68833    0.04961 -13.875  < 2e-16 ***
## forest      -0.28914    0.16725  -1.729   0.0838 .  
## shrub        1.13560    0.17299   6.565 5.21e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
  summary(fit.replicated) # replicated (dependent) data
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: indiv.data7.replicated
## 
##      AIC      BIC   logLik deviance df.resid 
##   4963.9   4990.7  -2478.0   4955.9     5996 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -2.33068    0.10605 -21.978   <2e-16 ***
## dist.dev    -0.69055    0.03768 -18.326   <2e-16 ***
## forest      -0.29470    0.12541  -2.350   0.0188 *  
## shrub        1.14786    0.12909   8.892   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Okay, so this is not so good. We can see the standard errors of the estimates from the model `fit.replicated` are smaller because we inappropriately doubled our sample size and treated dependent data as independent. This effect will translate into much too small confidence intervals for our estimates and predictions. It also means that our p-values are too small. Notice that our coefficient for forest is now statistically clearly different than zero (p-value of 0.0187789) while it was not so using only independent data (p-value of 0.0838469).

<br>

### The data generating model and the model being fit

Context only.

<br>

### The model being approximated

It is the responsibility of each researcher to make sure their modeling process is done such that they are approximating the true underlying Inhomogeneous Poisson Point process (IPP) model well. To demonstrate this, we will consider how estimated coefficients change with increasing numbers of available samples. We will also be able to see how the intercept changes, due its partial (mostly) interpretation as a measure of the ratio of used to available samples.

In practice, there are two components that need to be considered in the approximation. First, is how to sample the complexity of variation in the spatial covariates being considered. Simply, you want to make sure that the spatial variation is captured with these discrete locations (points placed within the defined spatial region of what is available that will be the available sample and respective covariate values) by spreading them out well. Second, is to have enough of these samples such that the approximation works well. Capturing the spatial variation is best done by making a systematic grid of locations within the area that is considered available and extracting the covariate values. Having enough available locations can be achieved by having a high density of points for the systematic grid and also by weighting the available locations within the model fitting process.

Since we are not explicitly spatially sampling here, we can't directly evaluate how to capture spatial variation, but we can demonstrate how to weight the available locations. Essentially, you want to tell the model that for every used location to weight this data by 1, indicating that there is one of them. For the available locations, you want to weight each location by a large number (e.g.,1000). This tells the model to consider this data as having 1000 of them. You can think of it as a shortcut to replicating this available location and its covariate values 1000 times. Either way works, but weighting is more computationally efficient. By default, weighting should be done in all analyses to help with reducing approximation error. Note that this type of replicating of the available sample is a good thing, while replicating the used data (as seen in Critical Assumption 2: Independence of location data) is not a good thing.


``` r
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
```


```{.r .fold-hide}
par(mfrow=c(2,1))
plot(coef.save$N.Avail, coef.save$beta1,lwd=3,type="l",col=2,ylim=c(-1.5,1.5),
     main="Slope Coefficients",
     xlab="Number of Available Samples",ylab="Coefficient Estimate")
lines(coef.save$N.Avail, coef.save$beta2,lwd=3,col=2)
lines(coef.save$N.Avail, coef.save$beta3,lwd=3,col=2)
plot(coef.save$N.Avail, coef.save$Intercept,lwd=3,type="l",col=1,main="Intercept",
     xlab="Number of Available Samples",ylab="Coefficient Estimate")
```

![](TraditionalHSF_files/figure-html/plotting.sensitivity-1.png)<!-- -->

We can see that our estimates of the slope coefficients are sensitive to the number of available locations when there are less than a few thousand available locations; there is still a bit of jumping around until the available samples are above 4000. In contrast, we can see the intercept just keeps getting smaller and smaller because our available sample size is getting larger and larger. The intercept is not biologically meaningful.

Another way to look at this sensitivity is by showing the variability of the coefficients within and across different sizes of available samples.


``` r
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
```


```{.r .fold-hide}
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
```

![](TraditionalHSF_files/figure-html/plot.sensitive.org2-1.png)<!-- -->

We see more easily in this plot that there is high variability in the estimated coefficients when the available sample is small. This variability decreases as the available sample grows, showing how the estimates are stabilizing (except the intercept). What this means is that you could grab many different but large available samples and would get almost the same estimated coefficients. This is not true when using a small available sample; your estimated coefficients will vary and we do not want that to happen. The intercept estimates decrease in the variability but also continues to decline. It is up to each researcher to conduct a similar sensitivity analysis to see how many available samples are necessary for the coefficient estimates to stabilize.

<br>

### Individuals

We should \textit{a priori} assume there is individual variation in habitat selection estimates. This variation is masked when data are pooled. Lets consider the implications of pooling all data vs fitting a model separately to each individual.


``` r
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
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: sims2
## 
##      AIC      BIC   logLik deviance df.resid 
##  58028.7  58067.1 -29010.4  58020.7   109996 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -3.26310    0.03319  -98.33  < 2e-16 ***
## dist.dev    -0.97998    0.01167  -83.96  < 2e-16 ***
## forest      -0.19943    0.03825   -5.21 1.85e-07 ***
## shrub        1.05867    0.03942   26.86  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

In the pooled data, we see that the standard errors of the coefficients and p-values are very small. We are using a lot of information, assuming no variation among individuals, and thus assuming every individual's effects can be characterized by these estimates. We are also forcing the intercept to be the same for each individual. That's fine here because we have the same number of used locations for each individual and thus the ratio of used to available is the same. But, if you have differing number of used sampled per individual, this would not be a good thing. Now let's, look at just two individuals' estimates when fitting separate models.


``` r
summary(indiv.fit[[2]])
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: x
## 
##      AIC      BIC   logLik deviance df.resid 
##   3009.6   3036.0  -1500.8   3001.6     5496 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -3.50477    0.15003 -23.360  < 2e-16 ***
## dist.dev    -0.81633    0.05057 -16.141  < 2e-16 ***
## forest       0.61806    0.17213   3.591  0.00033 ***
## shrub        0.97917    0.17399   5.628 1.83e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

First, notice that the estimated effects are a bit different than the pooled estimates. Importantly, also notice that the standard errors and p-values are larger.


``` r
summary(indiv.fit[[20]])
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub
## Data: x
## 
##      AIC      BIC   logLik deviance df.resid 
##   2817.6   2844.1  -1404.8   2809.6     5496 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -3.4354     0.1518 -22.634  < 2e-16 ***
## dist.dev     -1.1067     0.0542 -20.418  < 2e-16 ***
## forest        0.1244     0.1709   0.728    0.467    
## shrub         0.8778     0.1781   4.929 8.26e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Now lets look at all the estimates together.


``` r
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
```


```{.r .fold-hide}
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
```

![](TraditionalHSF_files/figure-html/plot.pool.separate-1.png)<!-- -->

The plot on the right are the pooled estimates for beta1 (black), beta2 (red), and beta3 (green). The plot on the left are each individual's estimates (colored the same). By ignoring individual variation, we are much too confident in our estimates of uncertainty and are ignoring a lot of clear variation by individual. The pooled estimate's certainty has to do with psuedoreplication - we are treating a subunit (each location) as the main unit of replication. Our true unit of replication is the individual.

Since our sample sizes for each individual are equal, we see that the pooled estimates generally relate to the average of each estimate across all individuals. When the number of used locations varies by individual this won't be the case. The individuals with more used locations will disproportionately influence the pooled estimates.

<br>

#### Sample size

A common question is how many used locations for a given individual (or pooled across individuals) are needed to provide statistical clarity about a habitat selection variable? We can determine this using the methods of [Street et al. 2021](https://doi.org/10.1111/2041-210X.13701); their data/code can be found at [figshare](https://figshare.com/articles/dataset/Datasets_and_Code_zip/11910831).

We will explicitly use spatial data in this section. First, lets grab one spatial covariates that was already loaded.


``` r
S = covs[[1]]
values(S) = scale(values(S))  
```


Now, we can define our intercept and slope/effect for this spatial covariate in determining the true habitat selection values.


``` r
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
```

![](TraditionalHSF_files/figure-html/power2-1.png)<!-- -->

You can think of this covariate as the Normalized Difference Vegetation Index (NDVI) if that makes it easier to visualize. To determine the minimum number of used samples for the covariate and habitat selection values (with coefficient of $\beta[2]$), we need to define our criteria for statistical clarity. Specifically, we are interested in finding N:

"Recall that $N_{\alpha, p}(\beta)$ was defined to be the minimum number of samples, N, required so that we expect to reject the null hypothesis of $\beta = 0$ (vs. alternative $\beta \neq 0$), at a significance level of $p$, in $100(1 − \alpha)\%$ of experiments". [Street et al. 2021 (Supplementary Appendix A)](https://doi.org/10.1111/2041-210X.13701)


``` r
# Define Inputs
  alpha = 0.05
  p_val = 0.05
```

We are now ready to determine the number of used locations required to achieve our goal.


``` r
N.used = sample.size.used.locs(alpha = alpha, 
                               p_val = p_val,
                               HSF.True = S.HSF,
                               S = S,
                               beta = beta[2]
                               )
N.used$Nalphapbetas
```

```
## [1] 353.5197
```

Given the spatial heterogeneity in our covariate, a coefficient of -0.2 and the above levels of statistical clarity, we only need 353.5196997 used locations.

Lets consider varying our Type I error rate ($\alpha$) to see how our sample size changes. We will hold all other variables constant.


``` r
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
```


```{.r .fold-hide}
  plot(alpha,N.used,col=2,type="l",lwd=3,xlab="Type I Error Rate",ylab="Sample Size of Used Locations Needed")
```

![](TraditionalHSF_files/figure-html/power.plot1-1.png)<!-- -->

We we might expect, as we relax the Type I error rate, it is easier to determine statistical clarity at lower sample sizes of used locations.

Next, lets reset $\alpha = 0.05$ and evaluate our sample size requirements as we change the effect size (i.e., coefficient) of our spatial variable


``` r
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
```


```{.r .fold-hide}
  plot(beta.combs,N.used,col=2,type="b",lwd=3,xlab="Coefficient",ylab="Sample Size of Used Locations Needed")
```

![](TraditionalHSF_files/figure-html/power.plot2-1.png)<!-- -->

We see that as the coefficient gets closer to zero, we require more samples to determine that it is statistically clearly different from zero. This should hopefully make some logical sense.

Last, we will examine how the heterogeneity of the spatial covariate affects our required sample size. We will specifically do this by changing the values of the spatial covariate in terms of the standard deviation but will measure it in terms of 'landscape complexity', as one in [Street et al. 2021](https://doi.org/10.1111/2041-210X.13701). Note that we will return to a slope coefficient of ($\beta = -0.2$).


``` r
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
```


```{.r .fold-hide}
  N.used.plot = data.frame(N.used = N.used,
                           N.variance = N.variance,
                           sd = rep(sd.combs,each=n.reps)
                          )
  N.used.plot = N.used.plot[order(N.used.plot$N.variance),]
plot(N.used.plot$N.variance,N.used.plot$N.used,lwd=3,col=2,type="l",xlab="Landscape Complexity",ylab="Sample Size of Used Locations Needed")
```

![](TraditionalHSF_files/figure-html/power.plot3-1.png)<!-- -->

We see how the number of used locations declines with increasing landscape complexity. As there is more variation on the landscape, we have more statistical power at lower sample sizes. Intuitively, as the landscape becomes more heterogeneous, it should be easier (i.e., fewer used samples needed) to identify their relative selection of each habitat type. 

<br>

### Population inference

If we are interested in obtaining inference to the individual- and population-level effects, we can consider bootstrapping the results from individual models or jointly fitting a single model with a random effect across individuals. We will first consider the bootstrapping method. Second, we will demonstrate the use of random intercepts and slopes, as well as how to fix the variance of the random intercept so there is no shrinkage.

<br>

#### Bootstrapping

Bootstrapping is a [resampling technique](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) used often in statistics to obtain a sampling distribution of possible values. We will resample the individual-level coefficents to estimate the population-level mean and variation.  

The core functions for bootstrapping are adopted from Dr. Guillaume Bastille-Rousseau's github repository for the R package [IndRSA](https://github.com/BastilleRousseau/IndRSA/).

We first need to get all estimated coefficients from each individual.


``` r
  coef.list = lapply(indiv.fit,FUN=function(x){fixef(x)[[1]]})
  
  coef.df=do.call(rbind.data.frame, coef.list)
  colnames(coef.df)=names(fixef(indiv.fit[[1]])[[1]])

# Remove intercept
  coef.df = coef.df[,-1]

# Add name (could be modified to keep track if animals have an ID)
  coef.df$name = 1:nrow(coef.df)

# Frequency could be greater than one if there were multiple estimates of a single individual, perhaps across years
  coef.df$freq = 1
```

We are ready to bootstrap the coefficient estimates.


``` r
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
```

We now want to summarize our bootstrapped results. The population mean and 95% lower and upper confidence intervals are outputted.


``` r
  boot.pop.esimates=mean_ci_boot(boot)
  knitr::kable(boot.pop.esimates)
```



|         |      Mean|        LCI|        UCI|
|:--------|---------:|----------:|----------:|
|dist.dev | -0.962506| -1.1033896| -0.8794034|
|forest   | -0.233231| -0.6817081|  0.2474493|
|shrub    |  1.044056|  0.9049891|  1.2042720|

<br>

#### Random effects across individuals

##### Random intercept only model

The most common use of random effects in habitat selection analysis is the use of a random intercept. As described in the manuscript, this does not account for variation in the slope coefficients, which is problematic. Here, we will allow the intercepts to vary, which accommodates different sample sizes and thus different ratios of used to available samples per each individual (not a major issue in this contrived setting of all individuals having the same sample size). However, if we allowed the random effect variance to be estimated we would get shrinkage of the individual-level intercept estimates towards the mean, which we do not want in this case. Therefore, we will fix the random effect variance to a large value to ensure there is no shrinkage. We can do this using a step wise approach, 1) setting up the model and not fitting it, 2) coding the parameter as `NA` and 3) then lastly fitting the model via optimization.  

While we are using common language in describing the models and model fitting process, note that there are very different meanings of `fixing a parameter` (i.e., random intercept variance) versus a `fixed effect`. The former refers to setting a parameter at a given value and is not estimated. The later refers to a description of a type of parameter that is estimated. 


``` r
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
```


``` r
# The fixed components - these are the population level (across individual) effects. 
  fixef(re_int)
```

```
## 
## Conditional model:
## (Intercept)     dist.dev       forest        shrub  
##  -7.418e-05   -9.820e-01   -2.011e-01    1.060e+00
```

``` r
# The random components- these are the effect differences for each individual from the population mean estimate (referred to as the fixed effect in the 'conditional model' statement)
  ranef(re_int)
```

```
## $ID
##    (Intercept)
## 1    -3.277075
## 2    -3.237661
## 3    -3.220824
## 4    -3.326315
## 5    -3.261974
## 6    -3.304170
## 7    -3.217700
## 8    -3.237493
## 9    -3.233313
## 10   -3.270460
## 11   -3.323575
## 12   -3.234716
## 13   -3.249997
## 14   -3.274132
## 15   -3.225749
## 16   -3.319067
## 17   -3.243254
## 18   -3.234961
## 19   -3.330002
## 20   -3.267983
```

Notice the lack of variability in these estimated individual differences from the overall mean intercept. The reason for this is because we simulated data with the same number of used and available locations, therefore the ratio of used to available is the same. These estimates are likely to be  more variable when you have different used sample sizes per individual (total individual animal locations) and are applying the same ratio of available samples. There are many reasons individuals end up having a different number of locations, such as battery life. 


``` r
# Summarize results  
  summary(re_int)
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ dist.dev + forest + shrub + (1 | ID)
## Data: sims2
## 
##      AIC      BIC   logLik deviance df.resid 
##  58414.3  58452.7 -29203.2  58406.3   109996 
## 
## Random effects:
## 
## Conditional model:
##  Groups Name        Variance Std.Dev.
##  ID     (Intercept) 1e+06    1000    
## Number of obs: 110000, groups:  ID, 20
## 
## Conditional model:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -7.418e-05  2.236e+02    0.00        1    
## dist.dev    -9.820e-01  1.169e-02  -83.99  < 2e-16 ***
## forest      -2.011e-01  3.829e-02   -5.25 1.52e-07 ***
## shrub        1.060e+00  3.944e-02   26.88  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**A random intercepts only model is not recommended.** This is because we are not dealing with the expected variation of responses to hypothesized spatial factors across individuals. 


<br>

##### Alternative to Random intercept only model

The use of random effects is really not particularly useful. We are really looking to estimate slopes across all individuals (pooling) while allowing the intercepts by individual to vary to accommodate different sample sizes. We can accomplish the same thing using a fixed-effect only model and thus not need to fix any random effect variance. Specifically, we can do that by setting up a design matrix using contrast sums, also know as effecting coding.


``` r
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
```

Now lets compare the estimated effects from `fix_int` and `re_int`. First lets look at the slopes for `re_int`... 


``` r
kable(summary(re_int)[[6]][[1]][2:4,])
```



|         |   Estimate| Std. Error|    z value| Pr(>&#124;z&#124;)|
|:--------|----------:|----------:|----------:|------------------:|
|dist.dev | -0.9820027|  0.0116923| -83.987436|              0e+00|
|forest   | -0.2010643|  0.0382935|  -5.250612|              2e-07|
|shrub    |  1.0600931|  0.0394368|  26.880794|              0e+00|

And here we have the same estimates for `fix_int`, 


``` r
kable(summary(fix_int)[[6]][[1]][(n.indiv+1):(n.indiv+1+2),])
```



|         |   Estimate| Std. Error|    z value| Pr(>&#124;z&#124;)|
|:--------|----------:|----------:|----------:|------------------:|
|dist.dev | -0.9817199|  0.0116898| -83.980850|              0e+00|
|forest   | -0.2010046|  0.0382881|  -5.249799|              2e-07|
|shrub    |  1.0597804|  0.0394311|  26.876733|              0e+00|

Notice that they are practically the same. The intercepts are the same as well, but we need to reorganize the results a tad to make our comparison...


``` r
int.re = ranef(re_int)[[1]][[1]]

int.fe = c(summary(fix_int)[[6]][[1]][1,1]+
           summary(fix_int)[[6]][[1]][2:20],
           summary(fix_int)[[6]][[1]][1,1]
           )
kable(round(data.frame(individual = 1:20,int.fe=int.fe,int.re=int.re),digits=2))
```



| individual| int.fe| X.Intercept.|
|----------:|------:|------------:|
|          1|  -3.28|        -3.28|
|          2|  -3.24|        -3.24|
|          3|  -3.22|        -3.22|
|          4|  -3.33|        -3.33|
|          5|  -3.26|        -3.26|
|          6|  -3.30|        -3.30|
|          7|  -3.22|        -3.22|
|          8|  -3.24|        -3.24|
|          9|  -3.23|        -3.23|
|         10|  -3.27|        -3.27|
|         11|  -3.32|        -3.32|
|         12|  -3.23|        -3.23|
|         13|  -3.25|        -3.25|
|         14|  -3.27|        -3.27|
|         15|  -3.23|        -3.23|
|         16|  -3.32|        -3.32|
|         17|  -3.24|        -3.24|
|         18|  -3.23|        -3.23|
|         19|  -3.33|        -3.33|
|         20|  -3.26|        -3.27|
We get the same estimates, again. The point is that the random effect is not being used in the manner we usually might use it. Choosing to use a random effect with a large fixed variance (i.e., do not estimate). or fixed effects can be chosen based on convenience. There may be some convenience to the random effect process when there are many individuals. 

<br>

#### Random intercept and slopes model

This is the recommended model structure using random effects.


``` r
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
```

Let's look at the results


``` r
  fixef(re_int_slopes)
```

```
## 
## Conditional model:
## (Intercept)     dist.dev       forest        shrub  
##   -0.001188    -0.983164    -0.194497     1.040308
```

``` r
  broom.mixed::tidy(re_int_slopes, effects = "fixed", conf.int = TRUE)[-1,c(3,4,8,9)]
```

```
## # A tibble: 3 × 4
##   term     estimate conf.low conf.high
##   <chr>       <dbl>    <dbl>     <dbl>
## 1 dist.dev   -0.983   -1.06     -0.903
## 2 forest     -0.194   -0.577     0.188
## 3 shrub       1.04     0.951     1.13
```

Here are estimated population-level coefficients (across individual-level effects). Compared to the bootstrapped results above, they are generally similar. We should not expect them to be the same, as the random effect model shares information across individuals and the bootstrapped estimates do not.

We can also extract the estimated difference by individual from the population level coefficients.


``` r
  ranef(re_int_slopes)
```

```
## $ID
##    (Intercept)    dist.dev      forest       shrub
## 1   -10.299348 -0.11130168 -0.02090383  0.06758999
## 2   -10.414969  0.15779918  0.76843245 -0.02645743
## 3    -9.882139  0.28789052 -0.23574983  0.12923279
## 4   -10.570038 -0.25249448  0.18065484  0.03798059
## 5   -10.648244  0.02743712  0.97919430 -0.01236556
## 6    -9.801698 -0.10886945 -1.06614415 -0.02520669
## 7    -9.897331  0.28606535 -0.09030067  0.03004482
## 8    -9.755803  0.11628282 -0.56285305 -0.04205768
## 9   -10.129800  0.13680787  0.08445766  0.11262520
## 10  -10.556768  0.08693897  0.96851949 -0.12290572
## 11  -10.399506 -0.30346128 -0.11944704 -0.09799560
## 12  -10.015957 -0.11198466 -0.35259407 -0.05684952
## 13  -11.130693  0.15244917  1.82792587  0.08317337
## 14   -9.529007 -0.03325338 -1.41168291 -0.08380647
## 15   -9.587807  0.08392103 -1.05863919 -0.01054412
## 16  -10.633384 -0.12491422  0.52343503  0.05682853
## 17   -9.570631  0.05667112 -1.27550215  0.04337108
## 18  -10.761501  0.06733219  1.12346753  0.12540716
## 19  -10.212692 -0.31274622 -0.56843852 -0.11791448
## 20  -10.361531 -0.10115549  0.31052304 -0.09381608
```

``` r
  broom.mixed::tidy(re_int_slopes, effects = "ran_vals", conf.int = TRUE)[-c(1:20),c(5,6,8,9)]
```

```
## # A tibble: 60 × 4
##    term     estimate conf.low conf.high
##    <chr>       <dbl>    <dbl>     <dbl>
##  1 dist.dev  -0.111  -0.227     0.00449
##  2 dist.dev   0.158   0.0453    0.270  
##  3 dist.dev   0.288   0.174     0.402  
##  4 dist.dev  -0.252  -0.365    -0.140  
##  5 dist.dev   0.0274 -0.0838    0.139  
##  6 dist.dev  -0.109  -0.222     0.00379
##  7 dist.dev   0.286   0.173     0.399  
##  8 dist.dev   0.116   0.00245   0.230  
##  9 dist.dev   0.137   0.0232    0.250  
## 10 dist.dev   0.0869 -0.0243    0.198  
## # ℹ 50 more rows
```

Lastly, we can see a full summary of the model.


``` r
summary(re_int_slopes)
```

```
##  Family: binomial  ( logit )
## Formula:          
## status ~ dist.dev + forest + shrub + (1 | ID) + (0 + dist.dev |  
##     ID) + (0 + forest | ID) + (0 + shrub | ID)
## Data: sims2
## Weights: weight
## 
##      AIC      BIC   logLik deviance df.resid 
## 193336.4 193403.7 -96661.2 193322.4   109993 
## 
## Random effects:
## 
## Conditional model:
##  Groups Name        Variance  Std.Dev. 
##  ID     (Intercept) 1.000e+06 1000.0000
##  ID.1   dist.dev    3.114e-02    0.1765
##  ID.2   forest      7.368e-01    0.8584
##  ID.3   shrub       1.619e-02    0.1273
## Number of obs: 110000, groups:  ID, 20
## 
## Conditional model:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -0.001188 223.621785   0.000    1.000    
## dist.dev     -0.983164   0.040690 -24.163   <2e-16 ***
## forest       -0.194497   0.195184  -0.996    0.319    
## shrub         1.040308   0.045784  22.722   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Note the variance and standard deviation estimates (fourth column of under the "Conditional Model" section), indicating how variable coefficients are across individuals. We see that `forest` is estimated as the most variable. This corresponds well to how the data were generated with forest coefficients having the highest standard deviation of `beta.sd` which was set as 1.

Another thing to notice is that the population level effect for forest is not statistically clearly different than zero. Thus, at the population level, percent forest is being selected in proportion to what is available to each individual. Let's dive into this a bit more. Let's plot the individual estimates along with the population-level estimate.


``` r
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
```


```{.r .fold-hide}
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
```

![](TraditionalHSF_files/figure-html/RE.plot.forest.indiv.pop-1.png)<!-- -->

Our plot shows the individual effects of `forest` along with the population-level. Compare this output and inference to the figure (right-side) for the pooled model fit at the end of Section 3.11.  In the pooled case we would conclude that at the population level, this species avoids forest (i.e., negative coefficient), and we would say that with some level of statistical certainty because the confidence intervals do not include zero. Here we get a more nuanced picture of what’s happening, and that leads to a very different conclusion.

What is clear from this plot is that the reason there is no statistically clear difference of the population-level effect from zero is because there is a wide range of effects across individuals. Some individuals have positive coefficients and some have negative. Since these are generally equal, they balance out to a population-level effect of zero. The story is more complicated!

The point of this model is that we have effects per each individual and at the population-level for each estimated effect/slope. We also are accommodating different sample sizes with the random intercepts. However, one could use the trick from above to allow intercepts to vary by individual using fixed effects and then only use random effects for the spatial variables. 

<br>

#### Sample size

An important question to ask is how many individuals should be tracked to be able to provide statistical clarity about a population-level coefficient for a habitat variable? We can determine this using the methods of [Street et al. 2021](https://doi.org/10.1111/2041-210X.13701); their data/code can be found at [figshare](https://figshare.com/articles/dataset/Datasets_and_Code_zip/11910831). For details see the equation 8.

To evaluate this question, we need to consider the population level mean of the coefficient, variation across individuals, landscape complexity, and the number of used locations per individual. For simplicity, we will use a single covariate (`covs[[1]]`) for each individual;.


``` r
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
```

```
## [1] "M (tracked individuals) should be greater than or equal to  85"
```

Lets consider how the number of individuals needed changes with different population-mean coefficients.


```{.r .fold-hide}
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
```


``` r
plot(beta_mu,min.indiv,ylab="Minimum Number of Individuals Tracked",lwd=3,col=4,type="b")
```

![](TraditionalHSF_files/figure-html/power.pop.plot-1.png)<!-- -->

We can see that as the population mean gets larger it becomes easier to have statistical clarity with fewer individuals tracked.

Next, lets consider how the number of individuals needed changes with the changing the variation of the coefficient across indivduals (i.e, `beta_s` (standard deviation) and `beta_s2` (variance). We will reset the population mean back to 0.2.


```{.r .fold-hide}
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
```


``` r
plot(beta_s2,min.indiv,ylab="Minimum Number of Individuals Tracked",lwd=3,col=4,type="b")
```

![](TraditionalHSF_files/figure-html/power.pop.plot2-1.png)<!-- -->

Our results agree with the statement from [Street et al. 2021](https://doi.org/10.1111/2041-210X.13701) that "Note that [Equation 8] is a non-decreasing function of $\sigma^2$, meaning that more variation among individuals is likely to mean one has to sample a higher number of animals.

<br>

### Considering context in habitat selection analyses

In this section, we discuss ways of considering behavior in habitat selection analyses. To demonstrate the effect of ignoring behavior, we will simulate data where selection is different for two sets of animal locations. Imagine an animal is highly selective of being in forest landcover when it is resting, but when foraging (and otherwise) it selects against forest in proportion to its availability.


```{.r .fold-hide}
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
```

Let's fit a model where we ignore the behavioral differences in selection and fit a single model with all the data.


``` r
  model.ignore.behav = glmmTMB::glmmTMB(status ~ forest, 
                                        family = binomial(), 
                                        data = data.ignore 
                                        )
  summary(model.ignore.behav)
```

```
##  Family: binomial  ( logit )
## Formula:          status ~ forest
## Data: data.ignore
## 
##      AIC      BIC   logLik deviance df.resid 
##   2684.0   2696.8  -1340.0   2680.0     4398 
## 
## 
## Conditional model:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -2.34189    0.06949  -33.70   <2e-16 ***
## forest       0.09334    0.10592    0.88    0.378    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

We see that the estimated forest coefficient is not near the true values of -2 or 2. It's somewhat in between - near zero. Essentially, when animals are selecting different habitat features because of behavior, mixing across behaviors can lead to a muddled inference.

<br>

### Interpreting coefficients and predicting

Interpreting coefficients and predicting is outlined above in subsection `3.7. What are habitat-selection functions used for?`.

<br>

### Model selection

If focused on inference, fitting a single model is not only okay, but desirable.

<br>

### Concluding remarks

We hope this vignette provided some utility. If so, let us know with an email ([brian.gerber\@colostate.edu](mailto:brian.gerber@colostate.edu){.email}).

Have a nice day.

<br>

## Software

This report was generated from the R Statistical Software (v4.4.2; R Core Team 2021) using the [Markdown language](https://www.markdownguide.org/) and [RStudio](https://posit.co/products/open-source/rstudio/). The R packages used are acknowledged below.


|Package           |Version |Citation                                       |
|:-----------------|:-------|:----------------------------------------------|
|amt               |0.3.0.0 |@amt                                           |
|base              |4.4.2   |@base                                          |
|broom.mixed       |0.2.9.6 |@broommixed                                    |
|circular          |0.5.1   |@circular                                      |
|geoR              |1.9.4   |@geoR                                          |
|glmmTMB           |1.1.10  |@glmmTMB                                       |
|knitr             |1.49    |@knitr2014; @knitr2015; @knitr2024             |
|plotrix           |3.8.4   |@plotrix                                       |
|raster            |3.6.30  |@raster                                        |
|remotes           |2.5.0   |@remotes                                       |
|ResourceSelection |0.3.6   |@ResourceSelection                             |
|Rfast             |2.1.0   |@Rfast                                         |
|rmarkdown         |2.29    |@rmarkdown2018; @rmarkdown2020; @rmarkdown2024 |
|sp                |2.1.4   |@sp2005; @sp2013                               |
|survival          |3.7.0   |@survival-book; @survival-package              |
|tidyverse         |2.0.0   |@tidyverse                                     |


