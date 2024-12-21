<style type="text/css">
body, td {
   font-size: 14px;
}
code.r{
  font-size: 16px;
}
pre {
  font-size: 16px
}
</style>

## Introduction

This vignette is associated with the manuscript titled, ‘A plain
language review and guidance for modeling animal habitat-selection’. We
will be demonstrating some of the points made in the manuscript on
fitting models to estimate parameters in a movement-based habitat
selection function (HSF) or step-selection functions (SSF)). The scale
of interest is the the fourth-order of selection (e.g. selection along a
track within a home range).

*Note:* some of the code-chunks are not automatically displayed. To show
the code, select ‘show’.

## Setup

### Environment

First, load required R packages, source two local functions, and load an
object containing spatial covariates.

``` r
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
  source("sim.ind.movement.hsf.r")

# Source bootstapping function
  source("mean_ci_boot.r")

# Load spatial covariates stored in a save R object
  load("Covs")
```

### Simulate data

We will consider a habitat selection analysis of individuals along a
movement track within a landscape (e.g., fourth-order of selection).
Individual’s effects are assumed to come from a distribution of effects
that can be characterized by a mean and standard deviation (i.e., random
effect).

``` r
# Number of Sampled individuals (e.g., tracked via GPS telemetry)
  n.indiv = 20

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
```

``` r
  par(mfrow=c(3,1))
  hist(betas[,1],xlab=bquote(beta[1]),xlim=c(-3,3),main="True Individual Coeficient Values",breaks=10,freq=FALSE)
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
```

![](MovementHSF_files/figure-markdown_github/visualize%20true%20coeficients-1.png)

We have pre-loaded spatial covariate data. Lets look at the three
covariates we will consider. The covariates included are the Euclidean
distance to nearest development (`dist.dev`), the percent forest cover
(`forest`), and the percent shrub cover (`shrub`). Each has been
standardized to a mean of 0 and standard deviation of 1; this is a
common procedure in linear modeling to help optimization algorithms
converge to maximum likelihood estimates and to put each covaraite on
the same scale (1 standard deviation of the covariate) so that
coefficient magnitudes can be compared.

``` r
# Combine into 1 raster stack
  covs=stack(covs[[1]], covs[[2]],covs[[3]]) 

# Change the extent to be larger to accommodate more realistic movements 
#  (not required but makes me feel better)
  extent(covs)=c(0,1000,0,1000)
  
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
```

![](MovementHSF_files/figure-markdown_github/covs.plot-1.png)

Now, we will use the spatial covariates (`covs`) and true
individual-level coefficients (`betas`) to create the linear
combinations to create the true individual HSF (as a raster).

``` r
  hsf.true = vector("list",n.indiv)
  for (i in 1:n.indiv){
    hsf.true[[i]] = (covs[[1]]*betas[i,1]+
                     covs[[2]]*betas[i,2]+
                     covs[[3]]*betas[i,3]
                     )
  }
```

We are now ready to simulate individual-level data. We will do this with
the `sim.ind.movement.hsf` function that is sourced from the file with
the same name.

``` r
# Set number of available samples per used
  n.avail = 100

# Number of movements
  n.time = 400  

# Number of possible steps to choose from at each iteration- these are not available locations, this is for the simulation of the movement path  
  n.step = 400  
  
# Population (across individual) turning angle parameters for von mises distribution
  angle.mean = 0
  angle.kappa = 0.00001  

# Step length rate of exponential distribution
  step.rate = 0.2
  
# Simulate individual-level movement-based habitat selection data
  sims =  sim.ind.movement.hsf(hsf.true=hsf.true,
                               n.time=n.time,
                               n.step=n.step,
                               n.avail=n.avail,
                               angle.mean=angle.mean,
                               angle.kappa=angle.kappa,
                               step.rate=step.rate
                               )

# There are two list elements
# There is one list element for each individual in each of the two elements 
    names(sims)
```

    ## [1] "indiv" "locs"

``` r
# Individual 1 data    
  head(sims$indiv[[1]])
```

    ##     status strata   dist.dev
    ## 1        1      1 -0.5773911
    ## 401      0      1 -0.5933174
    ## 402      0      1  0.4770586
    ## 403      0      1 -0.5773911
    ## 404      0      1 -0.5773911
    ## 405      0      1 -0.5773911
    ##       forest      shrub
    ## 1   1.218153 0.16764176
    ## 401 1.287160 0.01988518
    ## 402 1.546580 0.34669398
    ## 403 1.218153 0.16764176
    ## 404 1.218153 0.16764176
    ## 405 1.218153 0.16764176
    ##           step
    ## 1    0.8705149
    ## 401  3.6848288
    ## 402 10.4154689
    ## 403  2.5056502
    ## 404  0.2987808
    ## 405  0.1996507

``` r
  table(sims$indiv[[1]]$status)
```

    ## 
    ##     0     1 
    ## 40000   400

We see from this individual’s data frame that we have the columns
status, strata, dist.dev, forest, shrub, step. The column `status` is
our response variable. A **1** indicates an animal location (used point)
and a **0** indicates a potentially used location (available point).
Each used location is paired with 100 available locations and are
matched by the strata number (column `strata`). Each 1 and 0 have
corresponding spatial variables, which we have called `dist.dev`,
`forest`, and `shrub`. Each covariate has been standardized to have a
mean of 0 and standard deviation of 1 (also called Normalizing). Note
that we have a large available available sample matched with each used
location. This is correct. We need a large available sample (the 0’s)
for each used location to approximate the true underlying
spatio-temporal habitat selection model (specified via a weighted
distribution; see Equation 1 of [Michelot et
al. 2024](https://doi.org/10.1111/2041-210X.14248)).

Let’s look at one individual’s tracks on top of the true HSF.

``` r
#Plot the true HSF with locations for individual k
  k=10
  par(mfrow=c(1,1))
  plot(hsf.true[[k]],xlim=c(300,600),ylim=c(400,800))
  points(sims$locs[[k]]$use.x,
         sims$locs[[k]]$use.y,
         pch=16
         )
```

![](MovementHSF_files/figure-markdown_github/plot.tracks-1.png)

We now want to package are individual simulated data into a single data
frame to be able to fit all individuals together. We will use the
objects `sims` and `sims2` below when we want to fit models to each
individual separately and together, respectively. For `sims2`, we also
need a new set of strata identifiers so that they are unique for each
individual.

``` r
# Combine data into a single data.frame
  sims2 = do.call("rbind", sims$indiv)
  head(sims2)
```

    ##     status strata   dist.dev
    ## 1        1      1 -0.5773911
    ## 401      0      1 -0.5933174
    ## 402      0      1  0.4770586
    ## 403      0      1 -0.5773911
    ## 404      0      1 -0.5773911
    ## 405      0      1 -0.5773911
    ##       forest      shrub
    ## 1   1.218153 0.16764176
    ## 401 1.287160 0.01988518
    ## 402 1.546580 0.34669398
    ## 403 1.218153 0.16764176
    ## 404 1.218153 0.16764176
    ## 405 1.218153 0.16764176
    ##           step
    ## 1    0.8705149
    ## 401  3.6848288
    ## 402 10.4154689
    ## 403  2.5056502
    ## 404  0.2987808
    ## 405  0.1996507

``` r
  dim(sims2)
```

    ## [1] 808000      6

``` r
# Create ID vector for individual's data
  ID=rep(1:n.indiv,each=nrow(sims$indiv[[1]]))
  sims2$ID=ID

  
# Create new strata
  sims2$indiv.id.strata=unique(sims2$ID+sims2$strata*100)

# The number of unique identifiers needed
# n.step*n.indiv

  head(sims2)
```

    ##     status strata   dist.dev
    ## 1        1      1 -0.5773911
    ## 401      0      1 -0.5933174
    ## 402      0      1  0.4770586
    ## 403      0      1 -0.5773911
    ## 404      0      1 -0.5773911
    ## 405      0      1 -0.5773911
    ##       forest      shrub
    ## 1   1.218153 0.16764176
    ## 401 1.287160 0.01988518
    ## 402 1.546580 0.34669398
    ## 403 1.218153 0.16764176
    ## 404 1.218153 0.16764176
    ## 405 1.218153 0.16764176
    ##           step ID
    ## 1    0.8705149  1
    ## 401  3.6848288  1
    ## 402 10.4154689  1
    ## 403  2.5056502  1
    ## 404  0.2987808  1
    ## 405  0.1996507  1
    ##     indiv.id.strata
    ## 1               101
    ## 401             201
    ## 402             301
    ## 403             401
    ## 404             501
    ## 405             601

``` r
  dim(sims2)  
```

    ## [1] 808000      8

## Manuscript sections

We are now ready to fit models to estimate our movement-based habitat
selection parameters and demonstrate points made in the manuscript. We
have organized the sections below to match with the sections of the
manuscript.

### 1. What is habitat?

Definitions only.

### 2. What is habitat selection?

Definitions only.

### 3. What is a habitat-selection function?

Definitions only.

We want to remind the reader of the nuance between a habitat selection
function and the statistical model we will be fitting. The model we are
fitting (as coded) and the true underlying model (weighted distribution
of spatio-temporal model) is not the habitat selection function. The
habitat selection function is a component of the model. We will
approximate the triue model in two different ways: 1) using a
conditional logistic regression model and 2) using a Poisson regression
model with stratum-specific fixed intercepts ([Muff et
al. 2019](https://doi.org/10.1111/1365-2656.13087)). Each approach is
used to estimate parameters within the movement-based habitat selection
function. We will use the same data set to show the equivalence of the
estimated parameters.

### 4. What about scale?

Definitions only.

In our data setup and analysis we will be estimating habitat selection
at the fourth-order of the scales defined in Johnson (1980).

### 5. Why ‘habitat selection function’ and not ‘resource selection function’?

Definitions only.

### 6. Considering objectives and data collection

This is the hardest part of the whole habitat selection study. How do
you decide on the what, where, and how many of sampling to answer the
question at hand. Talk to many people. Ask questions to many people.

The modeling mentioned here is about differences between study goals,
such as inference and prediction. Modeling for inference versus
prediction is not a straightforward distinction. There are lots of
opinions from the statistical folks (which are great and everyone should
read them, e.g., [Shmueli, 2010](https://doi.org/10.1214/10-STS330) and
[Scholz and Burner, 2022](https://arxiv.org/abs/2210.06927)) and the
science philosophy folks. We reference the manuscript [Gerber and
Northrup,
2020](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2953)
in regards to when the study goal is preditiction. Associated with this
manuscript is code in the Supporting Information file
(ecy2953-sup-0004-DataS1.Zip) that pertains to optimizing for predicting
spatial distributions of habitat selection (i.e., making a map). This
process can jeopardize inference, e.g., make the interpretation of
estimated effects unreliable. In contrast, if inference is sought you
should think hard about a single model that includes the most important
variables for the species and for your question that you want to
consider so that estimated effects and measures of uncertainty are
reliable ([Bolker 2023](https://doi.org/10.32942/X2Z01P), [Tredennick et
al. 2021](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3336)).

The data used to fit a movement-based HSF model is often individual’s
location is sampled at a high-fix rate, relative to the animal’s speed.
We curate our habitat selection data based on the movement process
(e.g., step lengths and turning angles). In the above, we have simulated
data with spatial locations at a standardized rate of movement and unit
of time between movements. It is common for real-world telemetry data to
have inconsistent times between consecutive spatial locations. It is
important to standardize this time difference, otherwise the movement
process parameters won’t be meaningful. For example, a time difference
of 30 minutes and 60 minutes are likely to produce different step
lengths and thus summarizing the what is availabile on each used
location will be compromised.

One way to standardize the time differences is by designing the data
setup such that relocations are specified into tracks based on a
threshold tolerance of allowable time between consecutive spatial
locations (e.g., done in the R package amt; [Signer et
al. 2019](https://doi.org/10.1002/ece3.4823)). Another way is to fit a
continuous time movement model to then predict locations at a
standarized interval (e.g., [Northrup et
al. 2018](https://doi.org/10.1111/gcb.13037) and [Johnsone et
al. 2008](https://doi.org/10.1890/07-1032.1))

### 7. What are habitat-selection functions used for?

In this section, we discuss some mechanics of inference and prediction,
as in interpreting coefficients and predicting relative habitat
selection. We will walk through the basics here, but refer the reader to
full in depth coverage in [Fieberg et
al. 2021](https://doi.org/10.1111/1365-2656.13441) and [Avgar et
al. 2017](https://doi.org/10.1002/ece3.3122). We also mention the goal
of using habitat selection to predict animal abundance. Demonstrating
this is beyond our work. We suggest starting with the reading of section
1.6, “When is species density a reliable reflection of habitat
suitability?” of [Matthiopoulos, Fieberg, and Aarts,
2023](http://hdl.handle.net/11299/217469).

#### Inference

Let’s fit a model to one individual’s matched used and available data
and discuss the estimated movement-based HSF coefficients. Remember, we
have setup an available sample to be paired with a used location to
control for the issues of serial dependence or autocorrelation between
sequential relocations. This process has also hopefully better defined
what is truly accessible to the individual at each location based on
their movement dynamics (specifically, turning angle and step-length
within a standardized time period). By doing this, we have focused our
attention on the dynamic selection of habitat as the animal movemees
through the landscape.

``` r
# We will use data from individual #1
  indiv.data1 = sims$indiv[[1]]

# Let's look at the data to remind us what it is
  head(indiv.data1)
```

    ##     status strata   dist.dev
    ## 1        1      1 -0.5773911
    ## 401      0      1 -0.5933174
    ## 402      0      1  0.4770586
    ## 403      0      1 -0.5773911
    ## 404      0      1 -0.5773911
    ## 405      0      1 -0.5773911
    ##       forest      shrub
    ## 1   1.218153 0.16764176
    ## 401 1.287160 0.01988518
    ## 402 1.546580 0.34669398
    ## 403 1.218153 0.16764176
    ## 404 1.218153 0.16764176
    ## 405 1.218153 0.16764176
    ##           step
    ## 1    0.8705149
    ## 401  3.6848288
    ## 402 10.4154689
    ## 403  2.5056502
    ## 404  0.2987808
    ## 405  0.1996507

We are now ready to approximate our true model using two difference
approaches: the conditional logistic regression model and the Poisson
regression model.

First, we can accomplish the first apporach using the `clogit function`
in package `survival`. We have related our response variable of the used
and available sample (`status`) to our covariates and tell the model how
the 1’s and 0’s are matched with the variable `strata`, which is part of
the linear equation but wrapped by an internal function `strata`. Note
that the variable combinations are one of an additive model, as opposed
to having interactions.

``` r
  model1 = clogit(status ~ dist.dev + forest + shrub + strata(strata), 
                  data = indiv.data1
                  ) 

# Look at estimated coefficients
  summary(model1)
```

    ## Call:
    ## coxph(formula = Surv(rep(1, 40400L), status) ~ dist.dev + forest + 
    ##     shrub + strata(strata), data = indiv.data1, method = "exact")
    ## 
    ##   n= 40400, number of events= 400 
    ## 
    ##              coef exp(coef)
    ## dist.dev -0.32642   0.72150
    ## forest   -0.13658   0.87234
    ## shrub     0.56172   1.75369
    ##          se(coef)      z
    ## dist.dev  0.08044 -4.058
    ## forest    0.08917 -1.532
    ## shrub     0.08083  6.949
    ##          Pr(>|z|)    
    ## dist.dev 4.95e-05 ***
    ## forest      0.126    
    ## shrub    3.67e-12 ***
    ## ---
    ## Signif. codes:  
    ##   0 '***' 0.001 '**' 0.01
    ##   '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##          exp(coef) exp(-coef)
    ## dist.dev    0.7215     1.3860
    ## forest      0.8723     1.1463
    ## shrub       1.7537     0.5702
    ##          lower .95 upper .95
    ## dist.dev    0.6163    0.8447
    ## forest      0.7325    1.0389
    ## shrub       1.4968    2.0547
    ## 
    ## Concordance= 0.623  (se = 0.012 )
    ## Likelihood ratio test= 75.71  on 3 df,   p=3e-16
    ## Wald test            = 66.77  on 3 df,   p=2e-14
    ## Score (logrank) test = 63.25  on 3 df,   p=1e-13

Let’s first consider the `Intercept`. OH, there is no intercept. That is
correct in these models and for our purposes, the intercepts have no
biological meaningful interpretation and thus are not necessary.

Next, we have estimated the effect of `dist.dev` as -0.33. This estimate
is negative, indicating that as the value of `dist.dev` increases from
its mean (defined at 0) habitat selection decreases, given that the
values of `forest` and `shrub` (the other variables in the model) are
held at their mean (0). In other words, habitat use relative to what is
available to this individual (as we defined it dynamically on a movement
track) decreases the further from development. That’s a lot of
qualifiers to understand this estimate. These are important though. In
terms of evidence of an effect, we can say that at a defined Type I
error of 0.05, this is a [statistically
clear](https://doi.org/10.1111/2041-210X.13159) effect that is not
likely zero. The evidence of this is the very small p-value of
4.9481644^{-5}.

Next, we have estimated the effect of `forest` as -0.14. This estimate
is negative, indicating that as the value of `forest` increases from its
mean (defined at 0) habitat selection decreases relative to what is
available to this individual (assuming `dist.dev` and `shrub` are at
their mean). At our defined Type I error, this is a statistically clear
effect that it is not likely zero (p-value = 0.1256098).

Next, we have estimated the effect of `shrub` as 0.56. This estimate is
positive, indicating that as the value of `shrub` increases from its
mean (defined at 0) habitat selection increases relative to what is
available to this individual (assuming `dist.dev` and `forest` are at
their mean). At our defined Type I error, this is a statistically clear
effect that is not likely zero (p-value = 3.6682191^{-12}).

The above output also conveniently shows us the exponentiation of each
coefficient. This is helpful because `exp(coefficient)` “quantifies the
relative intensity of use of two locations that differ by 1 standard
deviation unit of \[the variable\] but are otherwise equivalent
(i.e. they are equally available and have the same values of all other
habitat covariates)” ([Fieberg et
al. 2021](https://doi.org/10.1111/1365-2656.13441)).

##### **Side-bar**

We can estimate the same parameters and get the very same estimates with
other functions that implement the conditional logistic regression
model. For example, the functions `fit_clogit` or `fit_ssf` in the `amt`
package. These are are wrapper functions for using the `survival clogit`
function. We can see what the function is doing:

``` r
amt::fit_ssf
```

    ## function (data, formula, more = NULL, summary_only = FALSE, ...) 
    ## {
    ##     if (!any(grepl(pattern = "^strata\\(.+\\)$", attr(terms(formula), 
    ##         "term.labels")))) {
    ##         stop("No strata is provided, please make sure the formula includes a strata.")
    ##     }
    ##     m <- survival::clogit(formula, data = data, ...)
    ##     if (summary_only) {
    ##         m <- list(model = broom::tidy(m), sl_ = attributes(data)$sl_, 
    ##             ta_ = attributes(data)$ta_, more = more)
    ##         class(m) <- c("fit_clogit_summary_only", "fit_clogit", 
    ##             class(m))
    ##     }
    ##     else {
    ##         m <- list(model = m, sl_ = attributes(data)$sl_, ta_ = attributes(data)$ta_, 
    ##             more = more)
    ##         class(m) <- c("fit_clogit", class(m))
    ##     }
    ##     m
    ## }
    ## <bytecode: 0x000002b2e68f49c0>
    ## <environment: namespace:amt>

Now lets use these functions to estimate the same coefficients as we did
with `clogit`:

``` r
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
```

    ##   dist.dev     forest 
    ## -0.3264248 -0.1365766 
    ##      shrub 
    ##  0.5617218

``` r
  coef(model1.amt1)
```

    ##   dist.dev     forest 
    ## -0.3264248 -0.1365766 
    ##      shrub 
    ##  0.5617218

``` r
  coef(model1.amt2)
```

    ##   dist.dev     forest 
    ## -0.3264248 -0.1365766 
    ##      shrub 
    ##  0.5617218

We can see the estimated coefficients are all the same.

##### Poisson model equivalnce to the conditional logistic model

[Muff et al. 2019](https://doi.org/10.1111/1365-2656.13087) demonstrated
the equivalence between the conditional logistic regression model and
Poisson regression model with stratum-specific fixed intercepts. We can
fit this model using the `glmmTMB` package, which has computational
advantages over fitting the conditional logistic regression model. The
key to this implementation is that we want to include stratum as part of
the linear combination of variables wrapped in the function `strata` to
estimate fixed effect intercepts, or do the same procedure in a random
effect implementation but without shrinkage by fixing the variance to a
large value.

``` r
# Fit the model with fixed effect stratum-specific intercepts
    model1.tmb=glmmTMB(status ~ dist.dev+forest+shrub+strata(strata), 
                       data=indiv.data1,
                       family=poisson
                       )

# Or using a random effect with fixed variance
    model1.tmb2=glmmTMB(status ~ dist.dev+forest+shrub + (1| strata), 
                        data=indiv.data1,
                        family=poisson,
                        doFit = FALSE
                       )

# Make the intercept large and fixed
  model1.tmb2$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the intercept  
  model1.tmb2$mapArg = list(theta=factor(c(NA)))
  
# Now ready to fit the model  
  model1.tmb2 = glmmTMB:::fitTMB(model1.tmb2)
```

    ## Warning in
    ## finalizeTMB(TMBStruc, obj,
    ## fit, h, data.tmb.old): Model
    ## convergence problem; function
    ## evaluation limit reached
    ## without convergence (9). See
    ## vignette('troubleshooting'),
    ## help('diagnose')

Lets compare the coefficient estimates from these two model fits with
one of the clogit fits from above.

**Poisson random intercept (fixed variance)**

``` r
knitr::kable(summary(model1.tmb2)[[6]]$cond[2:4,],digits=3)
```

|          | Estimate | Std. Error | z value | Pr(\>\|z\|) |
|:---------|---------:|-----------:|--------:|------------:|
| dist.dev |   -0.326 |      0.080 |  -4.058 |       0.000 |
| forest   |   -0.137 |      0.089 |  -1.532 |       0.126 |
| shrub    |    0.562 |      0.081 |   6.949 |       0.000 |

**Poisson Fixed effect**

``` r
knitr::kable(summary(model1.tmb)[[6]]$cond[2:4,],digits=3)
```

|          | Estimate | Std. Error | z value | Pr(\>\|z\|) |
|:---------|---------:|-----------:|--------:|------------:|
| dist.dev |   -0.326 |      0.080 |  -4.058 |       0.000 |
| forest   |   -0.137 |      0.089 |  -1.532 |       0.126 |
| shrub    |    0.562 |      0.081 |   6.949 |       0.000 |

**clogit model**

``` r
knitr::kable(summary(model1)[[7]][,-2],digits=3)
```

|          |   coef | se(coef) |      z | Pr(\>\|z\|) |
|:---------|-------:|---------:|-------:|------------:|
| dist.dev | -0.326 |    0.080 | -4.058 |       0.000 |
| forest   | -0.137 |    0.089 | -1.532 |       0.126 |
| shrub    |  0.562 |    0.081 |  6.949 |       0.000 |

We find the same exact estimates. Woohoo!

Note that in the fixed effect implementation (model1.tmb) we are
actually estimating a lot of partial intercepts (strata). However, we
can ignore them as they are not of interest, just a means to an end.
That is way in the above code, we are limiting the estimates from the
summary with the part `cond[2:4,]`

**Relative Selection**

How do we quantitatively evaluate two locations in terms of habitat
selection using our model results? We can do so using Equation 8 of
[Fieberg et al. 2021](https://doi.org/10.1111/1365-2656.13441).

Perhaps we want to compare a location that is very near development
(dist.dev = -2) at high forest and shrub cover (forest = 2, shrub = 2)
with that of a location far from development (dist.dev = 2) and also at
high forest but low shrub cover (forest =2, shrub = -2).

``` r
# Get estimates
  coef = coef(model1)

# Relative Selection
  RS = exp(-2*coef[1] + 2*coef[2] + 2*coef[3]) / exp(2*coef[1] + 2*coef[2] + -2*coef[3])
  RS
```

    ## dist.dev 
    ## 34.90348

If given equal availability, this individual would relatively select the
first location by a factor of 34.9.

**Relative Selection Strength**

[Avgar et al. 2017](https://doi.org/10.1002/ece3.3122) refers to the
relative selection strength (RSS) as the exponentiation of each
coefficient.

For example,

``` r
# Coefficient for forest
  exp(coef[2])
```

    ##    forest 
    ## 0.8723395

Given two locations that differ by 1 standard deviation of percent
forest (`forest`), but are otherwise equal, this individual would be
0.87 as likely to choose the location with higher `forest`, or
equivalently, 1.15 times more likely to choose the location with the
lower `forest.`

#### Prediction

For convenience, we can use the `log_rss` function of the `amt` package
to predict the log relative strength when comparing many locations with
differing covariates values. To do so, we need create two data frames of
the covariate combinations of interest: `s1` and `s2`.

``` r
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
```

``` r
  plot(lr2,lwd=3)
```

    ## Guessing x_vars:
    ##   x_var1 = dist.dev
    ##   x_var2 = forest

``` r
  lines(lr2$df$dist.dev_x1,lr2$df$lwr,lty=2)
  lines(lr2$df$dist.dev_x1,lr2$df$upr,lty=2)
  abline(h=0,lwd=2,col=2)
```

![](MovementHSF_files/figure-markdown_github/plot.rss-1.png)

We have predicted a mean (solid line) and 95% confidence intervals
(dotted) of the log-RSS between each value of s1 relative to s2. These
locations differ only in their values of `dist.dev`, and the log-RSS is
equal to 0 when the two locations are identical (i.e., dist.dev = 0).
Above the zero line the individiaul is selected more than available and
below the zero line they are selecting less than available.

### 8. Traditional HSF or SSF?

In fitting these data using a movement-based HSF, we are simultaneously
trying to account for two important assumptions. First, we are limiting
what is available to the animal at each used location by defining the
availability as around this area and informed by the movement parameters
of step-length (how far is the animal likely to do) and step-length
(what direction is the animal likely to go). Second, we are defining our
data setup to deal with the fine-scale sequential dependency of
consecutive relocations and minimize issues of autocorrelation. If there
is interest in using these same data at a lower selection-order, its
important to consider the issues of autocorrelation and availability. In
the manuscript, we provide two options for doing so.

### 9. The model being fit

Context only.

### 10. The model being approximated

It is the responsibility of the researcher to make sure their modeling
process is done such that they are approximating the true underlying
model well. To demonstrate this, we will consider how estimated
coefficients change with increasing numbers of available samples per
each used location.

``` r
# Grab individual one's data
  indiv.dat1 =sims$indiv[[1]]

#Size of available samples per used locaiton
  n.avail2=seq(2,100,by=2)
  
#Save coefficients
  coef.save=NULL
  
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
```

``` r
par(mfrow=c(1,1))
plot(coef.save$N.Avail, coef.save$beta1,lwd=3,type="l",col=2,ylim=c(-01,1),
     main="Slope Coeficients",
     xlab="Number of Available Samples",ylab="Coeficient Estimate")
lines(coef.save$N.Avail, coef.save$beta2,lwd=3,col=2)
lines(coef.save$N.Avail, coef.save$beta3,lwd=3,col=2)
```

![](MovementHSF_files/figure-markdown_github/plotting.sensitiivty-1.png)

``` r
#knitr::kable(coef.save,digits=3)
```

We can see that our estimates of the slope coefficients are sensitive to
the number of available locations. The estimates visually stabilize
after about 80 available locations per used location.

Another way to look at this sensitivity is by showing the variability of
the coefficients within and across different sizes of available samples.

``` r
# Grab individual one's data
  indiv.dat1 =sims$indiv[[1]]

# Size of available samples per used locaiton
  n.avail2=c(2,10,20,40,60,80,100)
  
# Number of sample replications
 n.sample = 20
  
# Save coefficients
  coef.save=NULL
  
#Loop across available sample sizes  
  for(i in 1:length(n.avail2)){
  
  #re-sample each available sample size 20 times  
  for (z in 1:n.sample){  
  #Loop through each used location and reduce the number of available samples
    ind.dat=NULL
  for(j in 1:n.step){
    index = which(indiv.dat1$strata==j)
    #Get the used locaiton and the sub-sampled availables
    rand.sample= sample(index[-1],n.avail2[i])
    ind.dat.temp = indiv.dat1[c(index[1],rand.sample),]
    ind.dat = rbind(ind.dat,ind.dat.temp)
  } # end j step loop
    
# fit model with data
    model.temp = clogit(status ~ dist.dev + forest + shrub + strata(strata), 
                      data = ind.dat
                      ) 
  coef.save= rbind(coef.save, 
                   coef(model.temp)
                  )
  }#end z sample loop
}#end i loop
```

``` r
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


colnames(plot.data) = c("N.Available.Sample","Sample","Coeficient.Estiamte","Name")

ggplot2::ggplot(plot.data, aes(N.Available.Sample, Coeficient.Estiamte, fill=factor(Name))) +
  theme_bw()+
  geom_boxplot()
```

![](MovementHSF_files/figure-markdown_github/plot.sensitive-1.png)

We see more easily in this plot that there is high variability in the
estimated coefficients when the available sample is small. This
variability decreases as the available sample grows, showing how the
estimates are stabilizing.

### 11. Individuals

We should *a priori* assume there is individual variation in habitat
selection estimates. This variation is masked when data are pooled. Lets
consider the implications of pooling all data vs fitting a model
separately to each individual. Note that the object `sims2` has all the
individuals combined into a single data frame. We will use the `glmmTMB`
function to fit the model using the Poisson model regression approach
with strata as a random effect, but that is not really estimated. This
is fast and stable implementation.

``` r
# Pooled model - no consideration of individual variability
  model.pool = glmmTMB(status ~ dist.dev + forest + shrub + (1|indiv.id.strata), 
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
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## status ~ dist.dev + forest + shrub + (1 | indiv.id.strata)
    ## Data: sims2
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## 199819.3 199865.7 -99905.7 199811.3   807996 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups          Name        Variance Std.Dev.
    ##  indiv.id.strata (Intercept) 1e+06    1000    
    ## Number of obs: 808000, groups:  indiv.id.strata, 8000
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.004228  11.180350   0.000     1.00    
    ## dist.dev    -0.205381   0.012184 -16.856   <2e-16 ***
    ## forest      -0.005022   0.011773  -0.427     0.67    
    ## shrub        0.222130   0.012099  18.359   <2e-16 ***
    ## ---
    ## Signif. codes:  
    ## 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
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
```

In the pooled data, we see that the standard errors of the coefficients
and p-values are very small. We are using a lot of information, assuming
no variation among individuals, and thus assuming every individual’s
effects can be characterized by these estimates. Now let’s look at just
two individual’s estimates when fitting separate models.

``` r
summary(indiv.fit[[2]])
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## status ~ dist.dev + forest + shrub + (1 | strata)
    ## Data: x
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   9995.6  10030.1  -4993.8   9987.6    40396 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  strata (Intercept) 1e+06    1000    
    ## Number of obs: 40400, groups:  strata, 400
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -5.15471   50.00004  -0.103   0.9179    
    ## dist.dev    -0.10743    0.09245  -1.162   0.2452    
    ## forest       0.23440    0.09505   2.466   0.0137 *  
    ## shrub        0.81190    0.15177   5.350 8.81e-08 ***
    ## ---
    ## Signif. codes:  
    ## 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(indiv.fit[[10]])
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## status ~ dist.dev + forest + shrub + (1 | strata)
    ## Data: x
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   9988.9  10023.3  -4990.5   9980.9    40396 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  strata (Intercept) 1e+06    1000    
    ## Number of obs: 40400, groups:  strata, 400
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.005454  50.000042   0.000 0.999913    
    ## dist.dev    -0.069121   0.098739  -0.700 0.483904    
    ## forest       0.378282   0.108093   3.500 0.000466 ***
    ## shrub        0.527515   0.109547   4.815 1.47e-06 ***
    ## ---
    ## Signif. codes:  
    ## 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

First, notice that the estimated effects are a bit different than the
pooled estimates. Importantly, also notice that the standard errors and
p-values are larger.

Now lets look at all the estimates together.

``` r
# Extract pooled estimates with Confidence intervals
 pool.effect = broom.mixed::tidy(model.pool, effects = "fixed", conf.int = TRUE)

# Extract separate individual estimates
 indiv.coef = lapply(indiv.fit,FUN=function(x){
                     temp=broom.mixed::tidy(x, effects = "fixed", conf.int = TRUE)
                     data.frame(temp[-1,c(4,8,9)])
 })

# These are all coefs (1:3) for each individual 1:n.indiv
  estimates = sapply(indiv.coef,"[[",1)
  LCI = sapply(indiv.coef,"[[",2)
  UCI = sapply(indiv.coef,"[[",3)
```

``` r
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE),
       width=c(2,1))
plotCI(1:n.indiv,y=estimates[1,],
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
```

![](MovementHSF_files/figure-markdown_github/plot.pool.separate-1.png)

The plot on the right are the pooled estimates for *β*<sub>1</sub>
(black; dist.dev), *β*<sub>2</sub> (red; forest), and *β*<sub>3</sub>
(green; shrub). The plot on the left are each individual’s estimates
(colored the same). By ignoring individual variation, we are a much too
confident in our estimates of uncertainty and are ignoring a lot of
clear variation by individual. The pooled estimates certainty is do to
pseudoreplication - the treating a sub-nuit (each location) as are main
unit of replication. Our true unit of replication is the individual.

Since our sample sizes for each individual are equal, we see that the
pooled estimates generally relate to the average of each estimate across
all individuals. When the number of used locations varies by individual
this won’t be the case. The individuals with more used locations will
disproportionately influence the pooled estimates.

#### Sample size

The code for implementing the methods of [Street et
al. 2021](https://doi.org/10.1111/2041-210X.13701) can be found at
[figshare](https://figshare.com/articles/dataset/Datasets_and_Code_zip/11910831).
We are working on a more user friendly implementation.

### 12. Population inference

If we are interested in obtaining inference to the individual- and
population-level effects, we can consider bootstrapping the results from
individual models or jointly fitting a single model with a random effect
across individuals. We will first consider the bootstrapping method.
Second, we will demonstrate the use of random intercepts and slopes, as
well as how to fix the variance of the random intercept so there is no
shrinkage.

#### Bootrapping

The core functions for bootrapping are adopted from
Dr. Bastille-Rousseau’s github repository for the R package
[IndRSA](https://github.com/BastilleRousseau/IndRSA/).

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

#Which columns have coeficients in coef.df
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

We now want to summarize our bootstrapped results. The population mean
and 95% lower and upper confidence intervals are outputted.

``` r
#Source summary function
  boot.pop.esimates=mean_ci_boot(boot)
  rownames(boot.pop.esimates)=c("dist.dev","forest","shrub")
  knitr::kable(boot.pop.esimates,digits=3)
```

|          |   Mean |    LCI |    UCI |
|:---------|-------:|-------:|-------:|
| dist.dev | -0.424 | -0.560 | -0.325 |
| forest   |  0.008 | -0.242 |  0.197 |
| shrub    |  0.593 |  0.527 |  0.678 |

#### Random effects across individuals

##### Random intercept only model

The most common use of random effects in habitat selection analysis is
the use of a random intercept. As described in the manuscript, this does
not account for variation in the slope coefficients, which is
problematic.

``` r
# Setup HSF with Poisson regression approximation - random intercept model
# indiv.id.strata indicates individuals and strata and there the variance is fixed so there is no shrinkage
  re_int = glmmTMB(status ~ dist.dev + forest + shrub  +
                             (1|indiv.id.strata), 
                    family=poisson ,
                    data = sims2, 
                    doFit=FALSE
                    )

# Make the intercept large and fixed
  re_int$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the intercept  
  re_int$mapArg = list(theta=factor(c(NA)))
  
# Now ready to fit the model  
  options(warn=-1)
  re_int = glmmTMB:::fitTMB(re_int)
  options(warn=0)
```

``` r
# The fixed components - these are the mean effects when pooling all individuals. 
  fixef(re_int)
```

    ## 
    ## Conditional model:
    ## (Intercept)     dist.dev       forest        shrub  
    ##   -0.004228    -0.205381    -0.005022     0.222130

``` r
# The random components- these are the effect differences for each variation from the mean estimated (referred to as the fixed effect in the 'conditional model' statement). For this model we will have many many intercepts. Specifically, we will have a unique value for the number of individuals times the number of strata within individual
  head(ranef(re_int)[[1]][[1]])
```

    ##     (Intercept)
    ## 101   -4.643958
    ## 102   -4.620841
    ## 103   -4.622376
    ## 104   -4.600858
    ## 105   -4.641572
    ## 106   -4.608517

``` r
# Summarize results  
  summary(re_int)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## status ~ dist.dev + forest + shrub + (1 | indiv.id.strata)
    ## Data: sims2
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## 199819.3 199865.7 -99905.7 199811.3   807996 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups          Name        Variance Std.Dev.
    ##  indiv.id.strata (Intercept) 1e+06    1000    
    ## Number of obs: 808000, groups:  indiv.id.strata, 8000
    ## 
    ## Conditional model:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.004228  11.180350   0.000     1.00    
    ## dist.dev    -0.205381   0.012184 -16.856   <2e-16 ***
    ## forest      -0.005022   0.011773  -0.427     0.67    
    ## shrub        0.222130   0.012099  18.359   <2e-16 ***
    ## ---
    ## Signif. codes:  
    ## 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

A random intercepts only model is not recommended.

#### Random intercept and slopes model

This is the recommended model structure using random effects.

``` r
#Fit RSF with intercept and slopes with random effect
  re_int_slopes =  glmmTMB(status ~ -1 + dist.dev + forest + shrub  +
                                    (1|indiv.id.strata) +
                                    (0+dist.dev|ID) + (0+forest|ID) + 
                                    (0+shrub|ID), 
                            family=poisson, 
                            data = sims2, 
                            doFit=FALSE
                            )
# make the intercept large and fixed
  re_int_slopes$parameters$theta[1] = log(1e3)
  
# Tell glmmTMB to not estimate the intercept, but 
# to do so for the other three variables
  re_int_slopes$mapArg = list(theta=factor(c(NA,1:3)))
  
# Now ready to fit the model 
  options(warn=-1)
  re_int_slopes = glmmTMB:::fitTMB(re_int_slopes)
  options(warn=0)
```

Let’s look at the results. Here, we have the population-level (across)
individual-level means. So, generally, we see that at this level the
estimated coefcieents are negative, near zero, and positive for
`dist.dev`, `forest`, and `shrub`, respecitively.

``` r
  fixef(re_int_slopes)
```

    ## 
    ## Conditional model:
    ##   dist.dev      forest       shrub  
    ## -0.2272639  -0.0009686   0.2629444

We can get estimates with 95% confidence intervals as well. We see that
these statements are supported statistically in the above statement
based on examining whether intervals include 0 or not.

``` r
  broom.mixed::tidy(re_int_slopes, 
                    effects = "fixed", 
                    conf.int = TRUE)[,c(3,4,8,9)]
```

    ## # A tibble: 3 × 4
    ##   term      estimate conf.low conf.high
    ##   <chr>        <dbl>    <dbl>     <dbl>
    ## 1 dist.dev -0.227     -0.276    -0.179 
    ## 2 forest   -0.000969  -0.0553    0.0533
    ## 3 shrub     0.263      0.212     0.313

Here are estimated population-level coefficients (across
individual-level effects). Compared to the bootstrapped results above,
they are generally similar. We should not expect them to be the same, as
the random effect model shares information across individuals and the
bootstrapped esimates do not. However, the interpretation of the two
approaches leads to the same conclusions about the estimated sign and
statistical clarity.

We can also extract the estimated difference by individual from the
population level coefficients, as well as estimated with confidence
intervals. These estimates are an indication of whether an individual
has a different estimate than the overall population mean.

``` r
  ranef.out=ranef(re_int_slopes)
  ranef.out$cond$ID
```

    ##        dist.dev       forest        shrub
    ## 1   0.046332120 -0.102293176  0.036038573
    ## 2   0.158094515  0.038932581 -0.131883251
    ## 3   0.115246652  0.008089206 -0.018456905
    ## 4  -0.147065624 -0.044692247  0.115266792
    ## 5   0.030006044  0.040091381  0.088706451
    ## 6  -0.017802913 -0.087790672 -0.014549868
    ## 7   0.058197331  0.128001064 -0.004171541
    ## 8   0.023313556 -0.160819544  0.043869932
    ## 9  -0.001101453  0.077174588  0.076163523
    ## 10  0.125579425  0.034168000 -0.113093199
    ## 11 -0.071428191  0.103552620 -0.090337947
    ## 12 -0.041381693 -0.023319485 -0.018770742
    ## 13 -0.010188894  0.027302646 -0.149271096
    ## 14 -0.037571177 -0.144947210  0.030409607
    ## 15 -0.123397840 -0.042922553  0.111464462
    ## 16 -0.043438760 -0.003859300  0.001424223
    ## 17 -0.013715912 -0.146229907 -0.009986730
    ## 18 -0.074085315  0.217068652  0.089269558
    ## 19  0.034987420 -0.013709299  0.042263254
    ## 20 -0.004418097  0.099888654 -0.090881279

``` r
  ranef.out.ci = broom.mixed::tidy(re_int_slopes, 
                                   effects = "ran_vals", 
                                   conf.int = TRUE)

#Non-intercept estimates per indivual  
  ranef.out.ci[-which(ranef.out.ci$term=="(Intercept)"),-c(1,2,3,4)]
```

    ## # A tibble: 60 × 5
    ##    term     estimate std.error conf.low conf.high
    ##    <chr>       <dbl>     <dbl>    <dbl>     <dbl>
    ##  1 dist.dev  0.0463     0.0521  -0.0558    0.148 
    ##  2 dist.dev  0.158      0.0517   0.0568    0.259 
    ##  3 dist.dev  0.115      0.0507   0.0159    0.215 
    ##  4 dist.dev -0.147      0.0636  -0.272    -0.0225
    ##  5 dist.dev  0.0300     0.0493  -0.0667    0.127 
    ##  6 dist.dev -0.0178     0.0568  -0.129     0.0936
    ##  7 dist.dev  0.0582     0.0509  -0.0415    0.158 
    ##  8 dist.dev  0.0233     0.0543  -0.0831    0.130 
    ##  9 dist.dev -0.00110    0.0503  -0.0996    0.0974
    ## 10 dist.dev  0.126      0.0523   0.0232    0.228 
    ## # ℹ 50 more rows

Lastly, we can a full summary of the model.

``` r
  summary(re_int_slopes)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## status ~ -1 + dist.dev + forest + shrub + (1 | indiv.id.strata) +  
    ##     (0 + dist.dev | ID) + (0 + forest | ID) + (0 + shrub | ID)
    ## Data: sims2
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## 199749.6 199819.2 -99868.8 199737.6   807994 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups          Name        Variance  Std.Dev. 
    ##  indiv.id.strata (Intercept) 1.000e+06 1.000e+03
    ##  ID              dist.dev    8.284e-03 9.102e-02
    ##  ID.1            forest      1.196e-02 1.093e-01
    ##  ID.2            shrub       8.833e-03 9.398e-02
    ## Number of obs: 808000, groups:  
    ## indiv.id.strata, 8000; ID, 20
    ## 
    ## Conditional model:
    ##            Estimate Std. Error z value Pr(>|z|)    
    ## dist.dev -0.2272639  0.0247632  -9.177   <2e-16 ***
    ## forest   -0.0009686  0.0277011  -0.035    0.972    
    ## shrub     0.2629444  0.0257609  10.207   <2e-16 ***
    ## ---
    ## Signif. codes:  
    ## 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note the variance and standard deviation estimates, indicating how
variable coefficients are across individuals. We see that `forest` is
estimated as the most variable. This corresponds well to how the data
were generated with forest coefficients having the highest standard
deviation of 1\`.

Another thing to notice is that the population level effect for forest
is not statistically clearly different than zero. Thus, at the
population level, percent forest is being selected in proportion to what
is available to each individual. Let’s dive into this a bit more. Let’s
plot the individual estimated along with the population-level estimate.

``` r
  indiv.coef.popModel = broom.mixed::tidy(re_int_slopes, effects = "ran_vals", conf.int = TRUE)
  index=which(indiv.coef.popModel$term=="forest")
  indiv.forest.popModel = data.frame(indiv.coef.popModel[index,c(5,6,8,9)])

# Add back mean to get individual estimates  
  indiv.forest.popModel[,2:4] = indiv.forest.popModel[,2:4]+rep(fixef(re_int_slopes)[[1]][2],each=20)
    
# Extract population mean and uncertainty for forest effect
  pop.coef.popModel = data.frame(broom.mixed::tidy(re_int_slopes, effects = "fixed", conf.int = TRUE))
  pop.coef.popModel=pop.coef.popModel[,c(4,8,9)]
```

``` r
#Plot
  plotCI(x=1:20,
         y=indiv.forest.popModel$estimate,
         ui=indiv.forest.popModel$conf.high,
         li=indiv.forest.popModel$conf.low,
         lwd=3,col=1,
         xlab="Individual",ylab="Forest Coeficient")
  abline(h=pop.coef.popModel$estimate[2],lwd=3,col=1,lty=4)
  abline(h=pop.coef.popModel$conf.low[2],lwd=3,col=2,lty=4)
  abline(h=pop.coef.popModel$conf.high[2],lwd=3,col=2,lty=4)
```

![](MovementHSF_files/figure-markdown_github/RE.plot.forest.indiv.pop-1.png)

Our plot shows the individual effects of `forest` (vertical lines) along
with the population-level mean and confidence intervals (horizontal
lines). What is clear is that the reason there is no statistically clear
difference of the population-level effect from zero is because there is
a wide range of effects across individuals. Some individuals have
slightly positive coefficients and some have slightly negative. Since
these are generally equal, they balance out to a population-level effect
of zero. The story is more complicated! To decide on statistical clarity
between an individidual’s effect and the population, you can go back to
the confidence intervals of the estimated effect differences.

#### Model fitting issues with random effects

### 13. Considering context in habitat selection analyses

In this section we discuss ways of considering behavior in habitat
selection analyses. To demonstrate the effect of ignoring behavior, we
will simulate data where selection is different for two sets of animal
locations. Imagine an individual is highly selective of forest cover
when it is resting, but when foraging (and otherwise) selects against
forest cover in proportion to its availability.

``` r
hsf.true.behav=vector("list",2)
hsf.true.behav[[1]]=covs[[2]]*2
hsf.true.behav[[2]]=covs[[2]]*-2

  sim.behav1 =  sim.ind.movement.hsf(hsf.true=hsf.true.behav,
                                     n.time=n.time,
                                     n.step=n.step,
                                     n.avail=n.avail,
                                     angle.mean=angle.mean,
                                     angle.kappa=angle.kappa,
                                     step.rate=step.rate
                                     )

# Combine the data  
  data.ignore= rbind(sim.behav1$indiv[[1]],
                     sim.behav1$indiv[[2]]
                     )  
  
#Create ID vector for the different behavioral data
  ID=rep(1:2,each=nrow(sim.behav1$indiv[[1]]))
  data.ignore$ID=ID
#Create ID vector for unique strata within individual
  data.ignore$indiv.id.strata=unique(data.ignore$ID+data.ignore$strata*100)
```

Lets fit a model where we ignore the behavioral differences in selection
and fit a single model with all the data.

``` r
model.ignore.behav =  clogit(status ~ dist.dev + strata(indiv.id.strata), 
                             data = data.ignore
                  ) 
summary(model.ignore.behav)
```

    ## Call:
    ## coxph(formula = Surv(rep(1, 80800L), status) ~ dist.dev + strata(indiv.id.strata), 
    ##     data = data.ignore, method = "exact")
    ## 
    ##   n= 80800, number of events= 800 
    ## 
    ##              coef exp(coef) se(coef)     z Pr(>|z|)
    ## dist.dev 0.009867  1.009916 0.035433 0.278    0.781
    ## 
    ##          exp(coef) exp(-coef) lower .95 upper .95
    ## dist.dev      1.01     0.9902    0.9422     1.083
    ## 
    ## Concordance= 0.504  (se = 0.01 )
    ## Likelihood ratio test= 0.08  on 1 df,   p=0.8
    ## Wald test            = 0.08  on 1 df,   p=0.8
    ## Score (logrank) test = 0.08  on 1 df,   p=0.8

We see that the estimated forest coefficient is not near the values of
the true values of -2 or 2. It’s somewhat in between near zero.
Essentially, when animals are selecting difference habitats because of
behavior, mixing across behaviors can lead to a muddle inference.

### 14. Interpreting coefficients and predicting

Interpreting coefficients and predicting is outlined above in subsection
`7. What are habitat-selection functions used for?`.

### 15. Model selection

If focused on inference, fitting a single model is not only okay, but
desirable.

### 16. Concluding remarks

We hope this vignette proided some utility. If so, let us know with an
email (<brian.gerber@colostate.edu>).

Have a nice day.

## Software

This report was generated from the R Statistical Software (v4.4.2; R
Core Team 2021) using the [Markdown
language](https://www.markdownguide.org/) and
[RStudio](https://posit.co/products/open-source/rstudio/). The R
packages used are acknowledged below.

| Package           | Version | Citation                                       |
|:-----------------|:--------|:--------------------------------------------|
| amt               | 0.3.0.0 | @amt                                           |
| base              | 4.4.2   | @base                                          |
| broom.mixed       | 0.2.9.6 | @broommixed                                    |
| circular          | 0.5.1   | @circular                                      |
| geoR              | 1.9.4   | @geoR                                          |
| glmmTMB           | 1.1.10  | @glmmTMB                                       |
| knitr             | 1.49    | @knitr2014; @knitr2015; @knitr2024             |
| plotrix           | 3.8.4   | @plotrix                                       |
| raster            | 3.6.30  | @raster                                        |
| remotes           | 2.5.0   | @remotes                                       |
| ResourceSelection | 0.3.6   | @ResourceSelection                             |
| Rfast             | 2.1.0   | @Rfast                                         |
| rmarkdown         | 2.29    | @rmarkdown2018; @rmarkdown2020; @rmarkdown2024 |
| sp                | 2.1.4   | @sp2005; @sp2013                               |
| survival          | 3.7.0   | @survival-book; @survival-package              |
| tidyverse         | 2.0.0   | @tidyverse                                     |
