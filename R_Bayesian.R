# Embed BUGS in the R code ---------------------------------------------------

# To add the BUGS model in the same code, you can add the BUGS model within 
# your R code. The following model would create a file osteo.model2.txt 
# within the R code.

cat("model
 {
  for (i in 1:N){
   tbbmc[i] ~ dnorm(mu[i],tau)
   mu[i] <- beta0+beta1*(bmi[i]-mean(bmi[]))
  }
   
  sigma2 <- 1/tau
  beta0 ~ dnorm(0,1.0E-6)
  beta1 ~ dnorm(0,1.0E-6)
  tau ~ dgamma(1.0E-3,1.0E-3)
 }", file = "osteo.model2.txt")

# example1_rjags - R code ---------------------------------------------------

# Download current version of JAGS
# install.packages("rjags")
# install.packages("coda")

# RUN JAGS FROM INSIDE OF R

library("rjags")
library(coda)

# DATA PREPARATION

N <- 1000
x <- rnorm(N, 0,5)

model.data <- list('x' = x, 'N' = N)

# DEFINE INITIAL VALUES

model.inits <- list(sigma = 1,mu = 0)

# MODEL SPECIFICATION 
# -> PUT MODEL SPECIFICATION IN A FILE CALLED example1.txt

file.show("example1.txt")

# SET UP MODEL
# specify model, data, number of parallel chains
jags <- jags.model('example1.txt',
                   data = model.data,
                   inits = model.inits,
                   n.chains = 2)

# Generate MCMC samples and save output for specified variables
out <- coda.samples(jags,
                    c('mu', 'sigma'),
                    n.iter = 10000, thin = 1)

# Posterior summary statistics
burnin <- 2000
summary(window(out,start = burnin))


# History plot & posterior distributions & autocorrelation plot
plot(out, trace = TRUE, density = TRUE)   
plot(window(out, start = burnin), trace = TRUE, density = TRUE)   

densplot(out)
HPDinterval(out)

# Checking correlations
autocorr.plot(out)

par(mfrow = c(1,1))
crosscorr.plot(out)


# Convergence tests

gelman.diag(out)
gelman.plot(out,ask = FALSE)

geweke.diag(out)
geweke.plot(out,ask = FALSE)

# obtain DIC
dic <- dic.samples(model = jags,
                   n.iter = 1500, 
                   thin = 1)

# example2_rjags - R code ---------------------------------------------------

# Download current version of JAGS
# install.packages("rjags")
# install.packages("coda")


# RUN JAGS FROM INSIDE OF R

library('rjags')


# DATA PREPARATION

N <- 10
weeks <- c(1,2,3,4,5,6,7,8,9,10)
calls <- c(0,2,1,8,5,17,24,23,19,17)

model.data <- list('calls' = calls, 'weeks' = weeks, 'N' = N)

# DEFINE INITIAL VALUES

model.inits <- list(alpha = 1,beta = 1,gamma = 1)

# MODEL SPECIFICATION 
# -> PUT MODEL SPECIFICATION IN A FILE CALLED example2.txt

cat("model
  {
    # Specification data model
    for (i in 1:N)
    {
      eta[i] <- alpha+beta*weeks[i]
      lambda[i] <- gamma/(1+exp(-eta[i]))
      calls[i] ~ dpois(lambda[i])
    }
    # Prior specification
    alpha ~ dnorm(-5,16)
    beta ~ dnorm(0.75,4)
    gamma ~dgamma(3.5,0.08333)
  }", file = "example2.txt")
file.show("example2.txt")

# SET UP MODEL OBJECT 
# specify model, data, number of parallel chains
jags <- jags.model('example2.txt',
                   data = model.data,
                   inits = model.inits,
                   n.chains = 2)

# Generate MCMC samples and save output for specified variables
out <- coda.samples(jags,
                    c('alpha', 'beta', 'gamma'),
                    n.iter = 10000, thin = 1)

# Posterior summary statistics
burnin <- 5000
summary(window(out, start = burnin))

# History plot & posterior distributions
plot(out, trace = TRUE, density = FALSE)   
plot(out, trace = FALSE, density = TRUE)  

# Varicella vaccine data analysis -------------------------------------------

## load the packages ---------------------------------------------------------

library("R2OpenBUGS")
library("coda")
library("readr")
library(rjags)
library(ggmcmc)

## load the data ------------------------------------------------------------

varicella_vaccine_coverage <- read_csv("varicella_vaccine_coverage.csv")

# convert first two columns to factors
varicella_vaccine_coverage$Geography <- as.factor(varicella_vaccine_coverage$Geography)
varicella_vaccine_coverage$Age <- as.factor(varicella_vaccine_coverage$Age)
# change Age to Age_former
varicella_vaccine_coverage$Age_former <- varicella_vaccine_coverage$Age

# Add Age column
Age <- rep(c(13, 19, 24, 35), 5)
varicella_vaccine_coverage$Age <- Age

## write the BUGS program and put it in a txt file ---------------------------

cat("model 
{
  for (i in 1:J) {
    for (j in 1:M) {
      Y[i,j] ~ dbin(p[i,j], N[i,j])
      p[i,j] <- alpha + beta / (1 + exp(-gamma * (Age[j] - delta)))
    }
  }
  
  # Priors
  alpha ~ dnorm(0, 1.0E-2)  # Non-informative prior for base level
  beta ~ dgamma(0.01, 0.01)  # Non-informative prior for positive increment
  gamma ~ dgamma(0.01, 0.01)  # Non-informative prior for positive growth rate
  delta ~ dnorm(mean(Age), 1.0E-2)  # Non-informative prior centered around mean age
}", file = "varicella_BUGS.txt")

## prepare the data and collect them into the object `my.data' ---------------

# Prepare the data
Y <- matrix(varicella_vaccine_coverage$Vaccinated, nrow = 5, ncol = 4, byrow = TRUE)
N <- matrix(varicella_vaccine_coverage$Sample_Size, nrow = 5, ncol = 4, byrow = TRUE)
Age <- unique(varicella_vaccine_coverage$Age)

my_data <- list(
  J = nrow(Y),
  M = ncol(Y),
  Y = Y,
  N = N,
  Age = Age
)

## set the initial values ----------------------------------------------------

# Initial values
# Initial values

my.inits <- list(
  list(alpha = 0.5,
       beta = 0.5,  # Adjusted to have mean 1
       gamma = 0.5,  # Adjusted to have mean 1
       delta = 20),
  list(alpha = 0.3,
       beta = 0.3,  # Adjusted to have mean 1
       gamma = 0.3,  # Adjusted to have mean 1
       delta = 30),
  list(alpha = 0.4,
       beta = 0.4,  # Adjusted to have mean 1
       gamma = 0.4,  # Adjusted to have mean 1
       delta = 25))


## collect the parameters to be monitored ------------------------------------

# Parameters to monitor
parameters <- c("alpha", "beta", "gamma", "delta")

## run the MCMC chain ------------------------------------------------------------

# Run the MCMC
jags_model <- jags.model(file = 'varicella_BUGS.txt',
                         data = my_data,
                         inits = my.inits,
                         n.chains = 3)

coverage.sim_04 <- coda.samples(jags_model, 
                                parameters,
                                n.iter = 4000000,
                                thin = 10)

coverage.sim_05 <- coda.samples(jags_model, # the best till now, from 6mil to 8 mil
                                parameters,
                                n.iter = 8000000,
                                thin = 10)

coverage.sim_06 <- coda.samples(jags_model, 
                                parameters,
                                n.iter = 4000000,
                                thin = 5)

coverage.sim_07 <- coda.samples(jags_model,
                                parameters,
                                n.iter = 16000000,
                                thin = 20)


# Take burn in into account - very important!
coverage.sim_04_2_4_mil <- window(coverage.sim_04, start = 2000000)
coverage.sim_04_2_3_mil <- window(coverage.sim_04, start = 2000000, end = 3000000)
coverage.sim_04_0_15_mil <- window(coverage.sim_04, end = 1500000)

coverage.sim_05_6_8_mil <- window(coverage.sim_05, start = 6000000)

# Posterior summary statistics
# With burnin
burin_04 <- 2000000
summary_04 <- summary(window(coverage.sim_04, start = burin_04, end = 3000000))
# Or, a simpler way
summary_04 <- summary(coverage.sim_04_burnin)

# summary Without burnin
summary(coverage.sim_04)

summary(coverage.sim_05_6_8_mil)

### check with ggmcmc --------------------------------------------------------

out.ggs_2_4_mil <- ggs(coverage.sim_04_2_4_mil) # take burnin into account
out.ggs_2_3_mil <- ggs(coverage.sim_04_2_3_mil) # take burnin into account

out.ggs_0_8_mil <- ggs(coverage.sim_05)
out.ggs_6_8_mil <- ggs(coverage.sim_05_6_8_mil)

out.ggs_0_8_mil_thin5 <- ggs(coverage.sim_06)

ggs_histogram(out.ggs_2_4_mil)
ggs_histogram(out.ggs_2_3_mil)

ggs_histogram(out.ggs_6_8_mil)

trace_plot_2_4_mil <- ggs_traceplot(out.ggs_2_4_mil)
trace_plot_2_3_mil <- ggs_traceplot(out.ggs_2_3_mil)
trace_plot_0_8_mil <- ggs_traceplot(out.ggs_0_8_mil)
trace_plot_6_8_mil <- ggs_traceplot(out.ggs_6_8_mil)

trace_plot_0_8_mil_thin5 <- ggs_traceplot(out.ggs_0_8_mil_thin5)

ggs_running(out.ggs_2_3_mil)
ggs_running(out.ggs_2_4_mil)

ggs_running(out.ggs_6_8_mil)

ggs_compare_partial(out.ggs)
ggs_geweke(out.ggs)
ggs_pairs(out.ggs)
ggs_autocorrelation(out.ggs, nLags = 100)
ggs_caterpillar(out.ggs)
ggs_crosscorrelation(out.ggs)
ggs_density(out.ggs)
ggs_diagnostics(out.ggs_6_8_mil)

## Convergence test ---------------------------------------------------------

gelman.diag(coverage.sim_04_2_4_mil, autoburnin = FALSE)
gelman.diag(coverage.sim_04_2_3_mil, autoburnin = FALSE)
gelman.diag(coverage.sim_04_0_15_mil, autoburnin = FALSE)

gelman.diag(coverage.sim_05_6_8_mil, autoburnin = FALSE)


gelman.plot(coverage.sim_04_2_4_mil, autoburnin = FALSE) #gamma looks better with this
gelman.plot(coverage.sim_04_2_3_mil, autoburnin = FALSE)

gelman.plot(coverage.sim_05_6_8_mil, autoburnin = FALSE)


gelman.diag(coverage.mcmc_04)
gelman.diag(coverage.mcmc_05)
gelman.diag(coverage.mcmc_06)

gelman.plot(coverage.mcmc_01, ask = FALSE)
gelman.plot(coverage.mcmc_02, ask = FALSE)
gelman.plot(coverage.mcmc_03, ask = FALSE)

geweke.diag(coverage.mcmc)
geweke.plot(coverage.mcmc, ask = FALSE)

## Saving plots -------------------------------------------------
setwd("./plots")
ggsave("trace_01.png", trace_04, dpi = 300, width = 19, height = 9, units = "cm")
# get back to the main directory
Path <- getwd()
setwd(dirname(Path))

## Conclusion ---------------------------------------------------------

# 10^6 iterations, in which burnin 5*10^5 iterations.

## Redundance ---------------------------------------------------------
# Draft code posterior summary statistics
burnin <- 1000000
summary_01 <- summary(window(coverage.sim_01, start = burnin))
summary_02 <- summary(window(coverage.sim_02, start = burnin))
summary_03 <- summary(window(coverage.sim_03, start = burnin))
summary_04_b <- summary(window(coverage.sim_04, start = burnin))
burin_05 <- 4000000
summary_05 <- summary(window(coverage.sim_05, start = burin_05))
summary_06 <- summary(window(coverage.sim_06, start = burin_04))


# draft code running the MCMC chain
coverage.sim_01 <- coda.samples(jags_model,
                                parameters,
                                n.iter = 2000000,
                                thin = 20)

coverage.sim_02 <- coda.samples(jags_model,
                                parameters,
                                n.iter = 2000000,
                                thin = 10)

coverage.sim_03 <- coda.samples(jags_model,
                                parameters,
                                n.iter = 3000000,
                                thin = 20)

# History plot & posterior distributions
plot(coverage.sim, trace = TRUE, density = FALSE)   
plot(coverage.sim, trace = FALSE, density = TRUE)
# will use these
trace_01 <- plot(window(coverage.sim_01, start = burnin), trace = TRUE, density = FALSE) # plot discarding burn-in iterations
trace_02 <- plot(window(coverage.sim_02, start = burnin), trace = TRUE, density = FALSE)
trace_03 <- plot(window(coverage.sim_03, start = burnin), trace = TRUE, density = FALSE)
trace_04 <- plot(window(coverage.sim_04, start = burin_04), trace = TRUE, density = FALSE) #this traceplot looks ugly
trace_04_b <- plot(window(coverage.sim_04, start = burnin), trace = TRUE, density = FALSE)
trace_05 <- plot(window(coverage.sim_05, start = burin_04), trace = TRUE, density = FALSE)

### Produce general summary of obtained MCMC sampling -------------------------

print(coverage.sim)
plot(coverage.sim)

### Convert osteo.sim into mcmc.list for processing with CODA -----------------

coverage.mcmc_01 <- as.mcmc.list(coverage.sim_01)
coverage.mcmc_02 <- as.mcmc.list(coverage.sim_02)
coverage.mcmc_03 <- as.mcmc.list(coverage.sim_03)
coverage.mcmc_05 <- as.mcmc.list(coverage.sim_05)
coverage.mcmc_06 <- as.mcmc.list(coverage.sim_06)

### Produce general summary of obtained MCMC sampling -------------------------

plot(coverage.mcmc)
summary(coverage.mcmc)

### Specific output obtained from CODA functions -------------------------

par(mfrow = c(2,2)) # plot figures in 2x2 format if function allows
traceplot(coverage.mcmc_04) # trace plots #produce ugly traceplots independently
cumuplot(coverage.mcmc,ask = FALSE) # running mean plots
acfplot(coverage.mcmc) # autocorrelation function plot
autocorr(coverage.mcmc) # autocorrelation values
crosscorr.plot(coverage.mcmc) # cross-correlation output
densplot(coverage.mcmc) # density plots of the marginal posteriors
effectiveSize(coverage.mcmc) # effective size
HPDinterval(coverage.mcmc) # HPD intervals of all parameters


