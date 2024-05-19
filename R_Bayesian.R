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

## load the data ------------------------------------------------------------

varicella_vaccine_coverage <- read_delim("varicella_vaccine_coverage.txt", " ", 
                                         escape_double = FALSE, trim_ws = TRUE)
# doesn't work, try csv file instead
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
      logit(p[i,j]) <- alpha + beta / (1 + exp(-gamma * (Age[j] - delta)))
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
my_inits <- function() {
  list(alpha = rnorm(1, 0, 0.1),
       beta = rgamma(1, 0.01, 0.01),
       gamma = rgamma(1, 0.01, 0.01),
       delta = rnorm(1, mean(Age), 0.1))
}

## collect the parameters to be monitored ------------------------------------

# Parameters to monitor
parameters <- c("alpha", "beta", "gamma", "delta")

## run the MCMC chain ------------------------------------------------------------

# Run the MCMC
library(rjags)
jags_model <- jags.model(file = 'varicella_BUGS.txt',
                         data = my_data,
                         inits = my_inits,
                         n.chains = 3)
coverage.sim <- coda.samples(jags_model,
                             parameters,
                             n.iter = 10000,
                             thin = 1)

# Posterior summary statistics
burnin <- 5000
summary(window(coverage.sim, start = burnin))

# History plot & posterior distributions
plot(coverage.sim, trace = TRUE, density = FALSE)   
plot(coverage.sim, trace = FALSE, density = TRUE)

## Produce general summary of obtained MCMC sampling -------------------------

print(coverage.sim)
plot(coverage.sim)

## Convert osteo.sim into mcmc.list for processing with CODA -----------------

coverage.mcmc <- as.mcmc.list(coverage.sim)

## Produce general summary of obtained MCMC sampling -------------------------

plot(coverage.mcmc)
summary(coverage.mcmc)

## Specific output obtained from CODA functions -------------------------

par(mfrow = c(2,2)) # plot figures in 2x2 format if function allows
traceplot(coverage.mcmc) # trace plots
cumuplot(coverage.mcmc,ask = FALSE) # running mean plots
acfplot(coverage.mcmc) # autocorrelation function plot
autocorr(coverage.mcmc) # autocorrelation values
crosscorr.plot(coverage.mcmc) # cross-correlation output
densplot(coverage.mcmc) # density plots of the marginal posteriors
effectiveSize(coverage.mcmc) # effective size
HPDinterval(coverage.mcmc) # HPD intervals of all parameters


# Test -------------------



# Burn-in
update(jags_model, 1000)

# Sample from the posterior
samples <- coda.samples(jags_model, variable.names = parameters, n.iter = 10000)

# Summarize results
summary(samples)
