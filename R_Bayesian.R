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

# 

