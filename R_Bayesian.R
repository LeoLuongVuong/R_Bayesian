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
  alpha ~ dbeta(1, 1)  # Non-informative prior for base level
  beta ~ dgamma(0.01, 0.01)  # Non-informative prior for positive increment
  gamma ~ dbeta(1, 1)  # Non-informative prior for positive growth rate
  delta ~ dnorm(0.01, 0.01)  # Non-informative prior 
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
       delta = 20,
       .RNG.name = "base::Wichmann-Hill",
       .RNG.seed = 1),
  list(alpha = 0.3,
       beta = 0.3,  # Adjusted to have mean 1
       gamma = 0.3,  # Adjusted to have mean 1
       delta = 30,
       .RNG.name = "base::Marsaglia-Multicarry",
       .RNG.seed = 2),
  list(alpha = 0.4,
       beta = 0.4,  # Adjusted to have mean 1
       gamma = 0.4,  # Adjusted to have mean 1
       delta = 25,
       .RNG.name = "base::Super-Duper",
       .RNG.seed = 3))

## collect the parameters to be monitored ------------------------------------

# Parameters to monitor
parameters <- c("alpha", "beta", "gamma", "delta")

## run the MCMC chain ------------------------------------------------------------

# Run the MCMC
jags_model <- jags.model(file = 'varicella_BUGS.txt',
                         data = my_data,
                         inits = my.inits,
                         n.chains = 3)

coverage.sim_02 <- coda.samples(jags_model, # this is the best I've got, everything is perfect!
                                parameters,
                                n.iter = 6000000,
                                thin = 1)

# Take burn in into account - very important!
coverage.sim_02_5_6_mil <- window(coverage.sim_02, start = 5000000)

# Posterior summary statistics
summary(coverage.sim_02_5_6_mil)

### check with ggmcmc --------------------------------------------------------

# convert from mcmc.list to a dataset
out.ggs_5_6_mil_thin1 <- ggs(coverage.sim_02_5_6_mil) 
# take burnin into account

# make histogram for each parameter
#ggs_histogram(out.ggs_5_6_mil_thin1)

# create traceplot object
trace_plot_5_6_mil_thin1 <- ggs_traceplot(out.ggs_5_6_mil_thin1)
#trace_plot_3_5_mil_thin1

# make a running mean plot
running_mean_5_6_mil_thin1 <- ggs_running(out.ggs_5_6_mil_thin1) +
  theme(axis.text = element_text(size = 6),
        panel.spacing = unit(1, "lines"))
#running_mean_5_6_mil_thin1

geweke.plot_5_6_mil_thin1 <- ggs_geweke(out.ggs_5_6_mil_thin1)

# try out these following (optional)
ggs_compare_partial(out.ggs_5_6_mil_thin1)

ggs_pairs(out.ggs_5_6_mil_thin1)

ggs_autocorrelation(out.ggs_5_6_mil_thin1, nLags = 100)

ggs_caterpillar(out.ggs_5_6_mil_thin1)

ggs_crosscorrelation(out.ggs_5_6_mil_thin1)

ggs_density(out.ggs_5_6_mil_thin1)

ggs_diagnostics(out.ggs_5_6_mil_thin1)

### Convergence test ---------------------------------------------------------

# shrunken factor
gelman.diag(coverage.sim_02_5_6_mil, autoburnin = FALSE, transform = TRUE)
# this test assumes normality. ggs_density shows that it's not the case, so 
# transformation is needed

# plot this diagnostic
gelman.plot(coverage.sim_02_5_6_mil, autoburnin = FALSE)
#ggsave("gelman_plot_5_6_mil_thin1.png", gelman.plot(coverage.sim_02_5_6_mil, autoburnin = FALSE), dpi = 300, width = 19, height = 19, units = "cm")

# compare the mean of first 10% to later 50%
geweke.diag(coverage.sim_02_5_6_mil)

# plot of this test
geweke.plot(coverage.sim_02_5_6_mil, ask = FALSE)

## Saving plots -------------------------------------------------

setwd("./plots")
ggsave("trace_plot_5_6_mil_thin1.png", trace_plot_5_6_mil_thin1, dpi = 300, width = 19, height = 19, units = "cm")
ggsave("running_mean_5_6_mil_thin1.png", running_mean_5_6_mil_thin1, dpi = 300, width = 19, height = 9, units = "cm")
ggsave("geweke.plot_5_6_mil_thin1.png", geweke.plot_5_6_mil_thin1, dpi = 300, width = 19, height = 19, units = "cm")

# get back to the main directory
Path <- getwd()
setwd(dirname(Path))

## Another toolkit for MCMC chains -------------------------

print(coverage.sim_02_25_4_mil)
plot(coverage.sim_02_25_4_mil)

effectiveSize(coverage.sim_02_25_4_mil) # effective size
HPDinterval(coverage.sim_02_25_4_mil) # HPD intervals of all parameters

## Question 7 -----------------------------------------------------------
# use out.ggs_5_6_mil_thin1 for that

# convert into a wide format
out.ggs_5_6_mil_thin1_wide <- out.ggs_5_6_mil_thin1 |> 
  pivot_wider(names_from = "Parameter", values_from = "value")

# rep out.ggs_5_6_mil_thin1_wide 4 times
out.ggs_5_6_mil_thin1_wide <- rbind(out.ggs_5_6_mil_thin1_wide) |> 
  rbind(out.ggs_5_6_mil_thin1_wide) |> 
  rbind(out.ggs_5_6_mil_thin1_wide) |> 
  rbind(out.ggs_5_6_mil_thin1_wide)

# add Age, which has 4 values (13, 19, 24, 35), each is repeated nrow times
out.ggs_5_6_mil_thin1_wide <- out.ggs_5_6_mil_thin1_wide |> 
  mutate(Age = rep(c(13, 19, 24, 35), nrow(out.ggs_5_6_mil_thin1_wide)/4))

# add coverage column, which equals: alpha + beta / (1 + exp(-gamma * (Age - delta)))
out.ggs_5_6_mil_thin1_wide <- out.ggs_5_6_mil_thin1_wide |> 
  mutate(coverage = alpha + beta / (1 + exp(-gamma * (Age - delta))))

# summary statistics
coverage_predict <- out.ggs_5_6_mil_thin1_wide |>
  group_by(Age) %>%
  summarize(median = median(coverage),
            min = quantile(coverage, 0.025),
            max = quantile(coverage, 0.975),
            .groups = 'drop')

# add coverage column to varicella_vaccine_coverage dataset, which is Vaccinated/Sample_size
varicella_vaccine_coverage <- varicella_vaccine_coverage |>
  mutate(coverage = Vaccinated/Sample_Size)

# add coverage and Geography columns from varicella_vaccine_coverage to coverage_predict dataset
# matching by Age
coverage_predict_observed <- coverage_predict |>
  left_join(varicella_vaccine_coverage |> select(Age, Geography, coverage), by = "Age")

# plot the observed and predicted coverage
coverage_by_age <- ggplot(data = coverage_predict_observed, aes(x = Age, y = median)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#21918c", alpha = 0.5) +
  geom_line(color = "#440154", size = 0.8) +
  geom_point(aes(y = coverage, color = Geography)) +
  ylab("Vaccination coverage") +
  scale_x_continuous(limits = c(12, 36), breaks = c(13, 19, 24, 35), expand = c(0, 0.01)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1),expand = c(0, 0.01)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 8, family = "sans"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(size = 7, family = "sans"),
        legend.text = element_text(size = 7, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(color = "Region") 
coverage_by_age

### export the plot -------------------------------------------------

setwd("./plots")
ggsave("coverage_by_age.png", coverage_by_age, dpi = 300, width = 19, height = 9, units = "cm")
# get back to the main directory
Path <- getwd()
setwd(dirname(Path))

## Question 8 ---------------------------------------------------------

# prediction at age 15

predict_15 <- out.ggs_5_6_mil_thin1_wide |>
  filter(Age == 13) |>
  mutate(Age = 15) |>
  mutate(coverage = alpha + beta / (1 + exp(-gamma * (Age - delta))))

predict_15_summary <- predict_15 |>
  summarize(median = median(coverage),
            min = quantile(coverage, 0.025),
            max = quantile(coverage, 0.975))

# add Geography with Sample_Size 
predict_15_summary_geo <- rbind(predict_15_summary, predict_15_summary, predict_15_summary,
                            predict_15_summary, predict_15_summary) |>
  mutate(Geography = unique(varicella_vaccine_coverage$Geography),
          Sample_Size = unique(varicella_vaccine_coverage$Sample_Size))

predict_15_summary_geo <- predict_15_summary_geo |>
  mutate(Median = round(median * Sample_Size, 0),
         Min = round(min * Sample_Size, 0),
         Max = round(max * Sample_Size, 0))

# Extract Geography, Min, Median, Max, and export to latex
predict_15_summary_geo |>
  select(Geography, Min, Median, Max) 
         
                  
         