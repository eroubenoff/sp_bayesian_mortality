#### Full proper CAR model with complete priors ####

library(nimble)
library(tigris)
library(rgeos)
library(dplyr)


#### Note ####
# This model assigns a phi parameter specific to [x, t, a].  So there will be 
# spatial effects within a given year and age.  I've smoothed year but am not 
# sure how to smooth across age-- I think this is getting too complex.

# Thus as far as spatial effects are concerned, age groups are independent.

#### Full Spatial Model ####
rm(list = ls())
CA_model_full_spatial <- nimbleCode({
  for (x in 1:X){         # age
    for (t in 1:T){       # year
      for (a in 1:A){     # county (area)
        y.xta[x, t, a] ~ dpois(mu.xta[x,t,a]) # deaths are poisson dist
        mu.xta[x, t, a] <- pop.xta[x, t, a] * mx.xta[x, t, a] # pop x mort rate
        mx.xta[x, t, a] <- exp(logmx.xta[x, t, a])
        
        # regression model
        logmx.xta[x, t, a] <- beta.ta[t, a, 1] * Yx[x, 1] + 
          beta.ta[t, a, 2] * Yx[x, 2] + 
          beta.ta[t, a, 3]*Yx[x, 3] + 
          phi[x, t, a] + 
          u.xta[x, t, a]
      }
      # prior on phi
      phi[x, t, 1:A] ~ dcar_proper(mu.phi[x, t, 1:A], adj[1:N],  num[1:A], tau.phi, gamma = gamma)
    }
  }
  
 
  
  # Priors on mu_phi:
  for (x in 1:X){
    for (a in 1:A){
      for (t in 1:2){
        mu.phi[x, t, a] ~ dnorm(0, tau.mu.phi)
      }
      for (t in 3:T){
        mu.phi[x, t, a] ~ dnorm(2*mu.phi[x, t-1, a] - mu.phi[x, t-2, a], tau.mu.phi)
      }
    }
  }
  
  # Other priors for mu.phi
  gamma ~ dunif(-1, 1)
  tau.phi <- pow(sigma.phi, -2)
  sigma.phi ~ dunif(0, 40)
  tau.mu.phi <- pow(sigma.mu.phi, -2)
  sigma.mu.phi ~ dunif(0, 40)
  
  # Priors on betas
  for (t in 1:T){
    for (p in 1:3){
      for (a in 1:A){
        beta.ta[t, a, p] ~ dnorm(mu.beta[t, p], tau.beta[t, p])
      }
      tau.beta[t, p] <- pow(sigma.beta[t, p], -2)
      sigma.beta[t, p] ~ dunif(0, 40)
    }
  }
  
  # Priors on u.xta (error)
  for (a in 1:A){
    for (x in 1:X){
      for (t in 1:T){
        u.xta[x, t, a] ~ dnorm(0, tau.u[x])
      }
    }
  }
  
  # Time smoothing and mu.betas
  for (t in 1:2){
    for (p in 1:3){
      mu.beta[t, p] ~ dnorm(0, tau.mu[p])
    }
  }
  for (t in 3:T){
    for (p in 1:3){
      mu.beta[t, p] ~ dnorm(2*mu.beta[t-1, p] - mu.beta[t-2, p], tau.mu[p])
    }
  }
  
  # priors on precision
  for (p in 1:3){
    tau.mu[p] <- pow(sigma.mu[p], -2)
    sigma.mu[p] ~ dunif(0, 40)
  }
  for (x in 1:X){
    tau.u[x] <- pow(sigma.u[x], -2)
    sigma.u[x] ~ dunif(0, 0.1)
  }
})