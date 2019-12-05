#### Simple, Base Model ####

library(nimble)
rm(list = ls())
setwd("/hdir/0/eroubenoff/R/HMD_ER/bayesian_spatial/CA_county_example")
load("data.RData")

constants <- list(
  X = data$X,
  T = data$T,
  A = data$n.a
)

# Variables:
# X: age
# T: year
# A: area

CA_model <- nimbleCode({
  for (x in 1:X){
    for (t in 1:T){
      for (a in 1:A){
        y.xtas[x, t, a] ~ dpois(mu.xta[x, t, a])
        mu.xta[x, t, a] <- pop.xtas[x, t, a] * mx.xtas[x, t, a]
        mx.xtas[x, t, a] <- exp(logmx.xta[x, t, a])
        logmx.xta[x, t, a] <- beta.tas[t, a, 1] * Yx[x, 1] + 
          beta.tas[t, a, 2] * Yx[x, 2] + 
          beta.tas[t, a, 3] * Yx[x, 3] + 
          u.xta[x, t, a]
        u.xta[x, t, a] ~ dnorm(0, 1)
      }
    }
  }
  
  # Priors
  for (t in 1:T){
    for (a in 1:A){
      for (p in 1:P){
        beta.tas[t,a,p] ~ dunif(0, 1)
      }
    }
  }
})


X <- constants$X
T <- constants$T
A <- constants$A
P <- 3

inits <- list(
  y.xtas = array(0, c(X, T, A)),
  beta.tas = array(0, c(T, A, P)),
  u.xta = array(0, c(X, T, A)),
  mx.xta = array(0, c(X, T, A)),
  mu.xta = array(0, c(X, T, A))
)


mcmc.base <- nimbleMCMC(code = CA_model, 
                        constants = constants,
                        inits = inits,
                        data = data,
                        nchains = 4, niter = 10000,
                        summary = TRUE, WAIC = TRUE,
                        monitors = c('y.xtas', 'mu.xta', 'mx.xtas', 'beta.tas', 'u.xta', 'logmx.xta'))







#### Simple Spatial Model ####
# Preparing adjacency list
library(tigris)
library(rgdal)
library(rgeos)
CA_shp <- counties("California", cb = TRUE)
CA_shp <- CA_shp[!(CA_shp@data$NAME %in% c("Alpine", "Sierra")),] # Drop Sierra and Alpine counties

CA_adj <- gTouches(CA_shp, byid = TRUE)

CA_list <- list()
CA_length <- c()
for (i in 1:nrow(CA_adj)) {
  l <- c()       # empty sublist
  c <- 1            # counter of for list l
  for (j in 1:ncol(CA_adj)) {   # for each row, if the value is true
    if (CA_adj[i, j] == TRUE){
      l[c] <- j     # set sublist l element c equal to j
      c <- c + 1
    }
  }
  
  CA_list[[i]] <- l
  CA_length <- c(CA_length, length(l))
}

CA_weights <- CA_list
for (i in 1:length(CA_list)) {
  for (j in 1:length(CA_list[[i]] )){
    CA_weights[[i]][[j]] <- 1
  }
}

CA_list
CA_weights
CA_length

CA_list <- unlist(CA_list)
CA_weights <- unlist(CA_weights)
CA_list
CA_weights

n.chains <-1 
load("data.RData")
data[["adj"]] <- as.integer(CA_list)
data[["weights"]] <- as.integer(CA_weights)
data[["num"]] <- as.integer(CA_length)
data$Yx <- as.matrix(data$Yx)
data$A <- data$n.a


CA_model_spatial <- nimbleCode({
  for (x in 1:X){
    for (t in 1:T){
      for (a in 1:A){
        y.xtas[x, t, a] ~ dpois(mu.xta[x, t, a])
        mu.xta[x, t, a] <- pop.xtas[x, t, a] * mx.xtas[x, t, a]
        mx.xtas[x, t, a] <- exp(logmx.xta[x, t, a])
        logmx.xta[x, t, a] <- beta.tas[t, a, 1] * Yx[x, 1] + 
          beta.tas[t, a, 2] * Yx[x, 2] + 
          beta.tas[t, a, 3] * Yx[x, 3] + 
          phi[a] + 
          u.xta[x, t, a]
        
        # prior on u 
        u.xta[x, t, a] ~ dnorm(0, 1)
      }
    }
  }
  
  # prior on phi
  phi[1:A] ~ dcar_normal(adj[1:260], weights[1:260], num[1:A], tau, zero_mean = 0)
  tau ~ dgamma(0.5, 0.0005)
  
  # Priors on beta
  for (t in 1:T){
    for (a in 1:A){
      for (p in 1:P){
        beta.tas[t,a,p] ~ dunif(0, 1)
      }
    }
  }
})

# Inits
inits <- list(
  phi = array(0, c(A)),
  y.xtas = array(0, c(X, T, A)),
  beta.tas = array(0, c(T, A, P)),
  u.xta = array(0, c(X, T, A)),
  mx.xta = array(0, c(X, T, A)),
  mu.xta = array(0, c(X, T, A))
)

mcmc.out <- nimbleMCMC(code = CA_model_spatial, 
                       constants = constants,
                       inits = inits,
                       data = data,
                       nchains = 4, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('y.xtas', 'mu.xta', 'mx.xtas', 'beta.tas', 'u.xta', 'logmx.xta', 'phi'))


save(mcmc.base$summary, file = "base.RData")
save(mcmc.out$summary, file = "spatial.RData")

#### Simple Model Comparison ####
# Analysis
library(tidyverse)
library(stringr)
library(tigris)
library(sf)
library(tmap)
load("base.RData")
load("spatial.RData")
base <- data.frame(base$all.chains)
spatial <- data.frame(spatial$all.chains)
base$index <- gsub("\\]*","",  gsub(".*\\[", "", rownames(base)))
base$parname <- gsub("\\[.*","", rownames(base))
spatial$index <- gsub("\\]*","",  gsub(".*\\[", "", rownames(spatial)))
spatial$parname <- gsub("\\[.*","", rownames(spatial))
mx_base <- base[base$parname == "mx.xtas",]
mx_spatial <- spatial[spatial$parname == "mx.xtas",]
mx_base<- mx_base %>% separate(index, into = c("age", "year", "county"), ",") %>%
  mutate(age = as.numeric(age), 
         year = as.numeric(year),
         county = as.numeric(county))
mx_spatial <- mx_spatial %>% separate(index, into = c("age", "year", "county"), ",") %>%
  mutate(age = as.numeric(age), 
         year = as.numeric(year),
         county = as.numeric(county))

# Determine significance for phi
phi <- spatial %>% filter(parname == "phi")
phi[(phi$X95.CI_low > 0 & phi$X95.CI_upp > 0 ) | (phi$X95.CI_low < 0 & phi$X95.CI_upp < 0 ),]
hist(phi$Mean)

# Plot mx for alameda county (county == 1)
ac_base <- mx_base %>% filter(county == 1)
ac_spatial <- mx_spatial %>% filter(county == 1)

ggplot() + 
  geom_point(data = ac_base, aes(x = age, y= log(Mean), group = year), color = "blue") +
  geom_point(data = ac_spatial, aes(x = age, y= log(Mean), group = year), color = "red")

# Infant mortality by state

imr_base <- mx_base %>% 
  filter(age == 1 ) %>% 
  group_by(county) %>% 
  summarize(IMR = mean(Mean))
imr_spatial <- mx_spatial %>% 
  filter(age == 1) %>% 
  group_by(county) %>% 
  summarise(IMR = mean(Mean))

CA_shp <- counties("California", cb = T) %>% st_as_sf() %>% filter(NAME != "Sierra") %>% filter(NAME != "Alpine")
imr_base$NAME <- CA_shp$NAME
imr_spatial$NAME <- CA_shp$NAME

CA_shp <- left_join(CA_shp, imr_base)
CA_shp <- left_join(CA_shp, imr_spatial)
tm_shape(CA_shp) + tm_polygons(col = "IMR")


