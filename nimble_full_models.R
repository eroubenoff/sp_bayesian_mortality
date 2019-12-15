library(nimble)
library(tigris)
library(rgeos)
library(dplyr)


#### Base Full Model ####
rm(list = ls())
CA_model_full <- nimbleCode({
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
          u.xta[x, t, a]
      }
    }
  }
  
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

load("data.RData")
data$Yx <- as.matrix(data$Yx)
A <- data$n.a
X <- data$X
P <- 3
S <- 1
T <- data$T
data$n.a <- NULL
data$n.amax <- NULL
data$X <- NULL
data$P <- NULL
data$S <- NULL
data$T <- NULL
constants <- list(A = A, X = X, P = P, S = S, T = T)
data$y.xta <- data$y.xtas
data$pop.xta <- data$pop.xtas
data$y.xtas <- NULL
data$pop.xtas <- NULL
str(data)

monitors <- c("y.xta", "mu.xta", "mx.xta", "logmx.xta", "beta.ta", "u.xta")

# Inits
inits <- list(
  y.xta = array(0, c(X, T, A)),
  mu.xta = array(0, c(X, T, A)),
  mx.xta = array(0, c(X, T, A)),
  logmx.xta = array(1, c(X, T, A)),
  beta.ta = array(0, c(T, A, P)),
  mu.beta = array(0, c(T, P)),
  sigma.beta = array(1, c(T, P)),
  tau.beta = array(1, c(T, P)),
  u.xta = array(0, c(X, T, A)),
  sigma.mu = array(20, c(P)),
  sigma.u = array(0.05, c(X))
)


mcmc.out <- nimbleMCMC(code = CA_model_full, 
                       constants = constants,
                       inits = inits,
                       data = data,
                       nchains = 4, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = monitors)

full_base <- mcmc.out$summary
save(full_base, file = "full_base.RData")





#### Full Spatial Model ####
rm(list = ls())
library(nimble)
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
      phi[x, t, 1:A] ~ dcar_proper(adj[1:N], weights[1:N], num[1:A], tau.p[x, t], zero_mean = 0)
      tau.p[x, t] ~ dgamma(0.001, 0.001)
    }
  }
  
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

load("data.RData")
data$Yx <- as.matrix(data$Yx)
A <- data$n.a
X <- data$X
P <- 3
S <- 1
T <- data$T
data$n.a <- NULL
data$n.amax <- NULL
data$X <- NULL
data$P <- NULL
data$S <- NULL
data$T <- NULL
constants <- list(A = A, X = X, P = P, S = S, T = T)
data$y.xta <- data$y.xtas
data$pop.xta <- data$pop.xtas
data$y.xtas <- NULL
data$pop.xtas <- NULL
str(data)

monitors <- c("y.xta", "mu.xta", "mx.xta", "logmx.xta", "beta.ta", "u.xta", "phi")

library(tigris)
library(rgeos)
CA_shp <- counties("California", cb = TRUE)
CA_shp <- CA_shp[!(CA_shp@data$NAME %in% c("Alpine", "Sierra")),] # Drop Sierra and Alpine counties


CA_shp <- CA_shp[order(CA_shp$NAME), ]

CA_adj <- gTouches(CA_shp, byid = TRUE, returnDense = FALSE)
# check: alameda is county 1 and borders  15 26 37 42 49
# 15: kings, 26: Monterey, 37:San Francisco, 42: Santa Clara, 49: Sutter
# 20, 25, 763, 1116, 3180 (contra costa, santa clara, stanislaus, SF, san joaquin)
CA_adj[1] 
CA_shp$NAME


CA_list <- CA_adj

CA_weights <- CA_list
for (i in 1:length(CA_list)) {
  for (j in 1:length(CA_list[[i]] )){
    CA_weights[[i]][[j]] <- 1
  }
}



CA_list <- unname(CA_list)
CA_weights <- unname(CA_weights)
CA_length <- sapply(CA_list, function(x) length(x))
CA_list
CA_weights
CA_length

CA_list <- unlist(CA_list)
CA_weights <- unlist(CA_weights)
CA_list
CA_weights
CA_length

n.chains <-1 
load("data.RData")
data$Yx <- as.matrix(data$Yx)
A <- data$n.a
X <- data$X
P <- 3
S <- 1
T <- data$T
N <- length(CA_list)
data$n.a <- NULL
data$n.amax <- NULL
data$X <- NULL
data$P <- NULL
data$S <- NULL
data$T <- NULL
constants <- list(A = A, X = X, P = P, S = S, T = T, N= N)
data$y.xta <- data$y.xtas
data$pop.xta <- data$pop.xtas
data$adj <- CA_list
data$weights <- CA_weights
data$num <- CA_length
data$y.xtas <- NULL
data$pop.xtas <- NULL
str(data)
str(constants)


# Inits
inits <- list(
  y.xta = array(0, c(X, T, A)),
  mu.xta = array(0, c(X, T, A)),
  mx.xta = array(0, c(X, T, A)),
  logmx.xta = array(1, c(X, T, A)),
  beta.ta = array(0, c(T, A, P)),
  mu.beta = array(0, c(T, P)),
  sigma.beta = array(1, c(T, P)),
  tau.beta = array(1, c(T, P)),
  u.xta = array(0, c(X, T, A)),
  sigma.mu = array(20, c(P)),
  sigma.u = array(0.05, c(X)),
  phi = array(0, c(X, T, A)),
  tau.p = array(0.01, c(X, T))
)


mcmc.out <- nimbleMCMC(code = CA_model_full_spatial, 
                       constants = constants,
                       inits = inits,
                       data = data,
                       nchains = 4, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = monitors)
full_spatial <- mcmc.out$summary
save(full_spatial, file = "full_spatial.RData")
#### end ####

# Took 4:16 to run

#### analysis ####
library(tidyverse)
library(stringr)
library(tigris)
library(sf)
library(tmap)
load("full_base.RData")
load("full_spatial.RData")
base <- data.frame(full_base$all.chains)
spatial <- data.frame(full_spatial$all.chains)
base$index <- gsub("\\]*","",  gsub(".*\\[", "", rownames(base)))
base$parname <- gsub("\\[.*","", rownames(base))
spatial$index <- gsub("\\]*","",  gsub(".*\\[", "", rownames(spatial)))
spatial$parname <- gsub("\\[.*","", rownames(spatial))


mx_base <- base[base$parname == "mx.xta",]
mx_spatial <- spatial[spatial$parname == "mx.xta",]
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

CA_shp <- CA_shp[order(CA_shp$NAME),]

# Plot mx for alameda county (county == 1)
alameda_county <- mx_base %>% 
  filter(county == 1, year %in% c(1, 8, 18))%>%
  mutate(model = "base") %>%
  bind_rows(
    mx_spatial %>% 
      filter(county == 1, year %in% c(1, 8, 18)) %>%
      mutate(model = "spatial")
  ) %>%
  mutate(year = if_else(year == 1, 1999, if_else(year == 8, 2006, if_else(year == 18, 2016, 0))))
ggplot(alameda_county) + 
  geom_line( aes(x = age, y= log(Mean), color = model)) +
  scale_x_continuous(name = "Age", breaks = c(1:13), labels =  c(0, 1, 5, 10, 15, 20, 
                                                                 25,35, 
                                                                 45, 55, 
                                                                 65, 75, 85)) +
  ylab("log(Mx)") + 
  ggtitle("Age-Specific Log Mortality Rates for Alameda County in Selected Years") + 
  facet_grid(~year)

# Alameda: 1, SF: 37, Contra Costa: 6, San Mateo: 40, Santa Clara: 42
ggplot() + 
  geom_line(data = mx_base %>% filter(county %in% c(1, 37, 6, 40, 42), year == 1), 
            aes(x = age, y= log(Mean)), color = "blue") +
  geom_line(data = mx_spatial %>% filter(county %in% c(1, 37, 6, 40, 42), year == 1), 
            aes(x = age, y= log(Mean)), color = "red") +
  facet_grid(~county)


# Infant mortality by state

imr_base <- mx_base %>% 
  filter(age == 1 ) %>% 
  group_by(county) %>% 
  summarize(IMR_base = mean(Mean))
imr_spatial <- mx_spatial %>% 
  filter(age == 1) %>% 
  group_by(county) %>% 
  summarise(IMR_spatial = mean(Mean))

CA_shp <- counties("California", cb = T) %>% st_as_sf() %>% filter(NAME != "Sierra") %>% filter(NAME != "Alpine")
imr_base$NAME <- CA_shp$NAME
imr_spatial$NAME <- CA_shp$NAME

CA_shp <- left_join(CA_shp, imr_base)
CA_shp <- left_join(CA_shp, imr_spatial)
breaks <- c(0, 0.010,0.020, 0.030, 0.040, 0.050)
base_map <- tm_shape(CA_shp) + tm_polygons(col = "IMR_base", style = "fixed", breaks = breaks)
spatial_map <- tm_shape(CA_shp) + tm_polygons(col = "IMR_spatial", style = "fixed", breaks = breaks)
tmap_arrange(base_map, spatial_map)


# Simple test of means difference
comparison <- mx_base %>% 
  left_join(mx_spatial, by = c("parname", "age", "county", "year")) %>%
  select(Mean.x, Mean.y)
t.test(comparison$Mean.x, comparison$Mean.y, paired = T)
hist(spatial %>% filter(parname == "phi") %>% select("Mean"))


# Test autocorrellated errors in base model
library(sf)
library(tmap)
library(tigris)
library(spdep)
CA_shp <- counties("California", cb = T) %>% 
  st_as_sf() %>% 
  filter(NAME != "Sierra") %>% 
  filter(NAME != "Alpine") %>%
  arrange(NAME)
base_errors <- base %>% 
  filter(parname == "u.xta") %>% 
  separate(index, into = c("age", "year", "county"), ",") %>%
  mutate(age = as.numeric(age), 
         year = as.numeric(year),
         county = as.numeric(county))
# Take mean total error 
base_errors <- base_errors %>%
  group_by(county) %>%
  summarize(Mean = mean(Mean))
# Join with shapefile
CA_shp$Mean_error <- base_errors$Mean
tm_shape(CA_shp) + 
  tm_polygons(col = "Mean_error")
# Gracious Moran's function from https://mgimond.github.io/Spatial/spatial-autocorrelation-in-r.html
CA_listw <- CA_shp %>% poly2nb() %>% nb2listw()
CA_x <- CA_shp$Mean_error
moran_I <- moran.test(x = CA_x,listw =  CA_listw)
geary_C <- geary.test(x = CA_x,listw =  CA_listw)
moran_I
geary_C


# For the report
library(stargazer)
library(reshape2)
library(magrittr)
y.xta_out <- data$y.xta
ages <- c("< 1 year" ,   "1-4 years",   "5-9 years",   "10-14 years", "15-19 years", "20-24 years", "25-34 years", "35-44 years", "45-54 years", "55-64 years", "65-74 years", "75-84 years", "85+ years" )
years <- c(1999:2016)
counties <- c(CA_shp$NAME)
dimnames(y.xta_out) <- list(ages, years, counties)
str(y.xta_out)
y.xta_out <- y.xta_out %>% melt() 
y.xta_out %<>% rename(Age = Var1, Year = Var2, County = Var3, Deaths = value)

pop.xta_out <- data$pop.xta 
dimnames(pop.xta_out) <- list(ages, years, counties)
pop.xta_out <- pop.xta_out %>% melt()
pop.xta_out %<>% rename(Age = Var1, Year = Var2, County = Var3, Pop= value)


df_out <- left_join(y.xta_out, pop.xta_out)
stargazer(rbind(head(df_out), tail(df_out)), summary = FALSE)

pcs_out <- data$Yx
dimnames(pcs_out) <- list(ages, c("PC1", "PC2", "PC3"))
stargazer(pcs_out)

pcs_out$Age <- ages
pcs_out %<>% as.data.frame(pcs_out)
ggplot(pcs_out) + 
  geom_path(aes(x= Age, y = PC1))




