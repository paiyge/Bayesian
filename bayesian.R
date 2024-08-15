## Hierarchical Movement Modeling for BRI Quantitative Ecologist Position
## 
## Data : https://www.kaggle.com/datasets/saurabhshahane/predicting-animal-behavior-using-gps/data
## 
## Sources:
## Bayesian Hierarchical Models in Ecology by Steve Midway https://bookdown.org/steve_midway/BHME/#
## https://github.com/andrewcparnell/jags_examples
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2894957/

library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(rgdal)
library(sp)
library(raster)
library(rjags)
library(ggmap)
library(bayesplot)
library(psych)
library(rstanarm)

install.packages("rstanarm")
yes# Load the data
data <- read.csv('/Users/paigenorris/bayeisan/gps/anon_gps_tracks_with_dive.csv')

# EDA

summary(data)

# Check out missing data
missing_data <- data %>%
  summarise_all(~ sum(is.na(.))) %>%
  gather(key = "column", value = "missing_count")

ggplot(missing_data, aes(x = reorder(column, -missing_count), y = missing_count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Missing Data by Column", x = "Columns", y = "Number of Missing Values") +
  theme_minimal()

# 1000 out of 260000 is not detrimental, will delete
data <- data[!is.na(data$coverage_ratio), ]

# Data Pre-processing
data <- data %>%
  mutate(date_time = ymd_hms(date_time),,
         year = as.numeric(gsub("t", "", year)),
         species = factor(species, levels = c("tRAZO", "tCOGU", "tEUSH")),
         is_dive = as.logical(is_dive))

# Define the model (note this should be done in an external .txt file but for
# the sake of simplicity, everything is included in a single script.)
# Originally modeled on each individual bird but realized during convergence that 
# was likely a bad idea as it was difficult to garner insight.
model_string <- "
model {
  for (i in 1:N) {
    ## Observation models for latitude and longitude, assuming normal distributions with mean 0 and variance 100
    # Nature is stochastic with error, variability, and chance
    # stochastic = indexing via y[i] ∼ dnorm(mu, tau)
    
    lat[i] ~ dnorm(mu_lat[i], tau)
    lon[i] ~ dnorm(mu_lon[i], tau)
    
    ## Process model
    mu_lat[i] <- alpha_lat[bird[i]] + beta_lat * habitat[i]
    mu_lon[i] <- alpha_lon[bird[i]] + beta_lon * habitat[i]
  }
  
  ## Priors
  
  # for each bird species
  for (j in 1:max(bird)) {
    alpha_lat[j] ~ dnorm(0, 0.01)
    alpha_lon[j] ~ dnorm(0, 0.01)
  }
  
  # for location
  beta_lat ~ dnorm(0, 0.01)
  beta_lon ~ dnorm(0, 0.01)
  tau ~ dgamma(1, 1) # inverse of variance from normal distribution -> try to reduce over-parameterization 
}
"

# Bundle data
jags_data <- list(
  lat = data$lat,
  lon = data$lon,
  colony = data$colony2,
  species = data$species,
  bird = as.numeric(data$bird), # JAGS handles indices better as numbers
  habitat = data$coverage_ratio,
  N = nrow(data)
)

# Initial values for each chain
n.chains <- 3
inits <- lapply(1:n.chains, function(i) {
  list(
    alpha_lat = rnorm(nlevels(factor(data$bird))), 
    alpha_lon = rnorm(nlevels(factor(data$bird))),
    beta_lat = rnorm(1),
    beta_lon = rnorm(1),
    tau = 1
  )
})

# Parameters
parameters <- c("alpha_lat",
                "alpha_lon", 
                "beta_lat", 
                "beta_lon", 
                "tau")

# Run the model
model <- jags.model(textConnection(model_string), data = jags_data, inits = inits, n.chains = n.chains)
update(model, 1000)  # Burn-in
samples <- coda.samples(model, parameters, n.iter = 10000)

###_________
# Checking convergence - this needs work, not yet certain how to visualize correctly
# Might have to utilize proper JAGS syntax, ie separate files

# Plot trace plots for all parameters
plot(samples)
posterior <- as.matrix(samples)
plot_title <- ggtitle("Posterior Distributions of Bird Movement Parameters",
                      "with medians and 80% intervals")
params_to_plot <- c("alpha_lat[i]", "alphs_lon[i]", "beta_lat", "beta_lon", "tau")
mcmc_areas(posterior,
           pars = params_to_plot,
           prob = 0.8) + plot_title

# Goal: 
# Check the output - are the true values inside the 95% CI?
# Also look at the R-hat values - they need to be close to 1 if convergence has been achieved
# Create a plot of the posterior mean regression line
# Maybe a boxplot?


################
## Mapping Visuals 
library(ggmap)
register_google(key = "AIzaSyC25ay9IXSvxGo5rognii8Niv0oKA7ClhE")


# Get a map of the area via Google API
g_map <- get_map(location = c(
  lon = mean(data$lon), 
  lat = mean(data$lat)), 
  zoom = 8, 
  maptype = "terrain")

# Tracks by individual birds
ggmap(map) +
  geom_path(data = data, aes(x = lon, y = lat, color = factor(bird)), size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Bird Movement", x = "Longitude", y = "Latitude")

# Tracks by colonies
ggmap(map) +
  geom_path(data = data, aes(x = lon, y = lat, color = (colony2)), size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Colony Movement", x = "Longitude", y = "Latitude")

#_______
