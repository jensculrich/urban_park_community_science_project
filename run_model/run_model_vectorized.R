# for community sampling events inferred by [taxonomic] order, source this file:
# source("./run_model/prep_data_iNat_w_community_sampling_events.R")

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data_vectorized.R")

min_species_detections <- 1 # binary park/year/species detections
min_species_per_family <- 2
min_species_for_community_sampling_event = 1
city_name <- "los_angeles"

my_data <- prep_data(min_species_detections,
                     min_species_per_family,
                     min_species_for_community_sampling_event,
                     city_name
                     )

#saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))
my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))

# prepare to fit occupncy model
library(rstan)

# data to feed to the model
V <- my_data$V_detections
R <- nrow(V)
V_NA <- my_data$V_NA

covariate_data <- my_data$covariate_data
n_sites <- my_data$n_sites # number of sites
site <- covariate_data$site
n_species <- my_data$n_species # number of species
species <- covariate_data$species_number
n_years <- my_data$n_years # number of surveys
year <- covariate_data$year
n_surveys <- my_data$n_surveys
survey <- sequence(n_surveys)
survey <- (survey - mean(survey)) / sd(survey)

# categorical year dummy variables
X_year <- model.matrix(~ as.factor(year), data = covariate_data)

## predictors
# species
feature_diversity <- covariate_data$featureDiversity_scaled
ease_of_id <- covariate_data$research_grade_proportion_scaled
wingspan <- covariate_data$aveWingspan_scaled
# site
park_size <- covariate_data$park_size_scaled

species_info <- my_data$species_info 
site_data <- my_data$site_data

stan_data <- c("V", "V_NA", "R",
               "n_species", "n_sites", "n_years", "n_surveys",
               "species", "site", "X_year", "survey", 
               "feature_diversity", "ease_of_id", "wingspan",
               "park_size"
) 

## Parameters monitored 
params <- c(#"psi_0", 
            "sigma_psi_species",
            "sigma_psi_site",
            "psi_year",
            "psi_wingspan",
            "psi_park_size",
            
            #"p0", 
            "sigma_p_species",
            "sigma_p_site",
            "p_year",
            "p_wingspan",
            "p_feature_diversity",
            "p_ease_of_id",
            "p_park_size",
            "mu_p_species_date",
            "sigma_p_species_date",
            "mu_p_species_date_sq",
            "sigma_p_species_date_sq",
            
            "psi_species", "p_species",
            "psi_site", "p_site",
            "p_date", "p_date_sq"
            )

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(#psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0, 1),
       sigma_psi_site = runif(1, 0, 1),
       #psi_year = runif(1, -1, 1),
       psi_wingspan = runif(1, -1, 1),
       psi_park_size = runif(1, -1, 1),
       #p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       p_park_size = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 0),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

## --------------------------------------------------
### Run model

stan_model <- "./models/static_occupancy_model.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 control=list(adapt_delta=delta),
                 open_progress = FALSE,
                 cores = n_cores)

saveRDS(stan_out, "./model_outputs/stan_out.rds")
#stan_out <- readRDS("./model_outputs/stan_out.rds")

print(stan_out, digits = 3, 
      pars = c(#"psi_0", 
               "sigma_psi_species",
               "sigma_psi_site",
               "psi_year",
               "psi_wingspan",
               "psi_park_size",
               
               #"p0", 
               "sigma_p_species",
               "sigma_p_site",
               "p_year",
               "p_wingspan",
               "p_feature_diversity",
               "p_ease_of_id",
               "p_park_size",
               "mu_p_species_date",
               "sigma_p_species_date",
               "mu_p_species_date_sq",
               "sigma_p_species_date_sq"
      ))

print(stan_out, digits = 3, 
      pars = c("p_family"
      ))

print(stan_out, digits = 3, 
      pars = c("p_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "sigma_psi_site"
))
