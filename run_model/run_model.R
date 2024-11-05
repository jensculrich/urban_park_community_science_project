# for community sampling events inferred by [taxonomic] order, source this file:
# source("./run_model/prep_data_iNat_w_community_sampling_events.R")

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data.R")

min_species_detections <- 1 # binary park/year/species detections
min_species_per_family <- 2
min_species_for_community_sampling_event = 1
city_name <- "los_angeles"

set.seed(1)
my_data <- prep_data(min_species_detections,
                     min_species_per_family,
                     min_species_for_community_sampling_event,
                     city_name
                     )

saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))
my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))

# prepare to fit occupncy model
library(rstan)

# data to feed to the model
V <- my_data$V 
V_NA <- my_data$V_NA

species_names <- my_data$species_names 
# n (from species_names) is the total number of detections (not unique site visit detections)

n_species <- my_data$n_species # number of species
n_families <- my_data$n_families # number of families
n_sites <- my_data$n_sites # number of sites
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
n_surveys <- my_data$n_surveys

species_character <- my_data$species
sites_character <- my_data$sites
species <- as.integer(as.factor(my_data$species))
family_lookup <- as.integer(as.factor(my_data$families))
sites <- as.integer(as.factor(my_data$sites))
years_full <- as.integer(my_data$years)
years <- seq(1, n_years_minus1)
surveys_raw <- as.integer(my_data$surveys)
surveys <- (surveys_raw - mean(surveys_raw)) / sd(surveys_raw)

# do this for now
n_families <- length(unique(family_lookup)) # number of families

# predictors
migratory <- my_data$migratory
park_size_scaled <- my_data$park_size_scaled

View(cbind(species_names, as.data.frame(family_lookup)))


stan_data <- c("V", "V_NA", "species", "sites", "years", "surveys", 
               #"n_families", "family_lookup",
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys",
               "migratory", "park_size_scaled"
) 

## Parameters monitored 
params <- c("psi1_0", 
            "sigma_psi1_species",
            "psi1_migratory",
            "psi1_park_size",
            
            "gamma0", 
            "sigma_gamma_species",
            "gamma_migratory",
            "gamma_park_size",
            
            "phi0", 
            "sigma_phi_species",
            "phi_migratory",
            "phi_park_size",
            
            "p0", 
            "sigma_p_species",
            #"sigma_p_family",
            "p_migratory",
            "p_park_size",
            "mu_p_species_date",
            "sigma_p_species_date",
            "mu_p_species_date_sq",
            "sigma_p_species_date_sq",
            
            "W_species_rep",
            "psi1_species", "gamma_species", "phi_species", "p_species"#, "p_family"
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
  
  list(psi1_0 = runif(1, -1, 1),
       
       gamma0 = runif(1, -1, 0),
       
       phi0 = runif(1, 1, 2),
       
       p0 = runif(1, -1, 1)
       
  )
)

## --------------------------------------------------
### Run model

stan_model <- "./models/dynamic_occupancy_model5.stan"

## Call Stan from R
set.seed(1)
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

saveRDS(stan_out, "./model_outputs/stan_out5.2_w_all_species.rds")
#stan_out <- readRDS("./model_outputs/stan_out4.rds")

print(stan_out, digits = 3, 
      pars = c("psi1_0", 
               "sigma_psi1_species",
               "psi1_migratory",
               "psi1_park_size",
               
               "gamma0", 
               "sigma_gamma_species",
               "gamma_migratory",
               "gamma_park_size",
               
               "phi0", 
               "sigma_phi_species",
               "phi_migratory",
               "phi_park_size",
               
               "p0", 
               "sigma_p_species",
               #"sigma_p_family",
               "p_migratory",
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
  "psi1_0","sigma_psi1_species",
  "gamma0", "sigma_gamma_species",
  "phi0", "sigma_phi_species",
  "p0", "sigma_p_species"
))
