library(rstan)

# for community sampling events inferred by [taxonomic] order, source this file:
# source("./run_model/prep_data_iNat_w_community_sampling_events.R")

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data_biodiversity_w_land_cover.R")

min_species_detections <- 100
min_species_for_community_sampling_event <- 1
grid_size <- 5000

my_data <- prep_data(min_species_detections,
                     min_species_for_community_sampling_event,
                     grid_size)

# save or load data for later
saveRDS(my_data, "./run_model/prepped_data/prepped_data.RDS")
my_data <- readRDS("./run_model/prepped_data/prepped_data.RDS")

# data to feed to the model
V <- my_data$V 
V_NA <- my_data$V_NA

# for species names, note that n == number of unique BINARY detections
# while our filter min_species_detections is for total number of detections.
species_names <- my_data$species_names 

n_species <- my_data$n_species # number of species
n_families <- my_data$n_families # number of families
n_sites <- my_data$n_sites # number of sites
n_cities <- my_data$n_cities # number of cities
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
n_surveys <- my_data$n_surveys

species_character <- my_data$species
sites_character <- my_data$sites
species <- as.integer(as.factor(my_data$species))
family_lookup <- as.integer(as.factor(my_data$families))
sites <- as.integer(as.factor(my_data$sites))
cities = as.integer(as.factor(my_data$cities))
years_full <- as.integer(my_data$years)
years <- seq(1, n_years_minus1)
surveys_raw <- as.integer(my_data$surveys)
surveys <- (surveys_raw - mean(surveys_raw)) / sd(surveys_raw)

# do this for now
n_families <- length(unique(family_lookup)) # number of families

# n (from species_names) is the total number of detections (not unique site visit detections)
View(cbind(species_names, as.data.frame(family_lookup)))

stan_data <- c("V", "V_NA", 
               "species", "sites", "years", "surveys", 
               "cities", "n_cities",
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys"
) 

## Parameters monitored 
params <- c("psi1_0", 
            "sigma_psi1_species",
            
            "gamma0", 
            "sigma_gamma_species",
            
            "phi0", 
            "sigma_phi_species",
            
            "p0", 
            "sigma_p_species",
            #"sigma_p_family",
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

stan_model <- "./models/dynamic_occupancy_model3_multicity.stan"

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

saveRDS(stan_out, "./model_outputs/stan_out3.2_w_all_species.rds")
#stan_out <- readRDS("./model_outputs/stan_out4.rds")

print(stan_out, digits = 3, 
      pars = c("psi1_0", "sigma_psi1_species",
               
               "gamma0", "sigma_gamma_species",
               
               "phi0", "sigma_phi_species",
               
               "p0", "sigma_p_species",
               
               "mu_p_species_date", "sigma_p_species_date",
               "mu_p_species_date_sq", "sigma_p_species_date_sq"
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
