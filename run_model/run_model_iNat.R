library(rstan)

source("./run_model/prep_data_iNat.R")

min_species_detections <- 30
min_park_size_acres <- 50 # acres
max_park_size_acres <- 400 # acres
buffer_distance <- 200 # meters

my_data <- prep_data(min_species_detections,
                     min_park_size_acres,
                     max_park_size_acres,
                     buffer_distance)

# data to feed to the model
V <- my_data$V 

species_names <- my_data$species_names

n_species <- my_data$n_species # number of species
n_sites <- my_data$n_sites # number of sites
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
n_surveys <- my_data$n_surveys

species_character <- my_data$species
sites_character <- my_data$sites
species <- as.integer(as.factor(my_data$species))
sites <- as.integer(as.factor(my_data$sites))
years_full <- as.integer(my_data$years)
years <- seq(1, n_years_minus1)
surveys_raw <- as.integer(my_data$surveys)
surveys <- (surveys_raw - mean(surveys_raw)) / sd(surveys_raw)

stan_data <- c("V", "species", "sites", "years", "surveys", 
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
            "p_date",
            "p_date_sq",
            #"mu_p_species_date",
            #"sigma_p_species_date",
            #"mu_p_species_date_sq",
            #"sigma_p_species_date_sq",
            
            "W_species_rep",
            "psi1_species", "gamma_species", "phi_species", "p_species")

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

stan_model <- "./models/dynamic_occupancy_model2.stan"

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

saveRDS(stan_out, "./model_outputs/stan_out.rds")

print(stan_out, digits = 3, 
      pars = c("psi1_0", "sigma_psi1_species",
               
               "gamma0", "sigma_gamma_species",
               
               "phi0", "sigma_phi_species",
               
               "p0", "sigma_p_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "psi1_0","sigma_psi1_species",
  "gamma0", "sigma_gamma_species",
  "phi0", "sigma_phi_species",
  "p0", "sigma_p_species"
))
