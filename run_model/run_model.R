
# select a region
regions <- c(
  "midwest",
  "northeast",
  "southeast",
  "southeast_atlantic",
  "southeast_texas",
  "southwest"
)

region <- regions[6]

# list of city names

# midwest
if(region == regions[1]){
  city_names <- c(
    "Chicago",
    "Denver",
    "Des_Moines",
    "Detroit", 
    "Minneapolis",
    "St_Louis"
  )
}

# northeast
if(region == regions[2]){
  city_names <- c(
    "Boston", 
    "DC",
    "NYC", 
    "Philadelphia"
  )
}

# southeast
if(region == regions[3]){
  city_names <- c(
    "Atlanta",
    "Charlotte",
    "Dallas",
    "Denton",
    "Houston",
    "Raleigh"
  )
}

# southeast_atlantic
if(region == regions[4]){
  city_names <- c(
    "Atlanta",
    "Charlotte",
    "Raleigh"
  )
}

# southeast_texas
if(region == regions[5]){
  city_names <- c(
    "Dallas",
    "Denton",
    "Houston"
  )
}

# southwest
if(region == regions[6]){
  city_names <- c(
    "LA",
    "Phoenix",
    "Riverside",
    "SD",
    "SF"
  )
} 

# or choose one city
# city_names <- "Philadelphia"

min_species_detections <- 2 # binary park/year/species detections
min_site_years_w_detection <- 2 # remove parks never surveyed across repeat years
min_species_for_community_sampling_event <- 1 # if 1 species detected, any other species in same fam could have been 
family_sampling <- TRUE # Should enter either TRUE or FALSE 
remove_outlier_parks <- TRUE # remove very small parks
# family_sampling:
# if false infer sampling event for all butterflies if any butterflies detected
# if true only infer sampling event for butterflies in same family as any butterflies detected

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data.R")

my_data <- prep_data(city_names,
                     min_species_detections,
                     min_site_years_w_detection,
                     min_species_for_community_sampling_event,
                     family_sampling,
                     remove_outlier_parks
)

saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))
my_data <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))

# prepare to fit occupncy model
library(rstan)

# data to feed to the model
# detection data
V <- my_data$V # detections (1==detected)
V_NA <- my_data$V_NA # NAs (0==no known survey effort made)

# covariate data
species_info <- my_data$species_info
n_species <- my_data$n_species # number of species
species <- as.integer(as.factor(species_info$species))
site_data <- my_data$site_data
n_sites <- my_data$n_sites # number of sites
sites <- site_data$multicity_site_id
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
years <- seq(1, n_years_minus1)
n_surveys <- my_data$n_surveys
surveys <- sequence(n_surveys)
surveys <- (surveys - mean(surveys)) / sd(surveys)

## predictors
# species
feature_diversity <- species_info$featureDiversity_scaled
ease_of_id <- species_info$research_grade_proportion_scaled
wingspan <- species_info$aveWingspan_scaled
# site
park_size <- site_data$log_total_green_space_area_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
isolation <- site_data$log_isolation_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
city <- as.integer(as.factor(unique(site_data$city)))
n_cities <- length(unique(city))

## ranges
ranges <- my_data$ranges

# extra stuff
R <- my_data$R
species_integer_vector <- my_data$species_integer_vector
city_integer_vector <- my_data$city_integer_vector
multicity_site_id_vector <- my_data$multicity_site_id_vector
city_id_vector <- my_data$city_id_vector
site_survey_year_vector <- my_data$site_survey_year_vector
prev_index_vector <- my_data$prev_index_vector
confirmed_occurrence <- my_data$confirmed_occurrence

# plot
ggplot(site_data, aes(
  x = log_total_green_space_area_scaled_across_all_cities, y = log_isolation_scaled_across_all_cities, colour = city)) +
  geom_point()


stan_data <- c("R", "n_surveys", "surveys", 
               "V", "V_NA", "ranges", "site_survey_year_vector",
               "n_species", "species", "species_integer_vector",
               "n_sites", "sites", "multicity_site_id_vector",
               "n_cities","city", "city_id_vector",
               "feature_diversity", "ease_of_id", "wingspan",
               "park_size", "isolation",
               "confirmed_occurrence", "prev_index_vector"
               
) 

## Parameters monitored 
params <- c(
  
  "psi1_0", 
  "sigma_psi1_species",
  "sigma_psi1_city",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation",

  "gamma0", 
  "sigma_gamma_species",
  "sigma_gamma_city",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation",

  "phi0", 
  "sigma_phi_species",
  "sigma_phi_city",
  "mu_phi_wingspan",
  "sigma_phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_isolation",
  "sigma_phi_isolation",
  
  "p0", 
  "sigma_p_species",
  "sigma_p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "mu_p_species_date",
  "sigma_p_species_date",
  "mu_p_species_date_sq",
  "sigma_p_species_date_sq",
  
  # city effects
  "psi1_city",
  "psi1_wingspan",
  "psi1_park_size",
  "psi1_isolation",
  "gamma_city",
  "gamma_wingspan",
  "gamma_park_size",
  "gamma_isolation",
  "phi_city",
  "phi_wingspan",
  "phi_park_size",
  "phi_isolation",
  "p_city",
  
  # species effects and PPC
  #"W_species_rep",
  "psi1_species",
  "gamma_species", "phi_species",
  "p_species"
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.97

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 1, 2),
       psi1_wingspan = runif(1, -1, 1),
       mu_psi1_park_size = runif(1, 0, 1),
       gamma0 = runif(1, -3, -2),
       gamma_wingspan = runif(1, 1, 2),
       mu_gamma_park_size = runif(1, 0, 1),
       #gamma_isolation = runif(1, -2, -1),
       phi0 = runif(1, 2, 3),
       phi_wingspan = runif(1, -1, 0),
       mu_phi_park_size = runif(1, 0, 1),
       #phi_isolation = runif(1, 1, 2),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 1, 2),
       sigma_p_city = runif(1, 0, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 0),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)


## --------------------------------------------------
### Run model

stan_model <- "./models/dynamic_occupancy_model.stan"

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

saveRDS(stan_out, paste0("./model_outputs/stan_out_", region, ".rds"))

stan_out <- readRDS( paste0("./model_outputs/stan_out_", region, ".rds"))

print(stan_out, digits = 3, 
      pars = c(
        "psi1_0", 
        "sigma_psi1_species",
        "sigma_psi1_city",
        "mu_psi1_wingspan",
        "sigma_psi1_wingspan",
        "mu_psi1_park_size",
        "sigma_psi1_park_size",
        "mu_psi1_isolation",
        "sigma_psi1_isolation",
        
        "gamma0", 
        "sigma_gamma_species",
        "sigma_gamma_city",
        "mu_gamma_wingspan",
        "sigma_gamma_wingspan",
        "mu_gamma_park_size",
        "sigma_gamma_park_size",
        "mu_gamma_isolation",
        "sigma_gamma_isolation",
        
        "phi0", 
        "sigma_phi_species",
        "sigma_phi_city",
        "mu_phi_wingspan",
        "sigma_phi_wingspan",
        "mu_phi_park_size",
        "sigma_phi_park_size",
        "mu_phi_isolation",
        "sigma_phi_isolation",
        
        "p0", 
        "sigma_p_species",
        "sigma_p_city",
        "p_wingspan",
        "p_feature_diversity",
        "p_ease_of_id",
        "mu_p_species_date",
        "sigma_p_species_date",
        "mu_p_species_date_sq",
        "sigma_p_species_date_sq",
        
        # city effects
        "psi1_city",
        "psi1_wingspan",
        "psi1_park_size",
        "psi1_isolation",
        "gamma_city",
        "gamma_wingspan",
        "gamma_park_size",
        "gamma_isolation",
        "phi_city",
        "phi_wingspan",
        "phi_park_size",
        "phi_isolation",
        "p_city"
      ))

print(stan_out, digits = 3, 
      pars = c("psi1_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "psi1_0", 
  "psi1_city",
  "sigma_psi1_city",
  "sigma_psi1_species",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation",
  "psi1_park_size",
  "psi1_isolation"
  
))

traceplot(stan_out, pars = c(
  "gamma0", 
  "gamma_city",
  "sigma_gamma_city",
  "sigma_gamma_species",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation",
  "gamma_park_size",
  "gamma_isolation"
))  

traceplot(stan_out, pars = c(
  "phi0", 
  "phi_city",
  "sigma_phi_city",
  "sigma_phi_species",
  "mu_phi_wingspan",
  "sigma_phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_isolation",
  "sigma_phi_isolation",
  "phi_park_size",
  "phi_isolation"
))

traceplot(stan_out, pars = c(
  "p0", 
  "sigma_p_species",
  "p_city",
  "sigma_p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "mu_p_species_date",
  "sigma_p_species_date",
  "mu_p_species_date_sq",
  "sigma_p_species_date_sq"
))

