# plot effects of site predictors on local species richness

library(tidyverse)
library(vegan)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

#size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools.csv")

n_cities <- length(city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Raleigh",
  "SD",
  "SF"
))

##------------------------------------------------------------------------------
# get model fit

## get param estimates from m2.1
stan_out_m2.1 <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_jan2.rds")

# summarise all variables with default and additional summary measures
estimates <- stan_out_m2.1$draws(
  variables = c(
    "psi_0", 
    "sigma_psi_species",
    "sigma_psi_city",
    "mu_psi_wingspan",
    "sigma_psi_wingspan",
    "mu_psi_park_size",
    "sigma_psi_park_size",
    "mu_psi_tree_cover",
    "sigma_psi_tree_cover",
    "mu_psi_plant_diversity",
    "sigma_psi_plant_diversity",
    "mu_psi_landscape_isolation",
    "sigma_psi_landscape_isolation",
    "mu_psi_landscape_grassherb",
    "sigma_psi_landscape_grassherb",
    "mu_psi_landscape_woody",
    "sigma_psi_landscape_woody",
    
    "p0", 
    "sigma_p_species",
    "sigma_p_city",
    "p_city_detections",
    "p_wingspan",
    "p_feature_diversity",
    "p_ease_of_id",
    "delta0",
    "delta_regional_cluster",
    "sigma_p_species_date",
    "epsilon0",
    "epsilon_regional_cluster",
    "sigma_p_species_date_sq",
    
    # city effects
    "psi_city",
    "psi_wingspan",
    "psi_park_size",
    "psi_tree_cover",
    "psi_plant_diversity",
    "psi_landscape_isolation",
    "psi_landscape_grassherb",
    "psi_landscape_woody",
    "psi_species"),
  
  format = "draws_matrix"
)

# clear space
rm(stan_out_m2.1)
gc()

## --------------------------------------------------
## predict species richness in real parks

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
# Here is what we want to do:

# 1
# get a list of all the parks in each city
# with corresponding info on site covs
# for each species that has an overlapping range
# predict whether or not it occurs

# 2 
# summarize average species richness per park / size regional species pool
# summarize beta diversity across parks
# summarize total diversity / size of regional species pool

## --------------------------------------------------
# First get the data from all sites and place on the same scale fed to the model

## get all the site data from all cities (for sites that were actually included in the model)
my_data <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_all.rds"))

site_data <- my_data$site_data
  
# source the prep function
source("./part2_local_landscape_predictors_of_occupancy/run_model/get_site_data.R")

site_data_all <- get_site_data(city_names)
gc()

# we need to get the data from ALL sites (not just those that were used to fit the model)
# onto the prediction scale that was used to fit the model parameter estimates

# park size
# scale the data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_size <- sd(site_data$log_total_green_space_area) # sd of size of parks used to create scale
mean_size <- mean(site_data$log_total_green_space_area) # mean size of parks used to create scale
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_size_pred <- (site_data_all$log_total_green_space_area - mean_size) / sd_size 

# tree cover
# scale the data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_tree_cover <- sd(site_data$tree_percent_cover) 
mean_tree_cover <- mean(site_data$tree_percent_cover) 
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_tree_cover_pred <- (site_data_all$tree_percent_cover - mean_tree_cover) / sd_tree_cover 

# plant diversity
# scale the data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_plant_diversity <- sd(log(site_data$n_plant_genera + 1) / log(site_data$total_area_sqm))
mean_plant_diversity <- mean(log(site_data$n_plant_genera + 1) / log(site_data$total_area_sqm)) 
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_plant_diversity_pred <- ( log(site_data_all$n_plant_genera + 1) / log(site_data_all$total_area_sqm)
                                            - mean_plant_diversity) / sd_plant_diversity 

# landscape isolation
# scale the isolation data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_isolation <- sd(log(site_data$isolation)) # sd of size of parks used to create scale
mean_isolation <- mean(log(site_data$isolation)) # mean size of parks used to create scale
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_isolation_pred <- (log(site_data_all$isolation) - mean_isolation) / sd_isolation 

# landscape grass / herb cover
# scale the data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_landscape_grassherb <- sd(site_data$proportion_landscape_grassherb) 
mean_landscape_grassherb <- mean(site_data$proportion_landscape_grassherb) 
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_landscape_grassherb_pred <- (site_data_all$proportion_landscape_grassherb - mean_landscape_grassherb) / sd_landscape_grassherb

# landscape woody cover
# scale the data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_landscape_woody <- sd(site_data$proportion_landscape_woody) 
mean_landscape_woody <- mean(site_data$proportion_landscape_woody) 
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_landscape_woody_pred <- (site_data_all$proportion_landscape_woody - mean_landscape_woody) / sd_landscape_woody 


## --------------------------------------------------
# Join the site data with all species that occur at each of those sites
# i.e. all species for which the city is "in range"

# join species list by in range and add species identifiers:
# species cluster integer vector identifier
# and the species integer vector identifier

# get range data so we properly omit species that can't occur
# source the prep function
#source("./part2_local_landscape_predictors_of_occupancy/run_model/get_species_ranges.R")
#range_data <- get_species_ranges(city_names)

species_info <- my_data$species_info
species_region_cluster_id <- my_data$species_region_cluster_id

cluster <-c( "southeast", # atlanta
             "northeast", # boston
             "southeast", # charlotte
             "midwest", # chicago
             "texas", # dallas
             "northeast", # dc
             "texas", # denton
             "texas", # houston
             "california", # LA
             "midwest", # minneapolis
             "northeast", # nyc
             "northeast", # philadelphia
             "southeast", # raleigh
             "california", # sd
             "san_francisco") # sf
x_name <- "city"
y_name <- "cluster"

city_cluster <- data.frame(city_names,cluster)
names(city_cluster) <- c(x_name,y_name)

site_data_all <- left_join(site_data_all, city_cluster)

# get whether the city is in range (may not be for all cities in the same cluster)
source("./part2_local_landscape_predictors_of_occupancy/run_model/get_species_ranges.R")
range_data <- get_species_ranges(city_names)

## --------------------------------------------------
# Now grab all the sites from a city, get the expected psi rate
# then simulate occurrence outcomes based on expected rate
# then summarize biodiversity metrics from the occurrence
# then repeat many times, resampling the params for expected psi 
# from the model's joint posterior

# get indices for species and city random effects distributions for particular region
psi_0 <- which( colnames(estimates)=="psi_0" )
first_psi_city <- which( colnames(estimates)=="psi_city[1]" )
first_psi_species <- which( colnames(estimates)=="psi_species[1]" )
first_psi_wingspan <- which( colnames(estimates)=="psi_wingspan[1]" )
first_psi_parksize <- which( colnames(estimates)=="psi_park_size[1]" )
first_psi_tree_cover <- which( colnames(estimates)=="psi_tree_cover[1]" )
first_psi_plant_diversity <- which( colnames(estimates)=="psi_plant_diversity[1]" )
first_psi_isolation <- which( colnames(estimates)=="psi_landscape_isolation[1]" )
first_psi_landscape_grassherb <- which( colnames(estimates)=="psi_landscape_grassherb[1]" )
first_psi_landscape_woody <- which( colnames(estimates)=="psi_landscape_woody[1]" )

# some random samples from the posterior
n_draws = 50 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors
random_draws_from_posterior = sample.int(nrow(estimates), n_draws) # use if not using the full posterior

mean_richness <- array(dim = c(n_cities, n_draws))
beta_diversity <- array(dim = c(n_cities, n_draws))
gamma_diversity <- array(dim = c(n_cities, n_draws))

for(city_number in 1:n_cities){
  for(draw in 1:n_draws){
    
    # select a random draw from the joint posterior
    rand <- random_draws_from_posterior[draw]
    
    # filter the site covariate data to the specific city
    temp <- filter(site_data_all, city == city_names[city_number])
    
    # get all of the species that could potentially occur at each site
    temp <- left_join(temp, species_region_cluster_id)
    temp <- left_join(temp, select(species_info, species, aveWingspan_scaled))
    
    # get the correct range data
    temp_ranges <- filter(range_data, city == city_names[city_number])
    # and filter to the species that are actually in range of the city
    temp <- temp %>% filter(species %in% temp_ranges$species)
    
    # construct expected and realized occurrence arrays that are of length n_sites*n_species
    psi <- array(dim = c(nrow(temp)))
    occurrence <- array(dim = c(nrow(temp)))
    
    #for(r in 1:nrow(temp)){
    psi <- 
      as.numeric(estimates[rand, first_psi_city - 1 + city_number]) +
      as.numeric(estimates[rand, first_psi_species - 1 + temp$species_cluster_integer_vector]) +
      estimates[rand, first_psi_wingspan - 1 + city_number] * temp$aveWingspan_scaled +
      estimates[rand, first_psi_parksize - 1 + city_number] * temp$site_size_pred +
      estimates[rand, first_psi_tree_cover - 1 + city_number] * temp$site_tree_cover_pred +
      estimates[rand, first_psi_plant_diversity - 1 + city_number] * temp$site_plant_diversity_pred +
      estimates[rand, first_psi_isolation - 1 + city_number] * temp$site_isolation_pred +
      estimates[rand, first_psi_landscape_grassherb - 1 + city_number] * temp$site_landscape_grassherb_pred +
      estimates[rand, first_psi_landscape_woody - 1 + city_number] * temp$site_landscape_woody_pred 
    
    occurrence <- rbinom(length(psi), 1, ilogit(psi)) 
    #} # get psi 
  
  # mean richness of dims n_cities, n_draws  
  mean_richness[city_number, draw] <- cbind(temp, occurrence) %>%
    group_by(new_id) %>%
    mutate(site_richness = sum(occurrence)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(city_mean_alpha_richness = mean(site_richness)) %>%
    slice(1) %>%
    pull(city_mean_alpha_richness)
  
  # mean dissimilarity among parks
  temp_wide <- cbind(temp, occurrence) %>%
    select(new_id, species, occurrence) %>%
    pivot_wider(names_from = species, values_from = occurrence) %>%
    as.matrix(.)
  temp_wide <- temp_wide[,-1] 
  # higher values indicate low similarity
  beta_diversity[city_number, draw] <- mean(vegan::vegdist(temp_wide, method = "jaccard", binary = TRUE))

  # overall diversity supported in greenspaces
  gamma_diversity[city_number, draw] <- nrow(cbind(temp, occurrence) %>%
    filter(occurrence == 1) %>%
    group_by(species) %>%
    slice(1))
    
  } # for each draw from the joint posterior
} # for each city


#
hist(mean_richness[,1])
hist(beta_diversity[,1])
hist(gamma_diversity[,1])

simmed_diversity <- list(mean_richness, beta_diversity, gamma_diversity)
saveRDS(simmed_diversity, "./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

#-------------------------------------------------------------------------------
# get the city covariate data
city_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")

city_data <- city_data %>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:n_cities))

city_data <- city_data %>%
  mutate(log_avg_park_size = log(average_park_size_sqm),
         log_park_size_scaled = center_scale(log_avg_park_size))


#-------------------------------------------------------------------------------
# summarize uncertainty and plot relationships

#-------------------------------------------------------------------------------
# mean species richness

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

mean_richness_quantiles_df <- cbind(city_data, mean_richness_quantiles_df)

a <- ggplot(mean_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Mean log(Park Size)") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# beta diversity (mean jaccard dissimilarity)

#  calculate Means and CI's for the diversity metrics for each city
mean_beta_quantiles <- apply(beta_diversity, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             na.rm=TRUE)
mean_beta_quantiles_df <- as.data.frame(t(mean_beta_quantiles))
colnames(mean_beta_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

mean_beta_quantiles_df <- cbind(city_data, mean_beta_quantiles_df)

b <- ggplot(mean_beta_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Dissimilarity\n(Jaccard Index)") +
  xlab("Mean log(Park Size)") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# gamma (city-wide) species richness

#  calculate Means and CI's for the diversity metrics for each city
gamma_richness_quantiles <- apply(gamma_diversity, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
gamma_richness_quantiles_df <- as.data.frame(t(gamma_richness_quantiles))
colnames(gamma_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

gamma_richness_quantiles_df <- cbind(city_data, gamma_richness_quantiles_df)

c <- ggplot(gamma_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("Mean log(Park Size)") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(a, b, c, ncol = 3)

#-------------------------------------------------------------------------------
# do it again with relative species richness

# alpha and gamma diversity could be relative to the regional species pools
size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools_BAMONA.csv")
size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools.csv")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# mean species richness

# standardize mean richness
mean_richness_relative <- mean_richness
for(i in 1:n_cities){
  mean_richness_relative[i,] = mean_richness_relative[i,] / size_of_regional_species_pools[i,1]
}

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness_relative, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

mean_richness_quantiles_df <- cbind(city_data, mean_richness_quantiles_df)


d <- ggplot(mean_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness (Relative)") +
  xlab("Mean log(Park Size)") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# gamma (city-wide) species richness

# standardize gamma richness
gamma_diversity_relative <- gamma_diversity
for(i in 1:n_cities){
  gamma_diversity_relative[i,] = gamma_diversity_relative[i,] / size_of_regional_species_pools[i,1]
}

#  calculate Means and CI's for the diversity metrics for each city
gamma_richness_quantiles <- apply(gamma_diversity_relative, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
gamma_richness_quantiles_df <- as.data.frame(t(gamma_richness_quantiles))
colnames(gamma_richness_quantiles_df) <- c("lower90", "lower50",
                                           "mean", "upper50", "upper90")

gamma_richness_quantiles_df <- cbind(city_data, gamma_richness_quantiles_df)

f <- ggplot(gamma_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Total Number of Species\nOccurring in City Parks (Relative)") +
  xlab("Mean log(Park Size)") +
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(d, b, f, ncol = 3)








