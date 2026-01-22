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

# scale the size data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_size <- sd(site_data$log_total_green_space_area) # sd of size of parks used to create scale
mean_size <- mean(site_data$log_total_green_space_area) # mean size of parks used to create scale
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_size_pred <- (site_data_all$log_total_green_space_area - mean_size) / sd_size 

# scale the isolation data (have to use the same mean and sd that were used to scale the data fed to the model)
sd_isolation <- sd(log(site_data$isolation)) # sd of size of parks used to create scale
mean_isolation <- mean(log(site_data$isolation)) # mean size of parks used to create scale
# how does the park (which may or may not have been modelled stck up to this scale)
site_data_all$site_isolation_pred <- (log(site_data_all$isolation) - mean_isolation) / sd_isolation 


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
first_psi_ <- which( colnames(estimates)=="psi_tree_cover[1]" )
first_psi_parksize <- which( colnames(estimates)=="psi_[1]" )
first_psi_isolation <- which( colnames(estimates)=="psi_landscape_isolation[1]" )

# some random samples from the posterior
n_draws = 1 # small number for testing bc it does take a few minutes to simulate results
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
    
    # construct expected and realized occurrence arrays that are of length n_sites*n_species
    psi <- array(dim = c(nrow(temp), n_draws))
    occurrence <- array(dim = c(nrow(temp), n_draws))
    
    for(r in 1:nrow(temp)){
    psi[r,draw] <- 
      estimates[rand, first_psi_city - 1 + city_number] +
      estimates[rand, first_psi_species - 1 + temp$species_cluster_integer_vector[r]] +
      estimates[rand, first_psi_wingspan - 1 + city_number] * temp$aveWingspan_scaled[r]+
      estimates[rand, first_psi_parksize - 1 + city_number] * temp$site_size_pred[r] +
      estimates[rand, first_psi_isolation - 1 + city_number] * temp$site_isolation_pred[r] 
    
    occurrence[r,draw] <- rbinom(1, 1, ilogit(psi[r,draw])) 
    } # get psi 
  
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



#-------------------------------------------------------------------------------
# get some prediction data for each city

## get all the site data from all cities
my_data <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_all.rds"))

site_data <- my_data$site_data

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_size_pred_data_list <- vector(mode='list', length=n_cities)
park_size_original_data_list <- vector(mode='list', length=n_cities)
park_isolation_pred_data_list <- vector(mode='list', length=n_cities)
park_isolation_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

# remember that the predictor values were scaled within each city, so we have to
# come up with some pred data specific to the real values within each city (not across all cities)
for(city_number in 1:n_cities){
  
  # filter the site covariate data to the specific city
  temp <- filter(site_data, city == city_names[city_number]) %>%
    arrange(., log_total_green_space_area)
  
  # get the scaled park size data
  park_size_pred_data <- temp$log_total_green_space_area_scaled_across_all_cities
  # get the real park size data
  original_scale_park_size_data <- temp$log_total_green_space_area
  
  # get the scaled isolation data
  park_isolation_pred_data <- temp$log_isolation_scaled_across_all_cities
  # get the real isolation data
  original_scale_park_isolation_data <- log(temp$isolation)
  
  # store the city-specific scaled and real values
  park_size_pred_data_list[[city_number]] <- park_size_pred_data
  park_size_original_data_list[[city_number]]  <- original_scale_park_size_data
  park_isolation_pred_data_list[[city_number]]  <- park_isolation_pred_data
  park_isolation_original_data_list[[city_number]]  <- original_scale_park_isolation_data
  
  # and figure out how many sites in the city
  pred_length[[city_number]]  <- length(park_size_pred_data)

}

# set the maximum pred length to the largest number of sites within a city
max_pred_length = max(pred_length)

# get species data
species_data <- my_data$species_info
n_species <- nrow(species_data)

R <- my_data$R
V <- my_data$V
confirmed_occurrence <- my_data$confirmed_occurrence
multicity_site_id_vector <- my_data$multicity_site_id_vector
species_integer_vector <- my_data$species_integer_vector
city_id_vector <- my_data$city_id_vector
site_survey_year_vector <- my_data$site_survey_year_vector
species_cluster_integer_vector <- my_data$species_cluster_integer_vector
region_cluster_integer_vector <- my_data$region_cluster_integer_vector

# get range data so we properly omit species that can't occur
# source the prep function
source("./part2_local_landscape_predictors_of_occupancy/run_model/get_species_ranges.R")
range_data <- get_species_ranges(city_names)

city_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")
city_data <- city_data %>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:length(n_cities)))

range_data <- left_join(range_data, city_data)
range_data <- left_join(as.data.frame(city_id_vector), range_data,
                        by = c("city_id_vector" = "city_factor"))



## --------------------------------------------------

# some random samples from the posterior
random_draws_from_posterior = sample.int(nrow(estimates), n_draws) # use if not using the full posterior

# take random draws for psi and predict occurrence, then sum across n species
richness <- array(data = NA, dim=c(n_cities, max_pred_length, n_draws))

# get indices for species and city random effects distributions for particular region
psi_0 <- which( colnames(estimates)=="psi_0" )
first_psi_city <- which( colnames(estimates)=="psi_city[1]" )
first_psi_species <- which( colnames(estimates)=="psi_species[1]" )
first_psi_wingspan <- which( colnames(estimates)=="psi_wingspan[1]" )
first_psi_parksize <- which( colnames(estimates)=="psi_park_size[1]" )
first_psi_isolation <- which( colnames(estimates)=="psi_landscape_isolation[1]" )

psi <- matrix(nrow=R, ncol=n_draws)
occurrence <- matrix(nrow=R, ncol=n_draws)

# work in the ranges thing in advance here
# could multiply by ranges[0,1]
for(r in 1:R){
  for(draw in 1:n_draws){
   
    rand <- random_draws_from_posterior[draw]
    
    psi[r, draw] = 
      as.numeric(estimates[rand, first_psi_city - 1 + city_id_vector[r]]) +
      as.numeric(estimates[rand, first_psi_species - 1 + species_cluster_integer_vector[r]]) +
      as.numeric(estimates[rand, first_psi_wingspan - 1 + city_id_vector[r]] * species_data$aveWingspan_scaled[species_integer_vector[r]]) + 
      as.numeric(estimates[rand, first_psi_parksize - 1 + city_id_vector[r]] * site_data$log_total_green_space_area_scaled_across_all_cities[multicity_site_id_vector[r]]) +
      as.numeric(estimates[rand, first_psi_isolation - 1 + city_id_vector[r]] * site_data$isolation_scaled_across_all_cities[multicity_site_id_vector[r]]) 
    
  }
  
  if(confirmed_occurrence[r] > 0){ 
    
    # if a species occurrence was observed, then the species occurs
    occurrence[r,draw] <- 1
    
  } else{
    
    # else the species occurrence will be simulated using model predictions
    # in all years after the first year we use a temporal autocorrelation component
    occurrence[r,draw] <- rbinom(1, 1, prob = ilogit(psi[r,draw]))
    
  } 
}

##


  
  for(city_number in 1:n_cities){
    
    # get the pred data for species/sites from the particular city
    city_name <-  city_names[city_number]
    # get detections, if a species was observed we know it was present
    # otherwise we predict based on our model based estimates
    detections <- my_data$V_detections
    
    pred_length_city <- pred_length[start_city_index[region_number] + city_number]
    
    # get the predictor data for the particular city
    park_size_pred_data <- park_size_pred_data_list[[start_city_index[region_number] + city_number]]
    park_isolation_pred_data <- park_isolation_pred_data_list[[start_city_index[region_number] + city_number]]
    
    # construct some arrays to hold the data
    psi_expected <- array(data = NA, dim=c(n_species, pred_length_city))
    occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length_city, n_years, n_draws))

    # prep to drop out species that were not modelled in city (based on inferred classified park range)
    min_site_index <- min(filter(site_data, city == city_name) %>% select(multicity_site_id))
    max_site_index <- max(filter(site_data, city == city_name) %>% select(multicity_site_id))
    ranges_city <- range_data[,min_site_index:max_site_index]
    
    for(draw in 1:n_draws){
      
      rand <- random_draws_from_posterior[draw]
      
      # expected occurrence in year 1
      for(i in 1:n_species){
        for(j in 1:pred_length_city){
            
            psi1_expected[i,j] =
              ilogit(
                # global initial occurrence intercept
                estimates[rand,psi1_0] + 
                  # plus a city specific random intercept effect   
                  estimates[rand,(first_psi1_city+(city_number-1))] + 
                  # a species specific intercept effect (the number here should be first column)
                  estimates[rand,(first_psi1_species+(i-1))] +
                  # effect of wingspan * wingspan of species i + 
                  estimates[rand,(first_psi1_wingspan+(city_number-1))] * species_data$aveWingspan_scaled[i] + 
                  # effect of parksize * parksize of site j + 
                  estimates[rand,(first_psi1_parksize+(city_number-1))] * park_size_pred_data[j] +
                  # effect of isolation * isolation of site j + 
                  estimates[rand,(first_psi1_isolation+(city_number-1))] * park_isolation_pred_data[j]
              )
            
        }
      }
      
      # simmed occurrence for each species
      for(i in 1:n_species){
        for(j in 1:pred_length_city){
          for(k in 1:n_years){
            
            if(sum(detections[i,j,k,1:12]) > 0){ 
              
              # if a species occurrence was observed, then the species occurs
              occurrence_simmed[i,j,1,rand] <- 1
              
            } else{
              
              # else the species occurrence will be simulated using model predictions
              # in all years after the first year we use a temporal autocorrelation component
                occurrence_simmed[i,j,1,rand] <- rbinom(1, 1, prob = psi1_expected[i,j])
              
            } # end if/else species occurrence was observed in real life
            
            occurrence_simmed[i,j,k,rand] <- ranges_city[i,j] * occurrence_simmed[i,j,k,rand]
            
          } # end for(k)
        } # end for(j)
      } # end for(i)
      
      for(j in 1:pred_length_city){
        for(k in 1:n_years){
          richness[start_city_index[region_number] + city_number,j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k,rand])
        }
      }
      
    } # end for random draw
    
  } # end for city

## --------------------------------------------------
# summarize the results and plot species richness against real park size or park isolation

# collapse across years (average richness 
# per rand draw from the posterior [array dimension 4] 
# per site [dim 2], 
# per city [dim 1],
# across all years [dim 3])
richness2 <- richness
#richness <- apply(richness,c(1,2,4),mean) # average across years
richness <- apply(richness,c(1,2,4),mean,na.rm=TRUE) # average across years

#test <- richness[1:2, 1:30, 1, 1:50]
#richness <- test

# make an empty df
df <- data.frame()

# the data is stored by region so we will have to access by region
for(region_number in 1:n_regions){
  
  region <- regions[region_number]
  
  cities_region <- eval(parse(text=paste0("city_names_", tolower(region))))
  n_cities_region <- length(cities_region)
  
  for(city_number in 1:n_cities_region){
    
    city_name <- cities_region[city_number]
    
    # prep sites to model in each city
    my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))
    site_data <- my_data$site_data
    min_site_index <- min(filter(site_data, city == city_name) %>% select(multicity_site_id))
    max_site_index <- max(filter(site_data, city == city_name) %>% select(multicity_site_id))
    n_sites_city <- (max_site_index - min_site_index + 1)
    
    mean = vector(length=n_sites_city)
    lower_50 = vector(length=n_sites_city)
    upper_50 = vector(length=n_sites_city)
    lower_95 = vector(length=n_sites_city)
    upper_95 = vector(length=n_sites_city)
    
    for(j in 1:n_sites_city){
      
      quants = as.vector(quantile(
        richness[start_city_index[region_number] + city_number,j,], 
        probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm=TRUE))
      
      mean[j] = quants[3]
      lower_50[j] = quants[2]
      upper_50[j] = quants[4]
      lower_95[j] = quants[1]
      upper_95[j] = quants[5]
      
    } # end for across sites to quantify quantiles
    
    pred_data <- park_size_pred_data_list[[start_city_index[region_number] + city_number]]
    original_scale_park_size_data <- park_size_original_data_list[[start_city_index[region_number] + city_number]]
    
    temp <- as.data.frame(cbind(pred_data, #original_scale_data,
                                mean,
                                lower_50, 
                                upper_50,
                                lower_95, 
                                upper_95,
                                original_scale_park_size_data
    )) %>%
      mutate(city = city_name) %>%
      left_join(., size_of_regional_species_pools, by = "city") %>%
      mutate(rel_mean = mean / size_of_regional_pool,
             rel_lower_50 = lower_50 / size_of_regional_pool,
             rel_upper_50 = upper_50 / size_of_regional_pool,
             rel_lower_95 = lower_95 / size_of_regional_pool,
             rel_upper_95 = upper_95 / size_of_regional_pool)
    
    df <- rbind(df, temp)
    
  }  # end for across cities within region
  
} # end for across regions


## --------------------------------------------------
## Draw species richness plot

df <- df %>%
  mutate(city = as.factor(city)) %>%
  mutate(city=fct_relevel(city,c("Boston", 
                                 "DC",
                                 "NYC", 
                                 "Philadelphia",
                                 "Atlanta",
                                 "Charlotte",
                                 "Dallas",
                                 "Denton",
                                 "Houston",
                                 "Raleigh"))) 

# create plot object with loess regression lines
# this will create a smooth line between the 
pre_p <- ggplot(df) + 
  stat_smooth(aes(x = original_scale_park_size_data, y = rel_mean, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = rel_lower_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = rel_upper_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = rel_lower_50, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = rel_upper_50, group=city), method = "loess", se = FALSE)
pre_p

# build plot object for rendering 
p <- ggplot_build(pre_p)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(original_scale_park_size_data = p$data[[1]]$x,
                  rel_mean = p$data[[1]]$y,
                  rel_lower_95 = p$data[[2]]$y,
                  rel_upper_95 = p$data[[3]]$y,
                  rel_lower_50 = p$data[[4]]$y,
                  rel_upper_50 = p$data[[5]]$y) 

df2 <- cbind(df2, rep(city_names, each = nrow(df2)/n_cities))

colnames(df2)[7] <- "City"

df2 <- df2 %>%
  mutate(City = as.factor(City))

# use the loess data to add the 'ribbon' to plot 
(p  <- ggplot(data = df2, aes(original_scale_park_size_data, group=City)) +
    geom_line(aes(y=rel_mean, colour=City), size = 2) +
    geom_ribbon(aes(ymin = rel_lower_95, ymax = rel_upper_95, fill=City), alpha = 0.2) +
    #geom_ribbon(aes(ymin = rel_lower_50, ymax = rel_upper_50, fill=City), alpha = 0.5) +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("log(Park Size in m^2)") +
    ylab("Species Richness /\nSize of Regional Species Pool") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)
    )
  
) 


# now plot isolation

# build plot object for rendering 
q <- ggplot_build(pre_q)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(original_scale_isolation_data = q$data[[1]]$x,
                  rel_mean = q$data[[1]]$y,
                  rel_lower_95 = q$data[[2]]$y,
                  rel_upper_95 = q$data[[3]]$y,
                  rel_lower_50 = q$data[[4]]$y,
                  rel_upper_50 = q$data[[5]]$y) 

df2 <- cbind(df2, rep(city_names, each = nrow(df2)/n_cities))

colnames(df2)[7] <- "City"

# use the loess data to add the 'ribbon' to plot 
(q  <- ggplot(data = df2, aes(original_scale_isolation_data, group=City)) +
    geom_line(aes(y=rel_mean, colour=City), size = 2) +
    geom_ribbon(aes(ymin = rel_lower_95, ymax = rel_upper_95, fill=City), alpha = 0.2) +
    geom_ribbon(aes(ymin = rel_lower_50, ymax = rel_upper_50, fill=City), alpha = 0.5) +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("Isolation\n(Avg. Distance to Parks Within 2km)") +
    ylab("Species Richness /\nSize of Regional Species Pool") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)
    )
  
) 