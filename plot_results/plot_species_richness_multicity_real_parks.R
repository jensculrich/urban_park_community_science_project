# plot effects of site predictors on local species richness

library(tidyverse)
library(rstan)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools.csv")

# list of city names
city_names <- c(
  "Boston", # 1
  "Dallas", # 2
  "Houston", # 3
  "LA", # 4
  "NYC", # 5
  #"Riverside", # 6
  "SF" # 7
)

n_cities <- length(city_names)

# create empty list of length n_cities
# each element of the list holds the posterior distributions for 
# all parameters for each individual city 
mega_list <- vector(mode='list', length=n_cities)

# read and summarize data across loop
for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- cbind(as.data.frame(
      readRDS(paste0(
        "./model_outputs/stan_out_", city, "_2km_connectivity_family_50buffers_simple.rds"))
  ), city)
  
  mega_list[[city_number]] <- temp

}

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_size_pred_data_list <- vector(mode='list', length=n_cities)
park_size_original_data_list <- vector(mode='list', length=n_cities)
park_connectivity_pred_data_list <- vector(mode='list', length=n_cities)
park_connectivity_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- readRDS( paste0("./run_model/prepped_data/prepped_data_", city, ".rds"))
   
  # get park size data
  park_size_pred_data <- temp$site_data$log_total_green_space_area_scaled
    
  mean_park_size <- mean(temp$site_data$log_total_green_space_area)
  sd_park_size <- sd(temp$site_data$log_total_green_space_area)
  # now do some algebra to get the scaled data back onto a real life m^2 scale
  original_scale_park_size_data <- (park_size_pred_data * sd_park_size) + mean_park_size
  
  park_size_pred_data_list[[city_number]] <- park_size_pred_data
  park_size_original_data_list[[city_number]] <- original_scale_park_size_data
  
  # get park connectivity data
  park_connectivity_pred_data <- temp$site_data$area_weighted_avg_dist_2000m_scaled
  
  mean_park_connectivity <- mean(temp$site_data$area_weighted_avg_dist_2000m)
  sd_park_connectivity <- sd(temp$site_data$area_weighted_avg_dist_2000m)
  # now do some algebra to get the scaled data back onto a real life m^2 scale
  original_scale_park_connectivity_data <- (park_connectivity_pred_data * sd_park_connectivity) + mean_park_connectivity
  
  park_connectivity_pred_data_list[[city_number]] <- park_connectivity_pred_data
  park_connectivity_original_data_list[[city_number]] <- original_scale_park_connectivity_data
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_size_pred_data)
  
}

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## plot species richness by park size

## --------------------------------------------------
## get prediction range

n_years = 5 # number of years of the study
n_years_minus1 = n_years - 1 # number of interannual transitions

n_draws = 20 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors

max_pred_length = max(pred_length)

## --------------------------------------------------

# some random samples from the posterior
random_draws_from_posterior = sample.int(n=n_draws) # use if not using the full posterior

# take random draws for psi and predict occurrence, then sum across n species
richness <- array(data = NA, dim=c(n_cities, max_pred_length, n_years, n_draws))
jaccard_site <- array(data = NA, dim=c(n_cities, max_pred_length, max_pred_length, n_years, n_draws))

for(city_number in 1:n_cities){
  
  # get the posterior distributions for a particular city
  city_estimates <- mega_list[[city_number]]
  
  # get indices for species random effects distributions for particular city
  first_psi1 <- which( colnames(city_estimates)=="psi1_species[1]" )
  first_gamma <- which( colnames(city_estimates)=="gamma_species[1]" )
  first_phi <- which( colnames(city_estimates)=="phi_species[1]" )
  
  # get the pred data for species/sites from the particular city
  city_name <- city_names[city_number]
  my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".rds"))
  species_data <- my_data$species_info
  n_species <- nrow(species_data)
  detections <- my_data$V_detections
  pred_length_city <- pred_length[city_number]
  
  # get the predictor data for the particular city
  park_size_pred_data <- park_size_pred_data_list[[city_number]]
  park_connectivity_pred_data <- park_connectivity_pred_data_list[[city_number]]

  # construct some arrays to hold the data
  psi1_expected <- array(data = NA, dim=c(n_species, pred_length_city))
  gamma_expected <- array(data = NA, dim=c(n_species, pred_length_city, n_years_minus1))
  phi_expected <- array(data = NA, dim=c(n_species, pred_length_city, n_years_minus1))
  occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length_city, n_years, n_draws))
  #psi <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
  
  for(draw in 1:n_draws){
  
    rand <- random_draws_from_posterior[draw]
    
    # expected occurrence in year 1
    for(i in 1:n_species){
      for(j in 1:pred_length_city){
        for(k in 2:n_years){
          
          psi1_expected[i,j] =
            ilogit(
              # YEAR 1 is the global intercept
              city_estimates[rand,1] + 
                # a species specific intercept effect (the number here should be first column)
                city_estimates[rand,(first_psi1+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,3] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,4] * park_size_pred_data[j] +
                # effect of connectivity * connectivity of site j + 
                city_estimates[rand,5] * park_connectivity_pred_data[j]
            )
          
          gamma_expected[i,j,k-1] = # gamma[,,k-1] yields gamma for the transition between year 1 and 2
            ilogit(#gamma0 +
              city_estimates[rand,6] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_gamma+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,8] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,9] * park_size_pred_data[j]  +
                # effect of connectivity * connectivity of site j + 
                city_estimates[rand,10] * park_connectivity_pred_data[j]
            )
          
          phi_expected[i,j,k-1] = # phi[,,k-1] yields phi for the transition between year 1 and 2
            ilogit(#phi0 +
              city_estimates[rand,11] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_phi+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,13] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,14] * park_size_pred_data[j]  +
                # effect of connectivity * connectivity of site j + 
                city_estimates[rand,15] * park_connectivity_pred_data[j]
            )
          
        } 
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
            if(k < 2){ # in year 1
              occurrence_simmed[i,j,1,rand] <- rbinom(1, 1, prob = psi1_expected[i,j])
            } else{ # contingent on previous year
              occurrence_simmed[i,j,k,rand] = occurrence_simmed[i,j,k-1,rand] * phi_expected[i,j,k-1] + 
                (1 - occurrence_simmed[i,j,k-1,rand]) * gamma_expected[i,j,k-1]
            } # end if/else simulate species occurrence
            
          } # end if/else species occurrence was observed in real life
          
        } # end for(k)
      } # end for(j)
    } # end for(i)
    
    for(j in 1:pred_length_city){
      for(k in 1:n_years){
        richness[city_number,j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k,rand])
      }
    }
    
    ## --------------------------------------------------
    ## calculate beta diversity of simmed communities
    
    # jaccard index needs a reference level
    #ref_site <- 1
    # now compute dissimilarity in the occurrence matrix in site i versus site 1
    #for(j in 1:pred_length_city){ # jaccard index for sites (in terms of shared species)
    #  for(k in 1:n_years){
    #    jaccard_site[city_number,j,k,draw] <- sum(occurrence_simmed[,ref_site,k,rand]*occurrence_simmed[,j,k,rand]) /
    #      (sum(occurrence_simmed[,ref_site,k,rand]) +
    #         sum(occurrence_simmed[,ref_site,k,rand]) - sum(occurrence_simmed[,ref_site,k,rand]*occurrence_simmed[,j,k,rand]))
    #  }
    #}
    
    # now compute dissimilarity in the occurrence matrix in site i versus site 1
    for(j in 1:pred_length_city){ # jaccard index for sites (in terms of shared species)
      for(j2 in 1:pred_length_city){ # jaccard index for sites (in terms of shared species)
        for(k in 1:n_years){
          jaccard_site[city_number,j,j2,k,draw] <- sum(occurrence_simmed[,j2,k,rand]*occurrence_simmed[,j,k,rand]) /
            (sum(occurrence_simmed[,j2,k,rand]) +
               sum(occurrence_simmed[,j2,k,rand]) - sum(occurrence_simmed[,j2,k,rand]*occurrence_simmed[,j,k,rand]))
        }
      }
    }
    
  } # end for random draw

} # end for city
  
## --------------------------------------------------
# summarize the results (species richness)

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
#richness <- apply(richness,c(1,2,4),mean) # average across years
richness <- apply(richness,c(1,2,4), mean, na.rm=TRUE) # average across years

#imaginary_city_covariate <- rnorm(n_cities, 0, 1)
#imaginary_city_covariate <- as.data.frame(cbind(city_names, imaginary_city_covariate))

# extract city wide covariates
city_wide_connectivity <- vector(length = n_cities)

for(city_number in 1:n_cities){

  city_name <- city_names[city_number]
  
  city_wide_connectivity[city_number] <- read.csv(paste0(
    "./data/city_shapefiles/", city_name, 
    "/", city_name, 
    "_landscape_metrics.csv"))[1,1] 

}

city_wide_connectivity_scaled <- center_scale(city_wide_connectivity)

city_wide_connectivity_df <- as.data.frame(cbind(
  city_names, city_wide_connectivity, city_wide_connectivity_scaled))

# initialize df with correct sites per city
df <- richness[,,1]

colnames(df) <- seq(1:max_pred_length)

df <- as.data.frame(df)

df <- df %>%
  cbind(., city_names) %>%
  pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(draw)) %>%
  filter(!is.na(.[,3])) %>%
  select(city_names, site)

# now add the richness per site per random draw of param estimates from the posterior distr
for(draw in 1:n_draws){
  
  temp <- richness[,,draw]
  
  colnames(temp) <- seq(1:max_pred_length)
  
  temp <- as.data.frame(temp)
  
  temp <- temp %>%
    cbind(., city_names) %>%
    pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(draw)) %>%
    filter(!is.na(.[,3]))
  
  df <- cbind(df, temp[,3])
  
}

# calculate city mean species richness
mean = matrix(nrow=n_cities)
lower_50 = matrix(nrow=n_cities)
upper_50 = matrix(nrow=n_cities)
lower_95 = matrix(nrow=n_cities)
upper_95 = matrix(nrow=n_cities)

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  
  temp <- filter(df, city_names == city_name) %>%
    select(-city_names, -site)
    
    quants = as.vector(quantile(as.matrix(temp), probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[city_number] = quants[3]
    lower_50[city_number] = quants[2]
    upper_50[city_number] = quants[4]
    lower_95[city_number] = quants[1]
    upper_95[city_number] = quants[5]
}

estimates <- as.data.frame(cbind(city_names,
                                 mean,
                                 lower_50, 
                                 upper_50,
                                 lower_95, 
                                 upper_95
                                 
  )) %>%
  rename("mean" = "V2",
         "lower_50" = "V3",
         "upper_50" = "V4",
         "lower_95" = "V5",
         "upper_95" = "V6") %>%
  mutate(city = city_names) %>%
  left_join(., size_of_regional_species_pools, by = "city") %>%
  mutate(mean =as.numeric(mean),
         lower_50 =as.numeric(lower_50),
         upper_50 =as.numeric(upper_50),
         lower_95 =as.numeric(lower_95),
         upper_95 =as.numeric(upper_95)) %>%
  mutate(rel_mean = mean / size_of_regional_pool,
         rel_lower_50 = lower_50 / size_of_regional_pool,
         rel_upper_50 = upper_50 / size_of_regional_pool,
         rel_lower_95 = lower_95 / size_of_regional_pool,
         rel_upper_95 = upper_95 / size_of_regional_pool) 

estimates <- estimates %>%
  left_join(., city_wide_connectivity_df) %>%
  mutate(city_wide_connectivity = as.numeric(city_wide_connectivity),
         city_wide_connectivity_scaled = as.numeric(city_wide_connectivity_scaled))

## --------------------------------------------------
## Draw species richness plot

# plot means and variation across city_wide_connectivity_scaled
(p <- ggplot(estimates, aes(x=city_wide_connectivity_scaled)) +
   theme_bw() +
   scale_y_continuous(str_wrap("Mean Alpha Diversity /\nSize of Regional Species Pool", width = 30),
                      limits = c(0, 1), breaks = c(0, 0.5, 1)) +
   scale_x_continuous(str_wrap("City-Wide Connectivity (Scaled)", width = 30),
                      limits = c(-2, 2.5), breaks = c(-2, -1, 0, 1, 2)) +
   ggtitle("") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

p <- p +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_95, ymax=rel_upper_95, colour=city, group=city),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_50, ymax=rel_upper_50, colour=city, group=city),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=city_wide_connectivity_scaled, y=mean, colour=city, group=city), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p

## --------------------------------------------------
## Quantify association between species richness and city-wide covariate


df_w_city_covs <- left_join(df, city_wide_connectivity_df)

all_posterior_samples <- matrix(nrow=100*n_draws, ncol=2)

for(draw in 1:n_draws){
  # make a for loop across draws
  means_by_city <- df_w_city_covs %>%
    select(city_names, city_wide_connectivity_scaled, as.character(draw))
  
  colnames(means_by_city)[3] <- "richness"
  
  means_by_city <- means_by_city %>%
    group_by(city_names) %>%
    mutate(richness = as.numeric(richness),
           mean_alpha_diversity = mean(richness)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(city = city_names) %>%
    left_join(., size_of_regional_species_pools, by = "city") %>%
    mutate(rel_mean_alpha_diversity = mean_alpha_diversity / size_of_regional_pool) %>%
    mutate(city_wide_connectivity_scaled = as.numeric(city_wide_connectivity_scaled))
  
  
  stan_fit <- rstanarm::stan_glm(rel_mean_alpha_diversity ~ city_wide_connectivity_scaled, data = means_by_city)
  
  posterior_samples <- as.matrix(sample_n(as.data.frame(stan_fit)[,1:2], size = 100))

  all_posterior_samples[(1+100*(draw - 1)):(100+100*(draw - 1)),] <- posterior_samples
}

all_posterior_samples_df <- as.data.frame(all_posterior_samples) %>%
  rename("intercept" = "V1",
         "slope" = "V2")

hist(all_posterior_samples_df$slope)

# 
## --------------------------------------------------
## Draw species richness plot

# plot means and variation across cities
(p <- ggplot(estimates, aes(x=city_wide_connectivity_scaled)) +
   theme_bw() +
   scale_y_continuous(str_wrap("Mean Alpha Diversity /\nSize of Regional Species Pool", width = 30),
                      limits = c(0, 1), breaks = c(0, 0.5, 1)) +
   scale_x_continuous(str_wrap("City-Wide Connectivity (Scaled)", width = 30),
                      limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
   ggtitle("") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

p <- p +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_95, ymax=rel_upper_95, colour=city, group=city),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_50, ymax=rel_upper_50, colour=city, group=city),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=city_wide_connectivity_scaled, y=mean, colour=city, group=city), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p

random_lines <- sample.int(nrow(all_posterior_samples), 100)

for(i in 1:100){
  p <- p + geom_abline(intercept = all_posterior_samples[random_lines[i],1], 
                  slope = all_posterior_samples[random_lines[i],2], 
                  alpha = 0.3, colour = "grey")
}

p


## --------------------------------------------------
# summarize the results (species dissimilarity)

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
#richness <- apply(richness,c(1,2,4),mean) # average across years
test <- apply(jaccard,c(1,2,4), mean, na.rm=TRUE) # average across years

#imaginary_city_covariate <- rnorm(n_cities, 0, 1)
#imaginary_city_covariate <- as.data.frame(cbind(city_names, imaginary_city_covariate))

# extract city wide covariates
city_wide_connectivity <- vector(length = n_cities)

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  
  city_wide_connectivity[city_number] <- read.csv(paste0(
    "./data/city_shapefiles/", city_name, 
    "/", city_name, 
    "_landscape_metrics.csv"))[1,1] 
  
}

city_wide_connectivity_scaled <- center_scale(city_wide_connectivity)

city_wide_connectivity_df <- as.data.frame(cbind(
  city_names, city_wide_connectivity, city_wide_connectivity_scaled))

# initialize df with correct sites per city
df <- richness[,,1]

colnames(df) <- seq(1:max_pred_length)

df <- as.data.frame(df)

df <- df %>%
  cbind(., city_names) %>%
  pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(draw)) %>%
  filter(!is.na(.[,3])) %>%
  select(city_names, site)

# now add the richness per site per random draw of param estimates from the posterior distr
for(draw in 1:n_draws){
  
  temp <- richness[,,draw]
  
  colnames(temp) <- seq(1:max_pred_length)
  
  temp <- as.data.frame(temp)
  
  temp <- temp %>%
    cbind(., city_names) %>%
    pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(draw)) %>%
    filter(!is.na(.[,3]))
  
  df <- cbind(df, temp[,3])
  
}

# calculate city mean species richness
mean = matrix(nrow=n_cities)
lower_50 = matrix(nrow=n_cities)
upper_50 = matrix(nrow=n_cities)
lower_95 = matrix(nrow=n_cities)
upper_95 = matrix(nrow=n_cities)

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  
  temp <- filter(df, city_names == city_name) %>%
    select(-city_names, -site)
  
  quants = as.vector(quantile(as.matrix(temp), probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
  
  mean[city_number] = quants[3]
  lower_50[city_number] = quants[2]
  upper_50[city_number] = quants[4]
  lower_95[city_number] = quants[1]
  upper_95[city_number] = quants[5]
}

estimates <- as.data.frame(cbind(city_names,
                                 mean,
                                 lower_50, 
                                 upper_50,
                                 lower_95, 
                                 upper_95
                                 
)) %>%
  rename("mean" = "V2",
         "lower_50" = "V3",
         "upper_50" = "V4",
         "lower_95" = "V5",
         "upper_95" = "V6") %>%
  mutate(city = city_names) %>%
  left_join(., size_of_regional_species_pools, by = "city") %>%
  mutate(mean =as.numeric(mean),
         lower_50 =as.numeric(lower_50),
         upper_50 =as.numeric(upper_50),
         lower_95 =as.numeric(lower_95),
         upper_95 =as.numeric(upper_95)) %>%
  mutate(rel_mean = mean / size_of_regional_pool,
         rel_lower_50 = lower_50 / size_of_regional_pool,
         rel_upper_50 = upper_50 / size_of_regional_pool,
         rel_lower_95 = lower_95 / size_of_regional_pool,
         rel_upper_95 = upper_95 / size_of_regional_pool) 

estimates <- estimates %>%
  left_join(., city_wide_connectivity_df) %>%
  mutate(city_wide_connectivity = as.numeric(city_wide_connectivity),
         city_wide_connectivity_scaled = as.numeric(city_wide_connectivity_scaled))

## --------------------------------------------------
## Draw species richness plot

# plot means and variation across city_wide_connectivity_scaled
(p <- ggplot(estimates, aes(x=city_wide_connectivity_scaled)) +
   theme_bw() +
   scale_y_continuous(str_wrap("Mean Alpha Diversity /\nSize of Regional Species Pool", width = 30),
                      limits = c(0, 1), breaks = c(0, 0.5, 1)) +
   scale_x_continuous(str_wrap("City-Wide Connectivity (Scaled)", width = 30),
                      limits = c(-2, 2.5), breaks = c(-2, -1, 0, 1, 2)) +
   ggtitle("") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

p <- p +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_95, ymax=rel_upper_95, colour=city, group=city),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_50, ymax=rel_upper_50, colour=city, group=city),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=city_wide_connectivity_scaled, y=mean, colour=city, group=city), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p

## --------------------------------------------------
## Quantify association between species richness and city-wide covariate


df_w_city_covs <- left_join(df, city_wide_connectivity_df)

all_posterior_samples <- matrix(nrow=100*n_draws, ncol=2)

for(draw in 1:n_draws){
  # make a for loop across draws
  means_by_city <- df_w_city_covs %>%
    select(city_names, city_wide_connectivity_scaled, as.character(draw))
  
  colnames(means_by_city)[3] <- "richness"
  
  means_by_city <- means_by_city %>%
    group_by(city_names) %>%
    mutate(richness = as.numeric(richness),
           mean_alpha_diversity = mean(richness)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(city = city_names) %>%
    left_join(., size_of_regional_species_pools, by = "city") %>%
    mutate(rel_mean_alpha_diversity = mean_alpha_diversity / size_of_regional_pool) %>%
    mutate(city_wide_connectivity_scaled = as.numeric(city_wide_connectivity_scaled))
  
  
  stan_fit <- rstanarm::stan_glm(rel_mean_alpha_diversity ~ city_wide_connectivity_scaled, data = means_by_city)
  
  posterior_samples <- as.matrix(sample_n(as.data.frame(stan_fit)[,1:2], size = 100))
  
  all_posterior_samples[(1+100*(draw - 1)):(100+100*(draw - 1)),] <- posterior_samples
}

all_posterior_samples_df <- as.data.frame(all_posterior_samples) %>%
  rename("intercept" = "V1",
         "slope" = "V2")

hist(all_posterior_samples_df$slope)

# 
## --------------------------------------------------
## Draw species richness plot

# plot means and variation across cities
(p <- ggplot(estimates, aes(x=city_wide_connectivity_scaled)) +
   theme_bw() +
   scale_y_continuous(str_wrap("Mean Alpha Diversity /\nSize of Regional Species Pool", width = 30),
                      limits = c(0, 1), breaks = c(0, 0.5, 1)) +
   scale_x_continuous(str_wrap("City-Wide Connectivity (Scaled)", width = 30),
                      limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
   ggtitle("") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

p <- p +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_95, ymax=rel_upper_95, colour=city, group=city),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=city_wide_connectivity_scaled, ymin=rel_lower_50, ymax=rel_upper_50, colour=city, group=city),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=city_wide_connectivity_scaled, y=mean, colour=city, group=city), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p

random_lines <- sample.int(nrow(all_posterior_samples), 100)

for(i in 1:100){
  p <- p + geom_abline(intercept = all_posterior_samples[random_lines[i],1], 
                       slope = all_posterior_samples[random_lines[i],2], 
                       alpha = 0.3, colour = "grey")
}

p




# draw a plot 
(p  <- ggplot(data = means_by_city, aes(city_wide_connectivity_scaled, 
                                        rel_mean_alpha_diversity,
                                        colour = city)) +
    geom_point(aes(), size = 2) +
    geom_line() +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("city_wide_connectivity_scaled") +
    ylab("Mean Alpha Diversity /\nSize of Regional Species Pool") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)
    )
  
) 

## --------------------------------------------------
## Draw by site covariate plot

for(city in 1:n_cities){
  
}

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

# use the loess data to add the 'ribbon' to plot 
(p  <- ggplot(data = df2, aes(original_scale_park_size_data, group=City)) +
   geom_line(aes(y=rel_mean, colour=City), size = 2) +
   geom_ribbon(aes(ymin = rel_lower_95, ymax = rel_upper_95, fill=City), alpha = 0.2) +
   geom_ribbon(aes(ymin = rel_lower_50, ymax = rel_upper_50, fill=City), alpha = 0.5) +
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

## --------------------------------------------------
## analyze beta diversity of simmed communities

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
jaccard_site[is.na(jaccard_site)] <- 1 # sometimes the no difference sites come back as NaN rather than 1
jaccard_site_avg <- apply(jaccard_site,c(1,2,4),mean) # average across years
#test <- richness[1:2, 1:30, 1, 1:50]
#richness <- test

# have to minus one from pred length because we don't want to compare first site to itself
mean = matrix(nrow=n_cities, ncol=pred_length-1)
lower_50 = matrix(nrow=n_cities, ncol=pred_length-1)
upper_50 = matrix(nrow=n_cities, ncol=pred_length-1)
lower_95 = matrix(nrow=n_cities, ncol=pred_length-1)
upper_95 = matrix(nrow=n_cities, ncol=pred_length-1)

for(city_number in 1:n_cities){
  for(j in 2:pred_length){
    
    quants = as.vector(quantile(jaccard_site_avg[city_number,j,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[city_number, j-1] = quants[3]
    lower_50[city_number, j-1] = quants[2]
    upper_50[city_number, j-1] = quants[4]
    lower_95[city_number, j-1] = quants[1]
    upper_95[city_number, j-1] = quants[5]
  }
}

# make an empty df
df <- data.frame()

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  pred_data <- park_size_pred_data_list[[city_number]]
  pred_data <- pred_data[-1]
  original_scale_park_size_data <- park_size_original_data_list[[city_number]]
  original_scale_park_size_data <- original_scale_park_size_data[-1]
  
  temp <- as.data.frame(cbind(pred_data, #original_scale_data,
                              mean[city_number,],
                              lower_50[city_number,], 
                              upper_50[city_number,],
                              lower_95[city_number,], 
                              upper_95[city_number,],
                              original_scale_park_size_data
  )) %>%
    rename("mean" = "V2",
           "lower_50" = "V3",
           "upper_50" = "V4",
           "lower_95" = "V5",
           "upper_95" = "V6") %>%
    mutate(city = city_name) 
  
  df <- rbind(df, temp)
  
}

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
# this will create a smooth line between the 
pre_p2 <- ggplot(df) + 
  stat_smooth(aes(x = original_scale_park_size_data, y = mean, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = lower_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = upper_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = lower_50, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_park_size_data, y = upper_50, group=city), method = "loess", se = FALSE)
pre_p2

# build plot object for rendering 
p2 <- ggplot_build(pre_p2)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(original_scale_park_size_data = p2$data[[1]]$x,
                  mean = p2$data[[1]]$y,
                  lower_95 = p2$data[[2]]$y,
                  upper_95 = p2$data[[3]]$y,
                  lower_50 = p2$data[[4]]$y,
                  upper_50 = p2$data[[5]]$y) 

df2 <- cbind(df2, rep(city_names, each = nrow(df2)/n_cities))

colnames(df2)[7] <- "City"

# use the loess data to add the 'ribbon' to plot 
(p2  <- ggplot(data = df2, aes(original_scale_park_size_data, group=City)) +
    geom_line(aes(y=mean, colour=City), size = 2) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill=City), alpha = 0.2) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, fill=City), alpha = 0.5) +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("log(Park Size in m^2)") +
    ylab("Species Dissimilarity\n(Jaccard Index)") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)
    )
  
)

## --------------------------------------------------
# site predictors


# create empty list of length n_cities
# each element of the list holds the posterior distributions for 
# all parameters for each individual city 
#site_data_list <- vector(mode='list', length=n_cities)
site_data <- data.frame()

# read and summarize data across loop
for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- readRDS(paste0(
      "./run_model/prepped_data/prepped_data_", city, ".rds"))
  
  city_site_data <- temp$site_data %>%
    mutate(City = city)
  
  site_data <- rbind(site_data, city_site_data)
  
}

# park size
p1 <- ggplot(site_data, aes(x = log_total_green_space_area, 
                                colour = City, fill = City)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  scale_x_continuous(breaks = c(6, 9, 12, 15, 18)) +
  theme_classic() +
  xlab("log(Park Size in m^2)") + 
  ylab("Number of Park Clusters\n ") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cowplot::plot_grid(p, p1, ncol = 1)

##

cowplot::plot_grid(p, p2, p1, ncol = 1)

##------------------------------------------------------------------------------
# now do it for connectivity

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_size_pred_data_list <- vector(mode='list', length=n_cities)
park_size_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- 30 # something short for testing

for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- readRDS( paste0("./run_model/prepped_data/prepped_data_", city, ".rds"))
  
  pred_data <- seq(from = min(temp$site_data$avg_dist_2000m_scaled), 
                   to = max(temp$site_data$avg_dist_2000m_scaled), 
                   length.out = pred_length)
  
  mean_isolation <- mean(temp$site_data$avg_dist_2000m)
  sd_isolation <- sd(temp$site_data$avg_dist_2000m)
  # now do some algebra to get the scaled data back onto a real life distance scale
  original_scale_isolation_data <- (pred_data * sd_isolation) + mean_isolation
  
  park_size_pred_data_list[[city_number]] <- pred_data
  park_size_original_data_list[[city_number]] <- original_scale_isolation_data
  
}

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## plot species richness by park size

## --------------------------------------------------
## get prediction range

n_years = 5 # number of years of the study
n_years_minus1 = n_years - 1 # number of interannual transitions

n_draws = 100 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors

#pred_length = 30
#pred_data <- seq(from = -2, to = 2, length.out = pred_length)

## --------------------------------------------------

# some random samples from the posterior
random_draws_from_posterior = sample.int(n=n_draws) # use if not using the full posterior

# take random draws for psi and predict occurrence, then sum across n species
richness <- array(data = NA, dim=c(n_cities, pred_length, n_years, n_draws))
jaccard_site <- array(data = NA, dim=c(n_cities, pred_length, n_years, n_draws))

for(city_number in 1:n_cities){
  
  # get the posterior distributions for a particular city
  city_estimates <- mega_list[[city_number]]
  
  # get indices for species random effects distributions for particular city
  first_psi1 <- which( colnames(city_estimates)=="psi1_species[1]" )
  first_gamma <- which( colnames(city_estimates)=="gamma_species[1]" )
  first_phi <- which( colnames(city_estimates)=="phi_species[1]" )
  
  # get the trait data for species from the particular city
  city_name <- city_names[city_number]
  my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".rds"))
  
  pred_data <- park_size_pred_data_list[[city_number]]
  
  species_data <- my_data$species_info
  n_species <- nrow(species_data)
  
  # construct some arrays to hold the data
  psi1_expected <- array(data = NA, dim=c(n_species, pred_length))
  gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
  phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
  occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years, n_draws))
  #psi <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))
  
  for(draw in 1:n_draws){
    
    rand <- random_draws_from_posterior[draw]
    
    # expected occurrence in year 1
    for(i in 1:n_species){
      for(j in 1:pred_length){
        for(k in 2:n_years){
          
          
          psi1_expected[i,j] =
            ilogit(
              # YEAR 1 is the global intercept
              city_estimates[rand,1] + 
                # a species specific intercept effect (the number here should be first column)
                city_estimates[rand,(first_psi1+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,3] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,5] * pred_data[j] 
              # not adding any other park predictors would keep park at same value
            )
          
          
          gamma_expected[i,j,k-1] = # gamma[,,k-1] yields gamma for the transition between year 1 and 2
            ilogit(#gamma0 +
              city_estimates[rand,7] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_gamma+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,9] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,11] * pred_data[j] 
              # not adding any other park predictors would keep park at same value
            )
          
          
          phi_expected[i,j,k-1] = # phi[,,k-1] yields phi for the transition between year 1 and 2
            ilogit(#phi0 +
              city_estimates[rand,13] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_phi+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,15] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,17] * pred_data[j] 
              # not adding any other park predictors would keep park at same value
            )
          
        } 
      }
    }
    
    # simmed occurrence for each species
    for(i in 1:n_species){
      for(j in 1:pred_length){
        for(k in 1:n_years){
          if(k < 2){ # in year 1
            occurrence_simmed[i,j,1,rand] <- rbinom(1, 1, prob = psi1_expected[i,j])
          } else{ # contingent on previous year
            occurrence_simmed[i,j,k,rand] = occurrence_simmed[i,j,k-1,rand] * phi_expected[i,j,k-1] + 
              (1 - occurrence_simmed[i,j,k-1,rand]) * gamma_expected[i,j,k-1]
          }
        }
      }
    }
    
    for(j in 1:pred_length){
      for(k in 1:n_years){
        richness[city_number,j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k,rand])
      }
    }
    
    ## --------------------------------------------------
    ## calculate beta diversity of simmed communities
    
    # jaccard index needs a reference level
    ref_site <- 1
    
    # now compute dissimilarity in the occurrence matrix in site i versus site 1
    for(j in 1:pred_length){ # jaccard index for sites (in terms of shared species)
      for(k in 1:n_years){
        jaccard_site[city_number,j,k,draw] <- sum(occurrence_simmed[,ref_site,k,rand]*occurrence_simmed[,j,k,rand]) /
          (sum(occurrence_simmed[,ref_site,k,rand]) +
             sum(occurrence_simmed[,ref_site,k,rand]) - sum(occurrence_simmed[,ref_site,k,rand]*occurrence_simmed[,j,k,rand]))
      }
    }
    
  } # end for random draw
  
} # end for city

## --------------------------------------------------
# summarize the results

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
richness <- apply(richness,c(1,2,4),mean) # average across years
#test <- richness[1:2, 1:30, 1, 1:50]
#richness <- test

mean = matrix(nrow=n_cities, ncol=pred_length)
lower_50 = matrix(nrow=n_cities, ncol=pred_length)
upper_50 = matrix(nrow=n_cities, ncol=pred_length)
lower_95 = matrix(nrow=n_cities, ncol=pred_length)
upper_95 = matrix(nrow=n_cities, ncol=pred_length)

for(city_number in 1:n_cities){
  for(j in 1:pred_length){
    
    quants = as.vector(quantile(richness[city_number,j,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[city_number, j] = quants[3]
    lower_50[city_number, j] = quants[2]
    upper_50[city_number, j] = quants[4]
    lower_95[city_number, j] = quants[1]
    upper_95[city_number, j] = quants[5]
  }
}

# make an empty df
df <- data.frame()

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  pred_data <- park_size_pred_data_list[[city_number]]
  original_scale_isolation_data <- park_size_original_data_list[[city_number]]
  
  temp <- as.data.frame(cbind(pred_data, #original_scale_data,
                              mean[city_number,],
                              lower_50[city_number,], 
                              upper_50[city_number,],
                              lower_95[city_number,], 
                              upper_95[city_number,],
                              original_scale_isolation_data
  )) %>%
    rename("mean" = "V2",
           "lower_50" = "V3",
           "upper_50" = "V4",
           "lower_95" = "V5",
           "upper_95" = "V6") %>%
    mutate(city = city_name) %>%
    left_join(., size_of_regional_species_pools, by = "city") %>%
    mutate(rel_mean = mean / size_of_regional_pool,
           rel_lower_50 = lower_50 / size_of_regional_pool,
           rel_upper_50 = upper_50 / size_of_regional_pool,
           rel_lower_95 = lower_95 / size_of_regional_pool,
           rel_upper_95 = upper_95 / size_of_regional_pool)
  
  df <- rbind(df, temp)
  
}

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
# this will create a smooth line between the 
pre_q <- ggplot(df) + 
  stat_smooth(aes(x = original_scale_isolation_data, y = rel_mean, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = rel_lower_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = rel_upper_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = rel_lower_50, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = rel_upper_50, group=city), method = "loess", se = FALSE)
pre_q

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

## --------------------------------------------------
## analyze beta diversity of simmed communities

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
jaccard_site[is.na(jaccard_site)] <- 1 # sometimes the no difference sites come back as NaN rather than 1
jaccard_site_avg <- apply(jaccard_site,c(1,2,4),mean) # average across years
#test <- richness[1:2, 1:30, 1, 1:50]
#richness <- test

# have to minus one from pred length because we don't want to compare first site to itself
mean = matrix(nrow=n_cities, ncol=pred_length-1)
lower_50 = matrix(nrow=n_cities, ncol=pred_length-1)
upper_50 = matrix(nrow=n_cities, ncol=pred_length-1)
lower_95 = matrix(nrow=n_cities, ncol=pred_length-1)
upper_95 = matrix(nrow=n_cities, ncol=pred_length-1)

for(city_number in 1:n_cities){
  for(j in 2:pred_length){
    
    quants = as.vector(quantile(jaccard_site_avg[city_number,j,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[city_number, j-1] = quants[3]
    lower_50[city_number, j-1] = quants[2]
    upper_50[city_number, j-1] = quants[4]
    lower_95[city_number, j-1] = quants[1]
    upper_95[city_number, j-1] = quants[5]
  }
}

# make an empty df
df <- data.frame()

for(city_number in 1:n_cities){
  
  city_name <- city_names[city_number]
  pred_data <- park_size_pred_data_list[[city_number]]
  pred_data <- pred_data[-1]
  original_scale_isolation_data <- park_size_original_data_list[[city_number]]
  original_scale_isolation_data <- original_scale_isolation_data[-1]
  
  temp <- as.data.frame(cbind(pred_data, #original_scale_data,
                              mean[city_number,],
                              lower_50[city_number,], 
                              upper_50[city_number,],
                              lower_95[city_number,], 
                              upper_95[city_number,],
                              original_scale_isolation_data
  )) %>%
    rename("mean" = "V2",
           "lower_50" = "V3",
           "upper_50" = "V4",
           "lower_95" = "V5",
           "upper_95" = "V6") %>%
    mutate(city = city_name) 
  
  df <- rbind(df, temp)
  
}

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
# this will create a smooth line between the 
pre_q2 <- ggplot(df) + 
  stat_smooth(aes(x = original_scale_isolation_data, y = mean, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = lower_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = upper_95, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = lower_50, group=city), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_isolation_data, y = upper_50, group=city), method = "loess", se = FALSE)
pre_q2

# build plot object for rendering 
q2 <- ggplot_build(pre_q2)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(original_scale_isolation_data = q2$data[[1]]$x,
                  mean = q2$data[[1]]$y,
                  lower_95 = q2$data[[2]]$y,
                  upper_95 = q2$data[[3]]$y,
                  lower_50 = q2$data[[4]]$y,
                  upper_50 = q2$data[[5]]$y) 

df2 <- cbind(df2, rep(city_names, each = nrow(df2)/n_cities))

colnames(df2)[7] <- "City"

# use the loess data to add the 'ribbon' to plot 
(q2  <- ggplot(data = df2, aes(original_scale_isolation_data, group=City)) +
    geom_line(aes(y=mean, colour=City), size = 2) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill=City), alpha = 0.2) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, fill=City), alpha = 0.5) +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("Isolation\n(Avg. Distance to Parks Within 2km)") +
    ylab("Species Dissimilarity\n(Jaccard Index)") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)
    )
  
)

## --------------------------------------------------
# site predictors


# create empty list of length n_cities
# each element of the list holds the posterior distributions for 
# all parameters for each individual city 
#site_data_list <- vector(mode='list', length=n_cities)
site_data <- data.frame()

# read and summarize data across loop
for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- readRDS(paste0(
    "./run_model/prepped_data/prepped_data_", city, ".rds"))
  
  city_site_data <- temp$site_data %>%
    mutate(City = city)
  
  site_data <- rbind(site_data, city_site_data)
  
}

# park size
q1 <- ggplot(site_data, aes(x = avg_dist_2000m, 
                            colour = City, fill = City)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500)) +
  theme_classic() +
  xlab("Isolation\n(Avg. Distance to Parks Within 2km)") + 
  ylab("Number of Park Clusters\n ") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##

cowplot::plot_grid(q, q2, q1, ncol = 1)

