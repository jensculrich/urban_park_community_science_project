# plot effects of predictors on occurrence and detection

# load libraries
library(tidyverse)
library(cmdstanr)

## get param estimates from the region
stan_out <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_apr9.rds")

estimates <- as.data.frame(stan_out$draws(variables = c("psi_0", 
                                                  "sigma_psi_species",
                                                  "sigma_psi_city",
                                                  "psi_wingspan",
                                                  "psi_migratory",
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
                                                  "p_migratory",
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
                                                  "p_city"),
                                    format = "draws_matrix"
                                    )) # take estimates from each HMC step as a df

rm(stan_out)
gc()

#n_samp <- 10 # how mann_chains#n_samp <- 10 # how many samples do we have from the HMC run?
n_samp <- length(estimates[,1]) # how many samples do we have from the HMC run?

## get data from region
df <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_all.rds"))$site_data

# length of number of parks in the cities
pred_length <- nrow(df)

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# get important column names
psi_0 <- which( colnames(estimates)=="psi_0" )
psi_parksize <- which( colnames(estimates)=="mu_psi_park_size" )
psi_tree_cover <- which( colnames(estimates)=="mu_psi_tree_cover" )
psi_plant_diversity <- which( colnames(estimates)=="mu_psi_plant_diversity" )
psi_isolation <- which( colnames(estimates)=="mu_psi_landscape_isolation" )
psi_landscape_grassherb <- which( colnames(estimates)=="mu_psi_landscape_grassherb" )
psi_landscape_woody <- which( colnames(estimates)=="mu_psi_landscape_woody" )

#-------------------------------------------------------------------------------
# get some prediction data

# these are a sequence of scaled values of park sizes (what the model sees)
scaled_pred <- seq(from = -2, to = 2, length.out = pred_length)

#-------------------------------------------------------------------------------
# occurrence - park size (psi)

occurrence <- vector(length = pred_length)

# get length of stan fit object (HMC iterations * n_chains)  
predMean_parksize <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predMean_plantdiversity <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predMean_treecover <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predMean_isolation <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predMean_landscapeherb <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predMean_landscapewoody <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean_parksize[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_parksize]*scaled_pred
  )
  # community means don't depend on city effects
  predMean_plantdiversity[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_plant_diversity ]*scaled_pred
  )
  # community means don't depend on city effects
  predMean_treecover[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_tree_cover]*scaled_pred
  )
  # community means don't depend on city effects
  predMean_isolation[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_isolation]*scaled_pred
  )
  # community means don't depend on city effects
  predMean_landscapeherb[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_landscape_grassherb]*scaled_pred
  )
  # community means don't depend on city effects
  predMean_landscapewoody[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    estimates[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      estimates[i,psi_landscape_woody]*scaled_pred
  )
    
}


# posterior means by community average 
criMean_parksize <- apply(predMean_parksize, c(1), function(x) quantile(x, 
              prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

criMean_plantdiversity <- apply(predMean_plantdiversity, c(1), function(x) quantile(x, 
                                                                        prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

criMean_treecover <- apply(predMean_treecover, c(1), function(x) quantile(x, 
                                                                        prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

criMean_isolation <- apply(predMean_isolation, c(1), function(x) quantile(x, 
                                                                        prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

criMean_landscapeherb <- apply(predMean_landscapeherb, c(1), function(x) quantile(x, 
                                                                        prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

criMean_landscapewoody <- apply(predMean_landscapewoody, c(1), function(x) quantile(x, 
                                                                        prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

#-------------------------------------------------------------------------------

# community plot - park size - psi1
# organize the mean, 50 and 90% BCIs into a data frame
df1 <- as.data.frame(cbind(scaled_pred, criMean_parksize[3,], 
                           criMean_parksize[1,], criMean_parksize[5,],
                           criMean_parksize[2,], criMean_parksize[4,]),
                         ) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Park Size")

df2 <- as.data.frame(cbind(scaled_pred, criMean_plantdiversity[3,], 
                           criMean_plantdiversity[1,], criMean_plantdiversity[5,],
                           criMean_plantdiversity[2,], criMean_plantdiversity[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Local Plant Diversity")

df3 <- as.data.frame(cbind(scaled_pred, criMean_treecover[3,], 
                           criMean_treecover[1,], criMean_treecover[5,],
                           criMean_treecover[2,], criMean_treecover[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Local Tree Cover")

df4 <- as.data.frame(cbind(scaled_pred, criMean_isolation[3,], 
                           criMean_isolation[1,], criMean_isolation[5,],
                           criMean_isolation[2,], criMean_isolation[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Landscape Isolation")

df5 <- as.data.frame(cbind(scaled_pred, criMean_landscapeherb[3,], 
                           criMean_landscapeherb[1,], criMean_landscapeherb[5,],
                           criMean_landscapeherb[2,], criMean_landscapeherb[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Landscape Herbaceous Cover")

df6 <- as.data.frame(cbind(scaled_pred, criMean_landscapewoody[3,], 
                           criMean_landscapewoody[1,], criMean_landscapewoody[5,],
                           criMean_landscapewoody[2,], criMean_landscapewoody[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "mean" = "V2",
         "lower95" = "V3",
         "upper95" = "V4",
         "lower50" = "V5",
         "upper50" = "V6") %>%
  mutate(Predictor = "Landscape Woody Cover")

df_combined <- rbind(df1, df2, df3, df4, df5, df6)
df_combined1 <- rbind(df1, df2, df3) # local
df_combined2 <- rbind(df4, df5, df6) # landscape

# order factor variables for consistency
df_combined1$Predictor <- fct_relevel(df_combined1$Predictor, 
      "Local Tree Cover", "Local Plant Diversity", "Park Size")

df_combined2$Predictor <- fct_relevel(df_combined2$Predictor, 
      "Landscape Woody Cover", "Landscape Herbaceous Cover", "Landscape Isolation")

# choose palette
my_palette <- viridis::viridis(n=7, option = "viridis")
my_palette1 <- my_palette[c(1, 3, 5)] # drop black colour
my_palette2 <- my_palette[c(2, 4, 6)] # drop black colour

# plot the estimated relationship for the average city in the region
p1 <- ggplot(data = df_combined1, aes(scaled_pred, mean, colour=Predictor, fill=Predictor)) +
  geom_ribbon(aes(
    ymin=lower95, 
    ymax=upper95), alpha=0.1, colour = NA) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.2, colour = NA) +
  geom_line(size=2, lty=1) +
  scale_colour_manual(name="", values=my_palette1) +  
  scale_fill_manual(name="", values=my_palette1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Park-Site Predictor Value\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  theme(legend.position = c(0.025, 0.9755), # x=1 (right), y=0 (bottom)
        legend.justification = c(0, 1), # Justify the top-left corner of the legend box to these coordinates
        legend.text = element_text(size=16),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1

# plot the estimated relationship for the average city in the region
p2 <- ggplot(data = df_combined2, aes(scaled_pred, mean, colour=Predictor, fill=Predictor)) +
  geom_ribbon(aes(
    ymin=lower95, 
    ymax=upper95), alpha=0.1, colour = NA) +
  geom_ribbon(aes(
    ymin=lower50, 
    ymax=upper50), alpha=0.2, colour = NA) +
  geom_line(size=2, lty=1) +
  scale_colour_manual(name="", values=my_palette2) +  
  scale_fill_manual(name="", values=my_palette2) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Park-Site Predictor Value\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  theme(legend.position = c(0.025, 0.9755), # x=1 (right), y=0 (bottom)
        legend.justification = c(0, 1), # Justify the top-left corner of the legend box to these coordinates
        legend.text = element_text(size=16),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2

cowplot::plot_grid(p1, p2, ncol = 2)

#-------------------------------------------------------------------------------
# now plot some city-specific effects

city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Denver",
  "Des_moines",
  "Detroit",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "SD",
  "SF",
  "St_louis",
  "Tampa"
)

city_names_labels <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "Washington D.C.",
  "Denton",
  "Denver",
  "Des Moines",
  "Detroit",
  "Houston",
  "Los Angeles",
  "Minneapolis",
  "New York City",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "San Diego",
  "San Fransisco",
  "St. Louis",
  "Tampa"
)

n_cities <- length(city_names)


#-------------------------------------------------------------------------------
# occupancy - local park size

#-------------------------------------------------------------------------------
# for real park sizes in each city

#-------------------------------------------------------------------------------
# for real park sizes in each city 
# (not scaled values which are hard to interpret in a practical sense)

#-------------------------------------------------------------------------------
# get some prediction data

## get data from region
df <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))$site_data

# length of number of parks in the cities
pred_length <- nrow(df)

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_size_pred_data_list <- vector(mode='list', length=n_cities)
park_size_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

# remember that the predictor values were scaled within each city, so we have to
# come up with some pred data specific to the real values within each city (not across all cities)
for(city_number in 1:n_cities){
  
  # filter the site covariate data to the specific city
  temp <- filter(df, city == city_names[city_number]) %>%
    arrange(., log_total_green_space_area)
  
  # get the scaled park size data
  park_size_pred_data <- temp$log_total_green_space_area_scaled_across_all_cities
  # get the real park size data
  original_scale_park_size_data <- temp$log_total_green_space_area
  
  # store the city-specific scaled and real values
  park_size_pred_data_list[[city_number]] <- park_size_pred_data
  park_size_original_data_list[[city_number]] <- original_scale_park_size_data
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_size_pred_data)
  
}

# set the maximum pred length to the largest number of sites within a city
max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# estimate occurrence rate

# prediction trends by city 
predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) 

# get indices for random effects distributions for any particular city
first_psi <- which( colnames(estimates)=="psi_city[1]" )
first_psi_size <- which( colnames(estimates)=="psi_park_size[1]" )

# loop across all sites in each city, and make a prediction about occurrence
# for each of the samples from the posterior probability distribution
for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # estimate city specific predictions about occurrence rate
      predCity[city_number,i,j] <- ilogit( 
        # psi1_0 + # global intercept
        #estimates[j,psi_0] + 
        # psi1_city + # city effect on the intercept
        # in our current model the city intercept includes both the mean intercept and the deviation
        # for the city, so don't also include the global mean
        estimates[j,(first_psi+(city_number-1))] +  
          # psi1_park_size + # city specific effect of park size given real park size data
          estimates[j,(first_psi_size+(city_number-1))] * park_size_pred_data_list[[city_number]][i]
      )
      
    }
  }
}

# now we have to transform all of this into a data frame to be able to plot in ggplot
# initialize df with correct sites per city
new_df <- predCity[,,1]

# col names are sites within cities
colnames(new_df) <- seq(1:max_pred_length)

# convert to df format
new_df <- as.data.frame(new_df)

# cast sites into long format (multiple rows per city rather than one row per city)
new_df <- new_df %>%
  cbind(., city_names) %>%
  pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(j)) %>%
  filter(!is.na(.[,3])) %>%
  select(city_names, site)

# now add the occurrence per random draw of param estimates from the posterior distr
# for every site. Then later we will summarize the quantiles of these predictions for each site
for(j in 1:n_samp){
  
  temp <- predCity[,,j]
  
  colnames(temp) <- seq(1:max_pred_length)
  
  temp <- as.data.frame(temp)
  
  temp <- temp %>%
    cbind(., city_names) %>%
    pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(j)) %>%
    filter(!is.na(.[,3]))
  
  new_df <- cbind(new_df, temp[,3])
  
}

# get the mean, 50 and 90% BCIs for initial occurrence for each site
quants <- array(NA, dim=c(nrow=nrow(new_df), ncol=5))

for(i in 1:nrow(new_df)){
  quants[i,] <- as.vector(quantile(as.numeric(new_df[i,3:n_samp]), probs = c(0.05, 0.25, 0.50, 0.75, 0.95))) 
}

# get site data in ordered form
park_size_original_ordered <- unlist(park_size_original_data_list)
park_size_pred_ordered <- unlist(park_size_pred_data_list)

# now replace the entire posterior distribution of draws in new_df with the quantile summaries
new_df <- as.data.frame(cbind(new_df$city_names, new_df$site, quants)) %>%
  rename("city" = "V1",
         "site" = "V2",
         "lower_90" = "V3",
         "lower_50" = "V4",
         "mean" = "V5",
         "upper_50" = "V6",
         "upper_90" = "V7") %>%
  mutate(mean =as.numeric(mean),
         lower_50 =as.numeric(lower_50),
         upper_50 =as.numeric(upper_50),
         lower_90 =as.numeric(lower_90),
         upper_90 =as.numeric(upper_90)) %>%
  # and join with the site data
  cbind(., park_size_original_ordered, park_size_pred_ordered)

## --------------------------------------------------
## Draw multicity plot

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

# replot q to remove the legend automatically
# plot the relationships
q <- ggplot(data = new_df, aes(x=park_size_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.05, colour = NA) +
  #geom_ribbon(aes(
  #ymin=lower_90, 
  #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=1, lty=1) +
  xlim(c(min(park_size_original_ordered), max(park_size_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate \n(City-Specific)") +
  xlab("log(Park Size m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette, labels=city_names_labels ) + 
  
  scale_fill_manual(values=my_palette, labels=city_names_labels) + 
  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
q


#-------------------------------------------------------------------------------
# occupancy - landscape isolation

#-------------------------------------------------------------------------------
# for real park sizes in each city

#-------------------------------------------------------------------------------
# get some prediction data

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_isolation_pred_data_list <- vector(mode='list', length=n_cities)
park_isolation_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

for(city_number in 1:n_cities){
  
  temp <- filter(df, city == city_names[city_number]) %>%
    arrange(., isolation)
  
  # get park size data
  park_isolation_pred_data <- temp$log_isolation_scaled_across_all_cities
  original_scale_park_isolation_data <- log(temp$isolation)
  
  park_isolation_pred_data_list[[city_number]] <- park_isolation_pred_data
  park_isolation_original_data_list[[city_number]] <- original_scale_park_isolation_data
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_isolation_pred_data)
  
}

max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# predict occupancy

predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) # trends by city 

# get indices for species random effects distributions for particular city
first_psi <- which( colnames(estimates)=="psi_city[1]" )
first_psi_isolation <- which( colnames(estimates)=="psi_landscape_isolation[1]" )

for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # community means don't depend on city effects
      predCity[city_number,i,j] <- ilogit( # park isolation trend
        # psi1_0 +
        #estimates[j,1] + 
          # psi1_city +
          estimates[j,(first_psi+(city_number-1))] + 
          # psi1_park_isolation +
          estimates[j,(first_psi_isolation+(city_number-1))] * park_isolation_pred_data_list[[city_number]][i]
      )
      
    }
  }
}

# initialize df with correct sites per city
new_df <- predCity[,,1]

colnames(new_df) <- seq(1:max_pred_length)

new_df <- as.data.frame(new_df)

new_df <- new_df %>%
  cbind(., city_names) %>%
  pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(j)) %>%
  filter(!is.na(.[,3])) %>%
  select(city_names, site)

# now add the initial occurrence per site per random draw of param estimates from the posterior distr
for(j in 1:n_samp){
  
  temp <- predCity[,,j]
  
  colnames(temp) <- seq(1:max_pred_length)
  
  temp <- as.data.frame(temp)
  
  temp <- temp %>%
    cbind(., city_names) %>%
    pivot_longer(1:max_pred_length, names_to = "site", values_to = as.character(j)) %>%
    filter(!is.na(.[,3]))
  
  new_df <- cbind(new_df, temp[,3])
  
}

quants <- array(NA, dim=c(nrow=nrow(new_df), ncol=5))

for(i in 1:nrow(new_df)){
  quants[i,] <- as.vector(quantile(as.numeric(new_df[i,3:n_samp]), probs = c(0.05, 0.25, 0.50, 0.75, 0.95))) 
}

# get site data in ordered form
park_isolation_original_ordered <- unlist(park_isolation_original_data_list)
park_isolation_pred_ordered <- unlist(park_isolation_pred_data_list)

new_df <- as.data.frame(cbind(new_df$city_names, new_df$site, quants)) %>%
  rename("city" = "V1",
         "site" = "V2",
         "lower_90" = "V3",
         "lower_50" = "V4",
         "mean" = "V5",
         "upper_50" = "V6",
         "upper_90" = "V7") %>%
  mutate(mean =as.numeric(mean),
         lower_50 =as.numeric(lower_50),
         upper_50 =as.numeric(upper_50),
         lower_90 =as.numeric(lower_90),
         upper_90 =as.numeric(upper_90)) %>%
  cbind(., park_isolation_original_ordered, park_isolation_pred_ordered)

## --------------------------------------------------
## Draw multicity plot

q2 <- ggplot(data = new_df, aes(x=park_isolation_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.05, colour = NA) +
  geom_line(size=1, lty=1) +
  xlim(c(min(park_isolation_original_ordered), max(park_isolation_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(City-Specific)") +
  xlab("log(Park Isolation)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette) + 
                                         
  scale_fill_manual(values=my_palette) + 
                                        
  theme(legend.position = "none",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
q2

q_legend <- ggplot(data = new_df, aes(x=park_isolation_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.05) +
  geom_line(size=1, lty=1) +
  xlim(c(min(park_isolation_original_ordered), max(park_isolation_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(City-Specific)") +
  xlab("log(Park Isolation)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(name = "City: ", labels=city_names_labels, values=my_palette) + 
  
  scale_fill_manual(name = "City: ", labels=city_names_labels, values=my_palette) + 
  
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

q_legend

legend <- cowplot::get_legend(q_legend)
plot(legend)

#-------------------------------------------------------------------------------
# cowplot

grid <- cowplot::plot_grid(p1, p2, q, q2, nrow = 2, 
                   rel_heights = c(1, 1, 0.5),
                   labels = c('a)', 'b)', 'c)', 'd)'),
                   label_size = 20)

grid_w_legend <- cowplot::plot_grid(grid, legend, 
                                    ncol = 1, 
                                    rel_heights = c(5,1))
grid_w_legend
# save as 