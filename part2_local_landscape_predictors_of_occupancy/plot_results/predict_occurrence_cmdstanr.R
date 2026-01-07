# plot effects of predictors on occurrence and detection

# load libraries
library(tidyverse)
library(cmdstanr)

# enter the region/regions you want to plot
# currently I think this will only work if you enter one single region,
# but I think eventually we want to plot multiple regionss simultaneously


# select a region
regions <- c(
  "northeast",
  "southeast",
  "texas",
  "california",
  "all"
)

region <- regions[5]

# all
if(region == regions[5]){
  city_names <- c(
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
  )
}


n_cities <- length(city_names)

## get param estimates from the region
stan_out <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_jan2.rds")

tmp <- as.data.frame(stan_out$draws(variables = c("psi_0", 
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
                                                  "p_city"),
                                    format = "draws_matrix"
                                    )) # take estimates from each HMC step as a df

rm(stan_out)
gc()

#n_samp <- 10 # how mann_chains#n_samp <- 10 # how many samples do we have from the HMC run?
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?

## get data from region
df <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_all.rds"))$site_data

# length of number of parks in the cities
pred_length <- nrow(df)

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

# get important column names
psi_0 <- which( colnames(tmp)=="psi_0" )
psi_parksize <- which( colnames(tmp)=="mu_psi_park_size" )
psi_tree_cover <- which( colnames(tmp)=="mu_psi_tree_cover" )
psi_plant_diversity <- which( colnames(tmp)=="mu_psi_plant_diversity" )
psi_isolation <- which( colnames(tmp)=="mu_psi_landscape_isolation" )
psi_landscape_grassherb <- which( colnames(tmp)=="mu_psi_landscape_grassherb" )
psi_landscape_woody <- which( colnames(tmp)=="mu_psi_landscape_woody" )

#-------------------------------------------------------------------------------
# get some prediction data

# these are a sequence of scaled values of park sizes (what the model sees)
scaled_pred <- seq(from = -2, to = 2, length.out = pred_length)

#-------------------------------------------------------------------------------
# occurrence - park size (psi)

occurrence <- vector(length = pred_length)

# get length of stan fit object (HMC iterations * n_chains)  
predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    tmp[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      tmp[i,psi_parksize]*scaled_pred
  )
    
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
              prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

#-------------------------------------------------------------------------------

# community plot - park size - psi1
# organize the mean, 50 and 90% BCIs into a data frame
size_df <- as.data.frame(cbind(scaled_pred, criMean[3,], 
                               criMean[1,], criMean[5,],
                               criMean[2,], criMean[4,]),
                         ) %>%
  rename("scaled_pred" = "scaled_pred",
         "psi_size_community_mean" = "V2",
         "psi_size_community_lower95" = "V3",
         "psi_size_community_upper95" = "V4",
         "psi_size_community_lower50" = "V5",
         "psi_size_community_upper50" = "V6")

# plot the estimated relationship for the average city in the region
p <- ggplot(data = size_df, aes(scaled_pred, psi_size_community_mean)) +
  geom_ribbon(aes(
    ymin=psi_size_community_lower50, 
    ymax=psi_size_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi_size_community_lower95, 
    ymax=psi_size_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Park Size\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#-------------------------------------------------------------------------------
# for real park sizes in each city 
# (not scaled values which are hard to interpret in a practical sense)

#-------------------------------------------------------------------------------
# get some prediction data

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
# initial occurrence (psi1)

# prediction trends by city 
predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) 

# get indices for random effects distributions for any particular city
first_psi <- which( colnames(tmp)=="psi_city[1]" )
first_psi_size <- which( colnames(tmp)=="psi_park_size[1]" )

# loop across all sites in each city, and make a prediction about initial occurrence
# for each of the samples from the posterior probability distribution
for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # estimate city specific predictions about initial occurrence rate
      predCity[city_number,i,j] <- ilogit( 
        # psi1_0 + # global intercept
        tmp[j,psi_0] + 
          # psi1_city + # city effect on the intercept
          tmp[j,(first_psi+(city_number-1))] + 
          # psi1_park_size + # city specific effect of park size given real park size data
          tmp[j,(first_psi_size+(city_number-1))] * park_size_pred_data_list[[city_number]][i]
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

# now add the initial occurrence per random draw of param estimates from the posterior distr
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

# plot the relationships
q <- ggplot(data = new_df, aes(x=park_size_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.3) +
  #geom_ribbon(aes(
    #ymin=lower_90, 
    #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_size_original_ordered), max(park_size_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate \n(City-Specific)") +
  xlab("log(Park Size m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette) + 
                                          
  scale_fill_manual(values=my_palette) + 
                                       
  theme(#legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
q


#-------------------------------------------------------------------------------
# cowplot

# cowplot::plot_grid(p, q, ncol = 2, rel_widths = c(1, 1.5))


#-------------------------------------------------------------------------------
# occurrence - plant diversity (psi)

# length of number of parks in the cities
pred_length <- nrow(df)

# get length of stan fit object (HMC iterations * n_chains)  
predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    tmp[i,psi_0] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      tmp[i,psi_plant_diversity]*scaled_pred
  )
  
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
                                                      prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

#-------------------------------------------------------------------------------

# community plot - park size - psi1
# organize the mean, 50 and 90% BCIs into a data frame
plant_diversity_df <- as.data.frame(cbind(scaled_pred, criMean[3,], 
                               criMean[1,], criMean[5,],
                               criMean[2,], criMean[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "psi_pd_community_mean" = "V2",
         "psi_pd_community_lower95" = "V3",
         "psi_pd_community_upper95" = "V4",
         "psi_pd_community_lower50" = "V5",
         "psi_pd_community_upper50" = "V6")

# plot the estimated relationship for the average city in the region
r <- ggplot(data = plant_diversity_df, aes(scaled_pred, psi_pd_community_mean)) +
  geom_ribbon(aes(
    ymin=psi_pd_community_lower50, 
    ymax=psi_pd_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi_pd_community_lower95, 
    ymax=psi_pd_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Park Size\n(Std. Deviations from Mean Park Size)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
r

#-------------------------------------------------------------------------------
# for real plant diversity observed in parks in each city 
# (not scaled values which are hard to interpret in a practical sense)

#-------------------------------------------------------------------------------
# get some prediction data

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_pd_pred_data_list <- vector(mode='list', length=n_cities)
park_pd_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

# remember that the predictor values were scaled within each city, so we have to
# come up with some pred data specific to the real values within each city (not across all cities)
for(city_number in 1:n_cities){
  
  # filter the site covariate data to the specific city
  temp <- filter(df, city == city_names[city_number]) %>%
    arrange(., plant_genera_density)
  
  # get the scaled park size data
  park_pd_pred_data <- temp$plant_genera_density_scaled
  # get the real park size data
  original_scale_park_pd_data <- temp$plant_genera_density
  
  # store the city-specific scaled and real values
  park_pd_pred_data_list[[city_number]] <- park_pd_pred_data
  park_pd_original_data_list[[city_number]] <- original_scale_park_pd_data
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_pd_pred_data)
  
}

# set the maximum pred length to the largest number of sites within a city
max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# and make predictions

# prediction trends by city 
predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) 

# get indices for random effects distributions for any particular city
first_psi <- which( colnames(tmp)=="psi_city[1]" )
first_psi_pd <- which( colnames(tmp)=="psi_plant_diversity[1]" )

# loop across all sites in each city, and make a prediction about initial occurrence
# for each of the samples from the posterior probability distribution
for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # estimate city specific predictions about initial occurrence rate
      predCity[city_number,i,j] <- ilogit( 
        # psi1_0 + # global intercept
        tmp[j,psi_0] + 
          # psi1_city + # city effect on the intercept
          tmp[j,(first_psi+(city_number-1))] + 
          # psi1_park_size + # city specific effect of park size given real park size data
          tmp[j,(first_psi_pd+(city_number-1))] * park_pd_pred_data_list[[city_number]][i]
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

# now add the initial occurrence per random draw of param estimates from the posterior distr
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
park_pd_original_ordered <- unlist(park_pd_original_data_list)
park_pd_pred_ordered <- unlist(park_pd_pred_data_list)

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
  cbind(., park_pd_original_ordered, park_pd_pred_ordered)

## --------------------------------------------------
## Draw multicity plot

# plot the relationships
s <- ggplot(data = new_df, aes(x=park_pd_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.3) +
  #geom_ribbon(aes(
  #ymin=lower_90, 
  #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_pd_original_ordered), max(park_pd_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate \n(City-Specific)") +
  xlab("log(n Plant Genera Detected) / Park Area") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette) + 
  
  scale_fill_manual(values=my_palette) + 
  
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
s



#-------------------------------------------------------------------------------
# occupancy - landscape isolation

# length of number of parks in the cities
pred_length <- nrow(df)

predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 +
    tmp[i,psi_0] + 
      # psi1_ +
      tmp[i,psi_isolation]*scaled_pred
  )
  
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
                                                      prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

#-------------------------------------------------------------------------------

# community plot - park isolation - psi
isolation_df <- as.data.frame(cbind(scaled_pred, criMean[3,], 
                               criMean[1,], criMean[5,],
                               criMean[2,], criMean[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "psi_isolation_community_mean" = "V2",
         "psi_isolation_community_lower95" = "V3",
         "psi_isolation_community_upper95" = "V4",
         "psi_isolation_community_lower50" = "V5",
         "psi_isolation_community_upper50" = "V6")

t <- ggplot(data = isolation_df, aes(scaled_pred, psi_isolation_community_mean)) +
  geom_ribbon(aes(
    ymin=psi_isolation_community_lower50, 
    ymax=psi_isolation_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi_isolation_community_lower95, 
    ymax=psi_isolation_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Park Isolation\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  #scale_fill_manual(values=my_palette) +
  #scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
t

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
first_psi <- which( colnames(tmp)=="psi_city[1]" )
first_psi_isolation <- which( colnames(tmp)=="psi_landscape_isolation[1]" )

for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # community means don't depend on city effects
      predCity[city_number,i,j] <- ilogit( # park isolation trend
        # psi1_0 +
        tmp[j,1] + 
          # psi1_city +
          tmp[j,(first_psi+(city_number-1))] + 
          # psi1_park_isolation +
          tmp[j,(first_psi_isolation+(city_number-1))] * park_isolation_pred_data_list[[city_number]][i]
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

u <- ggplot(data = new_df, aes(x=park_isolation_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.3) +
  #geom_ribbon(aes(
    #ymin=lower_90, 
    #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
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
                                        
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
u

#-------------------------------------------------------------------------------
# occupancy - landscape grassherb

# length of number of parks in the cities
pred_length <- nrow(df)

predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 +
    tmp[i,psi_0] + 
      # psi1_ +
      tmp[i,psi_landscape_grassherb]*scaled_pred
  )
  
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
                                                      prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

#-------------------------------------------------------------------------------

# community plot - landscape grass/herb cover - psi
landscape_grassherb_df <- as.data.frame(cbind(scaled_pred, criMean[3,], 
                                    criMean[1,], criMean[5,],
                                    criMean[2,], criMean[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "psi_landscape_grassherb_community_mean" = "V2",
         "psi_landscape_grassherb_community_lower95" = "V3",
         "psi_landscape_grassherb_community_upper95" = "V4",
         "psi_landscape_grassherb_community_lower50" = "V5",
         "psi_landscape_grassherb_community_upper50" = "V6")

v <- ggplot(data = landscape_grassherb_df, aes(scaled_pred, psi_landscape_grassherb_community_mean)) +
  geom_ribbon(aes(
    ymin=psi_landscape_grassherb_community_lower50, 
    ymax=psi_landscape_grassherb_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi_landscape_grassherb_community_lower95, 
    ymax=psi_landscape_grassherb_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Landscape Herbaceous Cover\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  #scale_fill_manual(values=my_palette) +
  #scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
v

#-------------------------------------------------------------------------------
# for real park sizes in each city

#-------------------------------------------------------------------------------
# get some prediction data

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_landscape_grassherb_pred_data_list <- vector(mode='list', length=n_cities)
park_landscape_grassherb_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

for(city_number in 1:n_cities){
  
  temp <- filter(df, city == city_names[city_number]) %>%
    arrange(., proportion_landscape_grassherb)
  
  # get park size data
  park_landscape_grassherb_pred_data <- temp$proportion_landscape_grassherb_scaled_across_all_cities
  original_scale_park_landscape_grassherb_data <- temp$proportion_landscape_grassherb
  
  park_landscape_grassherb_pred_data_list[[city_number]] <- park_landscape_grassherb_pred_data
  park_landscape_grassherb_original_data_list[[city_number]] <- log(original_scale_park_landscape_grassherb_data+0.01)
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_landscape_grassherb_pred_data)
  
}

max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# predict occupancy

predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) # trends by city 

# get indices for species random effects distributions for particular city
first_psi <- which( colnames(tmp)=="psi_city[1]" )
first_psi_landscape_grassherb <- which( colnames(tmp)=="psi_landscape_grassherb[1]" )

for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # community means don't depend on city effects
      predCity[city_number,i,j] <- ilogit( # park isolation trend
        # psi1_0 +
        tmp[j,1] + 
          # psi1_city +
          tmp[j,(first_psi+(city_number-1))] + 
          # psi1_park_isolation +
          tmp[j,(first_psi_landscape_grassherb+(city_number-1))] * park_landscape_grassherb_pred_data_list[[city_number]][i]
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
park_landscape_grassherb_original_ordered <- unlist(park_landscape_grassherb_original_data_list)
park_landscape_grassherb_pred_ordered <- unlist(park_landscape_grassherb_pred_data_list)

new_df <- as.data.frame(cbind(new_df$city_names, new_df$site, quants))%>%
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
  cbind(., park_landscape_grassherb_original_ordered, park_landscape_grassherb_pred_ordered)

## --------------------------------------------------
## Draw multicity plot

w <- ggplot(data = new_df, aes(x=park_landscape_grassherb_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.3) +
  #geom_ribbon(aes(
  #ymin=lower_90, 
  #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_landscape_grassherb_original_ordered), max(park_landscape_grassherb_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(City-Specific)") +
  xlab("log(Landscape % Grass/Herb Cover)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette) + 
  
  scale_fill_manual(values=my_palette) + 
  
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
w


#-------------------------------------------------------------------------------
# occupancy - landscape woody

# length of number of parks in the cities
pred_length <- nrow(df)

predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 +
    tmp[i,psi_0] + 
      # psi1_ +
      tmp[i,psi_landscape_woody]*scaled_pred
  )
  
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
                                                      prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

#-------------------------------------------------------------------------------

# community plot - landscape grass/herb cover - psi
landscape_woody_df <- as.data.frame(cbind(scaled_pred, criMean[3,], 
                                              criMean[1,], criMean[5,],
                                              criMean[2,], criMean[4,]),
) %>%
  rename("scaled_pred" = "scaled_pred",
         "psi_landscape_woody_community_mean" = "V2",
         "psi_landscape_woody_community_lower95" = "V3",
         "psi_landscape_woody_community_upper95" = "V4",
         "psi_landscape_woody_community_lower50" = "V5",
         "psi_landscape_woody_community_upper50" = "V6")

x <- ggplot(data = landscape_woody_df, aes(scaled_pred, psi_landscape_woody_community_mean)) +
  geom_ribbon(aes(
    ymin=psi_landscape_woody_community_lower50, 
    ymax=psi_landscape_woody_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi_landscape_woody_community_lower95, 
    ymax=psi_landscape_woody_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(scaled_pred), max(scaled_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(Mean Across All Cities)") +
  xlab("Landscape Woody Cover\n(Std. Deviations from Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  #scale_fill_manual(values=my_palette) +
  #scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
x

#-------------------------------------------------------------------------------
# for real park sizes in each city

#-------------------------------------------------------------------------------
# get some prediction data

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
park_landscape_woody_pred_data_list <- vector(mode='list', length=n_cities)
park_landscape_woody_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- vector(length=n_cities)

for(city_number in 1:n_cities){
  
  temp <- filter(df, city == city_names[city_number]) %>%
    arrange(., proportion_landscape_woody)
  
  # get park size data
  park_landscape_woody_pred_data <- temp$proportion_landscape_woody_scaled_across_all_cities
  original_scale_park_landscape_woody_data <- temp$proportion_landscape_woody
  
  park_landscape_woody_pred_data_list[[city_number]] <- park_landscape_woody_pred_data
  park_landscape_woody_original_data_list[[city_number]] <- log(original_scale_park_landscape_woody_data+0.01)
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_landscape_woody_pred_data)
  
}

max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# predict occupancy

predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) # trends by city 

# get indices for species random effects distributions for particular city
first_psi <- which( colnames(tmp)=="psi_city[1]" )
first_psi_landscape_woody <- which( colnames(tmp)=="psi_landscape_woody[1]" )

for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # community means don't depend on city effects
      predCity[city_number,i,j] <- ilogit( # park isolation trend
        # psi1_0 +
        tmp[j,1] + 
          # psi1_city +
          tmp[j,(first_psi+(city_number-1))] + 
          # psi1_park_isolation +
          tmp[j,(first_psi_landscape_woody+(city_number-1))] * park_landscape_woody_pred_data_list[[city_number]][i]
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
park_landscape_woody_original_ordered <- unlist(park_landscape_woody_original_data_list)
park_landscape_woody_pred_ordered <- unlist(park_landscape_woody_pred_data_list)

new_df <- as.data.frame(cbind(new_df$city_names, new_df$site, quants))%>%
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
  cbind(., park_landscape_woody_original_ordered, park_landscape_woody_pred_ordered)

## --------------------------------------------------
## Draw multicity plot

y <- ggplot(data = new_df, aes(x=park_landscape_woody_original_ordered, y=mean, colour=city)) +
  geom_ribbon(aes(
    ymin=lower_50, 
    ymax=upper_50, fill=city), alpha=0.3) +
  #geom_ribbon(aes(
  #ymin=lower_90, 
  #ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_landscape_woody_original_ordered), max(park_landscape_woody_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Occupancy Rate\n(City-Specific)") +
  xlab("log(Landscape % Woody Cover)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=my_palette) + 
  
  scale_fill_manual(values=my_palette) + 
  
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
y

#-------------------------------------------------------------------------------
# cowplot

cowplot::plot_grid(p, q, t, u, v, w, x, y, ncol = 2, rel_widths = c(1, 1.5),
                   labels = c('a)', 'b)', 'c)', 'd)',
                              'e)', 'f)', 'g)', 'h)'),
                   label_size = 20)

