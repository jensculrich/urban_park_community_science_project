# plot effects of predictors on occurrence and detection

# load libraries
library(tidyverse)
library(rstan)

region <- "northeast"

# list of city names
city_names <- c(
  "Boston", 
  "DC",
  "NYC", 
  "Philadelphia"
)

# number of cities
n_cities <- length(city_names)

## get param estimates from the region
stan_out <- readRDS(paste0(
  "./model_outputs/stan_out_", region, "_2km_isolation_0buffers_simple3.rds"))
tmp <- as.data.frame(stan_out) # take estimates from each HMC step as a df
#n_samp <- 10 # how many samples do we have from the HMC run?
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?

# handy for viewing column numbers
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

## get data from region
df <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))$site_data

# length of number of parks in the cities
pred_length <- nrow(df)

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

#-------------------------------------------------------------------------------
# get some prediction data

# these are a sequence of scaled values of park sizes (what the model sees)
size_pred <- seq(from = -2, to = 2, length.out = pred_length)

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

initial_occurrence <- vector(length = pred_length)

# get length of stan fit object (HMC iterations * n_chains)  
predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 + # the initial occurrence rate at a completely average park
    tmp[i,1] + 
      # psi1_park_size # the effect of park size on the initial occurrence rate
      tmp[i,6]*size_pred
  )
    
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
              prob = c(0.05, 0.25, 0.5, 0.75, 0.95))) # get 50 and 90% BCIs

#-------------------------------------------------------------------------------

# community plot - park size - psi1
# organize the mean, 50 and 90% BCIs into a data frame
size_df <- as.data.frame(cbind(size_pred, criMean[3,], 
                               criMean[1,], criMean[5,],
                               criMean[2,], criMean[4,]),
                         ) %>%
  rename("size_pred" = "size_pred",
         "psi1_size_community_mean" = "V2",
         "psi1_size_community_lower95" = "V3",
         "psi1_size_community_upper95" = "V4",
         "psi1_size_community_lower50" = "V5",
         "psi1_size_community_upper50" = "V6")

# plot the estimated relationship for the average city in the region
p <- ggplot(data = size_df, aes(size_pred, psi1_size_community_mean)) +
  geom_ribbon(aes(
    ymin=psi1_size_community_lower50, 
    ymax=psi1_size_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi1_size_community_lower95, 
    ymax=psi1_size_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(size_pred), max(size_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Initial Occurrence Rate \n(Regional Mean)") +
  xlab("Park Size (Std. Deviations from Within-City Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
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
  park_size_pred_data <- temp$log_total_green_space_area_scaled_2
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
first_psi1 <- which( colnames(tmp)=="psi1_city[1]" )
first_psi1_size <- which( colnames(tmp)=="psi1_park_size[1]" )

# loop across all sites in each city, and make a prediction about initial occurrence
# for each of the samples from the posterior probability distribution
for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # estimate city specific predictions about initial occurrence rate
      predCity[city_number,i,j] <- ilogit( 
        # psi1_0 + # global intercept
        tmp[j,1] + 
          # psi1_city + # city effect on the intercept
          tmp[j,(first_psi1+(city_number-1))] + 
          # psi1_park_size + # city specific effect of park size given real park size data
          tmp[j,(first_psi1_size+(city_number-1))] * park_size_pred_data_list[[city_number]][i]
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
  #geom_ribbon(aes(
    #ymin=lower_50, 
    #ymax=upper_50, fill=city), alpha=0.8) +
  geom_ribbon(aes(
    ymin=lower_90, 
    ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_size_original_ordered), max(park_size_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Initial Occurrence Rate \n(City-Specific)") +
  xlab("log(Park Size m^2)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=c("#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                        "#1a5acd", "#E69F90", "#FFFF00")) + 
                                          
  scale_fill_manual(values=c("#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                     "#1a5acd", "#E69F90", "#FFFF00")) + 
                                       
  theme(#legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
q

#-------------------------------------------------------------------------------
# cowplot

cowplot::plot_grid(p, q, ncol = 2, rel_widths = c(1, 1.5))

#-------------------------------------------------------------------------------
# now do it for isolation

#-------------------------------------------------------------------------------
# get some prediction data

pred_length <- nrow(df)
isolation_pred <- seq(from = -2, to = 2, length.out = pred_length)

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

initial_occurrence <- vector(length = pred_length)

predMean <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predCity <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predMean[,i] <- ilogit( # park size trend
    # psi1_0 +
    tmp[i,1] + 
      # psi1_ +
      tmp[i,8]*isolation_pred
  )
  
}

# posterior means by community average 
criMean <- apply(predMean, c(1), function(x) quantile(x, 
                                                      prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

#-------------------------------------------------------------------------------

# community plot - park size - psi1
isolation_df <- as.data.frame(cbind(isolation_pred, criMean[3,], 
                               criMean[1,], criMean[5,],
                               criMean[2,], criMean[4,]),
) %>%
  rename("isolation_pred" = "isolation_pred",
         "psi1_isolation_community_mean" = "V2",
         "psi1_isolation_community_lower95" = "V3",
         "psi1_isolation_community_upper95" = "V4",
         "psi1_isolation_community_lower50" = "V5",
         "psi1_isolation_community_upper50" = "V6")

r <- ggplot(data = isolation_df, aes(isolation_pred, psi1_isolation_community_mean)) +
  geom_ribbon(aes(
    ymin=psi1_isolation_community_lower50, 
    ymax=psi1_isolation_community_upper50), alpha=0.8) +
  geom_ribbon(aes(
    ymin=psi1_isolation_community_lower95, 
    ymax=psi1_isolation_community_upper95), alpha=0.4) +
  geom_line(size=2, lty=1) +
  xlim(c(min(isolation_pred), max(isolation_pred))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Initial Occurrence Rate \n(Regional Mean)") +
  xlab("Park Isolation (Std. Deviations from Within-City Mean)") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  #scale_fill_manual(values=my_palette) +
  #scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
r

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
  park_isolation_pred_data <- temp$isolation_scaled_2
  original_scale_park_isolation_data <- temp$isolation
  
  park_isolation_pred_data_list[[city_number]] <- park_isolation_pred_data
  park_isolation_original_data_list[[city_number]] <- original_scale_park_isolation_data
  
  # and figure out how many sites in the city
  pred_length[city_number] <- length(park_isolation_pred_data)
  
}

max_pred_length = max(pred_length)

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

predCity <- array(NA, dim=c(n_cities, max_pred_length, n_samp)) # trends by city 

# get indices for species random effects distributions for particular city
first_psi1 <- which( colnames(tmp)=="psi1_city[1]" )
first_psi1_isolation <- which( colnames(tmp)=="psi1_isolation[1]" )

for(city_number in 1:n_cities){
  for(i in 1:max_pred_length){
    for(j in 1:n_samp){
      
      # community means don't depend on city effects
      predCity[city_number,i,j] <- ilogit( # park isolation trend
        # psi1_0 +
        tmp[j,1] + 
          # psi1_city +
          tmp[j,(first_psi1+(city_number-1))] + 
          # psi1_park_isolation +
          tmp[j,(first_psi1_isolation+(city_number-1))] * park_isolation_pred_data_list[[city_number]][i]
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

s <- ggplot(data = new_df, aes(x=park_isolation_original_ordered, y=mean, colour=city)) +
  #geom_ribbon(aes(
  #ymin=lower_50, 
  #ymax=upper_50, fill=city), alpha=0.8) +
  geom_ribbon(aes(
    ymin=lower_90, 
    ymax=upper_90, fill=city), alpha=0.2) +
  geom_line(size=3, lty=1) +
  xlim(c(min(park_isolation_original_ordered), max(park_isolation_original_ordered))) +
  ylim(c(0, 1)) +
  theme_bw() +
  ylab("Initial Occurrence Rate \n(City-Specific)") +
  xlab("Park Isolation") +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent 
  ) +
  scale_color_manual(values=c("#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                       "#1a5acd", "#E69F90", "#FFFF00")) + 
                                         
  scale_fill_manual(values=c("#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                      "#1a5acd", "#E69F90", "#FFFF00")) + 
                                        
  theme(#legend.position = "none",
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
s

#-------------------------------------------------------------------------------
# cowplot

cowplot::plot_grid(p, q, r, s, ncol = 2, rel_widths = c(1, 1.5))
