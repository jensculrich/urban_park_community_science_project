# plot effects of site predictors on local species richness

library(tidyverse)
library(rstan)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools.csv")

# list of city names
city_names <- c(
  "LA", # 1
  "NYC" # 2
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
        "./model_outputs/stan_out_", city, "_2km_connectivity_family_100buffers.rds"))
  ), city)
  
  mega_list[[city_number]] <- temp

}

# read city specific prediction data across loop
# i.e. only predict diversity in a city across the range of parks actually observed in the specific city
site_pred_data_list <- vector(mode='list', length=n_cities)
site_original_data_list <- vector(mode='list', length=n_cities)
pred_length <- 30 # something short for testing

for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- readRDS( paste0("./run_model/prepped_data/prepped_data_", city, ".rds"))
   
  pred_data <- seq(from = min(temp$site_data$log_total_green_space_area_scaled), 
                   to = max(temp$site_data$log_total_green_space_area_scaled), 
                   length.out = pred_length)
  
  mean_park_size <- mean(temp$site_data$log_total_green_space_area)
  sd_park_size <- sd(temp$site_data$log_total_green_space_area)
  # now do some algebra to get the scaled data back onto a real life m^2 scale
  original_scale_park_size_data <- (pred_data * sd_park_size) + mean_park_size
  
  site_pred_data_list[[city_number]] <- pred_data
  site_original_data_list[[city_number]] <- original_scale_park_size_data
  
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

n_draws = 50 # small number for testing bc it does take a few minutes to simulate results
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
  first_psi1 <- which( colnames(city_estimates)=="species_intercepts[1,1]" )
  first_gamma <- which( colnames(city_estimates)=="gamma_species[1]" )
  first_phi <- which( colnames(city_estimates)=="phi_species[1]" )
  
  # get the trait data for species from the particular city
  city_name <- city_names[city_number]
  my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".rds"))
  
  pred_data <- site_pred_data_list[[city_number]]
  
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
              city_estimates[rand,4] + 
                # a species specific intercept effect (the number here should be first column)
                city_estimates[rand,(first_psi1+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,5] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,6] * pred_data[j] 
              # not adding any other park predictors would keep park at same value
            )
          
          
          gamma_expected[i,j,k-1] = # gamma[,,k-1] yields gamma for the transition between year 1 and 2
            ilogit(#gamma0 +
              city_estimates[rand,10] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_gamma+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,12] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,13] * pred_data[j] 
              # not adding any other park predictors would keep park at same value
            )
          
        
          phi_expected[i,j,k-1] = # phi[,,k-1] yields phi for the transition between year 1 and 2
            ilogit(#phi0 +
              city_estimates[rand,17] + 
                #species_effects[species[i],1] + // a species specific intercept
                # start at first row of species effects
                # then each next species will be + i
                city_estimates[rand,(first_phi+(i-1))] +
                # effect of wingspan * wingspan of species i + 
                city_estimates[rand,19] * species_data$aveWingspan_scaled[i] + 
                # effect of parksize * parksize of site j + 
                city_estimates[rand,20] * pred_data[j] 
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
  pred_data <- site_pred_data_list[[city_number]]
  original_scale_park_size_data <- site_original_data_list[[city_number]]
  
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
# site predictors

LA_data <- readRDS("./run_model/prepped_data/prepped_data_LA.rds")
NYC_data <- readRDS("./run_model/prepped_data/prepped_data_NYC.rds")
#SEA_data <- readRDS("./run_model/prepped_data/prepped_data_SEA.rds")

LA_site_data <- LA_data$site_data %>%
  mutate(City = "LA")

NYC_site_data <- NYC_data$site_data %>%
  mutate(City = "NYC")

#SEA_site_data <- SEA_data$site_data %>%
#mutate(city = "SEA")

all_site_data <- rbind(
  LA_site_data,
  NYC_site_data#,
  #SEA_site_data
)

# park size
p1 <- ggplot(all_site_data, aes(x = log_total_green_space_area, 
                                colour = City, fill = City)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18)) +
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
  pred_data <- site_pred_data_list[[city_number]]
  pred_data <- pred_data[-1]
  original_scale_park_size_data <- site_original_data_list[[city_number]]
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

cowplot::plot_grid(p, p2, p1, ncol = 1)

##------------------------------------------------------------------------------
# old stuff

# collapse across years (average richness by site [array dimension 4] 
# per rand draw from the posterior [dim 2], across all years [dim 3])
# per city [dim 1],
jaccard_site <- apply(jaccard_site,c(1,2,4),mean) # average across years
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




# collapse across years (average richness by site [array dimension 3] 
# per rand draw from the posterior [dim 1], across all years [dim 2])
occurrence_simmed_avg <- apply(occurrence_simmed,c(1,2,4),mean)

# for year one only
jaccard_site <- array(data = NA, dim=c(pred_length, n_draws))

# jaccard index needs a reference level
ref_site <- 1

# now compute dissimilarity in the occurrence matrix in site i versus site 1
for(k in 1:n_draws){
  for(i in 1:pred_length){ # jaccard index for sites (in terms of shared species)
    jaccard_site[i,k] <- sum(occurrence_simmed_avg[,ref_site,k]*occurrence_simmed_avg[,i,k]) /
      (sum(occurrence_simmed_avg[,ref_site,k]) +
      sum(occurrence_simmed_avg[,i,k]) - sum(occurrence_simmed_avg[,ref_site,k]*occurrence_simmed_avg[,i,k]))
  }
}

pm <- apply(jaccard_site, 1, mean, na.rm=TRUE)
cri <- apply(jaccard_site, 1, function(x) quantile(x, prob = c(0.05, 0.95)))
cbind(pm, "5%" = cri[1,], "90%" = cri[2,])

df <- as.data.frame(cbind(pred_data, original_scale_park_size_data,
                          pm,
                          cri[1,], cri[2,]))

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
pre_p <- ggplot(df) + 
  stat_smooth(aes(x = pred_data, y = pm), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = V4), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = V5), method = "loess", se = FALSE) 
pre_p

# build plot object for rendering 
p <- ggplot_build(pre_p)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(pred_data = p$data[[1]]$x,
                  mean = p$data[[1]]$y,
                  V4 = p$data[[2]]$y,
                  V5 = p$data[[3]]$y) 

# use the loess data to add the 'ribbon' to plot 
(p  <- ggplot(data = df2, aes(pred_data)) +
    geom_line(aes(y=mean), size = 2, colour = "lightskyblue4") +
    geom_ribbon(aes(ymin = V4, ymax = V5), fill = "lightskyblue3", alpha = 0.3) +
    #geom_ribbon(
    # aes(ymin = lower_95, ymax = upper_95), fill = "#DCBCBC") +
    #geom_ribbon(
    #aes(ymin = lower_50, ymax = upper_50), fill = "#B97C7C") +
    #geom_line(size = 2, colour = "#8F2727") +
    ylim(c(0, 1)) +
    theme_classic() +
    xlab("Park Size (log-transformed and scaled)") +
    ylab("Posterior Predictive Distribution\nof Jaccard Index\n(with respect to smallest site)") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)#,
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_blank(), axis.line = element_line(colour = "black")
    ) 
) 

(hist <- ggplot(site_data, aes(x=log_total_green_space_area_scaled)) +
  geom_histogram(binwidth=.5, colour="black", fill="grey") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size = 20))
)

cowplot::plot_grid(p, hist, ncol=1)

## --------------------------------------------------
## plot species richness by connectivity

## --------------------------------------------------
## get prediction range

n_years = 5
n_years_minus1 = n_years - 1 

n_draws = 50 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors

pred_length = 500

pred_data <- seq(from = min(connectivity), to = max(connectivity), length.out = pred_length)

## --------------------------------------------------

psi1_expected <- array(data = NA, dim=c(n_species, pred_length))

gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years))

random_draws_from_posterior = seq(length.out=n_draws) # use if not using the full posterior

psi <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

richness = array(data = NA, dim=c(pred_length,n_years, n_draws))

# take random draws for psi and predict occurrence, then sum across n species

for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 2:n_years){
        
        psi1_expected[i,j] =
          ilogit(
            # YEAR 1 is the global intercept
            list_of_draws[rand,1] + 
              # a species specific intercept effect (the number here should be first column)
              list_of_draws[rand,(93+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,3] * species_data$aveWingspan_scaled[i] + 
              # effect of connectivity * connectivity of site j + 
              list_of_draws[rand,9] * pred_data[j] 
            # not adding any other park predictors would keep park at same value
          )
        
        
        gamma_expected[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,7] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + i
              list_of_draws[rand,(158+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,9] * species_data$aveWingspan_scaled[i] + 
              # effect of connectivity * connectivity of site j + 
              list_of_draws[rand,11] * pred_data[j] 
            # not adding any other park predictors would keep park at same value
          )
        
        
        phi_expected[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,13] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + i
              list_of_draws[rand,(223+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,15] * species_data$aveWingspan_scaled[i] + 
              # effect of connectivity * connectivity of site j + 
              list_of_draws[rand,17] * pred_data[j] 
            # not adding any other park predictors would keep park at same value
          )
        
        if(k == 2){
          # calculate occurrence based on transition from first year
          psi[i,j,k-1] = psi1_expected[i,j] * phi_expected[i,j,k-1] + 
            (1 - psi1_expected[i,j]) * gamma_expected[i,j,k-1]
          
        } else{
          # calculate occurrence based on transition from first year
          psi[i,j,k-1] = psi[i,j,k-2] * phi_expected[i,j,k-1] + 
            (1 - psi[i,j,k-2]) * gamma_expected[i,j,k-1]
          
        }
      } 
    }
  }
  
  # simmed occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
      for(k in 1:n_years){
        if(k < 2){
          occurrence_simmed[i,j,1] <- rbinom(1, 1, prob = psi1_expected[i,j])
        } else{
          occurrence_simmed[i,j,k] = occurrence_simmed[i,j,k-1] * phi_expected[i,j,k-1] + 
            (1 - occurrence_simmed[i,j,k-1]) * gamma_expected[i,j,k-1]
        }
      }
    }
  }
  
  for(j in 1:pred_length){
    for(k in 1:n_years){
      richness[j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k])
    }
  }
  
} 

## --------------------------------------------------
# summarize the results

# collapse across years (average richness by site [array dimension 3] 
# per rand draw from the posterior [dim 1], across all years [dim 2])
richness <- apply(richness,c(1,3),mean)

mean = vector(length=pred_length)
lower_50 = vector(length=pred_length)
upper_50 = vector(length=pred_length)
lower_95 = vector(length=pred_length)
upper_95 = vector(length=pred_length)

for(j in 1:pred_length){
  quants = as.vector(quantile(richness[j,], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
  
  mean[j] = quants[3]
  lower_50[j] = quants[2]
  upper_50[j] = quants[4]
  lower_95[j] = quants[1]
  upper_95[j] = quants[5]
}

df <- as.data.frame(cbind(pred_data,
                          mean,
                          lower_50, upper_50,
                          lower_95, upper_95))

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
pre_q <- ggplot(df) + 
  stat_smooth(aes(x = pred_data, y = mean), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = lower_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = upper_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = lower_50), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = upper_50), method = "loess", se = FALSE) 
pre_q

# build plot object for rendering 
q <- ggplot_build(pre_q)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(pred_data = q$data[[1]]$x,
                  mean = q$data[[1]]$y,
                  lower_95 = q$data[[2]]$y,
                  upper_95 = q$data[[3]]$y,
                  lower_50 = q$data[[4]]$y,
                  upper_50 = q$data[[5]]$y) 

# use the loess data to add the 'ribbon' to plot 
(q  <- ggplot(data = df2, aes(pred_data)) +
    geom_line(aes(y=mean), size = 2, colour = "lightskyblue4") +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "lightskyblue3", alpha = 0.3) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "lightskyblue2", alpha = 0.3) +
    #geom_ribbon(
    # aes(ymin = lower_95, ymax = upper_95), fill = "#DCBCBC") +
    #geom_ribbon(
    #aes(ymin = lower_50, ymax = upper_50), fill = "#B97C7C") +
    #geom_line(size = 2, colour = "#8F2727") +
    ylim(c(0, 65)) +
    theme_classic() +
    xlab("Connectivity: Mean Distance to\nOther Greenspace within 2km (scaled)") +
    ylab("Predicted Local Species Richness  \n (averaged across years") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)#,
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
  
) 


## --------------------------------------------------
## cowplot

cowplot::plot_grid(p, q)
