# plot effects of site predictors on local species richness

library(tidyverse)
library(rstan)

stan_out <- readRDS("./model_outputs/stan_out_LA_2km_connectivity_family.rds")
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

list_of_draws <- as.data.frame(stan_out)

city_name <- "LA"
my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".rds"))

species_data <- my_data$species_info

n_species <- nrow(species_data)

site_data <- my_data$site_data

park_size <- site_data$log_total_green_space_area_scaled 

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## plot species richness by park size

## --------------------------------------------------
## get prediction range

n_years = 5
n_years_minus1 = n_years - 1 

n_draws = 100 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors

pred_length = 50

pred_data <- seq(from = min(park_size), to = max(park_size), length.out = pred_length)
#center_scale <- function(x) {
#  (x - mean(x)) / sd(x)
#}
mean_park_size <- mean(site_data$log_total_green_space_area)
sd_park_size <- sd(site_data$log_total_green_space_area)
# now do some algebra to get the scaled data back onto a real life m^2 scale
original_scale_data <- (pred_data * sd_park_size) + mean_park_size

## --------------------------------------------------

psi1_expected <- array(data = NA, dim=c(n_species, pred_length))

gamma_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

phi_expected <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length, n_years, n_draws))

random_draws_from_posterior = seq(length.out=n_draws) # use if not using the full posterior

psi <- array(data = NA, dim=c(n_species, pred_length, n_years_minus1))

richness <- array(data = NA, dim=c(pred_length,n_years, n_draws))

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
              list_of_draws[rand,(80+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,3] * species_data$aveWingspan_scaled[i] + 
              # effect of parksize * parksize of site j + 
              list_of_draws[rand,4] * pred_data[j] 
            # not adding any other park predictors would keep park at same value
          )
        
        
        gamma_expected[i,j,k-1] =
          ilogit(#gamma0 +
            list_of_draws[rand,7] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + i
              list_of_draws[rand,(132+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,9] * species_data$aveWingspan_scaled[i] + 
              # effect of parksize * parksize of site j + 
              list_of_draws[rand,10] * pred_data[j] 
            # not adding any other park predictors would keep park at same value
          )
        
      
        phi_expected[i,j,k-1] =
          ilogit(#phi0 +
            list_of_draws[rand,13] + 
              #species_effects[species[i],1] + // a species specific intercept
              # start at first row of species effects
              # then each next species will be + i
              list_of_draws[rand,(183+(i-1))] +
              # effect of wingspan * wingspan of species i + 
              list_of_draws[rand,15] * species_data$aveWingspan_scaled[i] + 
              # effect of parksize * parksize of site j + 
              list_of_draws[rand,16] * pred_data[j] 
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
          occurrence_simmed[i,j,1,rand] <- rbinom(1, 1, prob = psi1_expected[i,j])
        } else{
          occurrence_simmed[i,j,k,rand] = occurrence_simmed[i,j,k-1,rand] * phi_expected[i,j,k-1] + 
            (1 - occurrence_simmed[i,j,k-1,rand]) * gamma_expected[i,j,k-1]
        }
      }
    }
  }
  
  for(j in 1:pred_length){
    for(k in 1:n_years){
      richness[j,k,draw] <- sum(occurrence_simmed[1:n_species,j,k,rand])
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

df <- as.data.frame(cbind(pred_data, original_scale_data,
                           mean,
                           lower_50, upper_50,
                           lower_95, upper_95))

## --------------------------------------------------
## Draw species richness plot

# create plot object with loess regression lines
pre_p <- ggplot(df) + 
  stat_smooth(aes(x = pred_data, y = mean), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = lower_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = upper_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = lower_50), method = "loess", se = FALSE) +
  stat_smooth(aes(x = pred_data, y = upper_50), method = "loess", se = FALSE) 
pre_p

# build plot object for rendering 
p <- ggplot_build(pre_p)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(pred_data = p$data[[1]]$x,
                  mean = p$data[[1]]$y,
                  lower_95 = p$data[[2]]$y,
                  upper_95 = p$data[[3]]$y,
                  lower_50 = p$data[[4]]$y,
                  upper_50 = p$data[[5]]$y) 

# use the loess data to add the 'ribbon' to plot 
(p  <- ggplot(data = df2, aes(pred_data)) +
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
   xlab("Park Size (log-transformed and scaled)") +
   ylab("Posterior Predictive Distribution\nof Species Richness\n(averaged across years)") +
   theme(axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size = 20)#,
         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         #panel.background = element_blank(), axis.line = element_line(colour = "black")
         )
 
) 

## now do it on the original data scale
# create plot object with loess regression lines
pre_p <- ggplot(df) + 
  stat_smooth(aes(x = original_scale_data, y = mean), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_data, y = lower_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_data, y = upper_95), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_data, y = lower_50), method = "loess", se = FALSE) +
  stat_smooth(aes(x = original_scale_data, y = upper_50), method = "loess", se = FALSE) 
pre_p

# build plot object for rendering 
p <- ggplot_build(pre_p)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(original_scale_data = p$data[[1]]$x,
                  mean = p$data[[1]]$y,
                  lower_95 = p$data[[2]]$y,
                  upper_95 = p$data[[3]]$y,
                  lower_50 = p$data[[4]]$y,
                  upper_50 = p$data[[5]]$y) 

# use the loess data to add the 'ribbon' to plot 
(p  <- ggplot(data = df2, aes(original_scale_data)) +
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
    xlab("log(Park Size m^2)") +
    ylab("Posterior Predictive Distribution\nof Species Richness\n(averaged across years)") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size = 20)#,
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
  
) 

## --------------------------------------------------
## analyze beta diversity of simmed communities

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

df <- as.data.frame(cbind(pred_data, original_scale_data,
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
