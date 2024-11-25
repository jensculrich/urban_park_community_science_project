# plot effects of site predictors on local species richness

library(tidyverse)
library(rstan)

stan_out <- readRDS("./model_outputs/stan_out.rds")
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

list_of_draws <- as.data.frame(stan_out)

city_name <- "los_angeles"
my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))

n_species <- my_data$covariate_data %>%
  group_by(species) %>%
  slice(1) %>%
  nrow(.)

species_data <- my_data$species_info
  
#site_data <- my_data$site_data
park_size <- my_data$covariate_data %>%
  select(site, park_size_scaled) %>%
  group_by(site) %>%
  slice(1) %>%
  ungroup() %>%
  pull(park_size_scaled)

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

## --------------------------------------------------
## get prediction range

n_draws = 100 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors

pred_length = 1000

pred_data <- seq(from = min(park_size), to = max(park_size), length.out = pred_length)

## --------------------------------------------------

psi_expected <- array(data = NA, dim=c(n_species, pred_length))

occurrence_simmed <- array(data = NA, dim=c(n_species, pred_length))

random_draws_from_posterior = seq(length.out=n_draws) # use if not using the full posterior

richness = array(data = NA, dim=c(pred_length, n_draws))

# take random draws for psi and predict occurrence, then sum across n species
for(draw in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence
  for(i in 1:n_species){
    for(j in 1:pred_length){
     
      psi_expected[i,j] =
        ilogit(
          # YEAR 1 is the global intercept
          list_of_draws[rand,2] + 
            # a species specific intercept effect
            list_of_draws[rand,(19+(i-1))] +
            # effect of wingspan * wingspan of species i + 
            list_of_draws[rand,7] * species_data$aveWingspan_scaled[i] + 
            # effect of parksize * wingspan of species i + 
            list_of_draws[rand,8] * pred_data[j] 
        )
      
    }
  }
  
  # simmed occurrence in year 1
  for(i in 1:n_species){
    for(j in 1:pred_length){
        occurrence_simmed[i,j] <- rbinom(1, 1, prob = psi_expected[i,j])
    }
  }
  
  for(j in 1:pred_length){
      richness[j,draw] <- sum(occurrence_simmed[1:n_species,j])
  }
  
} 

## --------------------------------------------------
# summarize the results

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
   geom_line(aes(y=mean), size = 2, colour = "#8F2727") +
   geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#DCBCBC", alpha = 0.3) +
   geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "#B97C7C", alpha = 0.3) +
    #geom_ribbon(
    # aes(ymin = lower_95, ymax = upper_95), fill = "#DCBCBC") +
   #geom_ribbon(
     #aes(ymin = lower_50, ymax = upper_50), fill = "#B97C7C") +
   #geom_line(size = 2, colour = "#8F2727") +
   ylim(c(0, 60)) +
   theme_classic() +
   xlab("Park size (log-transformed and then scaled)") +
   ylab("Predicted species richness") +
   theme(axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size = 20)#,
         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         #panel.background = element_blank(), axis.line = element_line(colour = "black")
         )
 
) 
