# how many detections would you expect per species (W)?
# we will compare this to how many detections per species simulated in the 
# generated quantities block of our model (W_rep) for a visual PPC

library(rstan)
library(tidyverse)

#city_name <- "los_angeles"
#my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
# detection data
V <- my_data$V_detections # detections (1==detected)
V_NA <- my_data$V_NA # NAs (0==no known survey effort made)
R <- nrow(V)

# covariate data
covariate_data <- my_data$covariate_data
n_sites <- my_data$n_sites # number of sites
site <- covariate_data$new_id
n_species <- my_data$n_species # number of species
species <- covariate_data$species_number
n_years <- my_data$n_years # number of surveys
year <- covariate_data$year
n_surveys <- my_data$n_surveys
survey <- sequence(n_surveys)
survey <- (survey - mean(survey)) / sd(survey)
reverse_index <- covariate_data$reverse_index

# categorical year dummy variables
#X_year <- model.matrix(~ as.factor(year), data = covariate_data)

## predictors
# species
feature_diversity <- covariate_data$featureDiversity_scaled
ease_of_id <- covariate_data$research_grade_proportion_scaled
wingspan <- covariate_data$aveWingspan_scaled
# site
park_size <- covariate_data$log_classified_area_scaled
# add more site predictors here

species_info <- my_data$species_info 
site_data <- my_data$site_data
## --------------------------------------------------
### Calculate number of detections in the data (by species)

W_species <- as.data.frame(V) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  cbind(covariate_data$species) %>%
  rename("species" = "covariate_data$species") %>%
  group_by(species) %>%
  mutate(total_detections = sum(sum)) %>%
  slice(1) %>%
  select(species, total_detections) %>%
  ungroup()

## --------------------------------------------------
### Calculate number of detections predicted by the model (by species)

stan_out <- readRDS("./model_outputs/stan_out_new.rds")  
  
fit_summary <- rstan::summary(stan_out)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

list_of_draws <- as.data.frame(stan_out)
  
## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))  
  
## --------------------------------------------------

n_draws = 100 # small number for testing bc it does take a few minutes to simulate results
#n_draws = nrow(list_of_draws) # number of samples from the posteriors
random_draws_from_posterior = seq(length.out=n_draws) # use if not using the full posterior

psi_expected <- vector(length=R)
simmed_occurrence <- vector(length=R)

p_expected <- array(data = NA, dim=c(R, n_surveys))
simmed_detections <- array(data = NA, dim=c(R, n_surveys))

#total_detections_per_draw <- array(data = NA, dim=c(R, n_surveys))

for(i in 1:n_draws){
  
  rand <- random_draws_from_posterior[draw]
  
  # expected occurrence
  for(i in 1:R){
      
      species <- covariate_data$species_number[i] # which species random effects to add
      
      psi_expected[i] =
        ilogit(
          # YEAR 1 is the global intercept
          list_of_draws[rand,2] + 
          # a year effect
          list_of_draws[rand,3] * X_year[i,2] + 
          list_of_draws[rand,4] * X_year[i,3] + 
          list_of_draws[rand,5] * X_year[i,4] + 
          list_of_draws[rand,6] * X_year[i,5] + 
          # a species specific intercept effect
          list_of_draws[rand,(19+(species-1))] +
          # effect of wingspan * wingspan of species i + 
          list_of_draws[rand,7] * covariate_data$aveWingspan_scaled[i] + 
          # effect of parksize * wingspan of species i + 
          list_of_draws[rand,8] * covariate_data$park_size_scaled[i] 
        )
      
      simmed_occurrence[i] <- rbinom(1, 1, prob = psi_expected[i])
      
      for(j in 1:n_surveys){
        p_expected =
          ilogit(
            # YEAR 1 is the global intercept
            list_of_draws[rand,9] + 
            # a species specific intercept effect
            list_of_draws[rand,(79+(species-1))] +
            # effect of wingspan * wingspan of species i + 
            list_of_draws[rand,11] * covariate_data$aveWingspan_scaled[i] + 
            # effect of feature diversity * feature diversity of species i + 
            list_of_draws[rand,12] * covariate_data$featureDiversity_scaled[i] +
            # effect of ease of id * ease of id of species i + 
            list_of_draws[rand,13] * covariate_data$research_grade_proportion_scaled[i] +
            # effect of parksize * wingspan of species i + 
            list_of_draws[rand,14] * covariate_data$park_size_scaled[i] +
            # a species specific effect of survey month
            list_of_draws[rand,(79+(species-1))] * survey[j] +
            # a species specific effect of survey squared
            list_of_draws[rand,(79+(species-1))] * survey[j]^2
          )
      }
  }

  
}






## --------------------------------------------------
# old stuff
for(i in 1:n_draws){
  W_species[i] = sum(V[i,])
}

# for simulated data
W_df <- as.data.frame(cbind(species, W_species)) %>%
  mutate(W_species = as.numeric(W_species))

# for real data
W_df <- as.data.frame(cbind(species_names, W_species)) %>%
  mutate(W_species = as.numeric(W_species))

# get W distributions from model
#stan_out <- readRDS("./model_outputs/stan_out4.rds")
fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest


## --------------------------------------------------
# 

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

start = 1 # which species to start at (hard to see them all at once)
# start at 1, 37, and 73 is pretty good for visualization
n = 50 # how many species to plot (36 is a good number to look at the species in 3 slices)

stan_fit_first_W <- 14 # this changes depending on how many params you tracked

df_estimates <- data.frame(X = numeric(), 
                           Y = numeric(), 
                           lower_95 = numeric(),
                           upper_95 = numeric(),
                           lower_50 = numeric(),
                           upper_50 = numeric()
) 

for(i in 1:n){
  
  row <- c((i + start - 1), 
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),1],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),4],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),8],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),5],
           fit_summary$summary[(stan_fit_first_W+(start+i-2)),7])
  
  df_estimates[i,] <- row
  
}

labels=as.vector(c(W_df[start:(start + n - 1),1]))
ylims = c(0,(max(df_estimates$upper_95)+5))
end_point  = 0.5 + nrow(df_estimates) + nrow(df_estimates) - 1 #

par(mar = c(9,4,1,2))
plot(1, type="n",
     xlim=c(start, n + start - 1 + 0.5), 
     xlab="",
     xaxt = "n",
     ylim=ylims, 
     ylab="50% and 95% Marginal Posterior Quantiles",
     main="Real Detections vs. Model Expectations of Detections")

#axis(1, at=start:(start+n), labels=labels, las = 2, cex.axis=.75)
text(seq(start, start + n - 1, by = 1), par("usr")[3]-0.25, 
     srt = 60, adj = 1, xpd = TRUE,
     labels = labels, cex = 1)

for(i in 1:n){
  sliced <- df_estimates[i,]
  W_sliced <- W_df[i+start - 1, 5]
  
  rect(xleft = (sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_95, ybottom = sliced$upper_95,
       col = c_mid, border = NA
  )
  
  rect(xleft =(sliced$X-0.35), xright=(sliced$X+0.35), 
       ytop = sliced$lower_50, ybottom = sliced$upper_50,
       col = c_mid_highlight, border = NA
  )
  
  points(x=sliced$X, y=W_sliced, pch=1)
}

