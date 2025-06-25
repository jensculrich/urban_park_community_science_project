library(rstan) # to run analysis

## --------------------------------------------------
### Define simulation conditions

# choose sample sizes and 
n_species <- 50 # number of species
n_sites <- 50 # number of sites (must be an even number for simulation code)
n_years <- 5 # number of years
n_years_minus1 <- n_years - 1
n_surveys <- 12 # number of surveys per year

# set parameter values
psi1_0 <- 0 # prob of initial occupancy
sigma_psi1_species <- 2 # prob of initial occupancy
psi1_wingspan <- 1 # effect of wingspan on species intercept

gamma0 <- -2 # prob of initial occupancy
sigma_gamma_species <- 1 # prob of initial occupancy

phi0 <- 2 # prob of initial occupancy
sigma_phi_species <- 1 # prob of initial occupancy

p0 <- -1 # probability of detection (logit scaled)
sigma_p_species <- 1 # species-specific variation
mu_p_species_date <- 0
sigma_p_species_date <- 1
mu_p_species_date_sq <- -0.5  
sigma_p_species_date_sq <- 1

# simulate missing data
# for STAN will also need to make an NA indicator array
create_missing_data <- TRUE # create holes in the data? (MAR)
prob_missing <- 0.2 # if so, what proportion of data missing?

## --------------------------------------------------
### Define simulation function

simulate_data <- function(
    n_species, n_sites, n_years, n_years_minus1, n_surveys,
    psi1_0,
    sigma_psi1_species,
    psi1_wingspan,
    
    gamma0,
    sigma_gamma_species,
    
    phi0,
    sigma_phi_species,
    
    p0, # probability of detection (logit scaled)
    sigma_p_species, # species-specific variation
    mu_p_species_date,
    sigma_p_species_date,
    mu_p_species_date_sq,  
    sigma_p_species_date_sq,
    create_missing_data,
    prob_missing
){
  
  ## ilogit and logit functions
  ilogit <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) log(x/(1-x))
  
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  # prepare arrays for z and y
  z <- array(NA, dim = c(n_species, n_sites, n_years)) # latent presence/absence
  y <- array(NA, dim = c(n_species, n_sites, n_years, n_surveys)) # observed data
  
  ## --------------------------------------------------
  ## Create covariate data
  
  ## --------------------------------------------------
  ## species wingspan
  
  # in my real data most species tend to have pretty low specialization (skewed left)
  # the degree and d' are possibly positively correlated
  
  wingspan_scaled <- rnorm(n_species, mean=0, sd = 1)
  
  ## --------------------------------------------------
  ## site covariates
  
  ## --------------------------------------------------
  ## survey number (time of the year)
  surveys <- seq(1:n_surveys)
  
  ## let's scale the calendar day of year by the mean date (z-score scaled)
  surveys_scaled <- center_scale(surveys) 
  
  ## --------------------------------------------------
  ## Create random effects
  
  ## species-specific random intercepts
  psi1_species_mean <- vector(length = n_species)
  psi1_species <- vector(length = n_species)
  for(i in 1:n_species){
    psi1_species_mean[i] <-psi1_0 + psi1_wingspan*wingspan_scaled[i]
    psi1_species[i] <- rnorm(n=1, mean=psi1_species_mean[i], sd=sigma_psi1_species)
  }
  
  gamma_species <- rnorm(n=n_species, mean=0, sd=sigma_gamma_species)
  
  phi_species <- rnorm(n=n_species, mean=0, sd=sigma_phi_species)
  
  p_species <- rnorm(n=n_species, mean=0, sd=sigma_p_species)
  
  p_species_date <- rnorm(n=n_species, mean=mu_p_species_date, sd=sigma_p_species_date)
  
  p_species_date_sq <- rnorm(n=n_species, mean=mu_p_species_date_sq, sd=sigma_p_species_date_sq)
  
  ## --------------------------------------------------
  ## Create expected values
  
  # generate p with heterogeneity
  logit_p <- array(NA, dim = c(n_species, n_sites, n_years, n_surveys)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_surveys){
          
          logit_p[i,j,k,l] = 
            p0 +
            p_species[i] +
            p_species_date[i]*surveys_scaled[l] + # a spatiotemporally specific intercept
            p_species_date_sq[i]*(surveys_scaled[l])^2 # a spatiotemporally specific intercept
          
        }
      }
    }  
  }
  
  
  # generate ecological expected values with heterogeneity
  logit_psi1 <- array(NA, dim = c(n_species, n_sites)) 
  logit_gamma <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  logit_phi <- array(NA, dim = c(n_species, n_sites, n_years_minus1)) 
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years_minus1){
        
        logit_psi1[i,j] = 
          psi1_0 +
          psi1_species[i] 
        
        logit_gamma[i,j,k] = # gamma for transition (starting for between years 1 and 2)
          gamma0 +
          gamma_species[i] 
        
        logit_phi[i,j,k] = 
          phi0 +
          phi_species[i]
        
      }
    }
  }
  
  
  
  # generate initial presence/absence states
  for(i in 1:n_species){
    for(j in 1:n_sites){
      z[i,j,1] <- rbinom(n=1, size=1, prob=ilogit(logit_psi1[i,j])) 
    }
  }
  
  true_occupancy <- apply(z[,,],c(1,3),sum ) / n_sites # true occupancy proportion in year 1     
  
  # generate presence/absence in subsequent years
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 2:n_years){
        
        # use z as a switch so we are estimating 
        exp_z <- z[i,j,k-1] * ilogit(logit_phi[i,j,k-1]) + # survival if z=1
          (1 - z[i,j,k-1]) * ilogit(logit_gamma[i,j,k-1]) # or colonization if z=0
        
        # and then assign z stochastically
        # some sites may transition if they are colonized or local extinction occurs
        # but might otherwise retain their state across years
        # look at year starting with 2, for the first transition (phi or gamma[,,1])
        z[i,j,k] <- rbinom(n = 1, size = 1, prob = exp_z) 
        
      }
    }    
  }
  
  # detection / non-detection data
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        for(l in 1:n_surveys){
          y[i,j,k,l] <- rbinom(n = 1, size = 1, prob = z[i,j,k]*ilogit(logit_p[i,j,k,l]))
        }
      }
    }
  }
  
  y ; str(y)
  
  sum((y/n_surveys) / sum(z)) # proportion of times detection given presence
  
  # need to fix this if going to use
  if(create_missing_data == TRUE){
    # generate missing values: create simple version of unbalanced data set
    prob_missing <- prob_missing # constant NA probability
    y2 <- y # duplicate balanced dataset
    for(i in 1:n_species){
      for(j in 1:n_sites){
        for(k in 1:n_years){
          for(l in 1:n_surveys){
            turnNA <- rbinom(1,1,prob_missing)
            y2[i,j,k,l] <- ifelse(turnNA==1, NA, y2[i,j,k,l])
          }
        }
      }
    }
    y2 ; str(y2)
  } else {
    y2 <- y # duplicate dataset with no missing data
  }
  
  # create an NA indicator array
  y_NA <- y2
  # survey occurred (either 0 or 1 detection) == 1
  y_NA[y_NA == 0] <- 1
  # NA == 0
  y_NA[is.na(y_NA)] <- 0
  
  # which species never occurred
  species_not_occurring <- length(which(apply(true_occupancy, 1, sum) == 0))
  
  # which species were never detected (even though they occurred)
  species_not_observed <- length(which(apply(y2, 1, sum, na.rm=TRUE) == 0))
  
  # now turn the NAs in the detection data to zeroes
  # stan can't handle NAs but the numeric value here doesn't matter
  # because the NA indicator array will tell stan to ignore these detection events
  y2[is.na(y2)] <- 0
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    V = y2, # return detection data after potentially introducing NAs,
    V_NA = y_NA,
    
    wingspan = wingspan_scaled,
    surveys = surveys_scaled
  ))
  
} # end function


## --------------------------------------------------
### Simulate some data

#set.seed(1)
my_simulated_data <- simulate_data(  
  n_species, n_sites, n_years, n_years_minus1, n_surveys,
  
  psi1_0,
  sigma_psi1_species,
  psi1_wingspan,
  
  gamma0,
  sigma_gamma_species,
  
  phi0,
  sigma_phi_species,
  
  p0, # probability of detection (logit scaled)
  sigma_p_species, # species-specific variation
  mu_p_species_date,
  sigma_p_species_date,
  mu_p_species_date_sq,  
  sigma_p_species_date_sq,
  create_missing_data,
  prob_missing
)

## get data
# study dims
V <- my_simulated_data$V
V_NA <- my_simulated_data$V_NA
surveys <- my_simulated_data$surveys
species <- seq(1, n_species, by=1)
sites <- seq(1, n_sites, by=1)
years <- seq(1, n_years_minus1, by=1)
# covs
wingspan <- my_simulated_data$wingspan

## --------------------------------------------------
### Prep data and tweak model options

stan_data <- c("V", "V_NA", "species", "sites", "years", "surveys", 
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys",
               "wingspan"
) 

## Parameters monitored 
params <- c("psi1_0", 
            "sigma_psi1_species",
            "psi1_wingspan",
            
            "gamma0", 
            "sigma_gamma_species",
            
            "phi0", 
            "sigma_phi_species",
            
            "p0", 
            "sigma_p_species",
            "mu_p_species_date",
            "sigma_p_species_date",
            "mu_p_species_date_sq",
            "sigma_p_species_date_sq",
            
            "W_species_rep",
            "psi1_species", "gamma_species", "phi_species", "p_species"
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0, 1),
       gamma0 = runif(1, -2, 0),
       sigma_gamma_species = runif(1, 0, 1),
       phi0 = runif(1, 2, 3),
       sigma_phi_species = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 0),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)

# targets
parameter_values <-  c(
  psi1_0, 
  sigma_psi1_species,
  psi1_wingspan,
  
  gamma0, 
  sigma_gamma_species,
  
  phi0, 
  sigma_phi_species,
  
  p0, 
  sigma_p_species,
  mu_p_species_date, sigma_p_species_date, 
  mu_p_species_date_sq, sigma_p_species_date_sq, 
  NA, NA, NA, NA, NA # 5 generated quantities to track
  
)

targets <- as.data.frame(cbind(params, parameter_values))

## --------------------------------------------------
### Run model

stan_model <- "./models/dynamic_occupancy_model_simtest.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 control=list(adapt_delta=delta),
                 open_progress = FALSE,
                 cores = n_cores)

saveRDS(stan_out, paste0("./model_outputs/stan_out_", city, "_2km_connectivity_family.rds"))

print(stan_out, digits = 3, 
      pars = c("psi1_0", 
               "sigma_psi1_species",
               
               "gamma0", 
               "sigma_gamma_species",
               
               "phi0", 
               "sigma_phi_species",
               
               "p0", 
               "sigma_p_species",
               "mu_p_species_date",
               "sigma_p_species_date",
               "mu_p_species_date_sq",
               "sigma_p_species_date_sq"
      ))


## --------------------------------------------------
### Plot parameter estimates and targets
library(ggplot2)
library(tidyverse)

targets2 <- targets[1:29,]

fit_summary <- rstan::summary(stan_out_sim)
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

X <- as.factor(seq(1:nrow(targets2)))

estimates_lower <- c(
  fit_summary$summary[1,4], # psi1_0
  fit_summary$summary[2,4], # sigma_psi1_species
  fit_summary$summary[3,4], # psi1_herbaceous_flowers
  fit_summary$summary[4,4], # psi1_woody_flowers
  fit_summary$summary[5,4], # psi1_specialization
  fit_summary$summary[6,4], # psi1_interaction_1
  fit_summary$summary[7,4], # psi1_interaction_2
  fit_summary$summary[8,4], # gamma0
  fit_summary$summary[9,4], # sigma_gamma_species
  fit_summary$summary[10,4], # gamma_herbaceous_flowers
  fit_summary$summary[11,4], # gamma_woody_flowers
  fit_summary$summary[12,4], # gamma_specialization
  fit_summary$summary[13,4], # gamma_interaction_1
  fit_summary$summary[14,4], # gamma_interaction_2
  fit_summary$summary[15,4], # phi0
  fit_summary$summary[16,4], # sigma_phi_species
  fit_summary$summary[17,4], # phi_herbaceous_flowers
  fit_summary$summary[18,4], # phi_woody_flowers
  fit_summary$summary[19,4], # phi_specialization
  fit_summary$summary[20,4], # phi_interaction_1
  fit_summary$summary[21,4], # phi_interaction_2
  fit_summary$summary[22,4], # p0
  fit_summary$summary[23,4], # sigma_p_species
  fit_summary$summary[24,4], # p_specialization
  fit_summary$summary[25,4], # mu_p_species_date
  fit_summary$summary[26,4], # sigma_p_species_date
  fit_summary$summary[27,4], # mu_p_species_date_sq
  fit_summary$summary[28,4], # sigma_p_species_date_sq
  fit_summary$summary[29,4] # p_flower_abundance_any
)

estimates_upper <- c(
  fit_summary$summary[1,8], # psi1_0
  fit_summary$summary[2,8], # sigma_psi1_species
  fit_summary$summary[3,8], # psi1_herbaceous_flowers
  fit_summary$summary[4,8], # psi1_woody_flowers
  fit_summary$summary[5,8], # psi1_specialization
  fit_summary$summary[6,8], # psi1_interaction_1
  fit_summary$summary[7,8], # psi1_interaction_2
  fit_summary$summary[8,8], # gamma0
  fit_summary$summary[9,8], # sigma_gamma_species
  fit_summary$summary[10,8], # gamma_herbaceous_flowers
  fit_summary$summary[11,8], # gamma_woody_flowers
  fit_summary$summary[12,8], # gamma_specialization
  fit_summary$summary[13,8], # gamma_interaction_1
  fit_summary$summary[14,8], # gamma_interaction_2
  fit_summary$summary[15,8], # phi0
  fit_summary$summary[16,8], # sigma_phi_species
  fit_summary$summary[17,8], # phi_herbaceous_flowers
  fit_summary$summary[18,8], # phi_woody_flowers
  fit_summary$summary[19,8], # phi_specialization
  fit_summary$summary[20,8], # phi_interaction_1
  fit_summary$summary[21,8], # phi_interaction_2
  fit_summary$summary[22,8], # p0
  fit_summary$summary[23,8], # sigma_p_species
  fit_summary$summary[24,8], # p_specialization
  fit_summary$summary[25,8], # mu_p_species_date
  fit_summary$summary[26,8], # sigma_p_species_date
  fit_summary$summary[27,8], # mu_p_species_date_sq
  fit_summary$summary[28,8], # sigma_p_species_date_sq
  fit_summary$summary[29,8] # p_flower_abundance_any
)

df_estimates <- as.data.frame(cbind(X, targets2, estimates_lower, estimates_upper))
df_estimates$parameter_value <- as.numeric(df_estimates$parameter_value)

(p <- ggplot(df_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 
                                         9, 10, 11, 12, 13, 14, 15, 16, 
                                         17, 18, 19, 20, 21, 22, 23, 
                                         24, 25, 26, 27, 28, 29
    ),
    
    labels=c(bquote(psi[0]), 
             bquote(sigma[psi1["species"]]),
             bquote(psi["herb."]),
             bquote(psi["woody"]), 
             bquote(psi["specialization"]),
             bquote(psi["spec.*herb."]),
             bquote(psi["spec.*woody"]),
             
             bquote(gamma[0]), 
             bquote(sigma[gamma["species"]]),
             bquote(gamma["herb."]),
             bquote(gamma["woody"]), 
             bquote(gamma["specialization"]),
             bquote(gamma["spec.*herb."]),
             bquote(gamma["spec.*woody"]),
             
             bquote(phi[0]), 
             bquote(sigma[phi["species"]]),
             bquote(phi["herb."]),
             bquote(phi["woody"]), 
             bquote(phi["specialization"]),
             bquote(phi["spec.*herb."]),
             bquote(phi["spec.*woody"]),
             
             bquote("p"[0]),
             bquote(sigma["p"["species"]]),
             bquote("p"["specialization"]),
             bquote("p"["date"]),
             bquote(sigma["p"["date - species"]]),
             bquote("p"["date^2"]),
             bquote(sigma["p"["date^2 - species"]]),
             bquote("p"["flower abundance - survey"])
    )
    ) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-3.5, 3)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=0, vjust=0),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip()
)

p <- p +
  geom_errorbar(aes(x=X, ymin=estimates_lower, ymax=estimates_upper),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_point(aes(x=X, y=parameter_value),
             size = 5, alpha = 0.8, shape = 10, colour = "firebrick2") 

p
