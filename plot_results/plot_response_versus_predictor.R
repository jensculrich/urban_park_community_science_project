# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

# handy for viewing column numbers
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

region <- "northeast"

# list of city names
cities_northeast <- c(
  "Boston", 
  "DC",
  "NYC", 
  "Philadelphia"
)

n_cities <- length(cities_northeast)

## get param estimates from the region
stan_out <- readRDS(paste0(
  "./model_outputs/stan_out_", region, "_2km_connectivity_family_50buffers_simple.rds"))
tmp <- as.data.frame(stan_out) # take estimates from each HMC step as a df
n_samp <- 10 # how many samples do we have from the HMC run?
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?

## get data from region
df <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))$site_data
  
pred_length <- nrow(df)

## ilogit and logit functions
ilogit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

#-------------------------------------------------------------------------------
# get some prediction data

size_pred <- seq(from = -2, to = 2, length.out = pred_length)

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

initial_occurrence <- vector(length = pred_length)

predC <- array(NA, dim=c(pred_length, n_samp)) # community means (overall region mean)
predSpec <- array(NA, dim=c(pred_length, n_samp, 2, n_cities)) # trends by city 

for(i in 1:n_samp){
  
  # community means don't depend on city effects
  predC[,i] <- ilogit( # park size trend
    # psi1_0 +
    tmp[i,1] + 
      # psi1_ +
      tmp[i,6]*size_pred
  )
    
}

# posterior means by community average 
criC <- apply(predC, c(1), function(x) quantile(x, 
              prob = c(0.05, 0.25, 0.5, 0.75, 0.95)))

# posterior means for specialization specific 
criSpec <- apply(predSpec, c(1,3,4), function(x) quantile(x, prob = c(0.25, 0.5, 0.75)))

#-------------------------------------------------------------------------------

# community plot - park size - psi1
size_df <- as.data.frame(cbind(size_pred, criC[3,], 
                               criC[1,], criC[5,],
                               criC[2,], criC[4,])) %>%
  rename("size_pred" = "size_pred",
         "psi1_size_community_mean" = "V2",
         "psi1_size_community_lower95" = "V3",
         "psi1_size_community_upper95" = "V4",
         "psi1_size_community_lower50" = "V5",
         "psi1_size_community_upper50" = "V6")

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
  #scale_fill_manual(values=my_palette) +
  #scale_colour_manual(values=my_palette) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle=45, vjust=-0.5),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

