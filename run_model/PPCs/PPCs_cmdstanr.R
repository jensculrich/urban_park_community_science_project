# how many detections would you expect per species (W)?
# we will compare this to how many detections per species simulated in the 
# generated quantities block of our model (W_rep) for a visual PPC

library(rstan)
library(tidyverse)

#city_name <- "los_angeles"
#my_data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", city_name, ".RDS"))

## --------------------------------------------------
### Prepare data for model
# for community sampling events inferred by [taxonomic] family, source this file:

my_data <- readRDS( paste0("./run_model/prepped_data/prepped_data_all.rds"))

## --------------------------------------------------
### Prepare data for model

# data to feed to the model
V <- my_data$V # detection data
site_data <- my_data$site_data
city <- as.integer(as.factor(unique(site_data$city)))
n_cities <- length(unique(city))
city_names <- unique(site_data$city)

city_integer_vector <- my_data$city_integer_vector
W_df <- as.data.frame(cbind(city_integer_vector, V)) %>%
  mutate(row_sum_detections = rowSums(.[,2:ncol(.)])) %>%
  select(city_integer_vector, row_sum_detections) %>%
  group_by(city_integer_vector) %>%
  summarise(W_city = sum(row_sum_detections),
            .groups = 'drop') %>%
  cbind(city_names)


# get W distributions from model
## get param estimates from the region
stan_out <- readRDS(
  "./model_outputs/stan_out_dec19.rds")

tmp <- as.data.frame(stan_out$draws(variables = "W_city_rep",
                                    format = "draws_matrix"
)) # take estimates from each HMC step as a df

rm(stan_out)
gc()

#n_samp <- 10 # how mann_chains#n_samp <- 10 # how many samples do we have from the HMC run?
n_samp <- length(tmp[,1]) # how many samples do we have from the HMC run?

## --------------------------------------------------
# 

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# get the mean, 50 and 90% BCIs for initial occurrence for each site
quants <- apply(
  X = tmp,
  MARGIN = 2, # 2 indicates columns
  FUN = quantile,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95) # Optional: specify desired probabilities
)

df_estimates <- t(data.frame(quants))

df_estimates <- cbind(W_df, df_estimates) %>%
  mutate(city_names = as.factor(city_names)) %>%
  rename("lower_95" = "5%",
         "upper_95" = "95%",
         "lower_50" = "25%",
         "upper_50" = "75%") 

(ggplot(df_estimates) +
  geom_segment(aes(x = city_names, y = lower_95, yend = upper_95), 
               linewidth = 10, colour = c_light) +
  geom_segment(aes(x = city_names, y = lower_50, yend = upper_50), 
               linewidth = 10, colour = c_mid_highlight) +  
  geom_point(aes(x=city_names, y = W_city), size = 3, shape = 10) +
  theme_classic() +
  scale_y_continuous(name = "Number of Detections") +
  theme(axis.text.x = element_text(size = 18, angle=45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
)
