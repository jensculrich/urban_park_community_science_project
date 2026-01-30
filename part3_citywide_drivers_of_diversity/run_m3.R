library(tidyverse)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

n_cities <- length(city_names <- c(
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
))


simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_richness <- simmed_diversity[[1]]
mean_prop_disturbance_avoidant <- simmed_diversity[[2]]
beta_diversity <- simmed_diversity[[3]]
gamma_diversity <- simmed_diversity[[4]]

#-------------------------------------------------------------------------------
# get the city covariate data
city_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")

city_data <- city_data %>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:n_cities))

city_data <- city_data %>%
  mutate(log_avg_park_size = log(average_park_size_sqm),
         log_park_size_scaled = center_scale(log_avg_park_size))


#-------------------------------------------------------------------------------
# summarize uncertainty and plot relationships

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

#-------------------------------------------------------------------------------
# mean species richness

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

mean_richness_quantiles_df <- cbind(city_data, mean_richness_quantiles_df)

a <- ggplot(mean_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# mean proportion disturbance or edge avoidant
# are there more ruderal species in cities with smaller parks?

#  calculate Means and CI's for the diversity metrics for each city
mean_prop_disturbance_avoidant_quantiles <- apply(mean_prop_disturbance_avoidant, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_prop_disturbance_avoidant_quantiles_df <- as.data.frame(t(mean_prop_disturbance_avoidant_quantiles))
colnames(mean_prop_disturbance_avoidant_quantiles_df) <- c("lower90", "lower50",
                                                           "mean", "upper50", "upper90")

mean_prop_disturbance_avoidant_quantiles_df <- cbind(city_data, mean_prop_disturbance_avoidant_quantiles_df)

a2 <- ggplot(mean_prop_disturbance_avoidant_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Prop. Disturbance Avoidant") +
  xlab("Mean log(Park Size)") +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# beta diversity (mean jaccard dissimilarity)

#  calculate Means and CI's for the diversity metrics for each city
mean_beta_quantiles <- apply(beta_diversity, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                             na.rm=TRUE)
mean_beta_quantiles_df <- as.data.frame(t(mean_beta_quantiles))
colnames(mean_beta_quantiles_df) <- c("lower90", "lower50",
                                      "mean", "upper50", "upper90")

mean_beta_quantiles_df <- cbind(city_data, mean_beta_quantiles_df)

b <- ggplot(mean_beta_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Dissimilarity\n(Jaccard Index)") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# gamma (city-wide) species richness

#  calculate Means and CI's for the diversity metrics for each city
gamma_richness_quantiles <- apply(gamma_diversity, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
gamma_richness_quantiles_df <- as.data.frame(t(gamma_richness_quantiles))
colnames(gamma_richness_quantiles_df) <- c("lower90", "lower50",
                                           "mean", "upper50", "upper90")

gamma_richness_quantiles_df <- cbind(city_data, gamma_richness_quantiles_df)

c <- ggplot(gamma_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(a, a2, b, c, ncol = 3)

#-------------------------------------------------------------------------------
# do it again with relative species richness

# alpha and gamma diversity could be relative to the regional species pools
#size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools_BAMONA.csv")
size_of_regional_species_pools <- read.csv("./data/size_of_regional_species_pools.csv")

# get number of species actually modelled in each city
for(city_number in 1:n_cities){
  # get the correct range data
  temp_ranges <- filter(range_data, city == city_names[city_number])
  # only consider species that were modelled
  temp_ranges <- temp_ranges %>%
    filter(species %in% species_info$species)
  size_of_regional_species_pools[city_number,1] <- nrow(temp_ranges)
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# mean species richness

# standardize mean richness
mean_richness_relative <- mean_richness
for(i in 1:n_cities){
  mean_richness_relative[i,] = mean_richness_relative[i,] / size_of_regional_species_pools[i,1]
}

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness_relative, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

mean_richness_quantiles_df <- cbind(city_data, mean_richness_quantiles_df)


d <- ggplot(mean_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness (Relative)") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

#-------------------------------------------------------------------------------
# gamma (city-wide) species richness

# standardize gamma richness
gamma_diversity_relative <- gamma_diversity
for(i in 1:n_cities){
  gamma_diversity_relative[i,] = gamma_diversity_relative[i,] / size_of_regional_species_pools[i,1]
}

#  calculate Means and CI's for the diversity metrics for each city
gamma_richness_quantiles <- apply(gamma_diversity_relative, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
gamma_richness_quantiles_df <- as.data.frame(t(gamma_richness_quantiles))
colnames(gamma_richness_quantiles_df) <- c("lower90", "lower50",
                                           "mean", "upper50", "upper90")

gamma_richness_quantiles_df <- cbind(city_data, gamma_richness_quantiles_df)

f <- ggplot(gamma_richness_quantiles_df, aes(log_avg_park_size , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1)  +
  geom_point(aes(colour=city), size = 4) +
  ylab("Total Number of Species\nOccurring in City Parks (Relative)") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(d, b, f, ncol = 3)


