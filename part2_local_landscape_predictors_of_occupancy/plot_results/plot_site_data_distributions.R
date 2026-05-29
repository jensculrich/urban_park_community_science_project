library(tidyverse)
library(ggridges)

city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Denver",
  "Des_moines",
  "Detroit",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "SD",
  "SF",
  "St_louis",
  "Tampa"
)

city_names_labels <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "Washington D.C.",
  "Denton",
  "Denver",
  "Des Moines",
  "Detroit",
  "Houston",
  "Los Angeles",
  "Minneapolis",
  "New York City",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "San Diego",
  "San Francisco",
  "St. Louis",
  "Tampa"
)

n_cities <- length(city_names)

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

city_names_labels <- as.data.frame(cbind(city_names_labels, seq(1:length(city_names_labels)))) %>%
  rename("city_number" = "V2") %>%
  mutate(city_number = as.integer(city_number))

df <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))$site_data

df <- df %>%
  group_by(city) %>%
  mutate(city_number = cur_group_id()) %>%
  ungroup() %>%
  left_join(., city_names_labels)

df <- df %>%
  mutate(connectivity = -1 * log(isolation))

(p <- ggplot(df, aes(log_total_green_space_area, city_names_labels, fill=city_names_labels)) +
  geom_density_ridges() +
  theme_classic() +
  xlab("log(park size(m^2))") +
  ylab("") +
  scale_fill_manual(values=my_palette, labels=city_names_labels) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
)

(q <- ggplot(df, aes(connectivity, city_names_labels, fill=city_names_labels)) +
  geom_density_ridges() +
  theme_classic() +
  xlab("log(connectivity)") +
  ylab("") +
  scale_fill_manual(values=my_palette, labels=city_names_labels) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
)

cowplot::plot_grid(p, q, ncol = 2)

(r <- ggplot(df, aes(log_n_plant_genera, city_names_labels, fill=city_names_labels)) +
    geom_density_ridges() +
    theme_classic() +
    xlab("log(n plant genera)") +
    ylab("") +
    scale_fill_manual(values=my_palette, labels=city_names_labels) + 
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
)

(s <- ggplot(df, aes(proportion_landscape_vegetation, city_names_labels, fill=city_names_labels)) +
    geom_density_ridges() +
    theme_classic() +
    #xlim(c(0, 0.1)) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 0.5)) +
    xlab("landscape herbaceous or woody cover") +
    ylab("") +
    scale_fill_manual(values=my_palette, labels=city_names_labels) + 
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
)

cowplot::plot_grid(r, s, ncol = 2)
