library(tidyverse)

# select a region
regions <- c(
  "midwest",
  "northeast",
  "southeast",
  "southwest"
)

region <- regions[2]

data <- readRDS(paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))

## --------------------------------------------------
# site predictors
site_data <- data$site_data 

# park size
p1 <- ggplot(site_data, aes(x = log_total_green_space_area, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  #scale_x_continuous(limits = c(4, 20), breaks = c(4, 6, 8, 10, 12, 14, 16, 18, 20)) +
  theme_bw() +
  xlab("log(Park Size in m^2)") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~city)
  
#  connectivity
p2 <- ggplot(site_data, aes(x = isolation, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("log(size and distance weighted isolation\nfrom other parks within 2km)") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~city)

#  tree cover
p3 <- ggplot(site_data, aes(x = tree_percent_cover, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("percent tree cover") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~city)

#  plant genus density
p4 <- ggplot(site_data, aes(x = plant_genera_density, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("log(number of plant genera detected in iNat)\n/ log(park size (m^2))") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~city)

# cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)
p1
p2
p3
p4

## --------------------------------------------------
# species predictors

LA_species_data <- LA_data$species_info %>%
  mutate(city = "LA")

NYC_species_data <- NYC_data$species_info %>%
  mutate(city = "NYC")

SEA_species_data <- SEA_data$species_info %>%
  mutate(city = "SEA")

all_species_data <- rbind(LA_species_data,
                          NYC_species_data,
                          SEA_species_data)

# n binary detection
q1 <- ggplot(all_species_data, aes(x = n_binary_detections, 
                                colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("number of binary detections") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#  wingspan
q2 <- ggplot(all_species_data, aes(x = aveWingspan, 
                                colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("average wingspan (cm)") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#  feature diversity
q3 <- ggplot(all_species_data, aes(x = featureDiversity, 
                                colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("feature diversity") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#  ease of id
q4 <- ggplot(all_species_data, aes(x = research_grade_proportion, 
                                colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("proportion of detections for genus rated as research grade") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cowplot::plot_grid(q1, q2, q3, q4, ncol = 2)
