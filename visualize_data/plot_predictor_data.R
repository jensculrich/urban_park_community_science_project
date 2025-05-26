library(tidyverse)

LA_data <- readRDS("./run_model/prepped_data/prepped_data_LA.rds")
NYC_data <- readRDS("./run_model/prepped_data/prepped_data_NYC.rds")

## --------------------------------------------------
# site predictors
LA_site_data <- LA_data$site_data %>%
  mutate(city = "LA")

NYC_site_data <- NYC_data$site_data %>%
  mutate(city = "NYC")

all_site_data <- rbind(LA_site_data,
                  NYC_site_data)

# park size
p1 <- ggplot(all_site_data, aes(x = log_total_green_space_area, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("log(greenspace area (m^2)))") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
#  connectivity
p2 <- ggplot(all_site_data, aes(x = avg_dist_2000m, 
                          colour = city, fill = city)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  xlab("average distance to other parks within 2km") + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#  tree cover
p3 <- ggplot(all_site_data, aes(x = tree_percent_cover, 
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
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#  plant genus density
p4 <- ggplot(all_site_data, aes(x = plant_genera_density, 
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
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)


## --------------------------------------------------
# site predictors
LA_species_data <- LA_data$species_info %>%
  mutate(city = "LA")

NYC_species_data <- NYC_data$species_info %>%
  mutate(city = "NYC")

all_species_data <- rbind(LA_species_data,
                          NYC_species_data)

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
