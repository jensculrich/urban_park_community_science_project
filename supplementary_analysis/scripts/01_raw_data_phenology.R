# This script is for look at the active preiod of each species in each city's regional pool
library(tidyverse)
library(cowplot)
library(vegan)

#list all the files in the repo
data_path <- "~/Documents/research/urban_park_community_science_project/data/inat_clean_data"
file_list<-list.files(path = data_path, full.names = TRUE, pattern = "\\.csv$") 
data_list <- lapply(file_list, read.csv)

list_names <- sapply(file_list, function(x) str_match(x, paste0("50m_", "(.*?)\\", "_regional"))[,2])

data_list <-set_names(data_list, list_names)

plot_list <- list()

for (city in list_names){
  filtered_data<-data_list[[city]]%>%
    select(family, genus, species, scientificName, eventDate, day, month, year, speciesKey)%>%
    group_by(speciesKey, month, year)%>%
    summarise(no_obv_per_month=n())%>%
    ungroup()%>%
    group_by(speciesKey, month)%>%
    summarise(mean_obv_no=mean(no_obv_per_month))
  
  plot<-filtered_data%>%
    ggplot(aes(x= month, y= mean_obv_no, group = speciesKey, color = as.factor(speciesKey)))+
    geom_line()+
    scale_x_discrete(name = "Month",                  
                     limits = month.abb,          
                     labels = month.abb)+
    ylab("average number of observations")+
    labs(title=paste0(ifelse(nchar(city) > 2, str_to_sentence(city), city))
)+
    theme_cowplot()+
    theme(legend.position = "NONE")
   
  plot_list[[city]] <- plot
  
}

lapply(names(plot_list), function(city) {
  ggsave(paste0("~/Documents/research/urban_park_community_science_project/supplementary_analysis/plots/", city, ".pdf"), 
         plot = plot_list[[city]], 
         dpi = 300, 
         width = 8, 
         height = 5, 
         bg = "white")
})


city_matrix <- map2(data_list, names(data_list), ~mutate(.x, city = .y))%>%
  map(~select(.x, city, speciesKey, month, year))%>%
  map(~group_by(.x, city, month, year))%>%
  map(~summarise(.x, n_obv_month=n()))%>%
  map(~ungroup(.x))%>%
  map(~mutate(.x, obv_month_pa=ifelse(n_obv_month>0, 1, 0)))%>%
  map(~group_by(.x, city, month))%>%
  map(~summarise(.x, n_obv_month=sum(obv_month_pa)))%>%
  map(~ungroup(.x))%>%
  bind_rows()%>%
  pivot_wider(names_from = month, values_from = n_obv_month, values_fill = 0.0)


bray_dist <- vegdist(city_matrix[, 2:ncol(city_matrix)], method = "bray")

pcoa_results <- cmdscale(bray_dist, k = 2, eig = TRUE)

explained_var <- pcoa_results$eig[1:2] / sum(pcoa_results$eig) * 100

data.frame(pcoa_results$points, city = as.vector(list_names))%>%
  ggplot(aes(x = X1, y = X2))+
  geom_point() +
  geom_text(aes(label = city), vjust = -2)+
  labs(
    title = "PCoA with Bray-Curtis Dissimilarity based on average total number of observation per month",
    x = paste0("PCo1 (", round(explained_var[1], 2), "%)"),
    y = paste0("PCo2 (", round(explained_var[2], 2), "%)")
  ) +
  ylim(-0.25,0.6)+
  theme_cowplot()+
  theme(legend.position = "NONE")


