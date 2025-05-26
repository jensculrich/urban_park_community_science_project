# estimate regional species richness per city
library(tidyverse)

# define butterfly families to include
butterfly_families <- c("Hesperiidae", "Lycaenidae", "Nymphalidae", 
                        "Papilionidae", "Pieridae")

# list of city names
city_names <- c(
  "LA", # 1
  "NYC" # 2
)

all_city_df <- matrix(nrow=length(city_names), ncol=2)

# read and summarize data across loop
for(i in 1:length(city_names)){
  
  city <- city_names[i]
  
  all_city_df[i,] <- cbind(
    nrow(
      read.csv(paste0(
        "./data/", city, "/02_", city, "_regional_species_pool.csv")) %>%
        filter(species != "") %>%
        filter(family %in% butterfly_families) %>%
        group_by(species) %>%
        slice(1)),
    city)
  
}

all_city_df <- all_city_df %>%
  as.data.frame(.) %>%
  rename("size_of_regional_pool" = "V1",
         "city" = "V2")

write.csv(all_city_df, "./data/size_of_regional_species_pools.csv")  

