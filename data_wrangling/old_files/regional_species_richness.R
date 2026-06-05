# estimate regional species richness per city
library(tidyverse)

##-----------------------------------------------------------------------------
# using BAMONA

# list of city names
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

# list of county names
n_counties <- length(county_names <- c(
  "Fulton",
  "Suffolk", 
  "Mecklenburg",
  "Cook",
  "Dallas",
  "DC",
  "Denton",
  "Harris",
  "LA",
  "Hennepin",
  "Queens",     
  "Philadelphia",
  "Wake",
  "SD",
  "SF"
))

all_city_df <- matrix(nrow=length(city_names), ncol=2)

# read and summarize data across loop
for(i in 1:n_cities){
  
  city <- city_names[i]
  
  all_city_df[i,] <- cbind(
    nrow(
      read.csv(paste0(
        "./data/city_wide_data/BAMONA_species_lists/", city_names[i], "_", county_names[i], "_county.csv"))
        ),
    city)
  
}

all_city_df <- all_city_df %>%
  as.data.frame(.) %>%
  rename("size_of_regional_pool" = "V1",
         "city" = "V2") %>%
  mutate(size_of_regional_pool = as.integer(size_of_regional_pool))

write.csv(all_city_df, "./data/size_of_regional_species_pools_BAMONA.csv", row.names = FALSE)  





##-----------------------------------------------------------------------------
# using iNat

# define butterfly families to include
butterfly_families <- c("Hesperiidae", "Lycaenidae", "Nymphalidae", 
                        "Papilionidae", "Pieridae")

# list of city names
city_names <- c(
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
)

all_city_df <- matrix(nrow=length(city_names), ncol=2)

# read and summarize data across loop
for(i in 1:length(city_names)){
  
  city <- city_names[i]
  
  all_city_df[i,] <- cbind(
    nrow(
      read.csv(paste0(
        "./data/detections_by_city/", city, "/02_0m_", city, "_regional_species_pool.csv")) %>%
        filter(species != "") %>%
        filter(family %in% butterfly_families) %>%
        group_by(species) %>%
        slice(1)),
    city)
  
}

all_city_df <- all_city_df %>%
  as.data.frame(.) %>%
  rename("size_of_regional_pool" = "V1",
         "city" = "V2") %>%
  mutate(size_of_regional_pool = as.integer(size_of_regional_pool))

write.csv(all_city_df, "./data/size_of_regional_species_pools.csv", row.names = FALSE)  

