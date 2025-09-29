###Cleaning iNaturalist data in Los Angeles and Seattle
library(stringr)
library(dplyr)

#Read in the raw iNat data
leps_CA_data<-read.csv("data/inat_data/02_filtered_data/Leps/leps_data_CA.csv")
leps_WA_data<-read.csv("data/inat_data/02_filtered_data/Leps/leps_data_WA.csv")

flowers_CA_data<-read.csv("E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_CA.csv")


#Filter out data with coordinate uncertainty greated than that of the obscured location and observation that is not identified to species level
#This workflow filters GBIF observations by adjusting the uncertainty threshold based on latitude, accounting for the 0.2° x 0.2° grid cell used for obscured locations, while excluding data with excessively high uncertainty to ensure data quality.


# Function to calculate the appropriate uncertainty threshold based on latitude
calculate_threshold <- function(latitude) {
  lat_distance <- 111320 * 0.2  # Constant for latitude (0.2 degrees)
  lon_distance <- 111320 * cos(latitude * pi / 180) * 0.2  # Varies by latitude for longitude
  return(max(lat_distance, lon_distance))
}

leps_CA_data_filtered <- leps_CA_data%>%
    mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
    filter(!is.na(coordinateUncertaintyInMeters))%>%
    filter(coordinateUncertaintyInMeters < uncertainty_threshold)
    write.csv(leps_CA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_CA_data_coordUncertainty_cleaned.csv")

leps_WA_data_filtered <- leps_WA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_WA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_WA_data_coordUncertainty_cleaned.csv")

leps_NY_data_filtered <- leps_data_NY%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_NY_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_NY_data_coordUncertainty_cleaned.csv")

leps_TX_data_filtered <- leps_data_TX%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_TX_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_TX_data_coordUncertainty_cleaned.csv")

leps_CO_data_filtered <- leps_data_CO%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_CO_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_CO_data_coordUncertainty_cleaned.csv")

leps_FL_data_filtered <- leps_data_FL%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_FL_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_FL_data_coordUncertainty_cleaned.csv")

leps_MA_data_filtered <- leps_data_MA%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_MA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_MA_data_coordUncertainty_cleaned.csv")

leps_PA_data_filtered <- leps_data_PA%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_PA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_PA_data_coordUncertainty_cleaned.csv")

leps_SD_data_filtered <- leps_data_SD%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_SD_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_SD_data_coordUncertainty_cleaned.csv")

leps_DC_data_filtered <- leps_data_DC%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_DC_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_DC_data_coordUncertainty_cleaned.csv")

leps_NC_data_filtered <- leps_data_NC%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(leps_NC_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_NC_data_coordUncertainty_cleaned.csv")

leps_MN_data_filtered <- leps_data_MN%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(leps_MN_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_MN_data_coordUncertainty_cleaned.csv")

leps_Essex_data_filtered <- fread("data/inat_data/01_raw_data/iNat_leps_essex.csv")%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(leps_Essex_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_Essex_data_coordUncertainty_cleaned.csv")

leps_MI_data_filtered <- leps_data_MI%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(leps_MI_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_MI_data_coordUncertainty_cleaned.csv")

leps_data_GA<- read.csv("data/inat_data/02_filtered_data/Leps/leps_data_GA.csv")
leps_GA_data_filtered <- leps_data_GA%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_GA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_GA_data_coordUncertainty_cleaned.csv")

leps_AZ_data_filtered <- leps_data_AZ%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_AZ_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_AZ_data_coordUncertainty_cleaned.csv")



leps_IL_data_filtered <- leps_data_IL%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_IL_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_IL_data_coordUncertainty_cleaned.csv")


leps_IA_data_filtered <- leps_data_IA%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_IA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_IA_data_coordUncertainty_cleaned.csv")


leps_MO_data_filtered <- leps_data_MO%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)
write.csv(leps_MO_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Leps/leps_MO_data_coordUncertainty_cleaned.csv")

###Flowering Plants
flowers_CA_data_filtered <- flowers_CA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_CA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_CA_data_coordUncertainty_cleaned.csv")

flowers_NY_data_filtered <- flowers_NY_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_NY_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_NY_data_coordUncertainty_cleaned.csv")


flowers_WA_data_filtered <- flowers_WA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_WA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_WA_data_coordUncertainty_cleaned.csv")

flowers_TX_data_filtered <- flowers_TX_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_TX_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_TX_data_coordUncertainty_cleaned.csv")

flowers_CO_data_filtered <- flowers_CO_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_CO_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_CO_data_coordUncertainty_cleaned.csv")


flowers_FL_data_filtered <- flowers_FL_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_FL_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_FL_data_coordUncertainty_cleaned.csv")


flowers_MA_data_filtered <- flowers_MA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_MA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_MA_data_coordUncertainty_cleaned.csv")

flowers_PA_data_filtered <- flowers_PA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_PA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_PA_data_coordUncertainty_cleaned.csv")

flowers_SD_data_filtered <- flowers_SD_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_SD_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_SD_data_coordUncertainty_cleaned.csv")


flowers_DC_data_filtered <- flowers_DC_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_DC_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_DC_data_coordUncertainty_cleaned.csv")



flowers_NC_data_filtered <- flowers_NC_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_NC_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_NC_data_coordUncertainty_cleaned.csv")


flowers_MN_data_filtered <- flowers_MN_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_MN_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_MN_data_coordUncertainty_cleaned.csv")

flowers_Essex_data_filtered <- fread("data/inat_data/01_raw_data/iNat_flowering_plants_essex.csv")%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_Essex_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_Essex_data_coordUncertainty_cleaned.csv")

flowers_MI_data_filtered <- flowers_MI_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_MI_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_MI_data_coordUncertainty_cleaned.csv")

flowers_GA_data <- read.csv("data/inat_data/02_filtered_data/Plants/inat_flowers_GA.csv")

flowers_GA_data_filtered <- flowers_GA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_GA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_GA_data_coordUncertainty_cleaned.csv")

flowers_AZ_data_filtered <- flowers_AZ_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_AZ_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_AZ_data_coordUncertainty_cleaned.csv")

flowers_IL_data_filtered <- flowers_IL_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_IL_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_IL_data_coordUncertainty_cleaned.csv")

flowers_IA_data_filtered <- flowers_IA_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_IA_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_IA_data_coordUncertainty_cleaned.csv")

flowers_MO_data_filtered <- flowers_MO_data%>%
  mutate(uncertainty_threshold = calculate_threshold(decimalLatitude)) %>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  filter(coordinateUncertaintyInMeters < uncertainty_threshold)

write.csv(flowers_MO_data_filtered,"E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/flowers_MO_data_coordUncertainty_cleaned.csv")

# # Create comparison with subspecies information
# comparison_CA <- leps_CA_data_filtered %>%
#   mutate(
#     # Extract genus + species from verbatimScientificName
#     verbatim_simple = str_extract(verbatimScientificName, "^\\w+\\s+\\w+"),
#     
#     # Check if species names match
#     species_match = species == verbatim_simple)%>%
#   select(species, verbatimScientificName, verbatim_simple, species_match, issue)%>%
#   filter(species_match=="FALSE")
#     
# comparison_WA <- leps_WA_data_filtered %>%
#   mutate(
#     # Extract genus + species from verbatimScientificName
#     verbatim_simple = str_extract(verbatimScientificName, "^\\w+\\s+\\w+"),
#     
#     # Check if species names match
#     species_match = species == verbatim_simple)%>%
#   select(species, verbatimScientificName, verbatim_simple, species_match, issue)%>%
#   filter(species_match=="FALSE")
#    
# 
# # Summary of matches/mismatches
# summary_stats <- comparison %>%
#   summarize(
#     total_records = n(),
#     matching_species_names = sum(species_match, na.rm = TRUE),
#     genus_mismatches = sum(genus_mismatch, na.rm = TRUE),
#     epithet_mismatches = sum(epithet_mismatch, na.rm = TRUE),
#     records_with_gbif_subspecies = sum(has_gbif_subspecies, na.rm = TRUE),
#     records_with_inat_subspecies = sum(has_inat_subspecies, na.rm = TRUE),
#     matching_subspecies = sum(subspecies_match, na.rm = TRUE),
#     subspecies_differences = sum(subspecies_difference, na.rm = TRUE)
#   )
# 
# print(summary_stats)
# 
# # Get examples of different types of mismatches
# species_mismatch_examples <- comparison %>%
#   filter(!species_match) %>%
#   select(species, verbatimScientificName, verbatim_simple, gbif_genus, inat_genus, gbif_species, inat_species) %>%
#   head(10)
# 
# subspecies_mismatch_examples <- comparison %>%
#   filter(subspecies_difference) %>%
#   select(species, verbatimScientificName, infraspecificEpithet, verbatim_subspecies) %>%
#   head(10)
# 
# print("Species mismatch examples:")
# print(species_mismatch_examples)
# 
# print("Subspecies difference examples:")
# print(subspecies_mismatch_examples)