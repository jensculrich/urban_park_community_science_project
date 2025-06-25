# Load the required packages ----
library(rgbif)
library(dplyr)
library(data.table)

# Set up query parameters
 leps_taxon_key <- name_backbone(name = "Lepidoptera")$usageKey
 plant_kingdom_key <- name_backbone(name = "Plantae")$usageKey  # Plantae kingdom key
 
 flower_class_key<-name_backbone_checklist(c("Magnoliopsida", "Liliopsida"))$usageKey
 
 dataset_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"  # Dataset key for iNaturalist
 start_date <- "2020-01-01"
 end_date <- "2024-12-31"
 country_code <- "US"  # United States

 

# Define the download request using date predicates
download_key <- occ_download(
  user = "cralibe",
  pwd = "A9aafbbd4a*",
  email = "yanyinje@usc.edu",
  pred("taxonKey", leps_taxon_key),
  pred("datasetKey", dataset_key),  # iNaturalist Research-grade dataset
  pred("country", country_code),  # United States
  pred("hasCoordinate", TRUE),  # Only records with coordinates
  #pred("hasGeospatialIssue", TRUE),  # Include records with geospatial issues
  pred_gte("eventDate", start_date),  # Event date greater than or equal to start_date
  pred_lte("eventDate", end_date),  # Event date less than or equal to end_date
  pred("occurrenceStatus", "present"),  # Only records with occurrence status "present"
  format = "SIMPLE_CSV"  # Format as a simple CSV file
)

# Define the download request using date predicates
download_key <- occ_download(
  user = "cralibe",
  pwd = "A9aafbbd4a*",
  email = "yanyinje@usc.edu",
  pred_or(
    pred("taxonKey", flower_class_key[1]),
    pred("taxonKey", flower_class_key[2])
  ),
  pred("datasetKey", dataset_key),  # iNaturalist Research-grade dataset
  pred("country", country_code),  # United States
  pred("hasCoordinate", TRUE),  # Only records with coordinates
  #pred("hasGeospatialIssue", TRUE),  # Include records with geospatial issues
  pred_gte("eventDate", start_date),  # Event date greater than or equal to start_date
  pred_lte("eventDate", end_date),  # Event date less than or equal to end_date
  pred("occurrenceStatus", "present"),  # Only records with occurrence status "present"
  format = "SIMPLE_CSV"  # Format as a simple CSV file
)

# Monitor the download
occ_download_wait(download_key)
# Download the data
lep_download <- occ_download_get(download_key)
flower_download <- occ_download_get(download_key)

#Check the first few line
file_lines <- readLines("data/inat_data/01_raw_data/iNat_leps_US.csv", n = 10)

#Read in the downloaded data ----
leps_us_data <- fread("data/inat_data/01_raw_data/iNat_leps_US.csv",header = TRUE, sep = "\t", fill = TRUE, quote = "")
head(leps_us_data)


leps_data_CA <- leps_us_data%>%
  filter(stateProvince=="California")

write.csv(leps_data_CA, "data/inat_data/02_filtered_data/Leps/leps_data_CA.csv")

leps_data_WA <- leps_us_data%>%
  filter(stateProvince=="Washington")

write.csv(leps_data_WA, "data/inat_data/02_filtered_data/Leps/leps_data_WA.csv")


leps_data_NY <- leps_us_data%>%
  filter(stateProvince=="New York" | stateProvince=="New Jersey")

write.csv(leps_data_NY, "data/inat_data/02_filtered_data/Leps/leps_data_NY.csv")






####Flowering PLANTS####

#Check the first few line
file_lines <- readLines("E:/phd_study/urban_park_community_science_project/data/inat_data/01_raw_data/iNat_flowering_plants_US.csv", n = 10)

#Read in the downloaded data ----
flowers_us_data<- fread("E:/phd_study/urban_park_community_science_project/data/inat_data/01_raw_data/iNat_flowering_plants_US.csv", header = TRUE, sep = "\t", fill = TRUE,  quote = "")


flowers_CA_data<-flowers_us_data%>%
  filter(stateProvince=="California")


write.csv(flowers_CA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_CA.csv")


flowers_NY_data<-flowers_us_data%>%
  filter(stateProvince=="New York")


write.csv(flowers_NY_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_NY.csv")


flowers_WA_data<-flowers_us_data%>%
  filter(stateProvince=="Washington")


write.csv(flowers_WA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_WA.csv")


