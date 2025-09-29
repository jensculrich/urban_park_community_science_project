# Load the required packages ----
library(rgbif)
library(dplyr)
library(data.table)

# Set up query parameters
leps_taxon_key <- name_backbone(name = "Lepidoptera")$usageKey
plant_kingdom_key <- name_backbone(name = "Plantae")$usageKey  # Plantae kingdom key
flower_class_key<-name_backbone_checklist(c("Magnoliopsida", "Liliopsida"))$usageKey

 #common keys
 dataset_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"  # Dataset key for iNaturalist
 start_date <- "2020-01-01"
 end_date <- "2024-12-31"
 country_code <- "US"  # United States

#Leps
# Define the download request using date predicates
download_key <- occ_download(
  user = "cralibe",
  pwd = "A9aafbbd4a*",
  email = "yanyinje@usc.edu",
  pred("taxonKey", leps_taxon_key),
  pred("datasetKey", dataset_key),  # iNaturalist Research-grade dataset
  #pred("country", country_code),  # United States
  #pred("gadm","MEX.3_1"), #for San Deigo only (we need observations in Baja California)
  pred("gadm","CAN.9.9_1"), #for Detroit, we need essex county in Canada
  pred("hasCoordinate", TRUE),  # Only records with coordinates
  #pred("hasGeospatialIssue", TRUE),  # Include records with geospatial issues
  pred_gte("eventDate", start_date),  # Event date greater than or equal to start_date
  pred_lte("eventDate", end_date),  # Event date less than or equal to end_date
  pred("occurrenceStatus", "present"),  # Only records with occurrence status "present"
  format = "SIMPLE_CSV"  # Format as a simple CSV file
)

##Flowers
# Define the download request using date predicates
download_key <- occ_download(
  user = "cralibe",
  pwd = "A9aafbbd4a*",
  email = "yanyinje@usc.edu",
  pred_or(
    pred("taxonKey", flower_class_key[1]),
    pred("taxonKey", flower_class_key[2])),
  pred("datasetKey", dataset_key),  # iNaturalist Research-grade dataset
  #pred("country", country_code),  # United States
  #pred("gadm","MEX.3_1"), #for San Deigo only (we need observations in Baja California)
  pred("gadm","CAN.9.9_1"), #for Detroit, we need essex county in Canada
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


leps_data_FL <- leps_us_data%>%
  filter(stateProvince=="Florida")

write.csv(leps_data_FL, "data/inat_data/02_filtered_data/Leps/leps_data_FL.csv")


leps_data_FL <- leps_us_data%>%
  filter(stateProvince=="Florida")

write.csv(leps_data_FL, "data/inat_data/02_filtered_data/Leps/leps_data_FL.csv")

leps_data_GA <- leps_us_data%>%
  filter(stateProvince=="Georgia")

write.csv(leps_data_GA, "data/inat_data/02_filtered_data/Leps/leps_data_GA.csv")

leps_data_TX <- leps_us_data%>%
  filter(stateProvince=="Texas")
write.csv(leps_data_TX, "data/inat_data/02_filtered_data/Leps/leps_data_TX.csv")

leps_data_CO <- leps_us_data%>%
  filter(stateProvince=="Colorado")
write.csv(leps_data_CO, "data/inat_data/02_filtered_data/Leps/leps_data_CO.csv")

leps_data_MA <- leps_us_data%>%
  filter(stateProvince=="Massachusetts")
write.csv(leps_data_MA, "data/inat_data/02_filtered_data/Leps/leps_data_MA.csv")

#For PA, we also need NJ data because philadelphia is right next to NJ
leps_data_PA <- leps_us_data%>%
  filter(stateProvince=="Pennsylvania" | stateProvince=="New Jersey")
write.csv(leps_data_PA, "data/inat_data/02_filtered_data/Leps/leps_data_PA.csv")

#We need a specialized data for San Deigo because it is next to Baja California
leps_data_baja <-fread("data/inat_data/01_raw_data/iNat_leps_Baja_California.csv")

leps_data_SD <-  leps_us_data%>%
  filter(stateProvince=="California")%>%
  rbind(leps_data_baja)
write.csv(leps_data_SD, "data/inat_data/02_filtered_data/Leps/leps_data_SD.csv")

#We need a specialized data for DC because it is next to Virgina and Maryland
leps_data_DC <- leps_us_data%>%
  filter(stateProvince=="District of Columbia" | stateProvince=="Maryland" | stateProvince=="Virginia")
write.csv(leps_data_DC, "data/inat_data/02_filtered_data/Leps/leps_data_DC.csv")

leps_data_NC <- leps_us_data%>%
  filter(stateProvince=="North Carolina" | stateProvince=="South Carolina")
write.csv(leps_data_NC, "data/inat_data/02_filtered_data/Leps/leps_data_NC.csv")

leps_data_MN <- leps_us_data%>%
  filter(stateProvince=="Minnesota")
write.csv(leps_data_MN, "data/inat_data/02_filtered_data/Leps/leps_data_MN.csv")

leps_data_MI <- leps_us_data%>%
  filter(stateProvince=="Michigan")
write.csv(leps_data_MI, "data/inat_data/02_filtered_data/Leps/leps_data_MI.csv")

leps_data_AZ <- leps_us_data%>%
  filter(stateProvince=="Arizona")
write.csv(leps_data_AZ, "data/inat_data/02_filtered_data/Leps/leps_data_AZ.csv")

leps_data_IL <- leps_us_data%>%
  filter(stateProvince=="Illinois" | stateProvince=="Indiana")
table(leps_data_IL$stateProvince)
write.csv(leps_data_AZ, "data/inat_data/02_filtered_data/Leps/leps_data_AZ.csv")

leps_data_IA <- leps_us_data%>%
  filter(stateProvince=="Iowa")
write.csv(leps_data_IA, "data/inat_data/02_filtered_data/Leps/leps_data_IA.csv")

leps_data_MO <- leps_us_data%>%
  filter(stateProvince=="Illinois" | stateProvince=="Missouri")
write.csv(leps_data_MO, "data/inat_data/02_filtered_data/Leps/leps_data_MO.csv")



####Flowering PLANTS####

#Check the first few line
file_lines <- readLines("E:/phd_study/urban_park_community_science_project/data/inat_data/01_raw_data/iNat_flowering_plants_US.csv", n = 10)

#Read in the downloaded data ----
flowers_us_data<- fread("E:/phd_study/urban_park_community_science_project/data/inat_data/01_raw_data/iNat_flowering_plants_US.csv", header = TRUE, sep = "\t", fill = TRUE,  quote = "")


flowers_CA_data<-flowers_us_data%>%
  filter(stateProvince=="California"
write.csv(flowers_CA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_CA.csv")


flowers_NY_data<-flowers_us_data%>%
  filter(stateProvince=="New York")
write.csv(flowers_NY_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_NY.csv")


flowers_WA_data<-flowers_us_data%>%
  filter(stateProvince=="Washington")
write.csv(flowers_WA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_WA.csv")

flowers_FL_data<-flowers_us_data%>%
  filter(stateProvince=="Florida")
write.csv(flowers_FL_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_FL.csv")

flowers_GA_data<-flowers_us_data%>%
  filter(stateProvince=="Georgia")
write.csv(flowers_GA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_GA.csv")

flowers_TX_data<-flowers_us_data%>%
  filter(stateProvince=="Texas")
write.csv(flowers_TX_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_TX.csv")

flowers_CO_data<-flowers_us_data%>%
  filter(stateProvince=="Colorado")
write.csv(flowers_CO_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_CO.csv")

flowers_MA_data<-flowers_us_data%>%
  filter(stateProvince=="Massachusetts")
write.csv(flowers_MA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_MA.csv")

flowers_PA_data <- flowers_us_data%>%
  filter(stateProvince=="Pennsylvania" | stateProvince=="New Jersey")
write.csv(flowers_PA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_PA.csv")

flowers_baja_data <-fread("data/inat_data/01_raw_data/iNat_flowering_plants_Baja_California.csv")

flowers_SD_data <-  flowers_us_data%>%
  filter(stateProvince=="California")%>%
  rbind(flowers_baja_data)
write.csv(flowers_SD_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_SD.csv")

#
flowers_DC_data <- flowers_us_data%>%
  filter(stateProvince=="District of Columbia" | stateProvince=="Maryland" | stateProvince=="Virginia")
write.csv(flowers_DC_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_DC.csv")

flowers_NC_data <- flowers_us_data%>%
  filter(stateProvince=="North Carolina"| stateProvince=="South Carolina")
write.csv(flowers_NC_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_NC.csv")

flowers_MN_data <- flowers_us_data%>%
  filter(stateProvince=="Minnesota")
write.csv(flowers_MN_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_MN.csv")

flowers_MI_data <- flowers_us_data%>%
  filter(stateProvince=="Michigan")
write.csv(flowers_MI_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_MI.csv")

flowers_AZ_data <- flowers_us_data%>%
  filter(stateProvince=="Arizona")
write.csv(flowers_AZ_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_AZ.csv")


flowers_IL_data <- flowers_us_data%>%
  filter(stateProvince=="Illinois" | stateProvince=="Indiana")
write.csv(flowers_IL_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_AZ.csv")

flowers_IA_data <- flowers_us_data%>%
  filter(stateProvince=="Iowa")
write.csv(flowers_IA_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_IA.csv")

flowers_MO_data <- flowers_us_data%>%
  filter(stateProvince=="Illinois" | stateProvince=="Missouri")
write.csv(flowers_MO_data, "E:/phd_study/urban_park_community_science_project/data/inat_data/02_filtered_data/Plants/inat_flowers_MO.csv")








