# Load the required packages ----
library(rgbif)

# Set up query parameters
taxon_key <- name_backbone(name = "Lepidoptera")$usageKey
dataset_key <- "50c9509d-22c7-4a22-a47d-8c48425ef4a7"  # Dataset key for iNaturalist

start_date <- "2020-01-01"
end_date <- "2024-04-30"


# Define the download request using date predicates
download_key <- occ_download(
  pred("taxonKey", taxon_key),
  pred("datasetKey", dataset_key),
  pred_gte("eventDate", start_date),
  pred_lte("eventDate", end_date),
  format = "SIMPLE_CSV"
)
# Monitor the download
occ_download_wait(download_key)
# Download the data
lep_download <- occ_download_get(download_key)
#Check the first few line
file_lines <- readLines("/Users/jennycheung/Documents/phd_study/iNat_project/data/inat_data/leps_us_data.csv", n = 10)

#Read in the downloaded data ----
lep_data<-read.table("/Users/jennycheung/Documents/phd_study/iNat_project/data/0018976-240506114902167.csv", header = TRUE, sep = "\t", fill = TRUE, quote = "")
head(lep_data)

## Filter for the US only ----
leps_us_data<-lep_data%>%
  filter(countryCode=="US")
write.csv(leps_us_data, "~/Documents/phd_study/iNat_project/data/leps_us_data.csv")

### Filter for the necessary columns ----
leps_filtered_data<-leps_data%>%
  select(gbifID, kingdom:species, taxonRank, stateProvince, countryCode, occurrenceStatus, decimalLatitude, decimalLongitude, day:year)
write.csv(leps_filtered_data,"~/Documents/phd_study/iNat_project/data/leps_us_filtered_data.csv")

#### Filter for the Washington State ----
leps_filtered_data<-read.csv("~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_us_filtered_data.csv")
leps_data_washington<-leps_filtered_data%>%
  filter(stateProvince=="Washington")
write.csv(leps_data_washington, "~/Documents/phd_study/iNat_project/data/leps_data_washington.csv")

#### Filter for California ----
leps_data_cali<-leps_filtered_data%>%
  filter(stateProvince=="California")
write.csv(leps_data_cali, "~/Documents/phd_study/iNat_project/data/leps_data_cali.csv")

#### Filter for Arizona ----
leps_data_arizona<-leps_filtered_data%>%
  filter(stateProvince=="Arizona")
write.csv(leps_data_arizona, "~/Documents/phd_study/iNat_project/data/leps_data_arizona.csv")

#### Filter for Colorado ----
leps_data_colorado<-leps_filtered_data%>%
  filter(stateProvince=="Colorado")
write.csv(leps_data_colorado, "~/Documents/phd_study/iNat_project/data/leps_data_colorado.csv")

#### Filter for Texas ----
leps_data_texas<-leps_filtered_data%>%
  filter(stateProvince=="Texas")
write.csv(leps_data_texas, "~/Documents/phd_study/iNat_project/data/leps_data_texas.csv")

#### Filter for Iowa ----
leps_data_iowa<-leps_filtered_data%>%
  filter(stateProvince=="Iowa")
write.csv(leps_data_iowa, "~/Documents/phd_study/iNat_project/data/leps_data_iowa.csv")

#### Filter for Minnesota ----
leps_data_minnesota<-leps_filtered_data%>%
  filter(stateProvince=="Minnesota")
write.csv(leps_data_minnesota, "~/Documents/phd_study/iNat_project/data/leps_data_minnesota.csv")

#### Filter for Missouri ----
leps_data_missouri<-leps_filtered_data%>%
  filter(stateProvince=="Missouri")
write.csv(leps_data_missouri, "~/Documents/phd_study/iNat_project/data/leps_data_missouri.csv")

#### Filter for Illinois ----
leps_data_illinois<-leps_filtered_data%>%
  filter(stateProvince=="Illinois")
write.csv(leps_data_illinois, "~/Documents/phd_study/iNat_project/data/leps_data_illinois.csv")

#### Filter for Georgia ----
leps_data_georgia<-leps_filtered_data%>%
  filter(stateProvince=="Georgia")
write.csv(leps_data_georgia, "~/Documents/phd_study/iNat_project/data/leps_data_georgia.csv")

#### Filter for Michigan ----
leps_data_michigan<-leps_filtered_data%>%
  filter(stateProvince=="Michigan")
write.csv(leps_data_michigan, "~/Documents/phd_study/iNat_project/data/leps_data_michigan.csv")

#### Filter for Florida ----
leps_data_florida<-leps_filtered_data%>%
  filter(stateProvince=="Florida")
write.csv(leps_data_florida, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_florida.csv")

#### Filter for North Carolina ----
leps_data_nc<-leps_filtered_data%>%
  filter(stateProvince=="North Carolina")
write.csv(leps_data_nc, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_north_carolina.csv")

#### Filter for Washingon DC ----
leps_data_dc<-leps_filtered_data%>%
  filter(stateProvince=="District of Columbia")
write.csv(leps_data_dc, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_dc.csv")

#### Filter for Pennsylvania ----
leps_data_pa<-leps_filtered_data%>%
  filter(stateProvince=="Pennsylvania")
write.csv(leps_data_pa, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_pennsylvania.csv")

#### Filter for New York ----
leps_data_ny<-leps_filtered_data%>%
  filter(stateProvince=="New York")
write.csv(leps_data_ny, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_new_york.csv")

#### Filter for Massachusetts ----
leps_data_ma<-leps_filtered_data%>%
  filter(stateProvince=="Massachusetts")
write.csv(leps_data_ma, "~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_massachusetts.csv")





















