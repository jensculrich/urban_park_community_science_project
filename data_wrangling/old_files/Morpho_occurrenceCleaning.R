## NCEAS Morpho Initiative Working Group: "Plant Prioritizer" 
## Full Occurrence Data Cleaning Pipeline
## Compiled by Chris Cosma
## Started: August 8, 2024
## Last updated: November 26, 2024

# Adapted primarily from BeeBDC Package:
# https://cran.r-project.org/web/packages/BeeBDC/vignettes/BeeBDC_main.html

# Clear workspace
rm(list = ls(all = TRUE))

#wd
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow')

#### 1: Raw Occurrence Data ####

# Load packages
library(tidyverse)

# Sources for each taxa
# Plants: Calflora
# Insects: GBIF, SCAN, Cheshire (bees)
# Hummingbirds: GBIF 

# Download criteria
# Bounding Box (WGS 84): 42.01 (N), 32.53 (S), -124.42 (W), -114.13 (E)
# BBOX is used for initial download, further spatial cleaning below

#### 1.1: Plants ####

#### CalFlora ####

# Go to https://www.calflora.org/entry/wsearch.html 
# Need a CalFlora Account to use the downloader
# Change "Native Status" to  "Native"
# Change "Download Format" to "CSV"
# Change "Geometry" to "point"
# Click on each County from the list one at a time (download limit of 500,000 records is reached when doing "Any") then click "Search", then click "Download". Skip "- Bay Area -"
# Should be 58 total files
# Last downloaded Nov 19, 2024

# Read in each file (for the most part, one file = one county) and combine into a single df
calflora_files <- list.files(path = '/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/calflora', pattern = "*.csv", full.names = TRUE)
calflora_dataframes <- lapply(calflora_files, read.csv)
calflora_raw <- do.call(rbind, calflora_dataframes)

#### CCH2 #### 

# Download occurrence from https://www.cch2.org/ using bounding box, georeferenced only, exclude cultivated/captive
# 1 mil download limit, have to do it in batches
# last downloaded Nov 20, 2024

# Define the main folder path
main_folder <- "/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/cch2"

# Find all "occurrences.csv" files in the folder and its subfolders
file_paths <- list.files(main_folder, pattern = "occurrences\\.csv$", recursive = TRUE, full.names = TRUE)

# Read all CSV files and combine them into a single data frame

#below was taking too long, had to do one at a time
#cch2_raw <- do.call(rbind, lapply(file_paths, read.csv))

file_paths

cch2_1 = read.csv("/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/cch2/SymbOutput_2024-11-20_174713_DwC-A/occurrences.csv")

cch2_2 = read.csv("/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/cch2/SymbOutput_2024-11-20_181116_DwC-A/occurrences.csv")

cch2_3 = read.csv("/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/cch2/SymbOutput_2024-11-20_182546_DwC-A/occurrences.csv")

## I think every row of data in the CCH2 datasets are given a unique ID, so even if there were duplicates between the 3 datasets (which there are), you first need to get rid of that unique ID column to remove them. Skipping for now
cch2_raw = rbind(cch2_1, cch2_2, cch2_3)

#### 1.2 Pollinators ####

#### GBIF ####

####Option 1: Manual download#####
#download from https://www.gbif.org/occurrence/search 
#called the zip folder "gbif_manual.zip"
#Last downloaded Nov 20, 2024

#Have to use data.table to read it in, read.csv throws errors

library(data.table)
gbif_raw1 = data.table::fread('/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/gbif/0024811-241107131044228.csv')

#Cite: GBIF.org (21 November 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.gsdubr

#### Option 2: Auto-download ####

# Download using rgbif
#install.packages("rgbif")
library(rgbif)

#Supply GBIF user credentials
gbif_creds <- list(
  user = "ccosm001",
  pwd = "Wyl2smb?",
  email = "ccosm001@ucr.edu"
)

# make list of GBIF taxon keys for each desired pollinator group
name_backbone(name = "Lepidoptera") #797
name_backbone(name = "Syrphidae") #6920
name_backbone(name = "Trochilidae")#5289
taxon_keys <- c(797,6920,5289)

#Generate WKT (geometry) for CA. rgbif needs WKT format
library(sf)
library(tigris)

# Get map sf object for CA
california_sf <- states(cb = TRUE, resolution = "500k") %>%
  filter(STUSPS == "CA")

#fix invalid geometries
california_sf <- st_make_valid(california_sf)

# Convert sf to WKT
california_wkt <- st_as_text(st_geometry(california_sf))

# Initiate a download request in California
# This takes a very long time (~4 hours). Much faster to download manually from GBIF
# Last Downloaded Nov 19, 2024
key <- occ_download(
  pred_in("taxonKey",taxon_keys),  # Filter for desired groups
  pred("hasCoordinate", TRUE),        # Include only records with coordinates
  pred ("geometry", california_wkt),          # Filter to California
  format = "SIMPLE_CSV",     # Request simple CSV format
  user = gbif_creds$user,
  pwd = gbif_creds$pwd,
  email = gbif_creds$email               
)

# See status of download
occ_download_wait(key)

# Retrieve the download file and import to R (takes a few minutes; saves a zipped copy in the indicated folder path)
gbif_raw2 <- occ_download_get(key, path = "/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/gbif") %>%
  occ_download_import()

# Or unzip the file manually and read in in
gbif_raw2 = fread('/Users/chriscosma/Desktop/Morpho Final Workflow/input/occurrence/gbif/0022445-241107131044228.csv')

# For now, combine auto and manual and remove duplicates
gbif_raw = unique(rbind(gbif_raw1, gbif_raw2))

# Fix common issue with GBIF data
gbif_raw$scientificName <- ifelse(grepl("^BOLD", gbif_raw$scientificName), gbif_raw$species, gbif_raw$scientificName)

#### SCAN ####

# Download data from https://scan-bugs.org/
# Set to Occurrences only
# Search Higher Taxonomy -> Insecta
# Bounding Box: 42.01 (N), 32.53 (S), -124.42 (W), -114.13 (E)
scan_raw = read.csv("input/occurrence/scan/webreq_DwC-A/occurrences.csv")

#### MPG #####
#Records supplied by Steve Nanz of Moth Photographers' Group on Nov 13, 2024

mpg_raw = read.csv("input/occurrence/mpg/mpg_raw.csv")

#### Chesshire ####

# Download from https://figshare.com/projects/Completeness_analyses_for_over_3000_United_States_bee_species_identifies_persistent_data_gaps/138673 
chesshire_raw <- read.csv("input/occurrence/chesshire/contiguousRecords_high_Only.csv") 

#### Save Raw Occurrence Data #####
# Save everything again so it works in cleaning steps
setwd("~/Desktop/Morpho Final Workflow/intermediate/occurrence/formatted")
write.csv(gbif_raw, "gbif_formatted.csv", row.names = F)
write.csv(scan_raw, "scan_formatted.csv", row.names = F)
write.csv(chesshire_raw, "chesshire_formatted.csv", row.names = F)
write.csv(calflora_raw, "calflora_formatted.csv", row.names = F)
write.csv(cch2_raw, "cch2_formatted.csv", row.names = F)
write.csv(mpg_raw, "mpg_formatted.csv", row.names = F)

#### 2: Raw Occurrence Cleaning ####

# Clear workspace
rm(list = ls(all = TRUE))

# WD 
setwd("~/Desktop/Morpho Final Workflow/intermediate/occurrence/formatted")

# Load occurrence data
gbif = read.csv("gbif_formatted.csv", nrows = 100)
scan = read.csv("scan_formatted.csv", nrows = 100)
chesshire = read.csv("chesshire_formatted.csv", nrows = 100)
calflora = read.csv("calflora_formatted.csv", nrows = 100)
cch2 = read.csv("cch2_formatted.csv", nrows = 100)
mpg = read.csv("mpg_raw.csv", nrows = 100)

# Prepare bdc metadata file

# gbif, scan, chesshire, cch2 all in Darwincore
# Find matching column names between them
matching_columns <- t(data.frame(names(gbif)[names(gbif) %in% intersect(intersect(names(scan), names(cch2)), names(chesshire))]))

#save and manually prepare the file in excel
write.csv(matching_columns, "matching_columns.csv", row.names = F)

# Calflora and mpg are the weird ones, have to manually add the corresponding columns in
names(mpg)
names(calflora)

##### Dataset Harmonization #####

# Clear workspace
rm(list = ls(all = TRUE))

# WD 
setwd("~/Desktop/Morpho Final Workflow/intermediate/occurrence/formatted")

# Load packages
library(dplyr)
library(readr)
library(taxadb)
library(bdc)
library(sf)
library(BeeBDC)
library(lubridate)
library(tidyverse)

# Load the required metadata file
bdc_meta <- read_csv("bdc_meta.csv")

# Merge and standardize the datasets with bdc
database <-
  bdc_standardize_datasets(
    metadata = bdc_meta,
    format = "csv",
    overwrite = TRUE,
    save_database = FALSE) %>%
  #remove non UTF-8 characters
  mutate(scientificName = iconv(scientificName, from = "latin1", to = "UTF-8", sub = ""))

# Save a copy (still contains duplicated rows)
setwd("~/Desktop/Morpho Final Workflow/intermediate/occurrence/formatted")
write.csv(database, "database.csv", row.names = F)

#Below was a note from Teagan, but as of Nov 21, 2024, I don't think this is a problem anymore?

# !!!!!!!!!!!!!!!!!!!!!!!! NEED TO FIGURE OUT WHY CHESSHIRE DATA NOT READING IN FULLY (CERTAIN COLUMNS) !!!!!!!!!!!!!!!!!!!!!!!! 

#### Basic Cleaning ####

# Clear workspace
rm(list = ls(all = TRUE))

# WD 
setwd("~/Desktop/Morpho Final Workflow/intermediate/occurrence/formatted")

# Load packages
library(dplyr)
library(readr)
library(taxadb)
library(bdc)
library(sf)
library(BeeBDC)
library(lubridate)
library(tidyverse)

# Load database
database = read.csv("database.csv")

# Records missing species names
check_pf <-
  bdc_scientificName_empty(
    data = database,
    sci_name = "scientificName")

# Records lacking information on geographic coordinates
check_pf <- bdc_coordinates_empty(
  data = check_pf,
  lat = "decimalLatitude",
  lon = "decimalLongitude")

# Records with out-of-range coordinates
check_pf <- bdc_coordinates_outOfRange(
  data = check_pf,
  lat = "decimalLatitude",
  lon = "decimalLongitude")

# Records from poor sources
check_pf <- bdc_basisOfRecords_notStandard(
  data = check_pf,
  basisOfRecord = "basisOfRecord",
  names_to_keep = c("Event","HUMAN_OBSERVATION", "HumanObservation", 
                    "LIVING_SPECIMEN", "LivingSpecimen", "MACHINE_OBSERVATION", 
                    "MachineObservation", "MATERIAL_SAMPLE", "None", "O", 
                    "Occurrence", "MaterialSample", "OBSERVATION", 
                    "Pinned Specimen", "Photograph", "Preserved Specimen", 
                    "PRESERVED_SPECIMEN", "preservedspecimen Specimen", 
                    "Preservedspecimen", "PreservedSpecimen", 
                    "preservedspecimen", "S", "Specimen", "Taxon", 
                    "UNKNOWN", "", NA))

# Since all records in US, to see if any lat/lon transposed
# Ensure no negative numbers in latitude and no positives in longitude
check_pf$.coordinates_transposed = ifelse(check_pf$decimalLatitude<0, FALSE, TRUE)

#create California boundary shape
library(sf)
library(tigris)

# Get map sf object for CA
CA_boundary <- states(cb = TRUE, resolution = "500k") %>%
  filter(STUSPS == "CA") %>%
  st_transform(crs = st_crs(4326)) 

#fix invalid geometries
CA_boundary <- st_make_valid(CA_boundary)

plot(CA_boundary)

# Check if points are in California (region of interest)
check_pf <- st_as_sf(check_pf, 
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = st_crs(CA_boundary), 
                     remove = FALSE, 
                     na.fail = FALSE) %>%
  # Use st_intersects() to check for overlap between points and polygon
  # If no intersection, set .overlaps_California to FALSE
  mutate(.overlaps_California = lengths(st_intersects(., CA_boundary)) > 0) %>%
  # Remove geometry to turn back into standard dataframe
  st_drop_geometry()

# Create summary column
check_pf <- bdc_summary_col(data = check_pf)

# Create report
report <-
  bdc_create_report(data = check_pf,
                    database_id = "database_id",
                    workflow_step = "prefilter",
                    save_report = FALSE)
report

##### Filter out flagged records #####

# Save flagged records in a separate dataset
flagged_records <- check_pf %>%
  dplyr::filter(.summary == FALSE) %>%
  select(-.summary)

# Filter out all flagged records
merged_occurrences <-
  check_pf %>%
  dplyr::filter(.summary == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

##### Cleaning and Harmonizing Scientific Names #####

# Be sure to first install gnparser from https://github.com/gnames/gnparser
parse_names <-
  bdc_clean_names(sci_names = merged_occurrences$scientificName, save_outputs = FALSE)

# Examine errors
wrong_names <- filter(parse_names, is.na(names_clean) | quality == 0)

# You could spend a lot of time manually cleaning some of the remaining errors, but realistically can just accept that some will be wrong. More will also be corrected during taxonomic name harmonization

# Replace cleaned names in merged_occurrences file
parse_names <-
  parse_names %>%
  dplyr::select(.uncer_terms, names_clean)

merged_occurrences <- dplyr::bind_cols(merged_occurrences, parse_names)

names(merged_occurrences)

# # Add subspecies name to all darwincore style datasets
# merged_occurrences <- merged_occurrences %>%
#   mutate(num_parts_name = sapply(strsplit(trimws(names_clean), "\\s+"), length)) %>%
#   mutate(scientificNameUpdated = ifelse(taxonRank %in% c("SUBSPECIES", "VARIETY", "Subspec") & (num_parts_name == 2) & !(is.na(infraspecificEpithet)), paste(names_clean, infraspecificEpithet), names_clean)) %>%
#   select(-num_parts_name)

#Gnverifier
# Generate list of unique names for GNV
names <- data.frame(unique(merged_occurrences$names_clean))

# WD
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')

# Save file as tsv
write.table(names, "occurrence_names.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, fileEncoding = "UTF-8")

# Using global names verifier: https://verifier.globalnames.org/
# Make sure you have gnverifier downloaded: https://github.com/gnames/gnverifier?tab=readme-ov-file#as-a-restful-api
system('gnverifier "occurrence_names.tsv" > names_harmonized.tsv')

# Read in harmonized names
library(data.table)
names_harm <- fread("names_harmonized.tsv", sep = ",", header = TRUE, encoding = "UTF-8")

## Run CurrentName back through parser because it has authorship attached
parse_names <-
  bdc_clean_names(sci_names = names_harm$CurrentName, save_outputs = FALSE)

# Clean up harmonization outputs
names_harm <- names_harm %>%
  # Join back with harmonized name df
  mutate(names_harm = parse_names$names_clean) %>%
  # Replace ones with no CurrentName
  mutate(names_harm = ifelse(is.na(names_harm), MatchedCanonical, names_harm)) %>%
  # Replace any blanks with NA
  mutate(names_harm = ifelse(MatchType == "NoMatch", NA, names_harm)) %>%
  # Now remove some of the inaccuracies in the harmonized names, mostly just when we have subspecies
  mutate(names_harm = stringr::str_remove_all(names_harm, paste(c("subsp. ", "var. ", "f. "), collapse = "|")))

# Join back with other data higher taxonomy
merged_occurrences_harmonized <- left_join(merged_occurrences, names_harm, by = c("names_clean" = "ScientificName"))

##### Update Higher Taxonomy #####

names(merged_occurrences_harmonized)

# !!!!!!!!!!!!!!!!!!!!!!!! NEED TO CLEAN HIGHER TAXONOMY !!!!!!!!!!!!!!!!!!!!!!!! 

# Add a taxa column that lists what specific taxa we are dealing with
# examples: plants, leps, bees, hoverflies, hummingbirds, other
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  mutate(taxonGroup = case_when(
    grepl("calflora", database_id, ignore.case = TRUE)  ~ "plants",
    grepl("cch2", database_id, ignore.case = TRUE)  ~ "plants",
    grepl("trochilidae", family, ignore.case = TRUE)  ~ "hummingbirds",
    grepl("trochilidae", ClassificationPath, ignore.case = TRUE)  ~ "hummingbirds",
    grepl("chesshire", database_id, ignore.case = TRUE)  ~ "bees",
    grepl("Apidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Andrenidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Colletidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Halictidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Megachilidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Mellittidae", ClassificationPath, ignore.case = TRUE)  ~ "bees",
    grepl("Apidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Andrenidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Colletidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Halictidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Megachilidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Mellittidae", family, ignore.case = TRUE)  ~ "bees",
    grepl("Syrphidae", ClassificationPath, ignore.case = TRUE)  ~ "hoverflies",
    grepl("Syrphidae", family, ignore.case = TRUE)  ~ "hoverflies",
    grepl("Lepidoptera", ClassificationPath, ignore.case = TRUE)  ~ "leps",
    grepl("mpg", database_id, ignore.case = TRUE)  ~ "leps",
    grepl("Lepidoptera", order, ignore.case = TRUE)  ~ "leps"
  ))


##Add a column for butterfly or moth
butt_fams = c("Hedylidae", "Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae", "Riodinidae", "Danaidae")

merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  mutate(lepGroup = case_when(
    taxonGroup == "leps" & (
      grepl(paste(butt_fams, collapse = "|"), family, ignore.case = TRUE) |
        grepl(paste(butt_fams, collapse = "|"), ClassificationPath, ignore.case = TRUE)
    ) ~ "butterflies",
    taxonGroup == "leps" ~ "moths",  # Assign "moths" for leps that are not butterflies
    TRUE ~ NA_character_  # For rows not in taxonGroup "leps", assign NA or leave blank
  ))

# Add a column for bee family
bee_fams <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae", 
              "Megachilidae", "Melittidae", "Stenotritidae")  # Add or adjust as needed

merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  mutate(beeGroup = case_when(
    grepl(paste(bee_fams, collapse = "|"), family, ignore.case = TRUE) |
      grepl(paste(bee_fams, collapse = "|"), ClassificationPath, ignore.case = TRUE) ~ 
      str_extract(paste(family, ClassificationPath, sep = " "), 
                  paste(bee_fams, collapse = "|")),  # Extract the matched group
    TRUE ~ NA_character_  # Everything else gets NA
  ))


#examine errors
group_miss = filter(merged_occurrences_harmonized, is.na(taxonGroup))

lep_miss = filter(merged_occurrences_harmonized, taxonGroup == "leps" & is.na(lepGroup))

bee_miss = filter(merged_occurrences_harmonized, taxonGroup == "bees" & is.na(beeGroup))

#create a subset to view more easily
test <- merged_occurrences_harmonized %>%
  slice(seq(1, n(), by = 20))


# # Take a look at irregularities to see what is going on (as an example)
# leps_irregularities <- merged_occurrences_harmonized %>%
#   filter(order == "Lepidoptera" & !(grepl("Lepidoptera", ClassificationPath, ignore.case = TRUE)))
# # Seems like mostly it is just when ClassificationPath is empty due to source used. Not a big deal

##### Filter out flagged records and save checkpoint #####

# Add more flagged records to the flagged_records dataset
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  # Create flag column for taxon group. FALSE if taxon group not in scope of project
  mutate(.correct_taxonGroup = ifelse(!is.na(taxonGroup), TRUE, FALSE))

# Filter out records with flag
incorrect_group <- merged_occurrences_harmonized %>%
  filter(.correct_taxonGroup == FALSE)

# Add the records with incorrect taxonomy to flagged records
flagged_records <- bind_rows(flagged_records, incorrect_group) %>%
  select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove any occurrences that don't fall into one of our taxon groups of interest
# These are going to be arthropods that are not bees, leps, or hoverflies
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  filter(.correct_taxonGroup == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Save here as a csv as checkpoint

setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')

write.csv(merged_occurrences_harmonized, "merged_occurrences_harmonized_temp1.csv", row.names = F)
write.csv(flagged_records, "flagged_records_temp1.csv", row.names = F)

# Clean up working directory
rm(list = setdiff(ls(), c("merged_occurrences_harmonized", "flagged_records")))

#### Spatial Cleaning ######

#re-load files if needed
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')

merged_occurrences_harmonized = read.csv("merged_occurrences_harmonized_temp1.csv")

flagged_records = read.csv("flagged_records_temp1.csv")

#Re load packages if needed
library(dplyr)
library(readr)
library(taxadb)
library(bdc)
library(sf)
library(BeeBDC)
library(lubridate)
library(tidyverse)

# Make coordinates as type numeric (instead of character)
merged_occurrences_harmonized <- merged_occurrences_harmonized %>%
  # make coordinates numeric instead of character
  mutate(decimalLatitude = as.numeric(decimalLatitude),
         decimalLongitude = as.numeric(decimalLongitude))

# Check other common spatial issues
check_space <- BeeBDC::jbd_coordinates_precision(data = merged_occurrences_harmonized, 
                                                 lon = "decimalLongitude", 
                                                 lat = "decimalLatitude", 
                                                 ndec = 2) %>%
  # Rename flag column to make more explicit
  rename(.coordinates_Precise = .rou)
# ^ This flagged nearly 2 million records when ndec=2

# !!!!!!!!!!!!!!!!!!!!!!!! DIFFERENT PRECISION? !!!!!!!!!!!!!!!!!!!!!!!! 

# Force countrycode to US so it does not flag anything because of this below
# (lots of records have no country code due to dataset it is coming from)
check_space$countryCode = "US"

#subset the dataset to run in chunks

# Define chunk size
chunk_size <- 2e6  # 2 million rows

# Initialize final dataframe
check_space2 <- tibble()

# Total number of rows
total_rows <- nrow(check_space)

# Process in chunks
for (start_row in seq(1, total_rows, by = chunk_size)) {
  # Define the end row for the current chunk
  end_row <- min(start_row + chunk_size - 1, total_rows)
  
  # Subset the data for the current chunk
  chunk <- check_space[start_row:end_row, ]
  
  # Run the cleaning function on the current chunk
  chunk_cleaned <-
    CoordinateCleaner::clean_coordinates(
      x =  chunk,
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "names_harm",
      tests = c(
        "capitals",     # records within 0.5 km of capitals centroids
        "centroids",    # records within 1 km around country and province centroids
        "equal",      # records with equal coordinates
        "gbif",         # records within 1 km of GBIF headquarters. (says 1 degree in package, but code says 1000 m)
        "institutions", # records within 100m of zoo and herbaria
        "zeros"       # records with coordinates 0,0
        # "seas"        # Not flagged as this should be flagged by coordinate country inconsistent
      ),
      capitals_rad = 1000,
      centroids_rad = 500,
      centroids_detail = "both", # test both country and province centroids
      inst_rad = 100, # remove zoo and herbaria within 100m
      range_rad = 0,
      zeros_rad = 0.5,
      capitals_ref = NULL,
      centroids_ref = NULL,
      country_ref = NULL,
      country_refcol = "countryCode",
      inst_ref = NULL,
      range_ref = NULL,
      # seas_scale = 50,
      value = "spatialvalid" # result of tests are appended in separate columns
    ) %>%
    # Remove duplicate .summary column that can be replaced later and turn into a tibble
    dplyr::select(!tidyselect::starts_with(".summary")) %>%
    dplyr::tibble()
  
  # Append the cleaned chunk to the final dataframe
  check_space2 <- bind_rows(check_space2, chunk_cleaned)
  
  # Optional: print progress
  print(paste("Processed rows", start_row, "to", end_row))
}



##### Filter out flagged records and save checkpoint #####

# Summarize flags
check_space_flagSummary <- BeeBDC::summaryFun(data = check_space2, 
                                              dontFilterThese = c(".uncer_terms"), 
                                              removeFilterColumns = FALSE,
                                              filterClean = FALSE)

# Convert data type of coordinates back to character for consistency with other objects
check_space_flagSummary <- check_space_flagSummary %>%
  mutate(decimalLatitude = as.character(decimalLatitude),
         decimalLongitude = as.character(decimalLongitude))

flagged_records <- flagged_records %>%
  mutate(decimalLatitude = as.character(decimalLatitude),
         decimalLongitude = as.character(decimalLongitude))

# Filter out records with flag
spatial_issues <- check_space_flagSummary %>%
  filter(.summary == FALSE)

# Add the records with incorrect taxonomy to flagged records
flagged_records <- bind_rows(flagged_records, spatial_issues) %>%
  select(-starts_with("."), starts_with(".")) %>% # arrange the columns so flag columns are at the end!
  select(-.summary)

# Remove the flagged rows
space_clean <- check_space_flagSummary %>%
  filter(.summary == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Clean up working directory
rm(list = setdiff(ls(), c("space_clean", "flagged_records")))

# Save a copy as check point
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')
write.csv(space_clean, "space_clean.csv", row.names = F)
write.csv(flagged_records, "flagged_records_temp2.csv", row.names = F)

#### Temporal Cleaning #### 

# Reload data if needed 
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')
space_clean = read.csv("space_clean.csv")
flagged_records = read.csv("flagged_records_temp2.csv")

#Re load packages if needed
library(dplyr)
library(readr)
library(taxadb)
library(bdc)
library(sf)
library(BeeBDC)
library(lubridate)
library(tidyverse)

### I'm not sure why Teagan had the date column called "date"- mine was "eventDate". Maybe I did it wrong? Also have day, month, year columns

# Date cleaning - get date information from multiple columns and clean
occurrences_time <- space_clean %>%
  mutate( # First, add dates to clean columns using the date column, which has variable format depending on data source
    dateClean = case_when(
      grepl("\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z", eventDate) ~ as.Date(eventDate, format = "%Y-%m-%dT%H:%M:%OSZ"),
      grepl("\\d{4}-\\d{2}-\\d{2}", eventDate) ~ as.Date(eventDate, format = "%Y-%m-%d"),
      TRUE ~ as.Date(NA)
    ),
    yearClean = as.numeric(lubridate::year(dateClean)),
    monthClean = as.numeric(lubridate::month(dateClean)),
    dayClean = as.numeric(lubridate::day(dateClean))
  ) %>%
  mutate( # Then replace any NAs with values from the year/month/day column if values exist there
    yearClean = ifelse(is.na(yearClean) & !is.na(year) & !is.na(month) & !is.na(day), as.numeric(year), yearClean),
    monthClean = ifelse(is.na(monthClean) & !is.na(year) & !is.na(month) & !is.na(day), as.numeric(month), monthClean),
    dayClean = ifelse(is.na(dayClean) & !is.na(year) & !is.na(month) & !is.na(day), as.numeric(day), dayClean)
  ) %>%
  select(-dateClean) 

##### Dataset Cleaning #####

# Filter out based on cleaner dates
occurrences_time_cleaned <- occurrences_time %>%
  mutate(.valid_date = case_when(
    yearClean >= 1900 & yearClean <= 2024 & monthClean %in% 1:12 & dayClean %in% 1:31 ~ TRUE,
    TRUE ~ FALSE
  ))

# Filter out records with flag
temporal_issues <- occurrences_time_cleaned %>%
  filter(.valid_date == FALSE)

# Add the records with incorrect taxonomy to flagged records
flagged_records <- bind_rows(flagged_records, temporal_issues) %>%
  select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove the flagged rows
time_clean <- occurrences_time_cleaned %>%
  filter(.valid_date == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Clean up working directory
rm(list = setdiff(ls(), c("time_clean", "flagged_records")))

#####Save checkpoint#####

setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')
write.csv(time_clean, "time_clean.csv", row.names = F)
write.csv(flagged_records, "flagged_records_temp3.csv", row.names = F)

##### Identify Duplicates #####

# Reload data if needed 
setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaning')
time_clean = read.csv("time_clean.csv")
flagged_records = read.csv("flagged_records_temp3.csv")

#Re load packages if needed
library(dplyr)
library(readr)
library(taxadb)
library(bdc)
library(sf)
library(BeeBDC)
library(lubridate)
library(tidyverse)

# Good resource on identifying duplicates:
# https://discourse.gbif.org/t/duplicate-occurrence-records/3735

#3 general types

# exact duplicates have the same entries in all fields, including unique identifiers like occurrenceID.

# strict duplicates have the same entries in all fields except those which give the record a unique ID, like the occurrenceID field.

# relaxed duplicates have the same entries in key occurrence fields, such as taxon, coordinates, date and collector/observer.


# names(time_clean_cleaned)
# 
# #save a abbreviated copy to examine 
# test <- time_clean_cleaned %>%
#   slice(seq(1, n(), by = 100))
# write.csv(test, "test.csv", row.names = F)

#Here again Teagan had the ID column, but mine was called occurrenceID. There is also catalogNumber, recordNumber

#Chris's new duplciate code as of Nov 22, 2024

# Convert the dataset to a data.table object because of memory/computing limits

library(data.table)

time_clean <- as.data.table(time_clean)

# Determine the duplicate types

#Type 1: .dupExact = rows that are duplicates in everything (except "database_id", which is a column we added)
time_clean[, `:=`(
  .dupExact = !(
    duplicated(.SD, by = names(time_clean)[!names(time_clean) %in% "database_id"]) |
      duplicated(.SD, by = names(time_clean)[!names(time_clean) %in% "database_id"], fromLast = TRUE)
  ),
  
  #Type 2: .dupStrict = rows that are duplicates in everything except "database_id","occurrenceID", "catalogNumber", "recordNumber"
  .dupStrict = !(
    duplicated(.SD, by = names(time_clean)[!names(time_clean) %in% c("database_id", "occurrenceID", "catalogNumber", "recordNumber")]) |
      duplicated(.SD, by = names(time_clean)[!names(time_clean) %in% c("database_id", "occurrenceID", "catalogNumber", "recordNumber")], fromLast = TRUE)
  ),
  
  #Type 3: .dupRelaxed = rows that are duplciates based only on "names_harm", "decimalLatitude", "decimalLongitude", "yearClean", "monthClean", "dayClean"
  .dupRelaxed = !(
    duplicated(.SD, by = c("names_harm", "decimalLatitude", "decimalLongitude", "yearClean", "monthClean", "dayClean")) |
      duplicated(.SD, by = c("names_harm", "decimalLatitude", "decimalLongitude", "yearClean", "monthClean", "dayClean"), fromLast = TRUE)
  )
)]

## Teagan's old duplicate code
# # ID and remove duplicates
# duped_flagged <- time_clean_cleaned %>%
#   # are there duplicates based on ID column?
#   mutate(.duplicated = duplicated(select(., "occurrenceID")) | duplicated(select(., "occurrenceID"), fromLast = TRUE)) %>%
#   # Do not filter duplicates identified in Calflora data of chesshire data, only the duplicates from GBIF or SCAN
#   mutate(.duplicated = ifelse(!(.duplicated & (database_id %in% c("h_arthropoda_gbif", "h_arthropoda_scan", "h_trochilidae_gbif"))), TRUE, FALSE))

# !!!!!!!!!!!!!!!!!!!!!!!! UPDATE DE-DUPLICATION? !!!!!!!!!!!!!!!!!!!!!!!!

##### Filter out duplicates #####

# Filter out records with flag
duped_remove <- time_clean %>%
  filter(.dupExact == FALSE | .dupStrict == FALSE, .dupRelaxed == FALSE)



# Add the records with incorrect taxonomy to flagged records, with extra code to enable joining with different column classes
flagged_records <- bind_rows(flagged_records, duped_remove %>% mutate(across(intersect(names(flagged_records), names(.)), ~ as(., class(flagged_records[[cur_column()]]))))) %>%
  select(-starts_with("."), starts_with(".")) # arrange the columns so flag columns are at the end!

# Remove the flagged rows
clean_occurrences <- time_clean %>%
  filter(.dupExact == TRUE & .dupStrict == TRUE & .dupRelaxed == TRUE) %>%
  bdc_filter_out_flags(data = ., col_to_remove = "all")

# Clean up working directory
rm(list = setdiff(ls(), c("clean_occurrences", "flagged_records")))

#### Column filtering ####

#Teagan's old code
# # Only select what columns we need
# occurrences_clean_niceColumns <- clean_occurrences %>%
#   select(taxonGroup, decimalLatitude, decimalLongitude, names_harm, yearClean, monthClean, dayClean, sex, lifeStage, phenology, phenologyCode, occurrenceRemarks, DataSourceTitle, ClassificationPath, database_id) %>%
#   rename(scientificName = names_harm,
#          dataSourceTitle = DataSourceTitle,
#          classificationPath = ClassificationPath,
#          databaseID = database_id)

names(clean_occurrences)

## Chris's new code 11/25/2024

# Fix database_id column
clean_occurrences <- clean_occurrences %>%
  mutate(database_id = sub("_[^_]*$", "", database_id)) %>%
  mutate(database_id = sub("_formatted", "", database_id))

# Only select what columns we need
occurrences_clean_niceColumns <- clean_occurrences %>%
  select(database_id, taxonGroup, lepGroup, beeGroup, names_harm, decimalLatitude, decimalLongitude, yearClean, monthClean, dayClean) %>%
  rename(scientificName = names_harm,
         databaseID = database_id)


##### Save Outputs #####

setwd('/Users/chriscosma/Desktop/Morpho Final Workflow/intermediate/occurrence/cleaned')

# Save occurrences
write.csv(occurrences_clean_niceColumns, "occurrences_clean.csv", row.names = F)

# Save occurrences - all columns
write.csv(clean_occurrences, "occurrences_clean_allColumns.csv", row.names = F)

# Save flagged records
write.csv(flagged_records, "occurrences_flagged.csv", row.names = F)
