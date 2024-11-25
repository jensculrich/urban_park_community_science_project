# create a detection array from iNaturalist data.
# this will join detections to sites.
# non-detections are currently added whenever we don't see a species.
# sampling effort could then be modeled as a function of species/site/year/time,
# using the months within year as repeat "survey periods".
# but possibly we could look to identify "community sampling events"
# and only model non-detections as those species not observed on a given community sampling event?

library(tidyverse) # data organization
library(sf) # spatial data processing
library(wesanderson) # colour palettes

prep_data <- function(min_species_detections,
                      min_species_per_family,
                      min_species_for_community_sampling_event, 
                      city_name) {
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  #-----------------------------------------------------
  # load the detection data
  
  butterfly_families <- c("Hesperiidae", "Papilionidae", "Pieridae", "Lycaenidae",
                          "Riodinidae", "Nymphalidae")
  
  # first read the data 
  df <- rbind(
    cbind(city=city_name, read.csv(
    "./data/wrangled_detection_data/leps_data_cali.csv"))
    ) %>%
    
  # use verbatim names for Riodinidae (GBIF does not accept them)
  mutate(species = ifelse(family == "Riodinidae", 
                          word(verbatimScientificName, 1, 2), 
                          species)) %>%
    
  # remove records with no species name
  filter(species != "") %>%
  
  # filter down to the butterfly families
  filter(family %in% butterfly_families)
    
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326))
  
  #-----------------------------------------------------
  # map the data onto a spatial file
  
  #-----------------------------------------------------
  # map the data on a spatial file
  
  
  city_shp <- readRDS(paste0(
    "./data/merged_parks_by_city/parks_merged_", city_name, ".rds")) %>%
    cbind(city = city_name) %>%
    cbind(site = 1:nrow(.)) 

  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070
  
  # los angeles
  city_shp <- city_shp  %>% 
    st_transform(., crs) %>%
    mutate(area = st_area(.)) %>%
    mutate(park_size_scaled = as.numeric(word(center_scale(log(area)), 1)))
  
  # view the shapefile
  #sub <- city_shp %>%
    #mutate(big = as.factor(ifelse(site== "185", 1, 0)))
  ggplot() + 
    #geom_sf(data = sub, aes(fill = big), lwd = 0.05) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle(paste0(city_name, ": ", nrow(city_shp), 
                   ", following 100m buffer merge parks"))
  
  # and now match detections with spatial units (grid cells)
  df2 <- st_transform(df_sf, crs = crs) %>%
    st_join(city_shp, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(site)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, gbifID, decimalLatitude, decimalLongitude), by="gbifID")
  
  df2 <- df2 %>%
    # and remove families with very few species
    group_by(family) %>%
    mutate(n_species_per_family=n_distinct(species)) %>% 
    filter(n_species_per_family >= min_species_per_family) %>%
    ungroup() %>%
    
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(species) %>%
    mutate(n_binary_detections_per_species=n_distinct(site, year)) %>% 
    add_tally() %>%
    ungroup()
  
  unique(df2$family)  

  # make the detection data a spatial file for plotting
  (df_sf <- st_as_sf(df2,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326))
  
  # view the shapefile
  mypalette <- as.vector(wesanderson::wes_palette("FantasticFox1"))
  
  ggplot() + 
    geom_sf(data = city_shp, fill = 'grey80', lwd = 0.05) +
    geom_sf(data = df_sf, aes(colour=as.factor(family))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    scale_colour_manual(name = "Family", values = mypalette) +
    theme_void() +
    ggtitle(paste0(
      "butterfly detections in ", city_name, " parks (2020 - 2024)"))
  
  gc()
  
  genus <- unique(df2$genus)
  
  #write.csv(genus, "data/lepidoptera_trait_data/ease_of_id/identifiability_by_genus.csv")
  
  #-----------------------------------------------------
  # prep data for array format
  # currently this only considers sites where 1 or more species ever detected
  
  # we will have one big array, with different lengths of species, sites and years
  
  df2 <- df2 %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(factor(month)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, family, site, park_size_scaled, survey, year,
                  n_species_per_family, n_binary_detections_per_species, n) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, site, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  # get dimensions of surveys and years
  survey_vector <- as.vector(levels(as.factor(df2$survey)))
  n_surveys <- length(survey_vector)
  
  year_vector <- as.vector(levels(as.factor(df2$year)))
  n_years <- length(year_vector)
  
  # how many species were detected?
  # and how many binary detections occurred? species/site/year/visit
  n_species <- length(species_vector <- pull(df2 %>%
                      # group by species ID
                      group_by(species) %>%
                      # and take one record
                      slice(1) %>%
                      select(species)))
  
  # get vector of families for each species
  n_families <- unique(family_vector <- pull(df2 %>%
                                             # group by species ID
                                             group_by(species) %>%
                                             # grab one row per species ID
                                             slice(1) %>%
                                             # and then pull out it's family
                                             select(family)))
  
  n_detections <- nrow(df2)
  
  # start a "species_info" table
  species_info <- df2  %>%
    mutate(genus = word(species, 1)) %>%
    select(species, genus, family, n_species_per_family, n_binary_detections_per_species, n) %>%
    group_by(species) %>%
    slice(1) %>%
    rename("n_total_detections_species" = "n") %>%
    ungroup()
  
  # n sites and site vector
  n_sites <- (nrow(site_data <- city_shp %>%
                     group_by(site) %>%
                     slice(1) %>%
                     ungroup() %>%
                     select(site, park_size_scaled)))
  
  site_vector <- site_data %>%
    pull(site)
  
  ## --------------------------------------------------
  ## get species trait data
  
  # ease of identification
  trait_df <- read.csv(
    "./data/lepidoptera_trait_data/butterfly_engagement/SpeciesListForTraits.csv") %>%
    select(GBIF_species, eButterfly_species, featureDiversity, aveWingspan) 
  
  ## separate out all possible names for the species
  # first figure out how many possible names
  trait_df <- trait_df %>%
    mutate(words = NULL) 
  for(i in 1:nrow(trait_df)){
    trait_df$words[i] <- length(str_split(trait_df$GBIF_species, "\\s+")[[i]])
  }
  max_names <- max(trait_df$words) / 2
  
  # make a column for each possible name 1:max_names
  name_cols <- as.character(vector(length = max_names))
  for(i in 1:max_names){
    name_cols[i] <- paste0("GBIF_species", as.character(i))
  }
  # Split name column into GBIF names
  trait_df[name_cols] <- str_split_fixed(
    trait_df$GBIF_species, ', ', max_names)
  
  # now get the trait data by genus incase anyone is missing for any name
  genus_id <- trait_df %>%
    mutate(genus = word(GBIF_species1, 1)) %>%
    group_by(genus) %>%
    mutate(featureDiversityGenus = mean(featureDiversity),
           aveWingspanGenus = mean(aveWingspan, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    select(genus, featureDiversityGenus, aveWingspanGenus)
  
  # now get the trait data by family incase anyone is missing for any name
  family_id <- trait_df %>%
    left_join(., select(df2, species, family), by=c("GBIF_species" = "species")) %>%
    group_by(family) %>%
    mutate(featureDiversityFamily = mean(featureDiversity),
           aveWingspanFamily = mean(aveWingspan, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    select(family, featureDiversityFamily, aveWingspanFamily)
  
  species_info <- left_join(species_info, trait_df, by=c("species" = "GBIF_species1")) %>%
    left_join(., genus_id) %>%
    mutate(featureDiversity = ifelse(is.na(featureDiversity), featureDiversityGenus, featureDiversity),
           aveWingspan = ifelse(is.na(aveWingspan), aveWingspanGenus, aveWingspan)) %>%
    left_join(., family_id) %>%
    mutate(featureDiversity = ifelse(is.na(featureDiversity), featureDiversityFamily, featureDiversity),
           aveWingspan = ifelse(is.na(aveWingspan), aveWingspanFamily, aveWingspan))
  species_info[species_info == ""] <- NA
  
  ## now add ease of identification info
  ease_id_df <- read.csv("./data/lepidoptera_trait_data/ease_of_id/identifiability_by_genus.csv") %>%
    select(genus, research_grade_proportion)
  
  species_info <- left_join(species_info, ease_id_df)
  
  # now scale and select all the variables 
  species_info <- species_info %>%
    mutate(aveWingspan_scaled = center_scale(aveWingspan),
           research_grade_proportion_scaled = center_scale(research_grade_proportion),
           featureDiversity_scaled = center_scale(featureDiversity)) %>%
    select(species, genus, family, n_binary_detections_per_species, n_total_detections_species,
           aveWingspan_scaled, featureDiversity_scaled, research_grade_proportion_scaled)
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V.
  
  # make a 4 dimensional array
  V <- array(data = NA, dim = c(n_species, n_sites, n_years, n_surveys))
  
  for(k in 1:n_years){
    for(l in 1:n_surveys){
      
      # iterate this across surveys within years
      temp <- df2 %>%
        
        # don't need the family column for this
        dplyr::select(-family) %>%
        
        # filter to indices for year and survey
        filter(year == (k), 
               survey == (l)) %>%  
        
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(., dplyr::select(species_info, species), by="species") %>%
        dplyr::select(-park_size_scaled) %>% # if you don't deselect the variable it expands the temp matrix for every unique one
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        #full_join(., select(site_names, grid_id), by="grid_id") %>%
        full_join(., as.data.frame(site_vector), by=c("site" = "site_vector")) %>% 
        
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(site, row, fill = 0) %>%
        
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(7:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(7:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species %in% levels(as.factor(species_info$species)))
      
      # remove the NA column at the end (why is this popping up? happens if no species detected at any sites)
      #last_column = ncol(temp_matrix)
      #temp_matrix <- temp_matrix[,-(last_column)]
      temp <- as.data.frame(temp)
      temp <- select(temp, -"<NA>")
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      
      # replace NAs for the interval i and visit j with the matrix
      V[1:n_species, 1:n_sites,k,l] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V) <- "numeric"
  
  print((paste0("prop detections = ", 
           signif(sum(V[1:n_species, 1:n_sites, 1:n_years, 1:n_surveys]) / 
             length(V) * 100, 3),
           "%")))
  

  #-----------------------------------------------------
  # identify whether a city is "in range" of each species
  # ie. whether the species was found in the city one or more times
  
  # to do
  
  #-----------------------------------------------------
  # identify community sampling events which we will use to infer non-detections
  # identified by site*survey months with > 1 species detected. Another option
  # might be to group by iNat user name by date*park but this may be a very narrow scope?
  
  # -----
  # first we need to regather the data but not filter out species below min detections
  # we want all species from the families considered when determining whether a sampling
  # event had occurred. So we have to reload the original df
  
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326))
  
  # and now match detections with spatial units (grid cells)
  df_full <- st_transform(df_sf, crs = crs) %>%
    st_join(city_shp, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(site)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") %>%
    
    # to confirm: did Jenny already remove samples without detailed coordinate uncertainty?
    
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(species) %>%
    add_tally() %>%
    rename("species_detections" = "n") %>%
    ungroup()
  
  df_full <- df_full %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(factor(month)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, family, site, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, site, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  community_samples <- df_full %>%
    
    # determine whether a community sampling event occurred
    # using collector/observer name
    group_by(family, site, year, survey) %>%
    mutate(n_species_sampled = n_distinct(species)) %>%
    ungroup() %>%
    
    #group_by(PARK_, occ_interval, occ_year, visit) %>%
    # slice max in case there are multiple institutions collecting from a site in a year
    #slice_max(n_species_sampled) %>%
    # then take one per year
    #slice(1) %>%
    
    # if we detected more than min_species_for_community_sampling_event at the site
    # in the "survey" then we assume a community sampling event occured
    mutate(community_sampled = ifelse(
      n_species_sampled >= min_species_for_community_sampling_event,
      1, 0)) %>%
    dplyr::select(-species, -n_species_sampled) %>%
    
    # we just need one row per community sampling event (not a row for all species positively detected)
    group_by(site, family, year, survey) %>%
    slice(1) %>%
    ungroup()
  
  ## --------------------------------------------------
  # Generate an indicator array for whether or not a community-wide museum 
  # sampling event occurred at the site in a year
  
  # first make a df of all possible site visits * species visits
  all_species_site_visits <- as.data.frame(cbind(
    rep(species_vector, times=(n_sites*n_years*n_surveys)),
    rep(site_vector, each=n_species, times=n_years*n_surveys),
    rep(year_vector, each=n_species*n_sites, times=(n_surveys)),
    rep(survey_vector, each=n_species*n_sites*n_years),
    rep(family_vector, times=(n_sites*n_years*n_surveys))
  )) %>%
    rename("species" = "V1",
           "site" = "V2",
           "year" = "V3",
           "survey" = "V4",
           "family" = "V5") %>%
    mutate(year = as.integer(year),
           survey = as.integer(survey),
           site = as.integer(site))
  
  # now propagate nondetections to any species 
  # this assumes that ALL other species could have been detected
  # could revise to propagate nondetection only to other species in same suborder, family, genus etc.
  all_visits_joined <- left_join(
    all_species_site_visits, community_samples, 
    by=c("site","family", "year", "survey")) %>%
    
    # create an indicator if the site visit was a sample or not (0 == not sampled, 1 == sampled)
    mutate(community_sampled = replace_na(community_sampled, 0)) 
  
  # for < min_species_for_community_sampling_event with a sampled indicator
  # by definition, we know that someone observed these species and so they were observed
  # even if we aren't considering it a community sampling event.
  # i.e., keep the data for the singletons, while not inferring absence for the rest of the community
  df2 <- df2 %>%
    # if the species was sampled than at least that species was sampled (may also be a comm sample)
    mutate(non_comm_sample = 1) 
  
  # a species was sampled if either a community sample or a targeted sample event occurred
  all_visits_joined <- left_join(all_visits_joined, df2) %>%
    # replace_na with zero, a targeted sampling event did not occur 
    mutate(non_comm_sample = replace_na(non_comm_sample, 0)) %>%
    # if a targeted sampling event occurred then sampling occurred,
    # if a targeted sampling event DID NOT occur, then look to see whether
    # or not a community sampling event occurred and then fill in accordingly
    mutate(any_sampled = ifelse(non_comm_sample == 1, 1, community_sampled)) 
  
  #test <- all_visits_joined %>%
  #mutate(any_sampled = +(!any_sampled))
  
  # now spread into 4 dimensions
  V_NA <- array(data = all_visits_joined$any_sampled, dim = c(n_species, n_sites, n_years, n_surveys))
  check <- which(V>V_NA) # this will give you numerical value
  # this check should be empty (can never detect a species where sampling has not occurred)
  
  if(length(check) > 0){
    print("WARNING: there are detections for instances marked as NA")
  } else{
    print("GOOD TO GO! detections are consistent with NA array")
  }
  
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    V = V, # community science detection data
    V_NA = V_NA, # sampling event indicator array
    
    species_info = species_info, # species sci name, family name, num detections and predictors
    site_data = site_data, # site names and predictors
    
    n_species = n_species, # number of species
    n_families = n_families, # number of Lepidoptera families
    n_sites = n_sites, # number of sites
    n_years = n_years, # number of years 
    n_surveys = n_surveys, # [max] number of surveys within year

    species = species_vector,
    families = family_vector, 
    sites = site_vector,
    years = year_vector,
    surveys = survey_vector

  ))

}
