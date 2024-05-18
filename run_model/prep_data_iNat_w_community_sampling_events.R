# create a detection array from iNaturalist data.
# this will join detections to sites.
# non-detections are currently added whenever we don't see a species.
# sampling effort could then be modeled as a function of species/site/year/time,
# using the months within year as repeat "survey periods".
# but possibly we could look to identify "community sampling events"
# and only model non-detections as those species not observed on a given community sampling event?

library(tidyverse)
library(sf) # spatial data processing

prep_data <- function(min_species_detections,
                      min_park_size_acres,
                      max_park_size_acres,
                      buffer_distance,
                      min_species_for_community_sampling_event) {
  
  #-----------------------------------------------------
  # summarize the detection data
  
  # first read the data 
  df <- read.csv(
    "./data/all_inat_lep_records_2020-2023_los_angeles_county/all_inat_lep_records_2020-2023_los_angeles_county.csv")
          
  # and perform some initial filters
  df <- df %>%
    
    # for now, let's only look at detections with a species ID
    filter(str_count(species, "\\w+") == 2) %>%
    
    # for now, let's filter out geo obscured detections
    #filter(geoprivacy != "obscured") %>%
    
    # and additionally, let's get rid of any points with huge coordinate uncertainty
    # say > 1000m for now
    filter(coordinateUncertaintyInMeters <= 1000) %>%
    filter(!is.na(coordinateUncertaintyInMeters)) %>%
    
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(species) %>%
    add_tally() %>%
    filter(n >= min_species_detections) %>%
    dplyr::select(-n) %>%
    ungroup()
  
  #-----------------------------------------------------
  # map the data on a spatial file
  
  parks_shp <- read_sf("./data/los_angeles_county_parks_shapefile/Regional_Site_Inventory.shp")
  
  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070
  
  parks_shp <- parks_shp  %>% 
    st_transform(., crs) # USA_Contiguous_Albers_Equal_Area_Conic
  
  # for now, I filtered out some things that were clearly not urban parks 
  parks_shp <- parks_shp  %>% 
    filter(!str_detect(PARK_NAME, "Angeles National Forest")) %>% 
    filter(!str_detect(PARK_NAME, "Los Padres National Forest")) %>% 
    filter(!str_detect(PARK_NAME, "Edwards AFB")) %>%
    filter(!str_detect(PARK_NAME, "Hungry Valley")) %>%
    filter(!str_detect(PARK_NAME, "Air Force")) %>%
    filter(!str_detect(PARK_NAME, "State Recreation Area")) %>%
    filter(!str_detect(PARK_NAME, "National Recreation")) %>%
    filter(!str_detect(PARK_NAME, "State Park")) %>%
    filter(!str_detect(PARK_NAME, "State Beach")) %>%
    filter(!str_detect(PARK_NAME, "County Beach")) %>%
    filter(!str_detect(PARK_NAME, "Santa Catalina Island")) %>%
    
    # for now I also filtered out small parks just to speed up the estimation times (fewer sites)
    filter(RRE_ACRES > min_park_size_acres) %>%
    # and also the really really big ones 
    filter(RRE_ACRES < max_park_size_acres)
  
  # let's add a buffer around each park and then merge parks that are touching or overlapping
  parks_shp <- st_buffer(parks_shp, buffer_distance)
  
  # include the detection data on the map
  
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("decimalLongitude", "decimalLatitude"), 
                     crs = 4326))
  
  # and then transform it to the crs
  df <- st_transform(df_sf, crs = crs) %>%
    st_join(parks_shp, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(PARK_NAME)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 
  
  rm(parks_shp, df_sf)
  
  #-----------------------------------------------------
  # for now, let's also just include only grid cells with 
  # one or more detections. This may be important to correct 
  # later because it's quite likely people looked for  
  # butterflies in other places but just didn't see any.
  # could add in all the sites by rejoining the clipped grid?
  
  # with this approach we assume that any grid cell in which a detection has occurred
  # was sampled on every day in every year. I don't think this is a very realistic
  # assumption, but maybe an ok place to start.
  
  # how many species were detected?
  n_species <- nrow(species_names <- df %>%
         # group by species ID
         group_by(species) %>%
         add_tally() %>%
         # and take one record
         slice(1) %>%
         select(species, family, n))
  
  species_vector <- species_names %>%
    pull(species)
  
  n_detections <- nrow(df)
  
  n_sites <- (nrow(site_names <- df %>%
         group_by(PARK_NAME) %>%
         slice(1) %>%
         ungroup() %>%
         select(PARK_NAME)))
  
  site_vector <- site_names %>%
    pull(PARK_NAME)
  
  
  #-----------------------------------------------------
  # prep data for array format
  
  df <- df %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(factor(month)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, PARK_NAME, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, PARK_NAME, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  # get dimensions of surveys and years
  survey_vector <- as.vector(levels(as.factor(df$survey)))
  n_surveys <- length(survey_vector)
  
  year_vector <- as.vector(levels(as.factor(df$year)))
  n_years <- length(year_vector)
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V.
  
  # make a 4 dimensional array
  V <- array(data = NA, dim = c(n_species, n_sites, n_years, n_surveys))
  
  for(k in 1:n_years){
    for(l in 1:n_surveys){
      
      # iterate this across surveys within years
      temp <- df %>%
        
        # filter to indices for year and survey
        filter(year == (k), 
               survey == (l)) %>%  
        
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(., select(species_names, species), by="species") %>%
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        full_join(., select(site_names, PARK_NAME), by="PARK_NAME") %>%
        
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(PARK_NAME, row, fill = 0) %>%
        
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(4:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(4:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species %in% levels(as.factor(species_names$species)))
      
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      # remove species names
      temp_matrix <- temp_matrix[,-1]
      # remove the NA column at the end (why is this popping up?)
      last_column = ncol(temp_matrix)
      temp_matrix <- temp_matrix[,-(last_column)]
      
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
  # identify community sampling events which we will use to infer non-detections
  # identified by site*survey months with > 1 species detected. Another option
  # might be to group by iNat user name by date*park but this may be a very narrow scope?
  
  community_samples <- df %>%
    
    # determine whether a community sampling event occurred
    # using collector/observer name
    group_by(PARK_NAME, year, survey) %>%
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
    group_by(PARK_NAME, year, survey) %>%
    slice(1) %>%
    ungroup()

  
  print((paste0("prop community sampling events = ", 
                signif(sum(community_samples$community_sampled) / 
                         nrow(community_samples) * 100, 3),
                "%")))
  
  ## --------------------------------------------------
  # Generate an indicator array for whether or not a community-wide museum 
  # sampling event occurred at the site in a year
  
  # first make a df of all possible site visits * species visits
  all_species_site_visits <- as.data.frame(cbind(
    rep(species_vector, times=(n_sites*n_years*n_surveys)),
    rep(site_vector, each=n_species, times=n_years*n_surveys),
    rep(year_vector, each=n_species*n_sites, times=(n_surveys)),
    rep(survey_vector, each=n_species*n_sites*n_years)
  )) %>%
    rename("species" = "V1",
           "PARK_NAME" = "V2",
           "year" = "V3",
           "survey" = "V4") %>%
    mutate(year = as.integer(year),
           survey = as.integer(survey))
  
  # now propagate nondetections to any species 
  # this assumes that ALL other species could have been detected
  # could revise to propagate nondetection only to other species in same suborder, family, genus etc.
  all_visits_joined <- left_join(
      all_species_site_visits, community_samples, 
      by=c("PARK_NAME", "year", "survey")) %>%
    
    # create an indicator if the site visit was a sample or not (0 == not sampled, 1 == sampled)
    mutate(community_sampled = replace_na(community_sampled, 0)) 
  
  # for < min_species_for_community_sampling_event with a sampled indicator
  # i.e., keep the data for the singletons, while not inferring absence for the rest of the community
  df <- df %>%
    # if the species was sampled than at least that species was sampled (may also be a comm sample)
    mutate(non_comm_sample = 1) 
  
  # a species was sampled if either a community sample or a targeted sample event occurred
  all_visits_joined <- left_join(all_visits_joined, df) %>%
    # replace_na with zero, a targeted sampling event did not occur 
    mutate(non_comm_sample = replace_na(non_comm_sample, 0)) %>%
    # if a targeted sampling event occurred then sampling occurred,
    # if a targeted sampling event DID NOT occur, then look to see whether
    # or not a community sampling event occurred and then fill in accordingly
    mutate(any_sampled = ifelse(non_comm_sample == 1, 1, community_sampled)) 
  
  test <- all_visits_joined %>%
    mutate(any_sampled = +(!any_sampled))
  
  # now spread into 4 dimensions
  V_NA <- array(data = all_visits_joined$any_sampled, dim = c(n_species, n_sites, n_years, n_surveys))
  check <- which(V>V_NA) # this will give you numerical value
  # this check should be empty (can never detect a species where sampling has not occurred)
  
  if(length(check) > 0){
    print("WARNING: there are detections for instances marked as NA")
  } 
  
  print((paste0("prop community sampling events = ", 
                signif(sum(community_samples$community_sampled) / 
                         nrow(community_samples) * 100, 3),
                "%")))
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    V = V, # community science detection data
    V_NA = V_NA, # sampling event indicator array
    
    species_names = species_names, # species sci name, common name and number of detections
    n_species = n_species, # number of species
    n_sites = n_sites, # number of sites
    n_years = n_years, # number of surveys 
    n_surveys = n_surveys,
    
    species = species_vector,
    sites = site_vector,
    years = year_vector,
    surveys = survey_vector
    
  ))

}
