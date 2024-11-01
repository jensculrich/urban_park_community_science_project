# create a detection array from iNaturalist data.
# this will join detections to sites.
# non-detections are currently added whenever we don't see a species.
# sampling effort could then be modeled as a function of species/site/year/time,
# using the months within year as repeat "survey periods".
# but possibly we could look to identify "community sampling events"
# and only model non-detections as those species not observed on a given community sampling event?

library(tidyverse) # data organization
library(sf) # spatial data processing
library(geojsonsf) # geojson spatial data processing

prep_data <- function(min_species_detections,
                      min_species_per_family,
                      min_species_for_community_sampling_event, 
                      n_parks_sampled_per_city,
                      buffer_distance) {
  
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
    cbind(city="los_angeles", read.csv(
    "./data/02_wrangled_data/leps_data_cali.csv"))
    ) %>%
  filter(species != "") #%>%
  
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
  
  
  los_angeles_shp <- readRDS("./data/parks_by_city/parks_merged_LA.rds") %>%
    cbind(city = "los_angeles") %>%
    cbind(site = 1:nrow(los_angeles_shp))

  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070
  
  # los angeles
  los_angeles_shp <- los_angeles_shp  %>% 
    st_transform(., crs) 
  
  # view the shapefile
  ggplot() + 
    geom_sf(data = los_angeles_shp, fill = 'white', lwd = 0.05) +
    #geom_sf(data = parks2, fill = 'skyblue', alpha= 0.7, lwd = 0.3) +
    #geom_sf(data = filter(df_sf, city == "los_angeles"), aes(colour=as.factor(year))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle(paste0("Los Angeles: ", nrow(los_angeles_shp), 
                   ", following 100m buffer merge parks"))
  
  # and now match detections with spatial units (grid cells)
  df2 <- st_transform(df_sf, crs = crs) %>%
    st_join(los_angeles_shp, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(site)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, gbifID, decimalLatitude, decimalLongitude), by="gbifID")
  
  df3 <- df2 %>%
    # and remove families with very few species
    group_by(family) %>%
    mutate(n_species_per_family=n_distinct(species)) %>% 
    filter(n_species_per_family >= min_species_per_family) %>%
    ungroup() %>%
    
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(species) %>%
    mutate(n_binary_detections_per_species=n_distinct(site, year)) %>% 
    add_tally()
  
  unique(df2$family)  

  
  # make the detection data a spatial file for plotting
  (df_sf <- st_as_sf(df,
                     coords = c("longitude", "latitude"), 
                     crs = 4326))
  
  # view the shapefile
  ggplot() + 
    geom_sf(data = los_angeles_shp, fill = 'white', lwd = 0.05) +
    geom_sf(data = parks2, fill = 'skyblue', alpha= 0.7, lwd = 0.3) +
    geom_sf(data = df_sf, aes(colour=as.factor(year))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle(paste0("Los Angeles: ", n_parks_sampled_per_city, 
                   " randomly sampled parks\n w/ buffers of ", buffer_distance, "m"))
  #theme(legend.position="none") 
  
  rm(parks_shp, df_sf, los_angeles_shp)
  gc()
  
  #-----------------------------------------------------
  # prep data for array format
  # currently this only considers sites where 1 or more species ever detected
  
  # we will have one big array, with different lengths of species, sites and years
  
  df <- df %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(factor(month)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # scale the park size covariate
    mutate(park_size_scaled = center_scale(log(Park_Size_))) %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, family, ParkID, park_size_scaled, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, ParkID, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  # get dimensions of surveys and years
  survey_vector <- as.vector(levels(as.factor(df$survey)))
  n_surveys <- length(survey_vector)
  
  year_vector <- as.vector(levels(as.factor(df$year)))
  n_years <- length(year_vector)
  
  # how many species were detected?
  # and how many binary detections occurred? species/site/year/visit
  n_species <- nrow(species_names <- df %>%
                      # group by species ID
                      group_by(species) %>%
                      add_tally() %>%
                      # and take one record
                      slice(1) %>%
                      select(species, family, n))
  
  species_vector <- species_names %>%
    pull(species)
  
  # get vector of families for each species
  n_families <- nrow(family_vector <- pull(df %>%
                                             # group by species ID
                                             group_by(species) %>%
                                             # grab one row per species ID
                                             slice(1) %>%
                                             # and then pull out it's family
                                             select(family)))
  
  n_detections <- nrow(df)
  
  # n sites and site vector
  n_sites <- (nrow(site_names <- df %>%
                     group_by(ParkID) %>%
                     slice(1) %>%
                     ungroup() %>%
                     select(ParkID)))
  
  site_vector <- site_names %>%
    pull(ParkID)
  
  park_size_scaled <- df %>%
    group_by(ParkID) %>%
    slice(1) %>%
    ungroup() %>%
    select(park_size_scaled) %>%
    pull(park_size_scaled)
  
  ## --------------------------------------------------
  ## get species trait data
  
  # https://www.sciencebase.gov/catalog/item/62d59ae5d34e87fffb2dda99
  #invasive_df <- read.csv("./data/lepidoptera_trait_data/USRIISv2_invasive_species_list_lepidoptera.csv") %>%
    #select(scientificName) %>%
    #rename("species" = "scientificName") %>%
    #mutate(invasive = 1) 
  
  #df <- left_join(df, invasive_df) %>%
    #mutate(invasive = replace_na(invasive, 0))
  
  #n_invasive = df %>%
    #group_by(species) %>%
    #slice(1) %>%
    #ungroup() %>%
    #filter(invasive == 1) %>%
    #nrow()
  
  #invasiveness <- df %>%
    #group_by(species) %>%
    #slice(1) %>%
    #ungroup() %>%
    #select(invasive) %>%
    #pull(invasive)
  
  #species_names <- cbind(species_names, invasiveness) %>%
    #rename("invasiveness" = "...4")
  
  # Chowdhury 2021 Migratory Status
  migratory_df <- read.csv("./data/lepidoptera_trait_data/migratory.csv") %>%
    select(species, migratory) 
  
  df <- left_join(df, migratory_df) %>%
    mutate(migratory = replace_na(migratory, 0))
  
  n_invasive = df %>%
    group_by(species) %>%
    slice(1) %>%
    ungroup() %>%
    filter(migratory == 1) %>%
    nrow()
  
  migratory <- df %>%
    group_by(species) %>%
    slice(1) %>%
    ungroup() %>%
    select(migratory) %>%
    pull(migratory)
  
  species_names <- cbind(species_names, migratory) %>%
    rename("migratory" = "...4")
  
  ## --------------------------------------------------
  ## Now we are ready to create the detection matrix, V.
  
  # make a 4 dimensional array
  V <- array(data = NA, dim = c(n_species, n_sites, n_years, n_surveys))
  
  for(k in 1:n_years){
    for(l in 1:n_surveys){
      
      # iterate this across surveys within years
      temp <- df %>%
        
        # don't need the family column for this
        select(-family) %>%
        
        # filter to indices for year and survey
        filter(year == (k), 
               survey == (l)) %>%  
        
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(., select(species_names, species), by="species") %>%
        select(-park_size_scaled) %>% # if you don't deselect the variable it expands the temp matrix for every unique one
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        #full_join(., select(site_names, grid_id), by="grid_id") %>%
        full_join(., as.data.frame(site_vector), by=c("ParkID" = "site_vector")) %>% 
        
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(ParkID, row, fill = 0) %>%
        
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(5:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(5:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species %in% levels(as.factor(species_names$species)))
      
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
  # event had occurred.
  
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df_original,
                     coords = c("longitude", "latitude"), 
                     crs = 4326))
  
  # and now match detections with spatial units (grid cells)
  df_full <- st_transform(df_sf, crs = crs) %>%
    st_join(parks2, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(ParkID)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df_original, gbifID, latitude, longitude), by="gbifID") %>%
    
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
    dplyr::select(species, family, ParkID, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, ParkID, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  community_samples <- df_full %>%
    
    # determine whether a community sampling event occurred
    # using collector/observer name
    group_by(family, ParkID, year, survey) %>%
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
    group_by(ParkID, family, year, survey) %>%
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
           "ParkID" = "V2",
           "year" = "V3",
           "survey" = "V4",
           "family" = "V5") %>%
    mutate(year = as.integer(year),
           survey = as.integer(survey))
  
  # now propagate nondetections to any species 
  # this assumes that ALL other species could have been detected
  # could revise to propagate nondetection only to other species in same suborder, family, genus etc.
  all_visits_joined <- left_join(
    all_species_site_visits, community_samples, 
    by=c("ParkID","family", "year", "survey")) %>%
    
    # create an indicator if the site visit was a sample or not (0 == not sampled, 1 == sampled)
    mutate(community_sampled = replace_na(community_sampled, 0)) 
  
  # for < min_species_for_community_sampling_event with a sampled indicator
  # by definition, we know that someone observed these species and so they were observed
  # even if we aren't considering it a community sampling event.
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
    
    species_names = species_names, # species sci name, family common name and number of detections
    n_species = n_species, # number of species
    n_families = n_families, # number of Lepidoptera families
    n_sites = n_sites, # number of sites
    n_years = n_years, # number of years 
    n_surveys = n_surveys, # [max] number of surveys within year
    #n_cities = n_cities, # number of cities
    
    species = species_vector,
    families = family_vector, 
    sites = site_vector,
    years = year_vector,
    surveys = survey_vector,
    #cities = city_vector,
    
    # covariate data
    migratory = migratory,
    park_size_scaled = park_size_scaled
    
    #species_city_occurrences = species_city_occurrences

  ))

}
