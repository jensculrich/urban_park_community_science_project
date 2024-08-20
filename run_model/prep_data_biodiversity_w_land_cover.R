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
                      min_species_for_community_sampling_event) {
  
  #-----------------------------------------------------
  # load the detection data
  
  # first read the data 
  # cbind each city's data with the corresponding city name,
  # rbind the data from different cities
  df <- rbind(
    cbind(city="los_angeles", read.csv(
    "./data/biodiversity_data_with_land_classification/leps_la_data_with_land_classification.csv")),
    cbind(city="san_diego", read.csv(
      "./data/biodiversity_data_with_land_classification/leps_san_diego_data_with_land_classification.csv")),
    cbind(city="dallas", read.csv(
      "./data/biodiversity_data_with_land_classification/leps_dallas_data_with_land_classification.csv"))
    #,
    #cbind(city="san_fransisco", read.csv(
      #"./data/biodiversity_data_with_land_classification/leps_sf_data_with_land_classification.csv"))
  ) %>%
  filter(species != "")
  
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("longitude", "latitude"), 
                     crs = 4326))
  
  #-----------------------------------------------------
  # map the data onto a spatial file
  
  dallas_shp <- geojson_sf("./data/urbanwatch_data/city_outlines/dallas_boundary/dallas.geojson")
  los_angeles_shp <- read_sf("./data/urbanwatch_data/city_outlines/los_angeles_boundary/geo_export_f62c93b8-4e53-4645-8ef0-b08c094f4cd8.shp")
  san_diego_shp <- read_sf("./data/urbanwatch_data/city_outlines/san_diego_boundary_datasd/san_diego_boundary_datasd.shp")
  
  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070
  grid_size <- 2000 # grid size in metres (1000 = 1km)
  
  # dallas
  dallas_shp <- dallas_shp  %>% 
    st_transform(., crs) 
  dallas_shp <- sf::st_buffer(dallas_shp, dist = 0)
  dallas_grid <- st_make_grid(dallas_shp, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  dallas_grid <- st_intersection(dallas_grid, dallas_shp) %>%
    mutate(grid_id = paste0("dallas_", grid_id)) %>%
    select(grid_id)
  # los angeles
  los_angeles_shp <- los_angeles_shp  %>% 
    st_transform(., crs) 
  los_angeles_grid <- st_make_grid(los_angeles_shp, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  los_angeles_grid <- st_intersection(los_angeles_grid, los_angeles_shp) %>%
    mutate(grid_id = paste0("los_angeles_", grid_id)) %>%
    select(grid_id)
  # san diego
  san_diego_shp <- san_diego_shp  %>% 
    st_transform(., crs) 
  san_diego_grid <- st_make_grid(san_diego_shp, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  san_diego_grid <- st_intersection(san_diego_grid, san_diego_shp) %>%
    mutate(grid_id = paste0("san_diego_", grid_id)) %>%
    select(grid_id)
  
  
  
  # combine the grids from all cities to get all sites
  all_grids <- bind_rows(dallas_grid, los_angeles_grid, san_diego_grid)
  
  # and now match detections with spatial units (grid cells)
  df <- st_transform(df_sf, crs = crs) %>%
    st_join(all_grids, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the urban grid
    filter(!is.na(grid_id)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, gbifID, latitude, longitude), by="gbifID") %>%
    
    # to confirm: did Jenny already remove samples without detailed coordinate uncertainty?
    
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(species) %>%
    add_tally() %>%
    # filter out species below some minimum number of detections
    # right now, this is detections across ALL cities
    filter(n >= min_species_detections) %>% 
    rename("species_detections" = "n") %>%
    ungroup()
  
  # make the detection data a spatial file for plotting
  (df_sf <- st_as_sf(df,
                     coords = c("longitude", "latitude"), 
                     crs = 4326))
  
  # view the shapefile
  ggplot() + 
    geom_sf(data = san_diego_shp, fill = 'white', lwd = 0.05) +
    geom_sf(data = san_diego_grid, fill = 'transparent', lwd = 0.3) +
    geom_sf(data = filter(df_sf, city == "san_diego"), aes(colour=as.factor(year))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle("San Diego")
  #theme(legend.position="none") 
  
  # view the shapefile
  ggplot() + 
    geom_sf(data = los_angeles_shp, fill = 'white', lwd = 0.05) +
    geom_sf(data = los_angeles_grid, fill = 'transparent', lwd = 0.3) +
    geom_sf(data = filter(df_sf, city == "los_angeles"), aes(colour=as.factor(year))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle("Los Angeles")
  #theme(legend.position="none") 
  
  # view the shapefile
  ggplot() + 
    geom_sf(data = dallas_shp, fill = 'white', lwd = 0.05) +
    geom_sf(data = dallas_grid, fill = 'transparent', lwd = 0.3) +
    geom_sf(data = filter(df_sf, city == "dallas"), aes(colour=as.factor(year))) +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    ggtitle("Dallas")
  #theme(legend.position="none")
  
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
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, family, city, grid_id, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(species, grid_id, year, survey) %>% 
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
                     group_by(grid_id) %>%
                     slice(1) %>%
                     ungroup() %>%
                     select(grid_id)))
  
  site_vector <- site_names %>%
    pull(grid_id)
  
  # n cities and city vector
  # how many species were detected?
  n_cities <- nrow(city_names <- df %>%
                      # group by species ID
                      group_by(city) %>%
                      slice(1) %>%
                      select(city))
  
  city_vector <- city_names %>%
    pull(city)
  
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
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        full_join(., select(site_names, grid_id), by="grid_id") %>%
        #full_join(., as.data.frame(site_vector), by=c("PARK_NAME" = "site_vector")) %>% 
        
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(grid_id, row, fill = 0) %>%
        
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
  # identify whether a city is "in range" of each species
  # ie. whether the species was found in the city one or more times
  
  all_species_city_combos <- as.data.frame(cbind(
    rep(species_vector, each=n_cities),
    rep(city_vector, times=n_species)
  )) %>% 
    rename("species" = "V1",
           "city" = "V2")
  
  species_city_occurrences <- rbind(
    cbind(city="dallas", read.csv(
      "./data/biodiversity_data_with_land_classification/leps_dallas_data_with_land_classification.csv")),
    cbind(city="los_angeles", read.csv(
      "./data/biodiversity_data_with_land_classification/leps_la_data_with_land_classification.csv")),
    cbind(city="san_diego", read.csv(
      "./data/biodiversity_data_with_land_classification/leps_san_diego_data_with_land_classification.csv"))
    #,
    #cbind(city="san_fransisco", read.csv(
    #"./data/biodiversity_data_with_land_classification/leps_sf_data_with_land_classification.csv"))
  ) %>%
    filter(species != "") %>%
    filter(species %in% species_vector) %>%
    group_by(species, city) %>%
    slice(1) %>%
    ungroup() %>%
    select(species, city) %>%
    mutate(occurs = 1) %>%
    left_join(all_species_city_combos, .) %>%
    mutate(occurs = replace_na(occurs, 0))

  # now create a df of all possible species city occurrences
  city <- sub("^(.*)[_].*", "\\1",site_vector)
  all_species_sites_possible <- as.data.frame(cbind(site_vector, city))
  
  # and add the city name
  ranges <- left_join(all_species_sites_possible, species_city_occurrences, 
                    by=c("city")) %>%
    arrange(., species)
  
  ranges <- ranges %>% 
    select(-city) %>% 
    pivot_wider(names_from = site_vector, values_from = occurs) %>% 
    column_to_rownames('species') %>%
    as.matrix()
  
  ranges <- unname(ranges)
  
  #-----------------------------------------------------
  # identify community sampling events which we will use to infer non-detections
  # identified by site*survey months with > 1 species detected. Another option
  # might be to group by iNat user name by date*park but this may be a very narrow scope?
  
  community_samples <- df %>%
    
    # determine whether a community sampling event occurred
    # using collector/observer name
    group_by(family, grid_id, year, survey) %>%
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
    group_by(grid_id, family, year, survey) %>%
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
           "grid_id" = "V2",
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
      by=c("grid_id","family", "year", "survey")) %>%
    
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
    n_cities = n_cities, # number of cities
    
    species = species_vector,
    families = family_vector, 
    sites = site_vector,
    years = year_vector,
    surveys = survey_vector,
    cities = city_vector,
    
    species_city_occurrences = species_city_occurrences,
    ranges = ranges
    
  ))

}
