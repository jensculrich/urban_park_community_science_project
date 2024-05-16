library(tidyverse)
library(sf) # spatial data processing

prep_data <- function(grid_size,
                      min_species_detections) {
  #-----------------------------------------------------
  # summarize the detection data
  
  # first read the data 
  df <- rbind(
    cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2020-los-angeles-county.csv"), year = 1),
    cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2021-los-angeles-county.csv"), year = 2),
    cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2022-los-angeles-county.csv"), year = 3),
    cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2023-los-angeles-county.csv"), year = 4),
    cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2024-los-angeles-county.csv"), year = 5)
  )
          
  # and perform some initial filters
  df <- df %>%
    
    # for now, let's only look at detections with a species ID
    filter(str_count(scientific_name, "\\w+") == 2) %>%
    
    # for now, let's filter out geo obscured detections
    filter(geoprivacy != "obscured") %>%
    
    # and additionally, let's get rid of any points with huge coordinate uncertainty
    # say > 1000m for now
    filter(positional_accuracy <= 1000) %>%
  
    # for now, to speed up the inference, let's also filter to the more common species
    group_by(scientific_name) %>%
    add_tally() %>%
    filter(n >= min_species_detections) %>%
    dplyr::select(-n) %>%
    ungroup()
  
  #-----------------------------------------------------
  # map the data on a spatial file
  
  county_shp <- read_sf("./data/los_angeles_county_boundary_shapefile/County_Boundary.SHP")
  
  # USA_Contiguous_Albers_Equal_Area_Conic
  crs <- 5070

  county_shp <- county_shp  %>% 
    st_transform(., crs) %>% # USA_Contiguous_Albers_Equal_Area_Conic
    filter(!OBJECTID %in% c(1,2,3,5,6,7))
  
  # create "grid_size" km grid over the area
  grid <- st_make_grid(county_shp, cellsize = c(grid_size, grid_size)) %>% 
    st_sf(grid_id = 1:length(.))
  
  clipped_grid <- st_intersection(grid, county_shp)
  
  rm(grid)
  
  # make the detection data a spatial file
  (df_sf <- st_as_sf(df,
                     coords = c("longitude", "latitude"), 
                     crs = 4326))
  
  # and then transform it to the crs
  df <- st_transform(df_sf, crs = crs) %>%
    st_join(clipped_grid, join = st_intersects) %>% as.data.frame %>%
    # filter out records from outside of the grid
    filter(!is.na(grid_id)) %>%
    # now rejoin the lat/long data for each point
    left_join(., dplyr::select(
      df, id, latitude, longitude), by="id")
  
  rm(df_sf)
  
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
                      group_by(scientific_name) %>%
                      add_tally() %>%
                      # and take one record
                      slice(1) %>%
                      select(scientific_name, common_name, n))
  
  species_vector <- species_names %>%
    pull(scientific_name)
  
  n_detections <- nrow(df)
  
  n_sites <- (nrow(site_names <- df %>%
         group_by(grid_id) %>%
         slice(1) %>%
         ungroup() %>%
         select(grid_id)))
  
  site_vector <- site_names %>%
    pull(grid_id)
  
  
  #-----------------------------------------------------
  # prep data for array format
  
  df <- df %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.numeric(factor(observed_on))) %>%
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(scientific_name, grid_id, survey, year) %>%
    
    # turn into binary detections (for occupancy rather than abundance model)
    group_by(scientific_name, grid_id, year, survey) %>% 
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
        full_join(., select(species_names, scientific_name), by="scientific_name") %>%
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        full_join(., select(site_names, grid_id), by="grid_id") %>%
        
        # group by SPECIES
        group_by(scientific_name) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(grid_id, row, fill = 0) %>%
        
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(4:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        dplyr::select(4:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(scientific_name %in% levels(as.factor(species_names$scientific_name)))
      
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


  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    V = V, # community science detection data
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
