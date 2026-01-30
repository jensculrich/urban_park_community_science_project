# create a detection array from iNaturalist data.
# this will join detections to sites.
# non-detections are currently added whenever we don't see a species.
# sampling effort could then be modeled as a function of species/site/year/time,
# using the months within year as repeat "survey periods".
# but possibly we could look to identify "community sampling events"
# and only model non-detections as those species not observed on a given community sampling event?

library(tidyverse) # data organization

prep_data <- function(city_names,
                      min_species_detections,
                      min_site_years_w_detection,
                      min_species_for_community_sampling_event,
                      family_sampling,
                      remove_outlier_parks,
                      write_city_data_csv
                      ) {
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  #-----------------------------------------------------
  # load the detection data
  
  df <- as.data.frame(matrix(nrow = 0, ncol = 70))
  
  for(i in 1:length(city_names)){
    
    city <- city_names[i]
    
    # first read the data 
    temp <- cbind(city, read.csv(paste0(
      "./data/detections_by_city/", city, "/01_0m_", city,
      "_observations_parkID_2km_clipped.csv"))
    )
    
    df <- rbind(df, temp)
    
  }
  
  rm(temp)
  
  # define butterfly families to include
  butterfly_families <- c("Hesperiidae", "Lycaenidae", "Nymphalidae", 
                          "Papilionidae", "Pieridae")
  
  # The GBIF output doesn't recognize Rioidinidae or the common species Apodemia virgulti
  df <- df %>%
    mutate(family = ifelse(taxonRank=="Riodinidae", "Riodinidae", family)) %>%
    mutate(species = ifelse(verbatimScientificName=="Apodemia virgulti", 
                            "Apodemia virgulti", species)) %>%
  
  # use Riodinidae as Lycaenidae for comm survey purposes (to be discussed with group)
  mutate(family = ifelse(family == "Riodinidae", 
                         "Lycaenidae", 
                         family)) %>%
    
  # remove records with no species name
  filter(species != "") %>%
  
  # filter down to the butterfly families
  filter(family %in% butterfly_families) %>%
  
  filter(coordinateUncertaintyInMeters < 100) 
  
  # genus <- unique(df$genus)
  
  # write a file with list of genera that we need to gather data for (only need to do this once)
  #write.csv(genus, "data/lepidoptera_trait_data/ease_of_id/identifiability_by_genus.csv")
  
  #-----------------------------------------------------
  # now we are going to filter out all of the points from unclassified parks 
  # (parks outside city boundaries)
  # you'd need to do this later if you want to model responses of species occurring in
  # classified parks but never in unclassified parks
  
  df <- df %>%
    filter(type == "classified")
  
  #-----------------------------------------------------
  # sort out sites with low temporal coverage of detections
  # for example, we could only model sites that have detections in at least two
  # years, so that we are more confident in estimating a temporal trend across those sites
  
  df <- df %>%
    group_by(city, new_id) %>%
    mutate(new_id_unique = cur_group_id()) %>%
    group_by(new_id_unique) %>%
    mutate(years_w_detection_by_site = length(unique(year))) %>%
    
    # now filter out data from sites with < min_site_years_w_detection
    # so that we do not model responses in these sites
    filter(years_w_detection_by_site >= min_site_years_w_detection)
    
  #-----------------------------------------------------
  # get raw number of detections by city
  total_detections_by_city <- df %>%
    group_by(city) %>%
    add_tally() %>%
    slice(1) %>%
    ungroup() %>%
    select(city, n) %>%
    rename("total_detections" = "n") %>%
    mutate(total_detections = center_scale(total_detections))
  
  #-----------------------------------------------------
  # get median number of detections per recorder by city
  recorder_detections_by_city <- df %>%
    group_by(city, recordedBy) %>%
    add_tally() %>%
    select(city, recordedBy, n) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city) %>%
    mutate(mean_detections_by_contributor = mean(n), 
           max_detections_by_contributor = max(n)) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, mean_detections_by_contributor, max_detections_by_contributor) %>%
    mutate(mean_detections_by_contributor_scaled = center_scale(mean_detections_by_contributor),
           max_detections_by_contributor_scaled = center_scale(max_detections_by_contributor))
  
  #-----------------------------------------------------
  # prep data for VECTOR format
  
  # currently this only considers sites where 1 or more species ever detected
  
  # we will have one big array, with different lengths of species, sites and years
  # for now we are going to define surveys as two-month periods
  
  df <- df %>%
    
    # change date to ordinal day
    #mutate(survey = as.numeric(factor(observed_on))) %>%
    
    # add tally of detections per species
    group_by(species) %>%
    add_tally(name = 'n_detections') %>%
    ungroup() %>%
    
    # to get two-month surveys, we can subtract one from all even numbered months
    #mutate(survey = ifelse(month %% 2 == 0, month - 1, month)) %>% 
    # or use one moth surveys
    mutate(survey = month) %>%
    
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(survey),
    #mutate(survey = as.integer(factor(survey)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(city, species, family, new_id, 
                  survey, year, n_detections) %>%
    
    # turn into binary detections (for occupancy rather than abundance)
    group_by(city, species, new_id, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # add tally of binary detections per species
    group_by(species) %>%
    add_tally(name = 'n_binary_detections') %>%
    ungroup() %>%
    
    # and filter out species with too few detections
    filter(n_binary_detections >= min_species_detections) %>%
    
    # arrange by survey within year  
    arrange(city, year, survey) 
  
  # get dimensions of surveys and years
  survey_vector <- seq(1:12) # monthly organization of survey events
  n_surveys <- length(survey_vector)
  
  year_vector <- as.vector(levels(as.factor(df$year)))
  n_years <- length(year_vector)
  
  # start a "species_info" table
  species_info <- df  %>%
    mutate(genus = word(species, 1)) %>%
    select(species, genus, family, n_detections, n_binary_detections) %>%
    group_by(species) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(species_integer = as.integer(as.factor(species)))
  
  ## --------------------------------------------------
  ## get site covariate data
  
  # source the prep function
  source("./part2_local_landscape_predictors_of_occupancy/run_model/get_site_data.R")
  
  site_data_temp <- get_site_data(city_names)
  
  site_data <- df %>%
    group_by(city, new_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, new_id)
  
  site_data <- left_join(site_data, site_data_temp, by = c("city", "new_id"))
  
  # the original scaled values are scaled within city ,considering all parks in those cities
  # scaled_2 are scaled within a city, ONLY CONSIDERING sites that will be modeled
  # scaled_across_all_cities are scaled across all cities, ONLY CONSIDERING sites that will be modeled
  site_data <- site_data %>%
    group_by(city) %>%
    mutate(isolation_scaled_2 = center_scale(isolation),
           log_total_green_space_area_scaled_2 = center_scale(log_total_green_space_area),
           log_isolation_scaled_2 = center_scale(log(isolation))) %>%
    ungroup() %>%
    # these are the ones we actually use for the model
    mutate(log_total_green_space_area_scaled_across_all_cities = center_scale(log_total_green_space_area),
           plant_genera_density_scaled_across_all_cities = center_scale(plant_genera_density),
           tree_cover_scaled_across_all_cities = center_scale(tree_percent_cover),
           isolation_scaled_across_all_cities = center_scale(isolation),
           log_isolation_scaled_across_all_cities = center_scale(log(isolation)),
           landcover_type_diversity_scaled_across_all_cities = center_scale(shannon_diversity),
           proportion_landscape_vegetation_scaled_across_all_cities = center_scale(log(proportion_landscape_vegetation+0.01)),
           proportion_landscape_open_developed_scaled_across_all_cities = center_scale(log(proportion_landscape_open_developed+0.01)),
           proportion_landscape_medhigh_developed_scaled_across_all_cities = center_scale(log(proportion_landscape_medhigh_developed+0.01)),
           proportion_landscape_woody_scaled_across_all_cities = center_scale(log(proportion_landscape_woody+0.01)),
           proportion_landscape_grassherb_scaled_across_all_cities = center_scale(log(proportion_landscape_grassherb+0.01))
    ) 
    
  
  if(remove_outlier_parks == TRUE){
    site_data <- site_data %>%
      filter(total_green_space_area > 0.11) 
  }
  
  n_sites <- nrow(site_data <- site_data %>%
    mutate(multicity_site_id = row_number()))
    
  # get a vector of site names
  site_vector <- site_data %>%
    pull(multicity_site_id)
  
  # plot the spread of the site covariate data
  par(mfrow=c(1,4))  
  hist(site_data$proportion_landscape_woody_scaled_across_all_cities)
  hist(site_data$proportion_landscape_grassherb_scaled_across_all_cities)  
  hist(site_data$plant_genera_density_scaled_across_all_cities)
  hist(site_data$tree_cover_scaled_across_all_cities)
  #hist(site_data$perarea_idx_scaled)
  #hist(site_data$proximity_scaled)
  hist(site_data$log_total_green_space_area_scaled_across_all_cities)
  hist(site_data$log_isolation_scaled_across_all_cities)
  hist(site_data$landcover_type_diversity_scaled_across_all_cities)
  
  par(mfrow=c(1,1)) 
  
  # how correlated are the site covariate data
  cor(site_data$log_total_green_space_area_scaled, site_data$log_isolation_scaled)
  cor(site_data$proportion_landscape_woody_scaled_across_all_cities, site_data$landcover_type_diversity_scaled_across_all_cities)
  cor(site_data$plant_genera_density_scaled, site_data$log_total_green_space_area_scaled_2)
  cor(site_data$log_isolation_scaled_across_all_cities, site_data$log_total_green_space_area_scaled_across_all_cities)
  
  #plot(site_data$log_total_green_space_area_scaled_2, site_data$log_isolation_scaled_2)

  
  # plot
  ggplot(site_data, aes(
    x = log_total_green_space_area_scaled_across_all_cities, y = log_isolation_scaled_across_all_cities, colour = city)) +
    geom_point()
  
  # plot
  ggplot(site_data, aes(
    x = log_total_green_space_area_scaled_across_all_cities, y = plant_genera_density_scaled_across_all_cities, colour = city)) +
    geom_point()
  
  # plot
  ggplot(site_data, aes(
    x = proportion_landscape_vegetation_scaled_across_all_cities, y = landcover_type_diversity_scaled_across_all_cities, colour = city)) +
    geom_point()
  
  # add detections by city
  site_data <- site_data %>%
    left_join(., total_detections_by_city) %>%
    left_join(., recorder_detections_by_city)
  
  ## --------------------------------------------------
  # Get species ranges
  
  # source the prep function
  source("./part2_local_landscape_predictors_of_occupancy/run_model/get_species_ranges.R")
  
  ranges_raw <- get_species_ranges(city_names)
  
  ## --------------------------------------------------
  # identify community sampling events
  if(family_sampling == TRUE){
    
    # first identify community sampling events 
    # (events with > min_species_for_community_sampling_event species detected)
    # we will have to re add any individual species detections below the threshold later
    # (times when we know single species detected, but not inferring search for rest of community)
    community_sampling_events <- df %>%
      group_by(city, new_id, year, survey, family) %>%
      mutate(num_species_detected = length(unique(species))) %>%
      filter(num_species_detected >= min_species_for_community_sampling_event) %>%
      slice(1) %>%
      ungroup() %>%
      select(city, new_id, year, survey, family) %>%
      
      # join the multicity site id and filter out any sites that were filtered out
      # due to the site data, e.g. park was too small to include.
      left_join(., select(site_data, city, new_id, multicity_site_id)) %>%
      filter(!is.na(multicity_site_id)) %>%
      
      # now add a unique community sample event id
      mutate(community_sample_id = row_number())
    
    # get a vector of community_sample event id's
    community_sample_id <- community_sampling_events$community_sample_id
    
    ## --------------------------------------------------
    # get detection data for all community sampling events
    # detection data will be cast into vector of lenght R (total number of yearXspeciesXsite with possible detections)
    
    # create a df containing all possible species for each community sampling event
    temp1 <- select(species_info, family, species)
    # expand the community sampling events to include all of the species that weren't detected (but could have been) 
    temp2 <- full_join(community_sampling_events, temp1, relationship = "many-to-many")
    
    # now join range data to filter out all species that are not in range of the city
    temp2 <- left_join(temp2, ranges_raw) %>%
      filter(!is.na(in_range)) %>%
      select(-in_range)
    
    # get names of all surveys (in case a cluster doesn't have them all surveyed)
    # this will avoid having to rewrite the stan model code for different # max survey lengths
    temp <- as.data.frame(survey_vector) %>%
      rename("survey" = "survey_vector")
    
    # now prepare the detections into a df and prepare to add 1's or 0's based on whether species were detected
    detections_df <- df %>%
      # join the data about the community sampling event from that year
      left_join(., temp2) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      #filter(!is.na(community_sample_id)) %>%
      # assign detection to those rows
      mutate(detection = 1) %>%
      # join all survey events within a year
      full_join(., temp) %>%
      # arrange surveys from 1 - 12
      arrange(match(survey, survey_vector)) %>%
      # and cast a detection matrix for each year with surveys 1 - 12
      pivot_wider(., names_from = survey, values_from = detection) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      filter(!is.na(community_sample_id)) %>%
      # join all undetected species for each comm survey event (they could have been detected but weren't)
      full_join(., temp2) %>%
      # now we only want community_sample_id, species, and the matrix of presence absences
      select(community_sample_id, family, species, "1", "2", "3", "4", "5", "6", 
             "7", "8", "9", "10", "11", "12") %>%
      
      # now we need to join all the data about site number / year for all those undetected species
      left_join(., temp2, by = c("community_sample_id", "family", "species")) %>%
      arrange(community_sample_id, species) %>%
      select(-survey)
    
    ## --------------------------------------------------
    # now we need to add in sampling events for individual species that were missed by the
    # broader community sampling events structure
    # if min species for community sampling event = 0 (the default) none of this will actually apply
    
    # reverse the filter and look at
    # (events with *<* min_species_for_community_sampling_event species detected)
    targeted_sampling_events <- df %>%
      group_by(city, new_id, year, survey, family) %>%
      mutate(num_species_detected = length(unique(species))) %>%
      filter(num_species_detected < min_species_for_community_sampling_event) %>%
      ungroup() %>%
      select(species, city, new_id, year, survey, family) %>%
      
      # join the multicity site id and filter out any sites that were filtered out
      # due to the site data, e.g. park was too small to include.
      left_join(., select(site_data, city, new_id, multicity_site_id)) %>%
      filter(!is.na(multicity_site_id)) 
    
    targeted_sampling_detections_df <- df %>%
      # join the data about the targeted sampling event from that year
      left_join(., targeted_sampling_events) %>%
      # filter out is na multicity_site_id - removes rows already added that included comm sampling events    
      filter(!is.na(multicity_site_id)) %>%
      # assign detection to those rows
      mutate(detection = 1) %>%
      # join all survey events within a year
      full_join(., temp) %>%
      # arrange surveys from 1 - 12
      arrange(match(survey, survey_vector)) %>%
      # and cast a detection matrix for each year with surveys 1 - 12
      pivot_wider(., names_from = survey, values_from = detection) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      filter(!is.na(multicity_site_id)) %>%
      # add a "community sample id", I know it's not actually a community sample but need to combine with comm samples df
      mutate(community_sample_id = (row_number() + max(community_sample_id))) %>%
      # now we only want community_sample_id, species, and the matrix of presence absences
      select(community_sample_id, family, species, "1", "2", "3", "4", "5", "6", 
             "7", "8", "9", "10", "11", "12", city, new_id, year, multicity_site_id)
    
    ## --------------------------------------------------
    # now combine the community surveys and targeted species surveys
    
    detections_df <- as.data.frame(rbind(detections_df, targeted_sampling_detections_df)) 
    # and replace NAs with 0's (species not detected)
    detections_df[,4:15][is.na(detections_df[,4:15])] <- 0
    
    # if siteXyear row already existed in detections df, because some community sampling events
    # happened at other survey months, then we want to fold in the targeted detections with the 
    # community detections
    detections_df <- detections_df %>%
      # couldn't figure out how to get R to recognize column names that are numbers
      rename("a" = "1", "b" = "2", 
             "c" = "3", "d" = "4", 
             "e" = "5", "f" = "6", 
             "g" = "7", "h" = "8", 
             "i" = "9", "j" = "10", 
             "k" = "11", "l" = "12") %>%
      # within a site, year, species, for each column, 
      # add together any taregted detections so that repeat rows
      # all recognize that the species was detected, whether through a
      # community sample or targeted sample
      group_by(multicity_site_id, year, species) %>%
      mutate(a = sum(a), 
             b = sum(b),
             c = sum(c),
             d = sum(d),
             e = sum(e),
             f = sum(f),
             g = sum(g),
             h = sum(h),
             i = sum(i),
             j = sum(j),
             k = sum(k),
             l = sum(l)
      ) %>%
      # add a group id that keeps track of the species, site, year combo
      # there may be unique detections in different months but we don't want
      mutate(group_id = cur_group_id()) %>%
      group_by(multicity_site_id, year, species, group_id) %>%
      # slice out repeats for a row with a targeted sampling event that was 
      # in a year with community sampling events at different months
      slice(1) %>%
      ungroup() %>%
      
      # Now find minimum year of survey for a park and tag this year as first year
      group_by(multicity_site_id, species) %>%
      mutate(first_year = ifelse(year == min(year), 1, 0)) %>%
      # and a unique survey year per site
      mutate(site_survey_year = as.integer(as.factor(year))) %>%
      ungroup() %>%
      
      # need to create an index to access row number of previous year for same species at same site
      # only need it for years after the first year 
      # (although will have to assign something other than NA for first year data)
      mutate(row_index = row_number()) %>%
      group_by(multicity_site_id, species) %>% 
      mutate("prev_index" = lag(row_index)) %>%
      ungroup() %>%
      mutate(prev_index = replace_na(prev_index, 0)) 
    
    ## --------------------------------------------------
    # now we need to get NA data to pass over zero effort months in 
    # site x years with some sampling effort 
    
    # there;s probably a more efficient way to do this, look into later!
    V_NA <- detections_df %>%
      
      # group by community sample id (year, site, family search event)
      group_by(year, multicity_site_id, family) %>%
      
      mutate(a2 = sum(a)) %>%
      mutate(check = ifelse(a2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(a3 = ifelse(check, max(a), a)) %>%
      mutate(a = a3) %>%
      
      mutate(b2 = sum(b)) %>%
      mutate(check = ifelse(b2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(b3 = ifelse(check, max(b), b)) %>%
      mutate(b = b3) %>%
      
      mutate(c2 = sum(c)) %>%
      mutate(check = ifelse(c2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(c3 = ifelse(check, max(c), c)) %>%
      mutate(c = c3) %>%
      
      mutate(d2 = sum(d)) %>%
      mutate(check = ifelse(d2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(d3 = ifelse(check, max(d), d)) %>%
      mutate(d = d3) %>%
      
      mutate(e2 = sum(e)) %>%
      mutate(check = ifelse(e2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(e3 = ifelse(check, max(e), e)) %>%
      mutate(e = e3) %>%
      
      mutate(f2 = sum(f)) %>%
      mutate(check = ifelse(f2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(f3 = ifelse(check, max(f), f)) %>%
      mutate(f = f3) %>%
      
      mutate(g2 = sum(g)) %>%
      mutate(check = ifelse(g2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(g3 = ifelse(check, max(g), g)) %>%
      mutate(g = g3) %>%
      
      mutate(h2 = sum(h)) %>%
      mutate(check = ifelse(h2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(h3 = ifelse(check, max(h), h)) %>%
      mutate(h = h3) %>%
      
      mutate(i2 = sum(i)) %>%
      mutate(check = ifelse(i2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(i3 = ifelse(check, max(i), i)) %>%
      mutate(i = i3) %>%
      
      mutate(j2 = sum(j)) %>%
      mutate(check = ifelse(j2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(j3 = ifelse(check, max(j), j)) %>%
      mutate(j = j3) %>%
    
      mutate(k2 = sum(k)) %>%
      mutate(check = ifelse(k2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(k3 = ifelse(check, max(k), k)) %>%
      mutate(k = k3) %>%
      
      mutate(l2 = sum(l)) %>%
      mutate(check = ifelse(l2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(l3 = ifelse(check, max(l), l)) %>%
      mutate(l = l3) 
               
    # now also pull out our detections
    V_NA <- V_NA %>%
      rename("1" = "a", "2" = "b", 
             "3" = "c", "4" = "d", 
             "5" = "e", "6" = "f", 
             "7" = "g", "8" = "h", 
             "9" = "i", "10" = "j", 
             "11" = "k", "12" = "l")
    
    # and then pull out the NAs stuff as a simple matrix
    V_NA <- as.matrix(V_NA[4:15])
    
    # now also pull out our detections
    detections_df <- detections_df %>%
      rename("1" = "a", "2" = "b", 
             "3" = "c", "4" = "d", 
             "5" = "e", "6" = "f", 
             "7" = "g", "8" = "h", 
             "9" = "i", "10" = "j", 
             "11" = "k", "12" = "l")
    
    V <- as.matrix(detections_df[4:15])
    
    # print mismatches in NA array (Detections when the NA says we haven't surveyed)
    print(paste0("Number of NA / detection mismatches: ", length(which(V > V_NA))))
    
  } else{ # else order level sampling inferred
    
    # first identify community sampling events 
    # (events with > min_species_for_community_sampling_event species detected)
    # we will have to re add any individual species detections below the threshold later
    # (times when we know single species detected, but not inferring search for rest of community)
    community_sampling_events <- df %>%
      #group_by(city, new_id, year, survey, family) %>%
      group_by(city, new_id, year, survey) %>%
      mutate(num_species_detected = length(unique(species))) %>%
      filter(num_species_detected >= min_species_for_community_sampling_event) %>%
      slice(1) %>%
      ungroup() %>%
      select(city, new_id, year, survey) %>%
    
      # join the multicity site id and filter out any sites that were filtered out
      # due to the site data, e.g. park was too small to include.
      left_join(., select(site_data, city, new_id, multicity_site_id)) %>%
      filter(!is.na(multicity_site_id)) %>%
      
      # now add a unique community sample event id
      mutate(community_sample_id = row_number())
    
    # get a vector of community_sample event id's
    community_sample_id <- community_sampling_events$community_sample_id
  
    ## --------------------------------------------------
    # get detection data for all community sampling events
    # detection data will be cast into vector of lenght R (total number of yearXspeciesXsite with possible detections)
    
    # create a df containing all possible species for each community sampling event
    temp1 <- as.data.frame(cbind(
      rep(species_info$species, times = length(community_sampling_events$community_sample_id)),
      rep(community_sampling_events$community_sample_id, each = length(species_info$species))
    )) %>%
      rename("species" = "V1",
             "community_sample_id" = "V2") %>%
      mutate(community_sample_id = as.integer(community_sample_id))
    # expand the community sampling events to include all of the species that weren't detected (but could have been) 
    temp2 <- full_join(community_sampling_events, temp1)
    
    # get names of all surveys (in case a cluster doesn't have them all surveyed)
    # this will avoid having to rewrite the stan model code for different # max survey lengths
    temp <- as.data.frame(survey_vector) %>%
      rename("survey" = "survey_vector")
    
    # now prepare the detections into a df and prepare to add 1's or 0's based on whether species were detected
    detections_df <- df %>%
      # join the data about the community sampling event from that year
      left_join(., temp2) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      #filter(!is.na(community_sample_id)) %>%
      # assign detection to those rows
      mutate(detection = 1) %>%
      # join all survey events within a year
      full_join(., temp) %>%
      # arrange surveys from 1 - 12
      arrange(match(survey, survey_vector)) %>%
      # and cast a detection matrix for each year with surveys 1 - 12
      pivot_wider(., names_from = survey, values_from = detection) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      filter(!is.na(community_sample_id)) %>%
      # join all undetected species for each comm survey event (they could have been detected but weren't)
      full_join(., temp2) %>%
      # now we only want community_sample_id, species, and the matrix of presence absences
      select(community_sample_id, species, "1", "2", "3", "4", "5", "6", 
             "7", "8", "9", "10", "11", "12")  %>%
    
      # now we need to join all the data about site number / year for all those undetected species
      left_join(., temp2, by = c("community_sample_id", "species")) %>%
      arrange(community_sample_id, species) %>%
      select(-survey)
    
    ## --------------------------------------------------
    # now we need to add in sampling events for individual species that were missed by the
    # broader community sampling events structure
    
    # reverse the filter and look at
    # (events with *<* min_species_for_community_sampling_event species detected)
    targeted_sampling_events <- df %>%
      group_by(city, new_id, year, survey) %>%
      mutate(num_species_detected = length(unique(species))) %>%
      filter(num_species_detected < min_species_for_community_sampling_event) %>%
      ungroup() %>%
      select(species, city, new_id, year, survey) %>%
    
      # join the multicity site id and filter out any sites that were filtered out
      # due to the site data, e.g. park was too small to include.
      left_join(., select(site_data, city, new_id, multicity_site_id)) %>%
      filter(!is.na(multicity_site_id)) 
    
    targeted_sampling_detections_df <- df %>%
      # join the data about the targeted sampling event from that year
      left_join(., targeted_sampling_events) %>%
      # filter out is na multicity_site_id - removes rows already added that included comm sampling events    
      filter(!is.na(multicity_site_id)) %>%
      # assign detection to those rows
      mutate(detection = 1) %>%
      # join all survey events within a year
      full_join(., temp) %>%
      # arrange surveys from 1 - 12
      arrange(match(survey, survey_vector)) %>%
      # and cast a detection matrix for each year with surveys 1 - 12
      pivot_wider(., names_from = survey, values_from = detection) %>%
      # filter out is na - removes rows added in to cover all survey columns     
      filter(!is.na(multicity_site_id)) %>%
      # add a "community sample id", I know it's not actually a community sample but need to combine with comm samples df
      mutate(community_sample_id = (row_number() + max(community_sample_id))) %>%
      # now we only want community_sample_id, species, and the matrix of presence absences
      select(community_sample_id, species, "1", "2", "3", "4", "5", "6", 
             "7", "8", "9", "10", "11", "12", city, new_id, year, multicity_site_id)
    
    ## --------------------------------------------------
    # now combine the community surveys and targeted species surveys
    
    detections_df <- as.data.frame(rbind(detections_df, targeted_sampling_detections_df)) 
    # and replace NAs with 0's (species not detected)
    detections_df[,3:14][is.na(detections_df[,3:14])] <- 0
    
    # if siteXyear row already existed in detections df, because some community sampling events
    # happened at other survey months, then we want to fold in the targeted detections with the 
    # community detections
    detections_df <- detections_df %>%
      # couldn't figure out how to get R to recognize column names that are numbers
      rename("a" = "1", "b" = "2", 
             "c" = "3", "d" = "4", 
             "e" = "5", "f" = "6", 
             "g" = "7", "h" = "8", 
             "i" = "9", "j" = "10", 
             "k" = "11", "l" = "12") %>%
      # within a site, year, species, for each column, 
      # add together any taregted detections so that repeat rows
      # all recognize that the species was detected, whether through a
      # community sample or targeted sample
      group_by(multicity_site_id, year, species) %>%
      mutate(a = sum(a), 
             b = sum(b),
             c = sum(c),
             d = sum(d),
             e = sum(e),
             f = sum(f),
             g = sum(g),
             h = sum(h),
             i = sum(i),
             j = sum(j),
             k = sum(k),
             l = sum(l)
      ) %>%
      # add a group id that keeps track of the species, site, year combo
      # there may be unique detections in different months but we don't want
      mutate(group_id = cur_group_id()) %>%
      group_by(multicity_site_id, year, species, group_id) %>%
      # slice out repeats for a row with a targeted sampling event that was 
      # in a year with community sampling events at different months
      slice(1) %>%
      ungroup() %>%
    
      # Now find minimum year of survey for a park and tag this year as first year
      group_by(multicity_site_id, species) %>%
      mutate(first_year = ifelse(year == min(year), 1, 0)) %>%
      # and a unique survey year per site
      mutate(site_survey_year = as.integer(as.factor(year))) %>%
      ungroup() %>%
    
      # need to create an index to access row number of previous year for same species at same site
      # only need it for years after the first year 
      # (although will have to assign something other than NA for first year data)
      mutate(row_index = row_number()) %>%
      group_by(multicity_site_id, species) %>% 
      mutate("prev_index" = lag(row_index)) %>%
      ungroup() %>%
      mutate(prev_index = replace_na(prev_index, 0)) 
    
    ## --------------------------------------------------
    # now we need to get NA data to pass over zero effort months in 
    # site x years with some sampling effort 
    
    # there;s probably a more efficient way to do this, look into later!
    V_NA <- detections_df %>%
      
      # group by site X year
      group_by(year, multicity_site_id) %>%
      
      mutate(a2 = sum(a)) %>%
      mutate(check = ifelse(a2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(a3 = ifelse(check, max(a), a)) %>%
      mutate(a = a3) %>%
      
      mutate(b2 = sum(b)) %>%
      mutate(check = ifelse(b2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(b3 = ifelse(check, max(b), b)) %>%
      mutate(b = b3) %>%
      
      mutate(c2 = sum(c)) %>%
      mutate(check = ifelse(c2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(c3 = ifelse(check, max(c), c)) %>%
      mutate(c = c3) %>%
      
      mutate(d2 = sum(d)) %>%
      mutate(check = ifelse(d2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(d3 = ifelse(check, max(d), d)) %>%
      mutate(d = d3) %>%
      
      mutate(e2 = sum(e)) %>%
      mutate(check = ifelse(e2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(e3 = ifelse(check, max(e), e)) %>%
      mutate(e = e3) %>%
      
      mutate(f2 = sum(f)) %>%
      mutate(check = ifelse(f2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(f3 = ifelse(check, max(f), f)) %>%
      mutate(f = f3) %>%
      
      mutate(g2 = sum(g)) %>%
      mutate(check = ifelse(g2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(g3 = ifelse(check, max(g), g)) %>%
      mutate(g = g3) %>%
      
      mutate(h2 = sum(h)) %>%
      mutate(check = ifelse(h2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(h3 = ifelse(check, max(h), h)) %>%
      mutate(h = h3) %>%
      
      mutate(i2 = sum(i)) %>%
      mutate(check = ifelse(i2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(i3 = ifelse(check, max(i), i)) %>%
      mutate(i = i3) %>%
      
      mutate(j2 = sum(j)) %>%
      mutate(check = ifelse(j2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(j3 = ifelse(check, max(j), j)) %>%
      mutate(j = j3) %>%
      
      mutate(k2 = sum(k)) %>%
      mutate(check = ifelse(k2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(k3 = ifelse(check, max(k), k)) %>%
      mutate(k = k3) %>%
      
      mutate(l2 = sum(l)) %>%
      mutate(check = ifelse(l2 >= min_species_for_community_sampling_event, TRUE, FALSE)) %>%
      mutate(l3 = ifelse(check, max(l), l)) %>%
      mutate(l = l3) 
    
    # now also pull out our detections
    V_NA <- V_NA %>%
      rename("1" = "a", "2" = "b", 
             "3" = "c", "4" = "d", 
             "5" = "e", "6" = "f", 
             "7" = "g", "8" = "h", 
             "9" = "i", "10" = "j", 
             "11" = "k", "12" = "l")
    
    # and then pull out the NAs stuff as a simple matrix
    V_NA <- as.matrix(V_NA[3:14])
    
    # now also pull out our detections
    detections_df <- detections_df %>%
      rename("1" = "a", "2" = "b", 
             "3" = "c", "4" = "d", 
             "5" = "e", "6" = "f", 
             "7" = "g", "8" = "h", 
             "9" = "i", "10" = "j", 
             "11" = "k", "12" = "l")
    
    V <- as.matrix(detections_df[3:14])
    
    # print mismatches in NA array (Detections when the NA says we haven't surveyed)
    print(paste0("Number of NA / detection mismatches: ", length(which(V > V_NA))))
    
  } # end if/else family sampling
  
  # get an indicator of whether species was detected 1 or more times
  # at the site within the year. If this value is greater than 0 we know it must be present
  confirmed_occurrence <- rowSums(V) 
  
  # get data of which species is being considered on each row of V
  species_integer_vector <- detections_df %>%
    mutate(species_integer_vector = as.integer(as.factor(species))) %>%
    pull(species_integer_vector)
  
  city_integer_vector <- detections_df %>%
    mutate(city_integer_vector = as.integer(as.factor(city))) %>%
    pull(city_integer_vector)
  
  multicity_site_id_vector <- detections_df %>%
    mutate(multicity_site_id_vector = multicity_site_id) %>%
    pull(multicity_site_id_vector)
  
  city_id_vector <- detections_df %>%
    mutate(city_id_vector = as.integer(as.factor(city))) %>%
    pull(city_id_vector)
  
  site_survey_year_vector <- detections_df %>%
    mutate(site_survey_year_vector = site_survey_year) %>%
    pull(site_survey_year_vector)
  
  prev_index_vector <- detections_df %>%
    mutate(prev_index_vector = prev_index) %>%
    pull(prev_index_vector)
  
  R <- nrow(V)
  
  mean_species_per_event <- detections_df %>%
    mutate(detection_one_or_more = ifelse(rowSums(.[,4:15]) > 0, 1, 0)) %>%
    group_by(community_sample_id) %>%
    summarize(sum(detection_one_or_more)) %>%
    ungroup() %>%
    as.matrix(.)
  
  mean_species_per_comm_sampling_event <- mean(as.numeric(mean_species_per_event[,2]))
  
  ## --------------------------------------------------
  ## get species by cluster integer to allow a species - cluster random effect on detection and initial occurrence
  
  city <-c( "Atlanta",
            "Boston", 
            "Charlotte",
            "Chicago",
            "Dallas",
            "DC",
            "Denton",
            "Houston",
            "LA",
            "Minneapolis",
            "NYC",     
            "Philadelphia",
            "Raleigh",
            "SD",
            "SF")
  cluster <-c( "southeast", # atlanta
               "northeast", # boston
               "southeast", # charlotte
               "midwest", # chicago
               "texas", # dallas
               "northeast", # dc
               "texas", # denton
               "texas", # houston
               "california", # LA
               "midwest", # minneapolis
               "northeast", # nyc
               "northeast", # philadelphia
               "southeast", # raleigh
               "california", # sd
               "san_francisco") # sf
  x_name <- "city"
  y_name <- "cluster"
  
  temp <- data.frame(city,cluster)
  names(temp) <- c(x_name,y_name)
  
  species_cluster <- left_join(ranges_raw, temp) %>%
    select(city, species, cluster)
  
  detections_df <- left_join(detections_df, species_cluster) %>%
    group_by(species, cluster) %>%
    mutate(species_cluster = cur_group_id()) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(region_cluster = cur_group_id()) %>%
    ungroup()
  
  # get data of which species*cluster combination is being considered on each row of V
  species_cluster_integer_vector <- detections_df %>%
    mutate(species_cluster_integer_vector = as.integer(as.factor(species_cluster))) %>%
    pull(species_cluster_integer_vector)
  
  # also get the species*cluster data in a form that is linked to species names
  # for making predictive plots after fitting the models
  species_region_cluster_id <- detections_df %>%
    mutate(species_cluster_integer_vector = as.integer(as.factor(species_cluster))) %>%
    select(species, cluster, region_cluster, species_cluster_integer_vector) %>%
    group_by(species_cluster_integer_vector) %>%
    slice(1) %>% ungroup
  
  # get data of which cluster is being considered on each row of V
  region_cluster_integer_vector <- detections_df %>%
    mutate(region_cluster_integer_vector = as.integer(as.factor(region_cluster))) %>%
    pull(region_cluster_integer_vector)
  
  species_city_cluster <- detections_df %>%
    group_by(city, species) %>%
    mutate(species_city_cluster = cur_group_id()) %>%
    ungroup() %>%
    pull(species_city_cluster)
  
  n_species_city_clusters <- length(unique(species_city_cluster))
  
  ## --------------------------------------------------
  ## get species trait data
  
  # species morpho traits
  trait_df <- read.csv(
    "./data/lepidoptera_trait_data/lepidoptera_trait_data/SpeciesListForTraits.csv") %>%
    select(GBIF_species, eButterfly_species, LepTraits_name, featureDiversity, aveWingspan, is_migratory) %>%
    
    # get a binary version of migratory status.
    # assume that a species unknown or not considered (blank) is non-migratory
    mutate(migratory = if_else(is_migratory == TRUE, 1, 0))
  
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
    left_join(., select(df, species, family), by=c("GBIF_species" = "species")) %>%
    group_by(family) %>%
    mutate(featureDiversityFamily = mean(featureDiversity),
           aveWingspanFamily = mean(aveWingspan, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    select(family, featureDiversityFamily, aveWingspanFamily)
  
  # in the event that a species has no data, this will replace the NA trait data
  # from the highest resolution taxonomic group with available data.
  # e.g. no data for species wingspan, replace species wingspan with the average wingspan
  # for all species from the same genus.
  species_info <- left_join(species_info, trait_df, by=c("species" = "GBIF_species1")) %>%
    left_join(., genus_id) %>%
    mutate(featureDiversity = ifelse(is.na(featureDiversity), featureDiversityGenus, featureDiversity),
           aveWingspan = ifelse(is.na(aveWingspan), aveWingspanGenus, aveWingspan)) %>%
    left_join(., family_id) %>%
    mutate(featureDiversity = ifelse(is.na(featureDiversity), featureDiversityFamily, featureDiversity),
           aveWingspan = ifelse(is.na(aveWingspan), aveWingspanFamily, aveWingspan))
  species_info[species_info == ""] <- NA
  
  ## now add ease of identification info (prop research grade detections in iNaturalist)
  ease_id_df <- read.csv("./data/lepidoptera_trait_data/ease_of_id/identifiability_by_genus.csv") %>%
    select(genus, research_grade_proportion)
  
  species_info <- left_join(species_info, ease_id_df)
  
  ## now add habitat affinity from LepTraits
  consensus <- read.csv("./data/lepidoptera_trait_data/leptraits/consensus.csv") %>%
    select(Species, FlightDuration, DiapauseStage, Voltinism, OvipositionStyle, 
           CanopyAffinity, EdgeAffinity, MoistureAffinity, DisturbanceAffinity,
           NumberOfHostplantFamilies) %>%
    rename("LepTraits_name" = "Species") %>%
    # when there are other names for the species the row will be duplicated but it 
    # appears that the first row is always the species in the accepted sense, e.g.:
    # Apodemia virgulti row 1 is Apodemia virgulti
    # Apodemia virgulti row 2 is Apodemia mormo. 
    # we want the traits for the closest definition of the species so always take row 1 here.
    group_by(LepTraits_name) %>%
    slice(1) %>%
    ungroup()
  
  #species_info <- left_join(species_info, consensus, by="LepTraits_name")
  species_info <- species_info %>%
    mutate(trait_join_name = ifelse(is.na(LepTraits_name), species, LepTraits_name))
  species_info <- left_join(species_info, consensus, by=c("trait_join_name" = "LepTraits_name")) 

  # now scale and select all the variables 
  species_info <- species_info %>%
    mutate(aveWingspan_scaled = center_scale(aveWingspan),
           research_grade_proportion_scaled = center_scale(research_grade_proportion),
           featureDiversity_scaled = center_scale(featureDiversity)) %>%
    select(species, genus, family, n_binary_detections, n_detections,
           aveWingspan, aveWingspan_scaled, 
           featureDiversity, featureDiversity_scaled, 
           research_grade_proportion, research_grade_proportion_scaled,
           migratory,
           FlightDuration, DiapauseStage, Voltinism, OvipositionStyle, 
           CanopyAffinity, EdgeAffinity, MoistureAffinity, DisturbanceAffinity,
           NumberOfHostplantFamilies)
  
  species_info_plot <- species_info %>%
    mutate(cond1 = ifelse(aveWingspan_scaled > 0, 0, 1),
           cond2 = ifelse(featureDiversity_scaled > 0, 0, 1),
           cond3 = ifelse(research_grade_proportion_scaled > 0, 0, 1))
  
  # plot the species trait data
  ggplot(species_info_plot, aes(x=as.factor(species), y=aveWingspan_scaled, 
                                fill = as.factor(cond1))) +
    geom_col() +
    theme_classic() +
    ylab("Wingspan (scaled)") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x=element_blank()) 
  
  ggplot(species_info_plot, aes(x=as.factor(species), y=featureDiversity_scaled, 
                                fill = as.factor(cond2))) +
    geom_col() +
    theme_classic() +
    ylab("Feature diversity (scaled)") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x=element_blank()) 
  
  ggplot(species_info_plot, aes(x=as.factor(species), y=research_grade_proportion_scaled, 
                                fill = as.factor(cond3))) +
    geom_col() +
    theme_classic() +
    ylab("Ease of ID (scaled)") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x=element_blank()) 
  
  # how many species were detected?
  # and how many binary detections occurred? species/site/year/visit
  n_species <- length(species_names_vector <- pull(species_info %>%
                                               # group by species ID
                                               group_by(species) %>%
                                               # and take one record
                                               slice(1) %>%
                                               select(species)))
  
  # get vector of families for each species
  n_families <- length(unique((family_names_vector <- pull(species_info %>%
                                                       # group by species ID
                                                       group_by(species) %>%
                                                       # grab one row per species ID
                                                       slice(1) %>%
                                                       # and then pull out it's family
                                                       select(family)))))
  
  
  ## --------------------------------------------------
  # summarize some city-wide metrics
  
  # how many community sampling events
  city_data <- as.data.frame(
    cbind(V_NA, 
          detections_df$community_sample_id, detections_df$city
          )) %>%
    rename("community_sample_id" = "V13",
           "city" = "V14") %>%
    group_by(community_sample_id) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city) %>%
    add_tally() %>%
    slice(1) %>%
    ungroup() %>%
    select(city, n) %>% 
    rename("total_community_sampling_events" = "n")
  
  # how many single v repeat sampling events per city
  city_repeat_events <- detections_df %>%
    group_by(community_sample_id) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(multicity_site_id, family) %>%
    add_tally() %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city) %>%
    mutate(single_year_only_events = length(which((n == 1)))) %>%
    select(city, single_year_only_events) %>%
    slice(1) %>%
    ungroup()
  
  # how many community sampling events by family
  city_family_events <- as.data.frame(
    cbind(V_NA, 
          detections_df$community_sample_id, detections_df$city, detections_df$family
    )) %>%
    rename("community_sample_id" = "V13",
           "city" = "V14",
           "family" = "V15") %>%
    group_by(community_sample_id) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city, family) %>%
    add_tally() %>%
    slice(1) %>%
    ungroup() %>%
    select(city, family, n) %>%  
    mutate(family = paste0(family, "_community_sampling_events")) %>% 
    rename("family_sampling_events" = "n") %>%
    pivot_wider(names_from = family, values_from = family_sampling_events)
  
  city_data <- left_join(city_data, city_repeat_events) %>%
    mutate(events_in_colonization_extinction_sequences = total_community_sampling_events - single_year_only_events)
  
  # start constructing some city stats df
  city_data <- left_join(city_data, city_family_events)
  
  # how many times individual sites surveyed repeatedly for same family
  city_survey_events <- as.data.frame(
    cbind(V_NA, 
          detections_df$community_sample_id, detections_df$city
    )) %>%
    rename("community_sample_id" = "V13",
           "city" = "V14") %>%
    group_by(community_sample_id) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city) %>%
    add_tally() %>%
    slice(1) %>%
    ungroup() %>%
    select(city, n) %>% 
    rename("total_community_sampling_events" = "n")
  
  # how many species detections in classified parks per city
  # and how many species detected >1 times per city
  total_annual_detections <- vector(length = nrow(detections_df))
  for(i in 1:nrow(detections_df)){
    total_annual_detections[i] <- sum(detections_df[i,4:15])
  }
  city_species <- as.data.frame(cbind(detections_df, total_annual_detections)) %>%
    group_by(city) %>%
    mutate(total_binary_detections_city = sum(total_annual_detections)) %>%
    ungroup() %>%
    mutate(any_detections = ifelse(total_annual_detections > 0, 1, 0)) %>%
    group_by(city, species) %>%
    filter(any_detections > 0) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(city) %>%
    add_tally() %>%
    rename("n_species_detected_city" = "n") %>%
    slice(1) %>%
    ungroup() %>%
    select(city, total_binary_detections_city, n_species_detected_city)
    
  city_data <- left_join(city_data, city_species)
  
  # get n sites per city
  # and park size/isolation of city
  city_site_data <- site_data %>%
    group_by(city) %>%
    mutate(n_sites = length(unique(multicity_site_id)),
           mean_log_park_size = mean(log_total_green_space_area),
           sd_park_size = sd(log_total_green_space_area),
           mean_log_isolation = mean(log(isolation)),
           sd_isolation = sd(log(isolation)),
           mean_tree_cover = mean(tree_percent_cover),
           sd_tree_cover = sd(tree_percent_cover),
           mean_plant_genera_density = mean(plant_genera_density),
           sd_plant_genera_density = sd(plant_genera_density)) %>%
    group_by(city) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, n_sites, mean_log_park_size, sd_park_size,
           mean_log_isolation, sd_isolation, mean_tree_cover, sd_tree_cover,
           mean_plant_genera_density, sd_plant_genera_density)
  
  city_data <- left_join(city_data, city_site_data)
  
  ## --------------------------------------------------
  # Return stuff
  return(list(
    
    V = V, # community science detection data
    V_NA = V_NA,
    R = R,

    species_info = species_info, # species sci name, family name, num detections and predictors
    site_data = site_data, # site names and predictors
    
    n_species = n_species, # number of species
    n_families = n_families, # number of Lepidoptera families
    n_sites = n_sites, # number of sites
    n_years = n_years, # number of years 
    n_surveys = n_surveys, # [max] number of surveys within year

    species_names_vector = species_names_vector,
    family_names_vector = family_names_vector, 
    sites = site_vector,
    years = year_vector,
    surveys = survey_vector,
    
    confirmed_occurrence = confirmed_occurrence,

    species_integer_vector = species_integer_vector,
    
    city_integer_vector = city_integer_vector,
    
    multicity_site_id_vector = multicity_site_id_vector,
    
    city_id_vector = city_id_vector, 
    
    site_survey_year_vector = site_survey_year_vector,
    
    prev_index_vector = prev_index_vector,
    
    species_cluster_integer_vector = species_cluster_integer_vector,
    
    region_cluster_integer_vector = region_cluster_integer_vector,
    
    city_data = city_data,

    species_region_cluster_id = species_region_cluster_id,
    
    species_city_cluster = species_city_cluster,
    
    n_species_city_clusters = n_species_city_clusters

  ))

}
