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
                      remove_outlier_parks
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
  source("./run_model/get_site_data.R")
  
  site_data_temp <- get_site_data(city_names)
  
  site_data <- df %>%
    group_by(city, new_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, new_id)
  
  site_data <- left_join(site_data, site_data_temp, by = c("city", "new_id"))
  
  site_data <- site_data %>%
    group_by(city) %>%
    mutate(isolation_scaled_2 = center_scale(isolation),
           log_total_green_space_area_scaled_2 = center_scale(log_total_green_space_area),
           log_isolation_scaled_2 = center_scale(log(isolation))) %>%
    ungroup()
  
  if(remove_outlier_parks == TRUE){
    site_data <- site_data %>%
      filter(log_total_green_space_area_scaled_2 > -3) 
  }
  
  n_sites <- nrow(site_data <- site_data %>%
    mutate(multicity_site_id = row_number()))
    
  # get a vector of site names
  site_vector <- site_data %>%
    pull(multicity_site_id)
  
  # plot the spread of the site covariate data
  par(mfrow=c(1,4))  
  hist(site_data$log_total_green_space_area_scaled)
  hist(site_data$log_isolation_scaled)  
  hist(site_data$log_total_green_space_area_scaled_2)
  hist(site_data$log_isolation_scaled_2)
  #hist(site_data$perarea_idx_scaled)
  #hist(site_data$proximity_scaled)
  
  par(mfrow=c(1,1)) 
  
  # how correlated are the site covariate data
  cor(site_data$log_total_green_space_area_scaled, site_data$log_isolation_scaled)
  cor(site_data$log_total_green_space_area_scaled_2, site_data$log_isolation_scaled_2)
  
  plot(site_data$log_total_green_space_area_scaled_2, site_data$log_isolation_scaled_2)
  
  ## --------------------------------------------------
  # identify occassions where sampling occurred
  community_sampling_event_site_years <- df %>%
    group_by(city, new_id, year) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, new_id, year)
  
  # join the multicity site id and filter out any sites that were filtered out
  # due to the site data, e.g. park was too small to include.
  community_sampling_event_site_years <- community_sampling_event_site_years %>%
    left_join(., select(site_data, city, new_id, multicity_site_id)) %>%
    filter(!is.na(multicity_site_id)) %>%
    
  # Now we want to do two things that will allow us to model a temporal autocorrelation
  # in the data (i.e., persistence/colonization that depends on the previous state)
    
  # first we want to know if this is the first year of survey for a park site
    # and add an indicator that will allow us to model psi1 effects for that year
  # then, for any years that follow we want to know what row number corresponds 
    # to the previous survey. We will model persistence/colonization based on inferred
    # presence/absence at the previous surveyed time step. 
  # IMPORTANTLY, in this framework, the previous time step could be anywhere between
    # 1 - 4 years earlier.
    
  # so, first, find minimum year of survey for a park and tag this year as first year
    group_by(multicity_site_id) %>%
    mutate(first_year = ifelse(year == min(year), 1, 0)) %>%
  # and a unique survey year per site
    mutate(site_survey_year = as.integer(as.factor(year))) %>%
    ungroup() %>%
  
  # now add a unique community sample event id
    mutate(community_sample_id = row_number())
  
  ## --------------------------------------------------
  # get detection data for all community sampling events r:R
  # detection data will be cast into vector of lenght R (total number of detections)
  
  community_sample_id <- community_sampling_event_site_years$community_sample_id
  
  # in a moment we will need a list of all species for all survey events
  # representing species that could have been detected when other species were
  all_comm_surveys_all_species <- as.data.frame(cbind(rep(community_sample_id, each=length(species_info$species)),
                                                      rep(species_info$species, times=length(community_sample_id)))) %>%
    rename("community_sample_id" = "V1", "species" = "V2") %>%
    mutate(community_sample_id = as.integer(community_sample_id))
  
  # get names of all surveys (in case a cluster doesn't have them all surveyed)
  # this will avoid having to rewrite the stan model code for different # max survey lengths
  temp <- as.data.frame(survey_vector) %>%
    rename("survey" = "survey_vector")
  
  detections_df <- df %>%
    # join the data about the community sampling event from that year
    left_join(., community_sampling_event_site_years) %>%
    # filter out is na - removes rows added in to cover all survey columns     
    filter(!is.na(community_sample_id)) %>%
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
    full_join(., all_comm_surveys_all_species) %>%
    # now we only want community_sample_id, species, and the matrix of presence absences
    select(community_sample_id, species, "1", "2", "3", "4", "5", "6", 
           "7", "8", "9", "10", "11", "12")
   
  # now we need to join all the data about site number / year for all those undetected species
  detections_df <- detections_df %>%
    left_join(., community_sampling_event_site_years, by = "community_sample_id") %>%
    arrange(community_sample_id, species)
  
  # now we need to extract survey-specific NA info
  # i.e., was any survey effort conducted 
  detections_df[,3:14][is.na(detections_df[,3:14])] <- 0
  
  V_NA <- detections_df %>%
    group_by(community_sample_id) %>%
    pivot_longer(cols=c("1", "2", "3", "4", "5", "6", 
                        "7", "8", "9", "10", "11", "12"), 
                 names_to="survey", values_to="V_NA") %>%
    group_by(community_sample_id, survey) %>%
    filter(V_NA == 1) %>%
    # join all survey events within a year
    mutate(survey = as.integer(survey)) %>%
    full_join(., temp) %>%
    # arrange surveys from 1 - 12
    arrange(match(survey, survey_vector)) %>%
    # now spread the NA data wide again
    pivot_wider(., names_from = survey, values_from = V_NA) %>%
    # now just take one row for the community (don't need one for every species)
    group_by(community_sample_id) %>%
    slice(1) %>%
    ungroup() %>%
    # filter out is na - removes rows added in to cover all survey columns     
    filter(!is.na(community_sample_id)) %>%
    # and arrange by community sampling event id
    arrange(., community_sample_id)
    
  # now replace NAs in the NA matrix with 0 which will serve to represent no survey effort
  V_NA[,9:20][is.na(V_NA[,9:20])] <- 0
  # and then pull out the NAs stuff as a simple matrix
  V_NA <- as.matrix(V_NA[9:20])

  # now also pull out our detections
  V <- as.matrix(detections_df[3:14])
  
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
  
  community_sample_id_vector <- detections_df %>%
    mutate(community_sample_id_vector = community_sample_id) %>%
    pull(community_sample_id_vector)
  
  R <- nrow(V)
  
  ## --------------------------------------------------
  ## get species trait data
  
  # species morpho traits
  trait_df <- read.csv(
    "./data/lepidoptera_trait_data/lepidoptera_trait_data/SpeciesListForTraits.csv") %>%
    select(GBIF_species, eButterfly_species, LepTraits_name, featureDiversity, aveWingspan) 
  
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
  
  species_info <- left_join(species_info, consensus, by="LepTraits_name")
  
  # now scale and select all the variables 
  species_info <- species_info %>%
    mutate(aveWingspan_scaled = center_scale(aveWingspan),
           research_grade_proportion_scaled = center_scale(research_grade_proportion),
           featureDiversity_scaled = center_scale(featureDiversity)) %>%
    select(species, genus, family, n_binary_detections, n_detections,
           aveWingspan, aveWingspan_scaled, 
           featureDiversity, featureDiversity_scaled, 
           research_grade_proportion, research_grade_proportion_scaled,
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
  # Get species ranges
  
  # source the prep function
  source("./run_model/get_species_ranges.R")
  
  ranges_raw <- get_species_ranges(city_names)
  
  all_species_city_combos <- as.data.frame(cbind(
    rep(species_info$species, times = length(city_names)),
    rep(city_names, each = n_species)
  )) %>%
    rename("species" = "V1",
           "city" = "V2")
  
  temp <- select(species_info, species)
  temp <- left_join(temp, ranges_raw, by = "species") %>%
    mutate(in_range = 1)
  
  ranges <- full_join(temp, all_species_city_combos) %>%
    mutate(in_range = replace_na(in_range, 0))
  
  temp2 <- select(site_data, city, multicity_site_id)
  
  ranges <- left_join(temp2, ranges, by = "city")
  
  # now spread into 4 dimensions
  #sort data frame by multiple columns alphabetically
  ranges <- ranges[with(ranges, order(multicity_site_id, species)), ]
  ranges <- array(data = ranges$in_range, dim = c(n_species, n_sites ))
  
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
    
    ranges = ranges,
    
    confirmed_occurrence = confirmed_occurrence,

    species_integer_vector = species_integer_vector,
    
    city_integer_vector = city_integer_vector,
    
    multicity_site_id_vector = multicity_site_id_vector,
    
    city_id_vector = city_id_vector, 
    
    site_survey_year_vector = site_survey_year_vector,
    
    community_sample_id_vector = community_sample_id_vector

  ))

}
