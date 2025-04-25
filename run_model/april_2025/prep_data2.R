# create a detection array from iNaturalist data.
# this will join detections to sites.
# non-detections are currently added whenever we don't see a species.
# sampling effort could then be modeled as a function of species/site/year/time,
# using the months within year as repeat "survey periods".
# but possibly we could look to identify "community sampling events"
# and only model non-detections as those species not observed on a given community sampling event?

library(tidyverse) # data organization
# library(sf) # spatial data processing
# library(wesanderson) # colour palettes

prep_data <- function(min_species_detections,
                      min_species_for_community_sampling_event,
                      family_sampling
                      ) {
  
  ## --------------------------------------------------
  ## Operation Functions
  ## predictor center scaling function
  center_scale <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  #-----------------------------------------------------
  # load the detection data
  
  # first read the data 
  df <- read.csv(
    "./data/LA_county_observations_park_coord_clipped.csv") 
  
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
  filter(family %in% butterfly_families)
  
  # genus <- unique(df$genus)
  
  #write.csv(genus, "data/lepidoptera_trait_data/ease_of_id/identifiability_by_genus.csv")
  
  #-----------------------------------------------------
  # now we are going to filter out all of the points from unclassified parks 
  # (parks outside city boundaries)
  # you'd need to do this later if you want to model responses of species occurring in
  # classified parks but never in unclassified parks
  
  df <- df %>%
    filter(type == "classified")
  
  #-----------------------------------------------------
  # prep data for array format
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
    mutate(month = ifelse(month %% 2 == 0, month - 1, month)) %>% 
  
    # add survey date within year
    group_by(year) %>% 
    mutate(survey = as.integer(factor(month)),
           year = as.integer(year - 2019)) %>% # used (- 2019) to make 2020 == year 1
    ungroup() %>%
    
    # for now, reducing down to mandatory data columns
    dplyr::select(species, family, new_id, 
                  total_green_space_area, tree_percent_cover, grass_shrub__percent_cover,
                  survey, year, n_detections) %>%
    
    # turn into binary detections (for occupancy rather than abundance)
    group_by(species, new_id, year, survey) %>% 
    slice(1) %>%
    ungroup() %>%
    
    # add tally of binary detections per species
    group_by(species) %>%
    add_tally(name = 'n_binary_detections') %>%
    ungroup() %>%
    
    # arrange by survey within year  
    arrange(year, survey) 
  
  # get dimensions of surveys and years
  survey_vector <- as.vector(levels(as.factor(df$survey)))
  n_surveys <- length(survey_vector)
  
  year_vector <- as.vector(levels(as.factor(df$year)))
  n_years <- length(year_vector)
  
  # start a "species_info" table
  species_info <- df  %>%
    mutate(genus = word(species, 1)) %>%
    select(species, genus, family, n_detections, n_binary_detections) %>%
    group_by(species) %>%
    slice(1) %>%
    ungroup()
  
  ## --------------------------------------------------
  ## get site covariate data
  
  # n sites and site vector
  n_sites <- (nrow(site_data <- df %>%
                     group_by(new_id) %>%
                     slice(1) %>%
                     ungroup() %>%
                     select(new_id,
                            total_green_space_area, 
                            tree_percent_cover, 
                            grass_shrub__percent_cover)  
                    ))
  
  # get a vector of site names
  site_vector <- site_data %>%
    pull(new_id)
  
  # get scaled versions of the site covariate data
  site_data <- site_data %>%
    mutate(total_green_space_area = center_scale(log(total_green_space_area)),
           tree_cover_scaled = center_scale(tree_percent_cover),
           grass_shrub_cover_scaled = center_scale(grass_shrub__percent_cover))
  
  # add connectivity data
  connectivity <- read.csv("./data/connectivity.csv") 
  
  # plot the spread of the site covariate data
  par(mfrow=c(1,3))  
  hist(site_data$total_green_space_area)
  hist(site_data$tree_cover_scaled)
  hist(site_data$grass_shrub_cover_scaled)
  par(mfrow=c(1,1)) 
  
  # how correlated are the site covariate data
  cor(site_data$total_green_space_area, site_data$tree_cover_scaled)
  cor(site_data$grass_shrub_cover_scaled, site_data$tree_cover_scaled)
  cor(site_data$grass_shrub_cover_scaled, site_data$total_green_space_area)
  
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
    left_join(., select(df, species, family), by=c("GBIF_species" = "species")) %>%
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
    select(species, genus, family, n_binary_detections, n_detections,
           aveWingspan_scaled, featureDiversity_scaled, research_grade_proportion_scaled)
  
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
  
  species_info_filtered <- species_info %>%
    filter(n_binary_detections >= min_species_detections)
  
  # how many species were detected?
  # and how many binary detections occurred? species/site/year/visit
  n_species <- length(species_vector <- pull(species_info_filtered %>%
                                               # group by species ID
                                               group_by(species) %>%
                                               # and take one record
                                               slice(1) %>%
                                               select(species)))
  
  # get vector of families for each species
  n_families <- length(unique((family_vector <- pull(species_info_filtered %>%
                                                       # group by species ID
                                                       group_by(species) %>%
                                                       # grab one row per species ID
                                                       slice(1) %>%
                                                       # and then pull out it's family
                                                       select(family)))))
  
  ## --------------------------------------------------
  ## Now lets make a detection matrix 
  
  ## --------------------------------------------------
  # Generate an indicator array for whether or not a community-wide  
  # sampling event occurred at the site in a year
  
  # make a 4 dimensional array
  V <- array(data = NA, dim = c(n_species, n_sites, n_years, n_surveys))
  
  for(k in 1:n_years){
    for(l in 1:n_surveys){
      
      # iterate this across surveys within years
      temp <- df %>%
        
        # don't need the family column for this
        dplyr::select(-family) %>%
        
        # filter to indices for year and survey
        filter(year == (k), 
               survey == (l)) %>% 
        
        # now join with all species (so that we include species not captured during 
        # this interval*visit but which might actually be at some sites undetected)
        full_join(., dplyr::select(species_info, c(species, n_binary_detections)), by=c("species", "n_binary_detections")) %>%
        
        # but then filter out the species we never plan to model
        filter(n_binary_detections >= min_species_detections) %>%
        
        dplyr::select(species, new_id) %>% # if you don't deselect the other variables it expands the temp matrix for every unique one
        
        # now join with all sites columns (so that we include sites where no species captured during 
        # this interval*visit but which might actually have some species that went undetected)
        #full_join(., select(site_names, grid_id), by="grid_id") %>%
        full_join(., as.data.frame(site_vector), by=c("new_id" = "site_vector")) %>% 
        
        # group by SPECIES
        group_by(species) %>%
        mutate(row = row_number()) %>%
        # spread sites by species, and fill with 0 if species never captured this interval*visit
        spread(new_id, row, fill = 0) %>%
        
        # replace number of unique site captures of the species (if > 1) with 1.
        mutate_at(2:(ncol(.)), ~replace(., . > 1, 1)) %>%
        # if more columns are added these indices above^ might need to change
        # 5:(n_species+5) represent the columns of each site in the matrix
        # just need the matrix of 0's and 1's
        #dplyr::select(2:(ncol(.))) %>%
        # if some sites had no species, this workflow will construct a row for species = NA
        # we want to filter out this row ONLY if this happens and so need to filter out rows
        # for SPECIES not in SPECIES list
        filter(species %in% levels(as.factor(species_info$species)))
      
      # remove the NA column at the end (why is this popping up? happens if no species detected at any sites)
      #last_column = ncol(temp_matrix)
      #temp_matrix <- temp_matrix[,-(last_column)]
      temp <- as.data.frame(temp)
      temp <- select(temp, -"<NA>", -"species")
      # convert from dataframe to matrix
      temp_matrix <- as.matrix(temp)
      
      # replace NAs for the interval i and visit j with the matrix
      V[1:n_species, 1:n_sites,k,l] <- temp_matrix[1:n_species, 1:n_sites]
      
    }
  }
  
  class(V) <- "numeric"
  
  #-----------------------------------------------------
  # identify community sampling events which we will use to infer non-detections
  # identified by site*survey months with > 1 species detected. Another option
  # might be to group by iNat user name by date*park but this may be a very narrow scope?
  
  if(family_sampling == TRUE){
    # determine sampling based on sampling from same family 
    community_samples <- df %>%
      
      # determine whether a community sampling event occurred
      # using sampling information from other species in same survey interval 
      group_by(family, new_id, year, survey) %>%
      mutate(n_species_sampled = n_distinct(species)) %>%
      ungroup() %>%
      
      # if we detected more than min_species_for_community_sampling_event at the site
      # in the "survey" then we assume a community sampling event occured
      mutate(community_sampled = ifelse(
        n_species_sampled >= min_species_for_community_sampling_event,
        1, 0)) %>%
      dplyr::select(-species, -n_species_sampled) %>%
      
      # we just need one row per community sampling event (not a row for all species positively detected)
      group_by(new_id, family, year, survey) %>%
      slice(1) %>%
      ungroup()
  } else{
    # determine sampling based on sampling from any butterfly family
    community_samples <- df %>%
      
      # determine whether a community sampling event occurred
      # using sampling information from other species in same survey interval 
      group_by(new_id, year, survey) %>%
      mutate(n_species_sampled = n_distinct(species)) %>%
      ungroup() %>%
      
      # if we detected more than min_species_for_community_sampling_event at the site
      # in the "survey" then we assume a community sampling event occured
      mutate(community_sampled = ifelse(
        n_species_sampled >= min_species_for_community_sampling_event,
        1, 0)) %>%
      dplyr::select(-species, -n_species_sampled) %>%
      
      # we just need one row per community sampling event (not a row for all species positively detected)
      group_by(new_id, family, year, survey) %>%
      slice(1) %>%
      ungroup()
  }
  
  
  ## --------------------------------------------------
  # Generate an indicator array for whether or not a community-wide museum 
  # sampling event occurred at the site in a year
  
  # we will use all of the species (whether we model them or not) to infer sampling
  #total_n_species = nrow(species_info)
  #all_species_vector = unique(species_info$species)
  
  # first make a df of all possible site visits * species visits
  all_species_site_visits <- as.data.frame(cbind(
    rep(species_vector, times=(n_sites*n_years*n_surveys)),
    rep(site_vector, each=n_species, times=n_years*n_surveys),
    rep(year_vector, each=n_species*n_sites, times=(n_surveys)),
    rep(survey_vector, each=n_species*n_sites*n_years),
    rep(family_vector, times=(n_sites*n_years*n_surveys))
  )) %>%
    rename("species" = "V1",
           "new_id" = "V2",
           "year" = "V3",
           "survey" = "V4",
           "family" = "V5") %>%
    mutate(year = as.integer(year),
           survey = as.integer(survey),
           new_id = as.integer(new_id))
  
  # now propagate nondetections to any species 
  # this assumes that ALL other species could have been detected
  # could revise to propagate nondetection only to other species in same suborder, family, genus etc.
  all_visits_joined <- left_join(
    all_species_site_visits, community_samples, 
    by=c("new_id","family", "year", "survey")) %>%
    
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
    
    V_detections = V, # community science detection data
    V_NA = V_NA,

    species_info = species_info_filtered, # species sci name, family name, num detections and predictors
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
