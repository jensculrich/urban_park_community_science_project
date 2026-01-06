// dynamic multi-species occupancy model for pollinators in urban parks
// this is for the full iNat data
// I'm using months as survey blocks and assuming all parks are surveyed in all months in all years
// so long as one or more species from corresponding family was detected
// I've included a month-phenology covariate for detection

// I've also added park site and species trait predictors (contrast with model 3)

// jcu, started may, 2024.

// throughout I denote dimensions as
// species == i
// site == j
// year == k
// visit == l

data {
  // dimensions of the data and the detections themselves
  int R; // length of the dataset (number of speciesXsiteXyear combos)
  int n_surveys; // number of repeat surveys within years
  array[n_surveys] real surveys; // surveys (difference from the mean) used as detection cov
  array[R, n_surveys] int<lower=0, upper=1> V; // binary detection / non detection
  array[R, n_surveys] int<lower=0, upper=1> V_NA; // sampling indicator 1==non-detection, 0==no evidence the species was sampled
  array[R] int<lower=1> site_survey_year_vector; // which year is the survey referencing?
  // species, site, and city indicators
  int<lower=0> n_species; // number of species
  array[R] int<lower=1, upper=n_species> species_integer_vector; // vector indicating which species is being observed
  int<lower=0> n_sites; // number of sites
  array[R] int<lower=1, upper=n_sites> multicity_site_id_vector; // vector indicating which site is being observed
  int<lower=0> n_cities; // number of cities
  array[R] int<lower=1, upper=n_cities> city_id_vector; // vector indicating which city is being observed 
  int<lower=0> n_species_clusters; // number of speciesXregions clusters
  array[R] int<lower=1, upper=n_species_clusters> species_cluster_id_vector;
  int<lower=0> n_regional_clusters; // number of speciesXregions clusters
  array[R] int<lower=1, upper=n_regional_clusters> regional_cluster_id_vector;
  // species and site covariate data
  vector[n_species] feature_diversity;
  vector[n_species] ease_of_id;
  vector[n_species] wingspan;
  vector[n_sites] park_size;
  vector[n_sites] tree_cover;
  vector[n_sites] landscape_isolation;
  vector[n_sites] landscape_grassherb;
  vector[n_sites] landscape_woody;
  vector[n_sites] total_detections_by_city;
  // other stuff
  array[R] int<lower=0> confirmed_occurrence;
} // end data

parameters {

  // initial state
  real psi_0; // intercept
  vector[n_species_clusters] psi_species_raw;  
  real<lower=0> sigma_psi_species;
  vector[n_cities] psi_city_raw;
  real<lower=0> sigma_psi_city;
  real mu_psi_wingspan;
  vector[n_cities] psi_wingspan_raw;  
  real<lower=0> sigma_psi_wingspan;
  real mu_psi_park_size;
  vector[n_cities] psi_park_size_raw;  
  real<lower=0> sigma_psi_park_size;
  real mu_psi_tree_cover;
  vector[n_cities] psi_tree_cover_raw;  
  real<lower=0> sigma_psi_tree_cover;
  real mu_psi_landscape_isolation;
  vector[n_cities] psi_landscape_isolation_raw;  
  real<lower=0> sigma_psi_landscape_isolation;
  real mu_psi_landscape_grassherb;
  vector[n_cities] psi_landscape_grassherb_raw;  
  real<lower=0> sigma_psi_landscape_grassherb;
  real mu_psi_landscape_woody;
  vector[n_cities] psi_landscape_woody_raw;  
  real<lower=0> sigma_psi_landscape_woody;

  // detection
  real p0; // intercept
  vector[n_species_clusters] p_species_raw;
  real<lower=0> sigma_p_species;
  vector[n_cities] p_city_raw;
  real<lower=0> sigma_p_city;
  real p_city_detections;
  real p_wingspan;
  real p_feature_diversity;
  real p_ease_of_id;
  vector[n_species_clusters] p_date; // phenology peak
  real delta0;
  vector[n_regional_clusters] delta_regional_cluster;
  real<lower=0> sigma_p_species_date; // variation
  vector[n_species_clusters] p_date_sq; // decay pattern of phenology
  real epsilon0;
  vector[n_regional_clusters] epsilon_regional_cluster;
  real<lower=0> sigma_p_species_date_sq; // community mean

} // end parameters

transformed parameters {

  // logit scale psi, gamma, phi
  array[R] real psi; // odds of occurrence
  array[R, n_surveys] real p; // odds of detection
  
  vector[n_species_clusters] psi_species;
  vector[n_cities] psi_city;
  vector[n_cities] psi_wingspan;
  vector[n_cities] psi_park_size;
  vector[n_cities] psi_tree_cover;
  vector[n_cities] psi_landscape_isolation;
  vector[n_cities] psi_landscape_grassherb;
  vector[n_cities] psi_landscape_woody;
  vector[n_species_clusters] p_species;
  vector[n_cities] p_city;
  vector[n_species_clusters] mu_p_species_date; // community mean
  vector[n_species_clusters] mu_p_species_date_sq; // variation
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi_species = sigma_psi_species * psi_species_raw;
  psi_city = psi_0 + sigma_psi_city * psi_city_raw;
  psi_wingspan = mu_psi_wingspan + sigma_psi_wingspan * psi_wingspan_raw;
  psi_park_size = mu_psi_park_size + sigma_psi_park_size * psi_park_size_raw;
  psi_tree_cover = mu_psi_tree_cover + sigma_psi_tree_cover * psi_tree_cover_raw;
  psi_landscape_isolation = mu_psi_landscape_isolation + sigma_psi_landscape_isolation * psi_landscape_isolation_raw;
  psi_landscape_grassherb = mu_psi_landscape_grassherb + sigma_psi_landscape_grassherb * psi_landscape_grassherb_raw;
  psi_landscape_woody = mu_psi_landscape_woody + sigma_psi_landscape_woody * psi_landscape_woody_raw;
  p_species = sigma_p_species * p_species_raw;
  p_city = p0 + sigma_p_city * p_city_raw;
  
  
  for(r in 1:R){
  
    mu_p_species_date[species_cluster_id_vector[r]] = delta0 + delta_regional_cluster[regional_cluster_id_vector[r]];
    mu_p_species_date_sq[species_cluster_id_vector[r]] = epsilon0 + epsilon_regional_cluster[regional_cluster_id_vector[r]];

  
    psi[r] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
      psi_city[city_id_vector[r]] +
      psi_species[species_cluster_id_vector[r]] + // a species specific intercept
      psi_wingspan[city_id_vector[r]] * wingspan[species_integer_vector[r]] + // a species effect of migratory
      psi_park_size[city_id_vector[r]] * park_size[multicity_site_id_vector[r]] + // a site effect of park size
      psi_tree_cover[city_id_vector[r]] * tree_cover[multicity_site_id_vector[r]] +
      psi_landscape_isolation[city_id_vector[r]] * landscape_isolation[multicity_site_id_vector[r]] + // a site effect of park isolation
      psi_landscape_grassherb[city_id_vector[r]] * landscape_grassherb[multicity_site_id_vector[r]] +
      psi_landscape_woody[city_id_vector[r]] * landscape_woody[multicity_site_id_vector[r]]
      ); // end psi[r]
      
  } // end loop across all data
  
  // have to do another loop because phi and gamma have a shorter k index length than p
  for(r in 1:R){
    for(l in 1:n_surveys){ // loop across all surveys

      p[r,l] = inv_logit( // probability (0-1) of detection is equal to..
        p_city[city_id_vector[r]] + 
        p_species[species_cluster_id_vector[r]] + // a species specific intercept
        p_city_detections * total_detections_by_city[multicity_site_id_vector[r]] +
        p_wingspan * wingspan[species_integer_vector[r]] + // a species effect of wingspan
        p_feature_diversity * feature_diversity[species_integer_vector[r]] + // a species effect of feature diversity
        p_ease_of_id * ease_of_id[species_integer_vector[r]] + // a species effect of ease of identification
        p_date[species_cluster_id_vector[r]] * surveys[l] + // a species-specific phenological detection effect (peak)
        p_date_sq[species_cluster_id_vector[r]] * (surveys[l])^2 // a species-specific phenological detection effect (decay)
        ); // end p[j,k,l]
        
    } // end loop across all surveys
  } // end loop across all R
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  // initial state
  psi_0 ~ normal(0, 1); // initial occurrence intercept
  psi_city_raw ~ std_normal();
  sigma_psi_city ~ normal(0, 1);
  psi_species_raw ~ std_normal();
  sigma_psi_species ~ normal(0, 1);
  mu_psi_wingspan ~ normal(0, 2);
  psi_wingspan_raw ~ std_normal();
  sigma_psi_wingspan ~ normal(0, 0.5);
  mu_psi_park_size ~ normal(0, 2);
  psi_park_size_raw ~ std_normal();
  sigma_psi_park_size ~ normal(0, 0.5);
  mu_psi_tree_cover ~ normal(0, 2);
  psi_tree_cover_raw ~ std_normal();
  sigma_psi_tree_cover ~ normal(0, 0.5);
  mu_psi_landscape_isolation ~ normal(0, 2);
  psi_landscape_isolation_raw ~ std_normal();
  sigma_psi_landscape_isolation ~ normal(0, 0.5);
  mu_psi_landscape_grassherb ~ normal(0, 2);
  psi_landscape_grassherb_raw ~ std_normal();
  sigma_psi_landscape_grassherb ~ normal(0, 0.5);
  mu_psi_landscape_woody ~ normal(0, 2);
  psi_landscape_woody_raw ~ std_normal();
  sigma_psi_landscape_woody ~ normal(0, 0.5);

  // detection
  p0 ~ normal(0, 1); // global intercept
  p_city_raw ~ std_normal();
  sigma_p_city ~ normal(0, 1);
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 1);
  p_city_detections ~ normal(0, 2);
  p_wingspan ~ normal(0, 2);
  p_feature_diversity ~ normal(0, 2);
  p_ease_of_id ~ normal(0, 2);
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  delta0 ~ normal(0, 2); // community mean
  delta_regional_cluster ~ normal(0, 1); // effect of region
  sigma_p_species_date ~ normal(0, 2); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  epsilon0 ~ normal(0, 1); // community mean
  epsilon_regional_cluster ~ normal(0, 0.5); // effect of region
  sigma_p_species_date_sq ~ normal(0, 1); // variation

  // LIKELIHOOD
  for(r in 1:R){ // for each potential species X site X year detection

    // model the occurrence at site multicity_site_id[r]
    // for the year site_survey_year[r], which (if site survey year > 1) 
    // will be contigent on occurrence at the same site in year site_survey_year[r - 1]
    
    // and also only for the surveys (months) with some survey effort (V_NA indicates evidence of survey effort)
    
    // e.g., target += (log(psi[i,multicity_site_id[r],site_survey_year[r]]))
      if(confirmed_occurrence[r] > 0){ // the species was ever detected across any surveys within the year of a community sampling event it must be present
          target += (log(psi[r]) + // present
                     // successful or failed detections across repeat surveys (12 months of year)
                     // contigent on whether a survey effort was recorded in that month'
                     // if no survey effort recorded, treat the month as NA
                     bernoulli_lpmf(V[r,1]|p[r,1])*V_NA[r,1] + 
                     bernoulli_lpmf(V[r,2]|p[r,2])*V_NA[r,2] + 
                     bernoulli_lpmf(V[r,3]|p[r,3])*V_NA[r,1] + 
                     bernoulli_lpmf(V[r,4]|p[r,4])*V_NA[r,4] +
                     bernoulli_lpmf(V[r,5]|p[r,5])*V_NA[r,5] + 
                     bernoulli_lpmf(V[r,6]|p[r,6])*V_NA[r,6] + 
                     bernoulli_lpmf(V[r,7]|p[r,7])*V_NA[r,7] + 
                     bernoulli_lpmf(V[r,8]|p[r,8])*V_NA[r,8] + 
                     bernoulli_lpmf(V[r,9]|p[r,9])*V_NA[r,9] + 
                     bernoulli_lpmf(V[r,10]|p[r,10])*V_NA[r,10] + 
                     bernoulli_lpmf(V[r,11]|p[r,11])*V_NA[r,11] + 
                     bernoulli_lpmf(V[r,12]|p[r,12])*V_NA[r,12]
        );
        
      } else{ // else the species was not detected during the comm survey event and could either be present or absent
          target += (// failed detections of an occurring species across all 12 months
                       log_sum_exp(log(psi[r]) + 
                       log1m(p[r,1])*V_NA[r,1] + 
                       log1m(p[r,2])*V_NA[r,2] + 
                       log1m(p[r,3])*V_NA[r,3] + 
                       log1m(p[r,4])*V_NA[r,4] +
                       log1m(p[r,5])*V_NA[r,5] + 
                       log1m(p[r,6])*V_NA[r,6] + 
                       log1m(p[r,7])*V_NA[r,7] + 
                       log1m(p[r,8])*V_NA[r,8] + 
                       log1m(p[r,9])*V_NA[r,9] + 
                       log1m(p[r,10])*V_NA[r,10] + 
                       log1m(p[r,11])*V_NA[r,11] + 
                       log1m(p[r,12])*V_NA[r,12],
                      // or just simple no occurrence
                      log1m(psi[r])));
      } // end if/else
    
  }  // end loop across all community sampling events

} // end model

generated quantities{
  
  //
  // posterior predictive check (number of detections, binned by city)
  //
  array[n_cities] int<lower=0> W_city_rep; // sum of simulated detections

  array[R] int z_simmed; // simulate occurrence

  for(r in 1:R){
    z_simmed[r] = bernoulli_rng(psi[r]); 
  }
  
  // initialize repped detections at 0
  for(i in 1:n_cities){
    W_city_rep[i] = 0;
  }
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(r in 1:R) { // loop across all site/year/species identities
    for(l in 1:n_surveys){ // loop across surveys
          
          // detections in replicated data (us z_simmed from above)
          W_city_rep[city_id_vector[r]] = W_city_rep[city_id_vector[r]] + 
            // multiply occupancy state by a simulated detection, AND...
            // multiply by the NA indicator - if we didn't survey in real life
            // we don't survey in this simulation.
            (z_simmed[r] * bernoulli_rng(p[r,l]) * V_NA[r,l]);
           
    } // end loop across surveys
  } // end loop across site/year/species identities
  
} // end generated quantities
