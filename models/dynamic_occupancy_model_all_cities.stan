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
  //array[n_species] int<lower=1> species; // vector of each species identity
  array[R] int<lower=1, upper=n_species> species_integer_vector; // vector indicating which species is being observed
  int<lower=0> n_sites; // number of sites
  //array[n_sites] int<lower=1> sites; // vector of site identities
  array[R] int<lower=1, upper=n_sites> multicity_site_id_vector; // vector indicating which site is being observed
  int<lower=0> n_cities; // number of cities
  //array[n_cities] int<lower=1> city; // vector of city identities
  array[R] int<lower=1, upper=n_cities> city_id_vector; // vector indicating which city is being observed 
  int<lower=0> n_species_clusters; // number of speciesXregions clusters
  //array[n_species_clusters] int<lower=1> species_cluster;
  array[R] int<lower=1, upper=n_species_clusters> species_cluster_id_vector;
  // species and site covariate data
  vector[n_species] feature_diversity;
  vector[n_species] ease_of_id;
  vector[n_species] wingspan;
  vector[n_sites] park_size;
  vector[n_sites] isolation;
  vector[n_sites] latitude;
  // other stuff
  //array[n_species, n_sites] int<lower=0> ranges; // matrix to constrain analysis within species ranges
  array[R] int<lower=0> confirmed_occurrence;
  array[R] int<lower=0> prev_index_vector;
} // end data

parameters {

  // initial state
  real psi1_0; // intercept
  vector[n_species_clusters] psi1_species_raw;  
  real<lower=0> sigma_psi1_species;
  vector[n_cities] psi1_city_raw;
  real<lower=0> sigma_psi1_city;
  real mu_psi1_wingspan;
  vector[n_cities] psi1_wingspan_raw;  
  real<lower=0> sigma_psi1_wingspan;
  real mu_psi1_park_size;
  vector[n_cities] psi1_park_size_raw;  
  real<lower=0> sigma_psi1_park_size;
  real mu_psi1_isolation;
  vector[n_cities] psi1_isolation_raw;  
  real<lower=0> sigma_psi1_isolation;

  // colonization
  real gamma0;
  vector[n_species] gamma_species_raw;
  real<lower=0> sigma_gamma_species;
  vector[n_cities] gamma_city_raw;
  real<lower=0> sigma_gamma_city;
  real mu_gamma_wingspan;
  vector[n_cities] gamma_wingspan_raw;  
  real<lower=0> sigma_gamma_wingspan;
  real mu_gamma_park_size;
  vector[n_cities] gamma_park_size_raw;  
  real<lower=0> sigma_gamma_park_size;
  real mu_gamma_isolation;
  vector[n_cities] gamma_isolation_raw;  
  real<lower=0> sigma_gamma_isolation;

  // persistence
  real phi0;
  vector[n_species] phi_species_raw;
  real<lower=0> sigma_phi_species;
  vector[n_cities] phi_city_raw;
  real<lower=0> sigma_phi_city;
  real mu_phi_wingspan;
  vector[n_cities] phi_wingspan_raw;  
  real<lower=0> sigma_phi_wingspan;
  real mu_phi_park_size;
  vector[n_cities] phi_park_size_raw;  
  real<lower=0> sigma_phi_park_size;
  real mu_phi_isolation;
  vector[n_cities] phi_isolation_raw;  
  real<lower=0> sigma_phi_isolation;

  // detection
  real p0; // intercept
  vector[n_species_clusters] p_species_raw;
  real<lower=0> sigma_p_species;
  vector[n_cities] p_city_raw;
  real<lower=0> sigma_p_city;
  real p_wingspan;
  real p_feature_diversity;
  real p_ease_of_id;
  vector[n_species] p_date; // phenology peak
  real mu_p_species_date; // community mean
  real<lower=0> sigma_p_species_date; // variation
  real p_date_latitude; // an effect of latitude on the effect of date
  vector[n_species] p_date_sq; // decay pattern of phenology
  real mu_p_species_date_sq; // variation
  real<lower=0> sigma_p_species_date_sq; // community mean
  real p_date_sq_latitude; // an effect of latitude on the effect of date^2
  
} // end parameters

transformed parameters {

  // logit scale psi1, gamma, phi
  array[R] real psi1; // odds of occurrence year 1
  array[R] real gamma; // odds of colonization
  array[R] real phi; // odds of persistence
  array[R, n_surveys] real p; // odds of detection
  
  vector[n_species_clusters] psi1_species;
  vector[n_cities] psi1_city;
  vector[n_cities] psi1_wingspan;
  vector[n_cities] psi1_park_size;
  vector[n_cities] psi1_isolation;
  vector[n_species] gamma_species;
  vector[n_cities] gamma_city;
  vector[n_cities] gamma_wingspan;
  vector[n_cities] gamma_park_size;
  vector[n_cities] gamma_isolation;
  vector[n_species] phi_species;
  vector[n_cities] phi_city;
  vector[n_cities] phi_wingspan;
  vector[n_cities] phi_park_size;
  vector[n_cities] phi_isolation;
  vector[n_species_clusters] p_species;
  vector[n_cities] p_city;
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_species = sigma_psi1_species * psi1_species_raw;
  psi1_city = psi1_0 + sigma_psi1_city * psi1_city_raw;
  psi1_wingspan = mu_psi1_wingspan + sigma_psi1_wingspan * psi1_wingspan_raw;
  psi1_park_size = mu_psi1_park_size + sigma_psi1_park_size * psi1_park_size_raw;
  psi1_isolation = mu_psi1_isolation + sigma_psi1_isolation * psi1_isolation_raw;
  gamma_species = sigma_gamma_species * gamma_species_raw;
  gamma_city = gamma0 + sigma_gamma_city * gamma_city_raw;
  gamma_wingspan = mu_gamma_wingspan + sigma_gamma_wingspan * gamma_wingspan_raw;
  gamma_park_size = mu_gamma_park_size + sigma_gamma_park_size * gamma_park_size_raw;
  gamma_isolation = mu_gamma_isolation + sigma_gamma_isolation * gamma_isolation_raw;
  phi_species = sigma_phi_species * phi_species_raw;
  phi_city = phi0 + sigma_phi_city * phi_city_raw;
  phi_wingspan = mu_phi_wingspan + sigma_phi_wingspan * phi_wingspan_raw;
  phi_park_size = mu_phi_park_size + sigma_phi_park_size * phi_park_size_raw;
  phi_isolation = mu_phi_isolation + sigma_phi_isolation * phi_isolation_raw;
  p_species = sigma_p_species * p_species_raw;
  p_city = p0 + sigma_p_city * p_city_raw;
  
  for(r in 1:R){
  
        psi1[r] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          psi1_city[city_id_vector[r]] +
          psi1_species[species_cluster_id_vector[r]] + // a species specific intercept
          psi1_wingspan[city_id_vector[r]] * wingspan[species_integer_vector[r]] + // a species effect of migratory
          psi1_park_size[city_id_vector[r]] * park_size[multicity_site_id_vector[r]] + // a site effect of park size
          psi1_isolation[city_id_vector[r]] * isolation[multicity_site_id_vector[r]] // a site effect of park isolation
          ); // end psi1[r]
        
        gamma[r] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma_city[city_id_vector[r]] +
          gamma_species[species_integer_vector[r]] + // a species specific intercept
          gamma_wingspan[city_id_vector[r]] * wingspan[species_integer_vector[r]] + // a species effect of migratory
          gamma_park_size[city_id_vector[r]] * park_size[multicity_site_id_vector[r]] + // a site effect of park size
          gamma_isolation[city_id_vector[r]] * isolation[multicity_site_id_vector[r]] // a site effect of park isolation
          ); // end gamma[i,j,k]
                
        phi[r] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_city[city_id_vector[r]] +
          phi_species[species_integer_vector[r]] + // a species specific intercept
          phi_wingspan[city_id_vector[r]] * wingspan[species_integer_vector[r]] + // a species effect of migratory
          phi_park_size[city_id_vector[r]] * park_size[multicity_site_id_vector[r]] + // a site effect of park size
          phi_isolation[city_id_vector[r]] * isolation[multicity_site_id_vector[r]] // a site effect of park isolation
          ); // end phi[i,j,k]
          
  } // end loop across all data
  
  // have to do another loop because phi and gamma have a shorter k index length than p
  for(r in 1:R){
        for(l in 1:n_surveys){ // loop across all surveys

          p[r,l] = inv_logit( // probability (0-1) of detection is equal to..
            p_city[city_id_vector[r]] + 
            p_species[species_cluster_id_vector[r]] + // a species specific intercept
            p_wingspan * wingspan[species_integer_vector[r]] + // a species effect of wingspan
            p_feature_diversity * feature_diversity[species_integer_vector[r]] + // a species effect of feature diversity
            p_ease_of_id * ease_of_id[species_integer_vector[r]] + // a species effect of ease of identification
            p_date[species_integer_vector[r]] * surveys[l] + // a species-specific phenological detection effect (peak)
            p_date_latitude * latitude[multicity_site_id_vector[r]] * surveys[l] + // an interactive effect of latitude on phenological detection effect (peak)
            p_date_sq[species_integer_vector[r]] * (surveys[l])^2 + // a species-specific phenological detection effect (decay)
            p_date_sq_latitude * latitude[multicity_site_id_vector[r]] * (surveys[l])^2 // an interactive effect of latitude on phenological detection effect (decay)
            ); // end p[j,k,l]
            
        } // end loop across all surveys
      } // end loop across all R
   
  // construct an occurrence array
  array[R] real psi;
  
  for(r in 1:R){
    
        if(site_survey_year_vector[r] < 2){ // define initial state
          psi[r] = psi1[r]; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          // reduce 1 from k for phi and gamma because there are n_years - 1 transitions 
          // and so there are only n_years - 1 speciesXsite "stacks" of phi and gamma
          // but phi[,,k-1] for k = 2 will actually consider the effects of 
          // e.g. flower abundacnce in year 2 (since year 2 is the first year we estimate phi)
          psi[r] = psi[prev_index_vector[r]] * phi[prev_index_vector[r]] + 
                    (1 - psi[prev_index_vector[r]]) * gamma[prev_index_vector[r]]; 
        } // end if/else
        
  } // end loop across all speciesXsiteXyear combinations r:R
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 1); // initial occurrence intercept
  psi1_city_raw ~ std_normal();
  sigma_psi1_city ~ normal(0, 0.25);
  psi1_species_raw ~ std_normal();
  sigma_psi1_species ~ normal(0, 1);
  mu_psi1_wingspan ~ normal(0, 2);
  psi1_wingspan_raw ~ std_normal();
  sigma_psi1_wingspan ~ normal(0, 0.5);
  mu_psi1_park_size ~ normal(0, 2);
  psi1_park_size_raw ~ std_normal();
  sigma_psi1_park_size ~ normal(0, 0.5);
  mu_psi1_isolation ~ normal(0, 2);
  psi1_isolation_raw ~ std_normal();
  sigma_psi1_isolation ~ normal(0, 0.5);

  // colonization
  gamma0 ~ normal(0, 1); // colonization intercept
  gamma_city_raw ~ std_normal();
  sigma_gamma_city ~ normal(0, 0.25);
  gamma_species_raw ~ std_normal();
  sigma_gamma_species ~ normal(0, 0.5);
  mu_gamma_wingspan ~ normal(0, 2);
  gamma_wingspan_raw ~ std_normal();
  sigma_gamma_wingspan ~ normal(0, 0.5);
  mu_gamma_park_size ~ normal(0, 2);
  gamma_park_size_raw ~ std_normal();
  sigma_gamma_park_size ~ normal(0, 0.5);
  mu_gamma_isolation ~ normal(0, 2);
  gamma_isolation_raw ~ std_normal();
  sigma_gamma_isolation ~ normal(0, 0.5);

  // persistence
  phi0 ~ normal(0, 1); // global intercept
  phi_city_raw ~ std_normal();
  sigma_phi_city ~ normal(0, 0.25);
  phi_species_raw ~ std_normal();
  sigma_phi_species ~ normal(0, 0.5);
  mu_phi_wingspan ~ normal(0, 2);
  phi_wingspan_raw ~ std_normal();
  sigma_phi_wingspan ~ normal(0, 0.5);
  mu_phi_park_size ~ normal(0, 2);
  phi_park_size_raw ~ std_normal();
  sigma_phi_park_size ~ normal(0, 0.5);
  mu_phi_isolation ~ normal(0, 2);
  phi_isolation_raw ~ std_normal();
  sigma_phi_isolation ~ normal(0, 0.5);

  // detection
  p0 ~ normal(0, 1); // global intercept
  p_city_raw ~ std_normal();
  sigma_p_city ~ normal(0, 0.25);
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 1);
  p_wingspan ~ normal(0, 2);
  p_feature_diversity ~ normal(0, 2);
  p_ease_of_id ~ normal(0, 2);
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0, 2); // mean
  sigma_p_species_date ~ normal(0, 2); // variation
  p_date_latitude ~ normal(0, 2); // effect of latitude on peak date
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0, 1); // mean
  sigma_p_species_date_sq ~ normal(0, 1); // variation
  p_date_sq_latitude ~ normal(0, 2); // effect of latitude on date decay

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
