// dynamic multi-species occupancy model for pollinators in urban parks
// this is for the full iNat data
// I'm using months as survey blocks and assuming all parks are surveyed in all months in all years
// so long as one or more species from corresponding family was detected during that month
// I've included a month-phenology covariate for detection

// I've also added park site and species trait predictors 
// contrast with model 1.1, I additionally added species traits as predictors of species random effects
// including an effect of migratory status on detectability and psi1

// jcu, started may, 2024.

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
  array[R] int<lower=1, upper=n_sites> multicity_site_integer_vector; // vector indicating which site is being observed
  int<lower=0> n_cities; // number of cities
  array[R] int<lower=1, upper=n_cities> city_id_vector; // vector indicating which city is being observed 
  int<lower=0> n_species_clusters; // number of speciesXregions clusters
  array[R] int<lower=1, upper=n_species_clusters> species_cluster_id_vector;
  int<lower=0> n_species_city_clusters; // number of speciesXregions clusters
  array[R] int<lower=1, upper=n_species_city_clusters> species_city_id_vector;
  int<lower=0> n_regional_clusters; // number of speciesXregions clusters
  array[R] int<lower=1, upper=n_regional_clusters> regional_cluster_id_vector;
  // species and site covariate data
  vector[n_species] feature_diversity;
  vector[n_species] ease_of_id;
  vector[n_species] wingspan;
  vector[n_species] migratory;
  vector[n_sites] park_size;
  vector[n_sites] isolation;
  vector[n_sites] total_detections_by_city;
  // other stuff
  array[R] int<lower=0> confirmed_occurrence;
  array[R] int<lower=0> prev_index_vector;
} // end data

parameters {

  // initial state
  real psi1_0; // intercept
  vector[n_cities] psi1_city_raw;
  real<lower=0> sigma_psi1_city;
  real mu_psi1_park_size;
  vector[n_cities] psi1_park_size_raw;  
  real<lower=0> sigma_psi1_park_size;
  real mu_psi1_isolation;
  vector[n_cities] psi1_isolation_raw;  
  real<lower=0> sigma_psi1_isolation;
  vector[n_species_city_clusters] psi1_species;
  real<lower=0> sigma_psi1_species;
  real psi1_wingspan;
  real psi1_migratory;

  // colonization
  real gamma0;
  vector[n_cities] gamma_city_raw;
  real<lower=0> sigma_gamma_city;
  real mu_gamma_park_size;
  vector[n_cities] gamma_park_size_raw;  
  real<lower=0> sigma_gamma_park_size;
  real mu_gamma_isolation;
  vector[n_cities] gamma_isolation_raw;  
  real<lower=0> sigma_gamma_isolation;
  vector[n_species_clusters] gamma_species;
  real<lower=0> sigma_gamma_species;
  real gamma_wingspan;
  real gamma_migratory;

  // persistence
  real phi0;
  vector[n_cities] phi_city_raw;
  real<lower=0> sigma_phi_city;
  real mu_phi_park_size;
  vector[n_cities] phi_park_size_raw;  
  real<lower=0> sigma_phi_park_size;
  real mu_phi_isolation;
  vector[n_cities] phi_isolation_raw;  
  real<lower=0> sigma_phi_isolation;
  vector[n_species_clusters] phi_species;
  real<lower=0> sigma_phi_species;
  real phi_wingspan;
  real phi_migratory;
  
  // detection
  real p0; // intercept
  vector[n_cities] p_city_raw;
  real<lower=0> sigma_p_city;
  real p_city_detections;
  vector[n_species_clusters] p_species;
  real p_wingspan;
  real p_migratory;
  real p_feature_diversity;
  real p_ease_of_id;
  real<lower=0> sigma_p_species;
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

  // logit scale psi1, gamma, phi
  array[R] real psi1; // odds of occurrence year 1
  array[R] real gamma; // odds of colonization
  array[R] real phi; // odds of persistence
  array[R, n_surveys] real p; // odds of detection
  
  vector[n_cities] psi1_city;
  vector[n_cities] psi1_park_size;
  vector[n_cities] psi1_isolation;
  vector[n_cities] gamma_city;
  vector[n_cities] gamma_park_size;
  vector[n_cities] gamma_isolation;
  vector[n_cities] phi_city;
  vector[n_cities] phi_park_size;
  vector[n_cities] phi_isolation;
  vector[n_cities] p_city;
  vector[n_species_city_clusters] mu_psi1_species; // phenology peak
  vector[n_species_clusters] mu_gamma_species; // phenology peak
  vector[n_species_clusters] mu_phi_species; // phenology peak
  vector[n_species_clusters] mu_p_species; // phenology peak
  vector[n_species_clusters] mu_p_species_date; // phenology peak
  vector[n_species_clusters] mu_p_species_date_sq; // phenology curve
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_city = psi1_0 + sigma_psi1_city * psi1_city_raw;
  psi1_park_size = mu_psi1_park_size + sigma_psi1_park_size * psi1_park_size_raw;
  psi1_isolation = mu_psi1_isolation + sigma_psi1_isolation * psi1_isolation_raw;
  gamma_city = gamma0 + sigma_gamma_city * gamma_city_raw;
  gamma_park_size = mu_gamma_park_size + sigma_gamma_park_size * gamma_park_size_raw;
  gamma_isolation = mu_gamma_isolation + sigma_gamma_isolation * gamma_isolation_raw;
  phi_city = phi0 + sigma_phi_city * phi_city_raw;
  phi_park_size = mu_phi_park_size + sigma_phi_park_size * phi_park_size_raw;
  phi_isolation = mu_phi_isolation + sigma_phi_isolation * phi_isolation_raw;
  p_city = p0 + sigma_p_city * p_city_raw;
  
  for(r in 1:R){
  
    // phenology predictors
    mu_p_species_date[species_cluster_id_vector[r]] = delta0 + delta_regional_cluster[regional_cluster_id_vector[r]];
    mu_p_species_date_sq[species_cluster_id_vector[r]] = epsilon0 + epsilon_regional_cluster[regional_cluster_id_vector[r]];
    
    // species traits predict expected value of species random effect
    // species can get a different random effect depending on the geographical cluster 
    // where the detection takes place (species_cluster_id_vector[r]), but the species traits
    // for that species are always the same (species_integer_vector[r]) e.g. we don't
    // have different measures of a species wingspan in different regions just one value of mean wingspan
    mu_psi1_species[species_city_id_vector[r]] = psi1_wingspan*wingspan[species_integer_vector[r]] +
        psi1_migratory*migratory[species_integer_vector[r]];
    mu_gamma_species[species_cluster_id_vector[r]] = gamma_wingspan*wingspan[species_integer_vector[r]] +
        gamma_migratory*migratory[species_integer_vector[r]];
    mu_phi_species[species_cluster_id_vector[r]] = phi_wingspan*wingspan[species_integer_vector[r]] +
        phi_migratory*migratory[species_integer_vector[r]];
    mu_p_species[species_cluster_id_vector[r]] = p_wingspan*wingspan[species_integer_vector[r]] +
        p_migratory*migratory[species_integer_vector[r]] + 
        p_feature_diversity * feature_diversity[species_integer_vector[r]] +
        p_ease_of_id * ease_of_id[species_integer_vector[r]];
        
    // ecological processes
    psi1[r] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
      psi1_city[city_id_vector[r]] +
      psi1_species[species_city_id_vector[r]] + // a species specific intercept
      psi1_park_size[city_id_vector[r]] * park_size[multicity_site_integer_vector[r]] + // a site effect of park size
      psi1_isolation[city_id_vector[r]] * isolation[multicity_site_integer_vector[r]] // a site effect of park isolation
      ); // end psi1[r]
    
    gamma[r] = inv_logit( // probability (0-1) of colonization is equal to..
      gamma_city[city_id_vector[r]] +
      gamma_species[species_cluster_id_vector[r]] + // a species specific intercept
      gamma_park_size[city_id_vector[r]] * park_size[multicity_site_integer_vector[r]] + // a site effect of park size
      gamma_isolation[city_id_vector[r]] * isolation[multicity_site_integer_vector[r]] // a site effect of park isolation
      ); // end gamma[i,j,k]
            
    phi[r] = inv_logit( // probability (0-1) of persistence is equal to..
      phi_city[city_id_vector[r]] +
      phi_species[species_cluster_id_vector[r]] + // a species specific intercept
      phi_park_size[city_id_vector[r]] * park_size[multicity_site_integer_vector[r]] + // a site effect of park size
      phi_isolation[city_id_vector[r]] * isolation[multicity_site_integer_vector[r]] // a site effect of park isolation
      ); // end phi[i,j,k]
          
  } // end loop across all data
  
  // have to do another loop because phi and gamma have a shorter k index length than p
  for(r in 1:R){
        for(l in 1:n_surveys){ // loop across all surveys

          p[r,l] = inv_logit( // probability (0-1) of detection is equal to..
            p_city[city_id_vector[r]] + 
            p_species[species_cluster_id_vector[r]] + // a species specific intercept
            p_city_detections * total_detections_by_city[multicity_site_integer_vector[r]] +
            p_date[species_cluster_id_vector[r]] * surveys[l] + // a species-specific phenological detection effect (peak)
            p_date_sq[species_cluster_id_vector[r]] * (surveys[l])^2 // a species-specific phenological detection effect (decay)
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
  sigma_psi1_city ~ normal(0, 0.5);
  mu_psi1_park_size ~ normal(0, 2);
  psi1_park_size_raw ~ std_normal();
  sigma_psi1_park_size ~ normal(0, 0.5);
  mu_psi1_isolation ~ normal(0, 2);
  psi1_isolation_raw ~ std_normal();
  sigma_psi1_isolation ~ normal(0, 0.5);
  psi1_species ~ normal(mu_psi1_species, sigma_psi1_species); // species-specific random effect
  psi1_wingspan ~ normal(0, 0.5); 
  psi1_migratory ~ normal(0, 0.5);
  sigma_psi1_species ~ normal(0, 1);

  // colonization
  gamma0 ~ normal(0, 1); // colonization intercept
  gamma_city_raw ~ std_normal();
  sigma_gamma_city ~ normal(0, 0.5);
  mu_gamma_park_size ~ normal(0, 2);
  gamma_park_size_raw ~ std_normal();
  sigma_gamma_park_size ~ normal(0, 0.5);
  mu_gamma_isolation ~ normal(0, 2);
  gamma_isolation_raw ~ std_normal();
  sigma_gamma_isolation ~ normal(0, 0.5);
  gamma_species ~ normal(mu_gamma_species, sigma_gamma_species); // species-specific random effect
  gamma_wingspan ~ normal(0, 0.5); 
  gamma_migratory ~ normal(0, 0.5);
  sigma_gamma_species ~ normal(0, 1);

  // persistence
  phi0 ~ normal(0, 1); // global intercept
  phi_city_raw ~ std_normal();
  sigma_phi_city ~ normal(0, 0.5);
  mu_phi_park_size ~ normal(0, 2);
  phi_park_size_raw ~ std_normal();
  sigma_phi_park_size ~ normal(0, 0.5);
  mu_phi_isolation ~ normal(0, 2);
  phi_isolation_raw ~ std_normal();
  sigma_phi_isolation ~ normal(0, 0.5);
  phi_species ~ normal(mu_phi_species, sigma_phi_species); // species-specific random effect
  phi_wingspan ~ normal(0, 0.5); 
  phi_migratory ~ normal(0, 0.5);
  sigma_phi_species ~ normal(0, 1);

  // detection
  p0 ~ normal(0, 1); // global intercept
  p_city_raw ~ std_normal();
  sigma_p_city ~ normal(0, 0.5);
  p_city_detections ~ normal(0, 1);
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  delta0 ~ normal(0, 2); // community mean
  delta_regional_cluster ~ normal(0, 1); // effect of region
  sigma_p_species_date ~ normal(0, 2); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  epsilon0 ~ normal(0, 1); // community mean
  epsilon_regional_cluster ~ normal(0, 0.5); // effect of region
  sigma_p_species_date_sq ~ normal(0, 1); // variation
  p_species ~ normal(mu_p_species, sigma_p_species); // species-specific random effect
  p_wingspan ~ normal(0, 1); 
  p_migratory ~ normal(0, 1);
  p_feature_diversity ~ normal(0, 1);
  p_ease_of_id ~ normal(0, 1);
  sigma_p_species ~ normal(0, 1);

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

