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

functions {
  
  // covariance matrix for detection rates and baseline occurrence rates
  matrix custom_cov_matrix(vector sigma, real rho) {
    matrix[2,2] Sigma;
    Sigma[1,1] = square(sigma[1]); // species variation in occurrence1 rates
    Sigma[2,2] = square(sigma[2]); // species variation in detection rates
    Sigma[1,2] = sigma[1] * sigma[2] * rho; // correlation between species-specific detection and occurrence1 rates
    Sigma[2,1] = Sigma[1,2]; // correlation between species-specific detection and occurrence1 rates
    return Sigma;
  }
  
  // mean values for multivariate normal distribution (center species on global intercepts)
  vector mu(real psi1_0, real p0){
    vector[2] global_intercepts;
    global_intercepts[1] = psi1_0;
    global_intercepts[2] = p0;
    return global_intercepts;
  }
  
}

data {
  
  int<lower=0> n_species; // number of species
  int<lower=1> species[n_species]; // vector of species identities
  int<lower=0> n_sites; // number of sites
  int<lower=1> sites[n_sites]; // vector of site identities
  int<lower=0> n_years; // total years
  int<lower=0> n_years_minus1;
  int<lower=1> years[n_years_minus1]; // vector of year transition indexes (index 1 == transition between years 1 and 2)
  int<lower=0> n_surveys; // surveys per year
  real surveys[n_surveys]; // surveys (difference from the mean)
  int<lower=0,upper=1> V[n_species, n_sites, n_years, n_surveys]; // binary detection / non detection
  int<lower=0,upper=1> V_NA[n_species, n_sites, n_years, n_surveys]; // sampling indicator 1==non-detection, 0==no evidence the species was sampled
  vector[n_species] feature_diversity;
  vector[n_species] ease_of_id;
  vector[n_species] wingspan;
  vector[n_sites] park_size;
  vector[n_sites] connectivity;
  vector[n_sites] plant_genera;
  vector[n_sites] tree_cover;

} // end data

parameters {
  
    // Covararying Parameters
  real<lower=-1,upper=1> rho;  // correlation of (species occurrence and detection intercepts)
  vector<lower=0>[2] sigma_species; // variance in species-specific detection rates 
    // (sigma_species[1] == occurrence and sigma_species[2] == detection)
  vector[2] species_intercepts[n_species];// species-level occurrence and detection intercepts
    // (species_intercepts_raw[1] == psi1 intercepts and species_intercepts_raw[2] == detection intercepts)

  // initial state
  real psi1_0; // intercept
  real psi1_wingspan;
  real psi1_park_size;
  real psi1_connectivity;
  real psi1_tree_cover;

  // colonization
  real gamma0;
  vector[n_species] gamma_species;
  real<lower=0> sigma_gamma_species;
  real gamma_wingspan;
  real gamma_park_size;
  real gamma_connectivity;
  real gamma_tree_cover;

  // persistence
  real phi0;
  vector[n_species] phi_species;
  real<lower=0> sigma_phi_species;
  real phi_wingspan;
  real phi_park_size;
  real phi_connectivity;
  real phi_tree_cover;

  // detection
  real p0; // intercept
  real p_wingspan;
  real p_feature_diversity;
  real p_ease_of_id;
  vector[n_species] p_date; // phenology peak
  real mu_p_species_date; // community mean
  real<lower=0> sigma_p_species_date; // variation
  vector[n_species] p_date_sq; // decay pattern of phenology
  real mu_p_species_date_sq; // variation
  real<lower=0> sigma_p_species_date_sq; // community mean
  
} // end parameters

transformed parameters {

  // logit scale psi1, gamma, phi
  real psi1[n_species, n_sites]; // odds of occurrence year 1
  real gamma[n_species, n_sites, n_years_minus1]; // odds of colonization
  real phi[n_species, n_sites, n_years_minus1]; // odds of persistence
  real p[n_species, n_sites, n_years, n_surveys];  // odds of detection

  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years_minus1){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          species_intercepts[species[i],1] + // a species specific intercept
          psi1_wingspan * wingspan[i] + // a species effect of migratory
          psi1_park_size * park_size[j] + // a site effect of park size
          psi1_connectivity * connectivity[j] + // a site effect of park connectivity
          psi1_tree_cover * tree_cover[j] // a site effect of tree cover
          ); // end psi1[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma_species[species[i]] + // a species specific intercept
          gamma_wingspan * wingspan[i] + // a species effect of migratory
          gamma_park_size * park_size[j] + // a site effect of park size
          gamma_connectivity * connectivity[j] + // a site effect of park connectivity
          gamma_tree_cover * tree_cover[j] // a site effect of tree cover
          ); // end gamma[i,j,k]
                
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_species[species[i]] + // a species specific intercept
          phi_wingspan * wingspan[i] + // a species effect of migratory
          phi_park_size * park_size[j] + // a site effect of park size
          phi_connectivity * connectivity[j] + // a site effect of park connectivity
          phi_tree_cover * tree_cover[j] // a site effect of tree cover
          ); // end phi[i,j,k]
          
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
  
  // have to do another loop because phi and gamma have a shorter k index length than p
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_surveys){ // loop across all surveys

          p[i,j,k,l] = inv_logit( // probability (0-1) of detection is equal to..
            species_intercepts[species[i],2] + // a species specific intercept
            p_wingspan * wingspan[i] + // a species effect of wingspan
            p_feature_diversity * feature_diversity[i] + // a species effect of feature diversity
            p_ease_of_id * ease_of_id[i] + // a species effect of ease of identification
            p_date[species[i]] * surveys[l] + // a species-specific phenological detection effect (peak)
            p_date_sq[species[i]] * (surveys[l])^2 // a species-specific phenological detection effect (decay)
            ); // end p[j,k,l]
            
        } // end loop across all surveys
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species 
   
  // construct an occurrence array
  real psi[n_species, n_sites, n_years];
  
  for(i in 1:n_species){
    for(j in 1:n_sites){
      for(k in 1:n_years){
        
        if(k < 2){ // define initial state
          psi[i,j,k] = psi1[i,j]; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          // reduce 1 from k for phi and gamma because there are n_years - 1 transitions 
          // and so there are only n_years - 1 speciesXsite "stacks" of phi and gamma
          // but phi[,,k-1] for k = 2 will actually consider the effects of 
          // e.g. flower abundacnce in year 2 (since year 2 is the first year we estimate phi)
          psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k-1] + (1 - psi[i,j,k-1]) * gamma[i,j,k-1]; 
        } // end if/else
        
      } // end loop across all years
    } // end loop across all sites
  } // end loop across all species
   
} // end transformed parameters

model {
  
  // PRIORS
  
     // correlated params:
  // correlated species effects
  sigma_species[1] ~ normal(0, 2);
  sigma_species[2] ~ normal(0, 2);
  (rho + 1) / 2 ~ beta(2, 2);
  
  // correlated species-specific detection rates
  // will send the mean (mu), variance and correlation to the covariance matrix
  species_intercepts ~ multi_normal(mu(psi1_0, p0), 
    custom_cov_matrix(sigma_species, rho));

  // occupancy
  // initial state
  psi1_0 ~ normal(0, 2); // persistence intercept
  psi1_wingspan ~ normal(0, 2);
  psi1_park_size ~ normal(0, 2);
  psi1_connectivity ~ normal(0, 2);
  psi1_tree_cover ~ normal(0, 2);

  // colonization
  gamma0 ~ normal(0, 1); // persistence intercept
  gamma_species ~ normal(gamma0, sigma_gamma_species);
  sigma_gamma_species ~ normal(0, 1);
  gamma_wingspan ~ normal(0, 2);
  gamma_park_size ~ normal(0, 2);
  gamma_connectivity ~ normal(0, 2);
  gamma_tree_cover ~ normal(0, 2);

  // persistence
  phi0 ~ normal(0, 1); // global intercept
  phi_species ~ normal(phi0, sigma_phi_species);
  sigma_phi_species ~ normal(0, 1);
  phi_wingspan ~ normal(0, 2);
  phi_park_size ~ normal(0, 2);
  phi_connectivity ~ normal(0, 2);
  phi_tree_cover ~ normal(0, 2);

  // detection
  p0 ~ normal(0, 2); // global intercept
  p_wingspan ~ normal(0, 2);
  p_feature_diversity ~ normal(0, 2);
  p_ease_of_id ~ normal(0, 2);
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0, 2); // mean
  sigma_p_species_date ~ normal(0, 2); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0, 1); // mean
  sigma_p_species_date_sq ~ normal(0, 1); // variation

  // LIKELIHOOD
  for(i in 1:n_species){
    for (j in 1:n_sites){
      for (k in 1:n_years){
        
        if(sum(V_NA[i,j,k]) > 0){
          
          if (sum(V[i,j,k]) > 0){ // lp observed 
            // detection on each visit given detection rate on each visit
            // V_NA == 0 indicates that no survey occurred. Multiplying 0 by the 
            // lpmf() statement should remove it from the target
            target += (log(psi[i,j,k]) + // present
                           bernoulli_lpmf(V[i,j,k,1]|p[i,j,k,1])*V_NA[i,j,k,1] + 
                           bernoulli_lpmf(V[i,j,k,2]|p[i,j,k,2])*V_NA[i,j,k,2] + 
                           bernoulli_lpmf(V[i,j,k,3]|p[i,j,k,3])*V_NA[i,j,k,3] + 
                           bernoulli_lpmf(V[i,j,k,4]|p[i,j,k,4])*V_NA[i,j,k,4] +
                           bernoulli_lpmf(V[i,j,k,5]|p[i,j,k,5])*V_NA[i,j,k,5] + 
                           bernoulli_lpmf(V[i,j,k,6]|p[i,j,k,6])*V_NA[i,j,k,6] + 
                           bernoulli_lpmf(V[i,j,k,7]|p[i,j,k,7])*V_NA[i,j,k,7] + 
                           bernoulli_lpmf(V[i,j,k,8]|p[i,j,k,8])*V_NA[i,j,k,8] + 
                           bernoulli_lpmf(V[i,j,k,9]|p[i,j,k,9])*V_NA[i,j,k,9] + 
                           bernoulli_lpmf(V[i,j,k,10]|p[i,j,k,10])*V_NA[i,j,k,10] + 
                           bernoulli_lpmf(V[i,j,k,11]|p[i,j,k,11])*V_NA[i,j,k,11] + 
                           bernoulli_lpmf(V[i,j,k,12]|p[i,j,k,12])*V_NA[i,j,k,12]
                           
                       );
          } else { // lp unobserved (set up for 6 annual surveys)
            // marginal likelihood of:
            // occurrence, yet...
            // non-detection on each visit given detection rate on each visit given occurrence
            target += (log_sum_exp(log(psi[i,j,k]) + log1m(p[i,j,k,1])*V_NA[i,j,k,1] + 
                                                     log1m(p[i,j,k,2])*V_NA[i,j,k,2] + 
                                                     log1m(p[i,j,k,3])*V_NA[i,j,k,3] + 
                                                     log1m(p[i,j,k,4])*V_NA[i,j,k,4] +
                                                     log1m(p[i,j,k,5])*V_NA[i,j,k,5] + 
                                                     log1m(p[i,j,k,6])*V_NA[i,j,k,6] + 
                                                     log1m(p[i,j,k,7])*V_NA[i,j,k,7] + 
                                                     log1m(p[i,j,k,8])*V_NA[i,j,k,8] + 
                                                     log1m(p[i,j,k,9])*V_NA[i,j,k,9] + 
                                                     log1m(p[i,j,k,10])*V_NA[i,j,k,10] + 
                                                     log1m(p[i,j,k,11])*V_NA[i,j,k,11] + 
                                                     log1m(p[i,j,k,12])*V_NA[i,j,k,12],
                                    // or just simple no occurrence
                                    log1m(psi[i,j,k])));
          } // end if/else
          
        } // end if (any in year not NA)
  
      } // end loop across all years
    } // end loop across all sites   
  } // end loop across all species

} // end model

generated quantities{
  
  //
  // posterior predictive check (number of detections, binned by species)
  //
  int<lower=0> W_species_rep[n_species]; // sum of simulated detections
  
  int z_simmed[n_species, n_sites, n_years]; // simulate occurrence

  for(i in 1:n_species){
   for(j in 1:n_sites){
     for(k in 1:n_years){
          z_simmed[i,j,k] = bernoulli_rng(psi[i,j,k]); 
      }    
    }
  }
  
  // initialize repped detections at 0
  for(i in 1:n_species){
    W_species_rep[i] = 0;
  }
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:n_species) { // loop across all species
    for(j in 1:n_sites) { // loop across all sites
      for(k in 1:n_years){ // loop across all years
        for(l in 1:n_surveys){
          
          // detections in replicated data (us z_simmed from above)
          W_species_rep[i] = W_species_rep[i] + 
            // multiply by the NA indicator - if we didn't survey in real life
            // we don't survey in this simulation.
            (z_simmed[i,j,k] * bernoulli_rng(p[i,j,k,l]) * V_NA[i,j,k,l]);
           
        } // end loop across surveys
      } // end loop across years
    } // end loop across sites
  } // end loop across species
  
} // end generated quantities
