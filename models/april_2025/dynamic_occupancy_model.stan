// dynamic multi-species occupancy model for pollinators in urban parks
// this is for the full iNat data
// I'm using months as survey blocks and assuming all parks are surveyed in all months in all years
// so long as one or more species from corresponding family was detected
// I've included a month-phenology covariate for detection

// I've also added park site and species trait predictors (contrast with model 3)

// jens ulrich, started april, 2025.

data {
  int<lower=1> R; // the number of rows in the data
  int<lower=0> n_surveys; // surveys per year
  int<lower=0,upper=1> V[R, n_surveys]; // binary detection / non detection
  int<lower=0,upper=1> V_NA[R, n_surveys]; // sampling indicator 1==non-detection, 0==no evidence the species was sampled
  int<lower=0> reverse_index[R]; // index of row number for speciesXsite at previous timestep
  int<lower=1> year[R]; // vector of year
  int<lower=1> species[R]; // vector of species identities
  int<lower=1> n_species; // total number of species
  int<lower=1> site[R]; // vector of site identities
  real survey[n_surveys]; // surveys (difference from the mean)
  vector[R] feature_diversity;
  vector[R] ease_of_id;
  vector[R] wingspan;
  vector[R] park_size;
} // end data

parameters {

  // initial state
  real psi1_0; // intercept
  vector[n_species] psi1_species_raw;  
  real<lower=0> sigma_psi1_species;
  real psi1_wingspan;
  real psi1_park_size;

  // colonization
  real gamma0;
  vector[n_species] gamma_species_raw;
  real<lower=0> sigma_gamma_species;
  real gamma_wingspan;
  real gamma_park_size;
  
  // persistence
  real phi0;
  vector[n_species] phi_species_raw;
  real<lower=0> sigma_phi_species;
  real phi_wingspan;
  real phi_park_size;

  // detection
  real p0; // intercept
  vector[n_species] p_species_raw;
  real<lower=0> sigma_p_species;
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
  real psi1[R]; // odds of occurrence year 1
  real gamma[R]; // odds of colonization
  real phi[R]; // odds of persistence
  real p[R, n_surveys];  // odds of detection
  
  vector[n_species] psi1_species;
  vector[n_species] gamma_species;
  vector[n_species] phi_species;
  vector[n_species] p_species;
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_species = psi1_0 + sigma_psi1_species * psi1_species_raw;
  gamma_species = gamma0 + sigma_gamma_species * gamma_species_raw;
  phi_species = phi0 + sigma_phi_species * phi_species_raw;
  p_species = p0 + sigma_p_species * p_species_raw;
  
  for(i in 1:R){
  
        psi1[i] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          psi1_species[species[i]] + // a species specific intercept
          psi1_wingspan * wingspan[i] + // a species effect of migratory
          psi1_park_size * park_size[i] // a site effect of park size
          ); // end psi1[j]
        
        gamma[i] = inv_logit( // probability (0-1) of colonization is equal to..
          gamma_species[species[i]] + // a species specific intercept
          gamma_wingspan * wingspan[i] + // a species effect of migratory
          gamma_park_size * park_size[i] // a site effect of park size
          ); // end gamma[i]
                
        phi[i] = inv_logit( // probability (0-1) of persistence is equal to..
          phi_species[species[i]] + // a species specific intercept
          phi_wingspan * wingspan[i] + // a species effect of migratory
          phi_park_size * park_size[i] // a site effect of park size
          ); // end phi[i]
          
    for(j in 1:n_surveys){ // loop across all surveys

          p[i,j] = inv_logit( // probability (0-1) of detection is equal to..
            p_species[species[i]] + // a species specific intercept
            p_wingspan * wingspan[i] + // a species effect of wingspan
            p_feature_diversity * feature_diversity[i] + // a species effect of feature diversity
            p_ease_of_id * ease_of_id[i] + // a species effect of ease of identification
            p_date[species[i]] * survey[j] + // a species-specific phenological detection effect (peak)
            p_date_sq[species[i]] * (survey[j])^2 // a species-specific phenological detection effect (decay)
            ); // end p[i,j]
            
        } // end loop across all surveys      
  } // end loop across all data units 
   
  // construct an occurrence array
  real psi[R];
  
  for(i in 1:R){
        
        if(year[i] < 2){ // define initial state
          psi[i] = psi1[i]; 
        } else { // describe temporally autocorrelated system dynamics
          // As psi approaches 1, there's a weighted switch on phi (survival)
          // As psi approaches 0, there's a weighted switch on gamma (colonization)
          // reduce 1 from k for phi and gamma because there are n_years - 1 transitions 
          // and so there are only n_years - 1 speciesXsite "stacks" of phi and gamma
          // but phi[,,k-1] for k = 2 will actually consider the effects of 
          // e.g. flower abundacnce in year 2 (since year 2 is the first year we estimate phi)
          
          // find row for previous timestep
          psi[i] = psi[reverse_index[i]] * phi[reverse_index[i]] + 
                      (1 - psi[reverse_index[i]]) * gamma[reverse_index[i]]; 
          // old: psi[i,j,k] = psi[i,j,k-1] * phi[i,j,k-1] + (1 - psi[i,j,k-1]) * gamma[i,j,k-1]; 
        } // end if/else

  } // end loop across all speciesXsiteXyears
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 2); // persistence intercept
  psi1_species_raw ~ std_normal();
  sigma_psi1_species ~ normal(0, 1);
  psi1_wingspan ~ normal(0, 2);
  psi1_park_size ~ normal(0, 2);

  // colonization
  gamma0 ~ normal(0, 2); // persistence intercept
  gamma_species_raw ~ std_normal();
  sigma_gamma_species ~ normal(0, 1);
  gamma_wingspan ~ normal(0, 2);
  gamma_park_size ~ normal(0, 2); 

  // persistence
  phi0 ~ normal(0, 2); // global intercept
  phi_species_raw ~ std_normal();
  sigma_phi_species ~ normal(0, 1);
  phi_wingspan ~ normal(0, 2);
  phi_park_size ~ normal(0, 2);

  // detection
  p0 ~ normal(0, 2); // global intercept
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 2);
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
  for(i in 1:R){
    // (set up for 6 annual sampling events / surveys)
    if (sum(V[i]) > 0){ // lp observed (species observed at site least once in a year)
      // detection on each visit given detection rate on each visit
      // V_NA == 0 indicates that no survey occurred. Multiplying 0 by the 
      // lpmf() statement should remove it from the target
      target += (log(psi[i]) + // present
                     bernoulli_lpmf(V[i,1]|p[i,1])*V_NA[i,1] + 
                     bernoulli_lpmf(V[i,2]|p[i,2])*V_NA[i,2] + 
                     bernoulli_lpmf(V[i,3]|p[i,3])*V_NA[i,3] + 
                     bernoulli_lpmf(V[i,4]|p[i,4])*V_NA[i,4] +
                     bernoulli_lpmf(V[i,5]|p[i,5])*V_NA[i,5] + 
                     bernoulli_lpmf(V[i,6]|p[i,6])*V_NA[i,6] 
                 );
    } else { // lp unobserved 
      // marginal likelihood of:
      // occurrence, yet...
      // non-detection on each visit given detection rate on each visit given occurrence
      target += (log_sum_exp(log(psi[i]) + log1m(p[i,1])*V_NA[i,1] + 
                                           log1m(p[i,2])*V_NA[i,2] + 
                                           log1m(p[i,3])*V_NA[i,3] + 
                                           log1m(p[i,4])*V_NA[i,4] +
                                           log1m(p[i,5])*V_NA[i,5] + 
                                           log1m(p[i,6])*V_NA[i,6],
                              // or just simple no occurrence
                              log1m(psi[i])));
    } // end if/else

  } // end loop across all speciesXsiteXyear sampling units

} // end model

generated quantities{
  
  //
  // posterior predictive check (number of detections, binned by species)
  //
  int<lower=0> W_R_rep[R]; // sum of simulated detections
  
  int z_simmed[R]; // simulate occurrence

  for(i in 1:R){
    z_simmed[i] = bernoulli_rng(psi[i]); 
  }
  
  // initialize repped detections at 0
  for(i in 1:R){
    W_R_rep[i] = 0;
  }
      
  // generating posterior predictive distribution
  // Predict Z at sites
  for(i in 1:R) { // loop across all species
    for(j in 1:n_surveys){
          
          // detections in replicated data (us z_simmed from above)
          W_R_rep[i] = W_R_rep[i] + 
            // multiply by the NA indicator - if we didn't survey in real life
            // we don't survey in this simulation.
            (z_simmed[i] * bernoulli_rng(p[i,j]) * V_NA[i,j]);
           
    } // end loop across surveys
  } // end loop across speciesXsiteXyear
  
} // end generated quantities
