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
  vector[n_sites] isolation;
  int<lower=0> n_cities; // number of cities
  int<lower=1> city[n_sites]; // vector of city identities
  int<lower=0> ranges[n_species, n_sites]; // matrix to constrain analysis within species ranges

} // end data

parameters {

  // initial state
  real psi1_0; // intercept
  vector[n_species] psi1_species_raw;  
  real<lower=0> sigma_psi1_species;
  vector[n_cities] psi1_city;
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
  vector[n_cities] gamma_city;
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
  vector[n_cities] phi_city;
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
  vector[n_species] p_species_raw;
  real<lower=0> sigma_p_species;
  vector[n_cities] p_city;
  real<lower=0> sigma_p_city;
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
  
  vector[n_species] psi1_species;
  //vector[n_cities] psi1_city;
  vector[n_cities] psi1_wingspan;
  vector[n_cities] psi1_park_size;
  vector[n_cities] psi1_isolation;
  vector[n_species] gamma_species;
  //vector[n_cities] gamma_city;
  vector[n_cities] gamma_wingspan;
  vector[n_cities] gamma_park_size;
  vector[n_cities] gamma_isolation;
  vector[n_species] phi_species;
  //vector[n_cities] phi_city;
  vector[n_cities] phi_wingspan;
  vector[n_cities] phi_park_size;
  vector[n_cities] phi_isolation;
  vector[n_species] p_species;
  //vector[n_cities] p_city;
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi1_species = sigma_psi1_species * psi1_species_raw;
  //psi1_city = sigma_psi1_city * psi1_city_raw;
  psi1_wingspan = mu_psi1_wingspan + sigma_psi1_wingspan * psi1_wingspan_raw;
  psi1_park_size = mu_psi1_park_size + sigma_psi1_park_size * psi1_park_size_raw;
  psi1_isolation = mu_psi1_isolation + sigma_psi1_isolation * psi1_isolation_raw;
  gamma_species = sigma_gamma_species * gamma_species_raw;
  //gamma_city = sigma_gamma_city * gamma_city_raw;
  gamma_wingspan = mu_gamma_wingspan + sigma_gamma_wingspan * gamma_wingspan_raw;
  gamma_park_size = mu_gamma_park_size + sigma_gamma_park_size * gamma_park_size_raw;
  gamma_isolation = mu_gamma_isolation + sigma_gamma_isolation * gamma_isolation_raw;
  phi_species = sigma_phi_species * phi_species_raw;
  //phi_city = sigma_phi_city * phi_city_raw;
  phi_wingspan = mu_phi_wingspan + sigma_phi_wingspan * phi_wingspan_raw;
  phi_park_size = mu_phi_park_size + sigma_phi_park_size * phi_park_size_raw;
  phi_isolation = mu_phi_isolation + sigma_phi_isolation * phi_isolation_raw;
  p_species = sigma_p_species * p_species_raw;
  //p_city = sigma_p_city * p_city_raw;
  
  for(i in 1:n_species){
    for(j in 1:n_sites){    // loop across all sites
      for(k in 1:n_years_minus1){ // loop across all years
  
        psi1[i,j] = inv_logit( // probability (0-1) of occurrence in year 1 is equal to..
          //psi1_0 + 
          psi1_city[city[j]] +
          psi1_species[species[i]] + // a species specific intercept
          psi1_wingspan[city[j]] * wingspan[i] + // a species effect of migratory
          psi1_park_size[city[j]] * park_size[j] + // a site effect of park size
          psi1_isolation[city[j]] * isolation[j] // a site effect of park isolation
          ); // end psi1[j,k]
        
        gamma[i,j,k] = inv_logit( // probability (0-1) of colonization is equal to..
          //gamma0 +
          gamma_city[city[j]] + 
          gamma_species[species[i]] + // a species specific intercept
          gamma_wingspan[city[j]] * wingspan[i] + // a species effect of migratory
          gamma_park_size[city[j]] * park_size[j] + // a site effect of park size
          gamma_isolation[city[j]] * isolation[j] // a site effect of park isolation
          ); // end gamma[i,j,k]
                
        phi[i,j,k] = inv_logit( // probability (0-1) of persistence is equal to..
          //phi0 +
          phi_city[city[j]] + 
          phi_species[species[i]] + // a species specific intercept
          phi_wingspan[city[j]] * wingspan[i] + // a species effect of migratory
          phi_park_size[city[j]] * park_size[j] + // a site effect of park size
          phi_isolation[city[j]] * isolation[j] // a site effect of park isolation
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
            //p0 +
            p_city[city[j]] + 
            p_species[species[i]] + // a species specific intercept
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
  
  // occupancy
  // initial state
  psi1_0 ~ normal(0, 1); // initial occurrence intercept
  psi1_city ~ normal(psi1_0, sigma_psi1_city);
  sigma_psi1_city ~ normal(0, 0.5);
  psi1_species_raw ~ std_normal();
  sigma_psi1_species ~ normal(0, 2);
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
  gamma_city ~ normal(gamma0, sigma_gamma_city);
  sigma_gamma_city ~ normal(0, 0.5);
  gamma_species_raw ~ std_normal();
  sigma_gamma_species ~ normal(0, 1);
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
  phi_city ~ normal(phi0, sigma_phi_city);
  sigma_phi_city ~ normal(0, 0.5);
  phi_species_raw ~ std_normal();
  sigma_phi_species ~ normal(0, 1);
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
  p0 ~ normal(0, 2); // global intercept
  p_city ~ normal(p0, sigma_p_city);
  sigma_p_city ~ normal(0, 0.5);
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
  for(i in 1:n_species){
    for (j in 1:n_sites){
      for (k in 1:n_years){
        
        if(ranges[i,j] > 0){ // if the site [determined by city] is in range of species...
          
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
        
        } // end if(site is in species range)
  
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
