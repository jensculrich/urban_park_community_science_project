// static multi-species occupancy model for pollinators in urban parks
// this is for the full iNat data
// I'm using months as survey blocks and assuming all parks are surveyed in all months in all years
// so long as one or more species from corresponding family was detected
// I've included a month-phenology covariate for detection

// I've also added park site and species trait predictors (contrast with model 3)
// I've made this static so that I could keep only the community sampling events 
// that actually occurred as a vector, with corresponding predictors of equivalent vector length
// Wasn't sure how to incorporate the temporal autocorrelation of a dynamic occ model
// into this structure. I did it this way in attempt to speed up the run times.

// jcu, started november, 2024.


data {
  
  int <lower=1> R; // total number of rows (species/site/year events in which comm sample occurs)
  int <lower=1> n_surveys; // total number of surveys per year (12)
  int<lower=0,upper=1> V[R, n_surveys]; // binary detection / non detection
  int<lower=0,upper=1> V_NA[R, n_surveys]; // NA matrix (0 == no community sample)
  int <lower=1> n_species;
  int<lower=1> species[R]; // vector of species identities
  int <lower=1> n_sites;
  int<lower=1> site[R]; // vector of site identities
  int <lower=1> n_years;
  matrix[R, n_years] X_year; // categorical year dummy matrix
  vector[n_surveys] survey; // surveys (difference from the mean)
  vector[R] feature_diversity;
  vector[R] ease_of_id;
  vector[R] wingspan;
  vector[R] park_size;
  
} // end data

parameters {

  // occurrence state
  vector[n_species] psi_species_raw;  
  real<lower=0> sigma_psi_species;
  vector[n_sites] psi_site_raw;  
  real<lower=0> sigma_psi_site;
  vector[n_years] psi_year;
  real psi_wingspan;
  real psi_park_size;

  // detection
  vector[n_species] p_species_raw;
  real<lower=0> sigma_p_species;
  vector[n_sites] p_site_raw;  
  real<lower=0> sigma_p_site;
  vector[n_years] p_year;
  real p_wingspan;
  real p_feature_diversity;
  real p_ease_of_id;
  real p_park_size;
  vector[n_species] p_date; // phenology peak
  real mu_p_species_date; // community mean
  real<lower=0> sigma_p_species_date; // variation
  vector[n_species] p_date_sq; // decay pattern of phenology
  real mu_p_species_date_sq; // variation
  real<lower=0> sigma_p_species_date_sq; // community mean
  
} // end parameters

transformed parameters {

  // logit scale psi and p
  real psi[R]; // odds of occurrence
  real p[R, n_surveys];  // odds of detection
  
  vector[n_species] psi_species;
  vector[n_sites] psi_site;
  vector[n_species] p_species;
  vector[n_sites] p_site;
  
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  psi_species = sigma_psi_species * psi_species_raw;
  psi_site = sigma_psi_site * psi_site_raw;
  p_species = sigma_p_species * p_species_raw;
  p_site = sigma_p_site * p_site_raw;

  for(i in 1:R){

    psi[i] = inv_logit( // probability (0-1) of occurrence is equal to..
      psi_species[species[i]] + // a species specific intercept
      psi_site[site[i]] + // a species specific intercept
      X_year[i] * psi_year  + // a year effect
      psi_wingspan * wingspan[i] + // a species effect of migratory
      psi_park_size * park_size[i] // a site effect of park size
      ); // end psi1[j,k]
      
    for(j in 1:n_surveys){  
      p[i, j] = inv_logit( // probability (0-1) of detection is equal to..
        p_species[species[i]] + // a species specific intercept
        p_site[site[i]] + // a species specific intercept
        X_year[i] * p_year  + // a year effect
        p_wingspan * wingspan[i] + // a species effect of wingspan
        p_feature_diversity * feature_diversity[i] + // a species effect of feature diversity
        p_ease_of_id * ease_of_id[i] + // a species effect of ease of identification
        p_park_size * park_size[i] + // a site effect of park size
        p_date[species[i]] * survey[j] + // a species-specific phenological detection effect (peak)
        p_date_sq[species[i]] * (survey[j])^2 // a species-specific phenological detection effect (decay)
        ); // end p[j,k,l]
    } // end loop across j
    
  } // end loop across R
   
} // end transformed parameters

model {
  
  // PRIORS
  
  // occupancy
  psi_species_raw ~ std_normal();
  sigma_psi_species ~ normal(0, 1);
  psi_site_raw ~ std_normal();
  sigma_psi_site ~ normal(0, 1);
  psi_year[1] ~ normal(0, 2);
  psi_year[2] ~ normal(0, 0.25);
  psi_year[3] ~ normal(0, 0.25);
  psi_year[4] ~ normal(0, 0.25);
  psi_year[5] ~ normal(0, 0.25);
  psi_wingspan ~ normal(0, 2);
  psi_park_size ~ normal(0, 2);

  // detection
  p_species_raw ~ std_normal();
  sigma_p_species ~ normal(0, 1);
  p_site_raw ~ std_normal();
  sigma_p_site ~ normal(0, 1);
  p_year[1] ~ normal(0, 2);
  p_year[2] ~ normal(0, 0.25);
  p_year[3] ~ normal(0, 0.25);
  p_year[4] ~ normal(0, 0.25);
  p_year[5] ~ normal(0, 0.25);
  p_wingspan ~ normal(0, 2);
  p_feature_diversity ~ normal(0, 2);
  p_ease_of_id ~ normal(0, 2);
  p_park_size ~ normal(0, 2);
  p_date ~ normal(mu_p_species_date, sigma_p_species_date); // species-specific phenology (peak)
  mu_p_species_date ~ normal(0, 2); // mean
  sigma_p_species_date ~ normal(0, 1); // variation
  p_date_sq ~ normal(mu_p_species_date_sq, sigma_p_species_date_sq); // species-specific phenology (decay)
  mu_p_species_date_sq ~ normal(0, 1); // mean
  sigma_p_species_date_sq ~ normal(0, 1); // variation

  // LIKELIHOOD
  for(i in 1:R){
    if (sum(V[i]) > 0){ // lp observed 
      // detection on each visit given detection rate on each visit
      // V_NA == 0 indicates that no survey occurred. Multiplying 0 by the 
      // lpmf() statement should remove it from the target
      target += (log(psi[i]) + // present
                     bernoulli_lpmf(V[i,1]|p[i,1])*V_NA[i,1] + 
                     bernoulli_lpmf(V[i,2]|p[i,2])*V_NA[i,2] + 
                     bernoulli_lpmf(V[i,3]|p[i,3])*V_NA[i,3] + 
                     bernoulli_lpmf(V[i,4]|p[i,4])*V_NA[i,4] +
                     bernoulli_lpmf(V[i,5]|p[i,5])*V_NA[i,5] + 
                     bernoulli_lpmf(V[i,6]|p[i,6])*V_NA[i,6] +
                     bernoulli_lpmf(V[i,7]|p[i,7])*V_NA[i,7] + 
                     bernoulli_lpmf(V[i,8]|p[i,8])*V_NA[i,8] +
                     bernoulli_lpmf(V[i,9]|p[i,9])*V_NA[i,9] + 
                     bernoulli_lpmf(V[i,10]|p[i,10])*V_NA[i,10] +
                     bernoulli_lpmf(V[i,11]|p[i,11])*V_NA[i,11] + 
                     bernoulli_lpmf(V[i,12]|p[i,12])*V_NA[i,12]
                 );
    } else { // lp unobserved (set up for 6 annual surveys)
      // marginal likelihood of:
      // occurrence, yet... non-detection on each visit given detection rate on each visit
      target += (log_sum_exp(log(psi[i]) + 
                 log1m(p[i,1])*V_NA[i,1] + log1m(p[i,2])*V_NA[i,2] + 
                 log1m(p[i,3])*V_NA[i,3] + log1m(p[i,4])*V_NA[i,4] +
                 log1m(p[i,5])*V_NA[i,5] + log1m(p[i,6])*V_NA[i,6] +
                 log1m(p[i,7])*V_NA[i,7] + log1m(p[i,8])*V_NA[i,8] +
                 log1m(p[i,9])*V_NA[i,9] + log1m(p[i,10])*V_NA[i,10] +
                 log1m(p[i,11])*V_NA[i,11] + log1m(p[i,12])*V_NA[i,12],
                 // or just simple no occurrence
                 log1m(psi[i])));
                                  
        } // end if/else
  } // end loop across R
} // end model
