# plot effects of predictors on occurrence and detection

library(tidyverse)
library(cmdstanr)

# enter the region/regions you want to plot
# currently I think this will only work if you enter one single region,
# but I think eventually we want to plot multiple regionss simultaneously


# select a region
regions <- c(
  "Mean"
)

region <- regions[1]

# list of city names

# all
if(region == regions[1]){
  city_names <- c(
    "Atlanta",
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
    "SF"
  )
}


n_cities <- length(city_names)
n_regions <- length(region)

## get param estimates from the region
stan_out <- readRDS(
  "./part1_urban_butterfly_community_dynamics/model_outputs/stan_out_jan29.rds")

# summarise all variables with default and additional summary measures
estimates <- as.data.frame(stan_out$summary(
  variables = c(
    "psi1_0", 
    "sigma_psi1_species",
    "sigma_psi1_city",
    "psi1_wingspan",
    "mu_psi1_park_size",
    "sigma_psi1_park_size",
    "mu_psi1_isolation",
    "sigma_psi1_isolation",
    "psi1_wingspan",
    "psi1_migratory",
    
    "gamma0", 
    "sigma_gamma_species",
    "sigma_gamma_city",
    "mu_gamma_park_size",
    "sigma_gamma_park_size",
    "mu_gamma_isolation",
    "sigma_gamma_isolation",
    "gamma_wingspan",
    "gamma_migratory",
    
    "phi0", 
    "sigma_phi_species",
    "sigma_phi_city",
    "mu_phi_park_size",
    "sigma_phi_park_size",
    "mu_phi_isolation",
    "sigma_phi_isolation",
    "phi_wingspan",
    "phi_migratory",
    
    "p0", 
    "sigma_p_species",
    "sigma_p_city",
    "p_city_detections",
    "p_wingspan",
    "p_migratory",
    "p_feature_diversity",
    "p_ease_of_id",
    "delta0",
    "delta_regional_cluster",
    "sigma_p_species_date",
    "epsilon0",
    "epsilon_regional_cluster",
    "sigma_p_species_date_sq",
    
    # city effects
    "psi1_city",
    "psi1_park_size",
    "psi1_isolation",
    "gamma_city",
    "gamma_park_size",
    "gamma_isolation",
    "phi_city",
    "phi_park_size",
    "phi_isolation",
    "p_city"),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates) <- estimates[, 1]

# handy for viewing column numbers
# this line of code won't work until you've actually read in a stan fit object
#View(cbind(1:nrow(estimates), estimates)) # View to see which row corresponds to the parameter of interest

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours
my_palette <- c("black", my_palette) # add black for the all cities mean

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  # cbind all cities, so we can look at city specific estimates
  # plus region name, so we can look at mean regional estimates
  cities <- c(cities, region_name)
  # length should equal the number of the cities in the region plus 1
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params_re <- 3 # intercept, park size, isolation, 
  params_fixed <- 2 #  wingspan, migratory
  params <- params_re + params_fixed
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  # create a vector of city names. 
  # We will want a df with city name and then param of interest, and then estimates of the param
  # repeated for each param of interest. and repeated for all cities in region
  city_name <- rep(cities, each=(params)) 
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (psi1_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  psi1_0 <- which( rownames(estimates)=="psi1_0" )
  psi1_parksize <- which( rownames(estimates)=="mu_psi1_park_size" )
  psi1_isolation <- which( rownames(estimates)=="mu_psi1_isolation" )
  psi1_wingspan <- which( rownames(estimates)=="psi1_wingspan" )
  psi1_migratory <- which( rownames(estimates)=="psi1_migratory" )
  first_psi1_city <- which( rownames(estimates)=="psi1_city[1]" )
  first_psi1_parksize <- which( rownames(estimates)=="psi1_park_size[1]" )
  first_psi1_isolation <- which( rownames(estimates)=="psi1_isolation[1]" )
  
  # now pull out param estimates for each city within the region
  # random effects
  for(j in 1:(n_cities-1)){
    
    city <- cities[j] # city name
    
    # index tells the loop where to store the param values
    # the index references by which city we are looping across
    # and which param we are interested in.
    index_lower <- 1 + ((j-1) * params_re) + ((j-1) * params_fixed)
    index_upper <- 1 + ((j-1) * params_re) + (params_re - 1) + ((j-1) * params_fixed)
    
    # Y is the mean response, we also plot out the 50 and 95% BCIs
    # using 50 and 95% because these are default values in stanfit summary
    # If you wanted to customize BCIs, transform the fit into a df and calculate desired quantiles.
    Y[index_lower:index_upper] <- c(
      estimates[first_psi1_city+(j-1),2], # psi1 - intercept
      estimates[first_psi1_parksize+(j-1),2],  # psi1 - park size
      estimates[first_psi1_isolation+(j-1),2]  # psi1 - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[first_psi1_city+(j-1),6], # psi1 - intercept
      estimates[first_psi1_parksize+(j-1),6],  # psi1 - park size
      estimates[first_psi1_isolation+(j-1),6]  # psi1 - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[first_psi1_city+(j-1),7], # psi1 - intercept
      estimates[first_psi1_parksize+(j-1),7],  # psi1 - park size
      estimates[first_psi1_isolation+(j-1),7]  # psi1 - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[first_psi1_city+(j-1),8], # psi1 - intercept
      estimates[first_psi1_parksize+(j-1),8],  # psi1 - park size
      estimates[first_psi1_isolation+(j-1),8]  # psi1 - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[first_psi1_city+(j-1),9], # psi1 - intercept
      estimates[first_psi1_parksize+(j-1),9],  # psi1 - park size
      estimates[first_psi1_isolation+(j-1),9]  # psi1 - isolation
    )
    
  }
  
  # fixed effects and means of random effects
  # j (last city in n_cities) actually represents the "mean" city or regional average
  # for this regional average, just don't add any city specific random effects
  # that would tell us how much a city deviates from the regional mean
  for(j in n_cities){ # will only index one unit (the regional average)
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params) 
    index_upper <- length(Y)
    
    Y[index_lower:index_upper] <- c(
      estimates[psi1_0,2], 
      estimates[psi1_parksize,2],
      estimates[psi1_isolation,2],
      estimates[psi1_wingspan,2],
      estimates[psi1_migratory,2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[psi1_0,6],
      estimates[psi1_parksize,6], 
      estimates[psi1_isolation,6], 
      estimates[psi1_wingspan,6], 
      estimates[psi1_migratory,6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[psi1_0,7],
      estimates[psi1_parksize,7],
      estimates[psi1_isolation,7],
      estimates[psi1_wingspan,7],
      estimates[psi1_migratory,7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[psi1_0,8],
      estimates[psi1_parksize,8],
      estimates[psi1_isolation,8],
      estimates[psi1_wingspan,8],
      estimates[psi1_migratory,8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[psi1_0,9],
      estimates[psi1_parksize,9],
      estimates[psi1_isolation,9],
      estimates[psi1_wingspan,9],
      estimates[psi1_migratory,9]
    )
    
  }
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$city_name <- as.factor(df_estimates$city_name)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  
  
}

# plot cities in alphabetical order
df_estimates$city_name <- fct_relevel(df_estimates$city_name, region_name)
                                                     
df_estimates1 <- df_estimates %>%
  filter(city_name == region_name)

df_estimates2 <- df_estimates %>%
  filter(city_name != region_name) %>%
  # filter to the params that have city-specific estimates to display
  filter(X %in% c(1, 2, 3))

temp <- df_estimates1 %>%
  filter(X %in% c(1, 2, 3))

df_estimates2 <- left_join(df_estimates2, temp, by = c("X", "city_name")) %>%
  mutate(Y = Y.x,
         lower_95 = lower_95.x,
         upper_95 = upper_95.x,
         lower_50 = lower_50.x,
         upper_50 =  upper_50.x) %>%
  select(X, city_name, Y, lower_95, upper_95, lower_50, upper_50)

df_estimates <- rbind(df_estimates1, df_estimates2)

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
p <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(psi["intercept"]),
                             bquote(psi["park size"]),
                             bquote(psi["isolation"]),
                             bquote(psi["wingspan"]),
                             bquote(psi["migratory"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3, 4), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=my_palette) + 
   geom_hline(yintercept = 0, lty = "dashed") +
   ggtitle("Initial Occurrence") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() 

# add estimates
p <- p +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 

# plot the plot
p

#-------------------------------------------------------------------------------
# everything below replicates the above but for the other processes (colonization, persistence, detection)
# I may not have commented everything as well yet.

#-------------------------------------------------------------------------------
# colonization (gamma)

for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  # cbind all cities, so we can look at city specific estimates
  # plus region name, so we can look at mean regional estimates
  cities <- c(cities, region_name)
  # length should equal the number of the cities in the region plus 1
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params_re <- 3 # intercept, park size, isolation, 
  params_fixed <- 2 #  wingspan, migratory
  params <- params_re + params_fixed
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  # create a vector of city names. 
  # We will want a df with city name and then param of interest, and then estimates of the param
  # repeated for each param of interest. and repeated for all cities in region
  city_name <- rep(cities, each=(params)) 
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (gamma_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  gamma_0 <- which( rownames(estimates)=="gamma0" )
  gamma_parksize <- which( rownames(estimates)=="mu_gamma_park_size" )
  gamma_isolation <- which( rownames(estimates)=="mu_gamma_isolation" )
  gamma_wingspan <- which( rownames(estimates)=="gamma_wingspan" )
  gamma_migratory <- which( rownames(estimates)=="gamma_migratory" )
  first_gamma_city <- which( rownames(estimates)=="gamma_city[1]" )
  first_gamma_parksize <- which( rownames(estimates)=="gamma_park_size[1]" )
  first_gamma_isolation <- which( rownames(estimates)=="gamma_isolation[1]" )
  
  # now pull out param estimates for each city within the region
  # random effects
  for(j in 1:(n_cities-1)){
    
    city <- cities[j] # city name
    
    # index tells the loop where to store the param values
    # the index references by which city we are looping across
    # and which param we are interested in.
    index_lower <- 1 + ((j-1) * params_re) + ((j-1) * params_fixed)
    index_upper <- 1 + ((j-1) * params_re) + (params_re - 1) + ((j-1) * params_fixed)
    
    # Y is the mean response, we also plot out the 50 and 95% BCIs
    # using 50 and 95% because these are default values in stanfit summary
    # If you wanted to customize BCIs, transform the fit into a df and calculate desired quantiles.
    Y[index_lower:index_upper] <- c(
      estimates[first_gamma_city+(j-1),2], # gamma - intercept
      estimates[first_gamma_parksize+(j-1),2],  # gamma - park size
      estimates[first_gamma_isolation+(j-1),2]  # gamma - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[first_gamma_city+(j-1),6], # gamma - intercept
      estimates[first_gamma_parksize+(j-1),6],  # gamma - park size
      estimates[first_gamma_isolation+(j-1),6]  # gamma - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[first_gamma_city+(j-1),7], # gamma - intercept
      estimates[first_gamma_parksize+(j-1),7],  # gamma - park size
      estimates[first_gamma_isolation+(j-1),7]  # gamma - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[first_gamma_city+(j-1),8], # gamma - intercept
      estimates[first_gamma_parksize+(j-1),8],  # gamma - park size
      estimates[first_gamma_isolation+(j-1),8]  # gamma - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[first_gamma_city+(j-1),9], # gamma - intercept
      estimates[first_gamma_parksize+(j-1),9],  # gamma - park size
      estimates[first_gamma_isolation+(j-1),9]  # gamma - isolation
    )
    
  }
  
  # fixed effects and means of random effects
  # j (last city in n_cities) actually represents the "mean" city or regional average
  # for this regional average, just don't add any city specific random effects
  # that would tell us how much a city deviates from the regional mean
  for(j in n_cities){ # will only index one unit (the regional average)
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params) 
    index_upper <- length(Y)
    
    Y[index_lower:index_upper] <- c(
      estimates[gamma_0,2], 
      estimates[gamma_parksize,2],
      estimates[gamma_isolation,2],
      estimates[gamma_wingspan,2],
      estimates[gamma_migratory,2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[gamma_0,6],
      estimates[gamma_parksize,6], 
      estimates[gamma_isolation,6], 
      estimates[gamma_wingspan,6], 
      estimates[gamma_migratory,6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[gamma_0,7],
      estimates[gamma_parksize,7],
      estimates[gamma_isolation,7],
      estimates[gamma_wingspan,7],
      estimates[gamma_migratory,7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[gamma_0,8],
      estimates[gamma_parksize,8],
      estimates[gamma_isolation,8],
      estimates[gamma_wingspan,8],
      estimates[gamma_migratory,8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[gamma_0,9],
      estimates[gamma_parksize,9],
      estimates[gamma_isolation,9],
      estimates[gamma_wingspan,9],
      estimates[gamma_migratory,9]
    )
    
  }
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$city_name <- as.factor(df_estimates$city_name)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  
  
}

# plot cities in alphabetical order
df_estimates$city_name <- fct_relevel(df_estimates$city_name, region_name)

df_estimates1 <- df_estimates %>%
  filter(city_name == region_name)

df_estimates2 <- df_estimates %>%
  filter(city_name != region_name) %>%
  # filter to the params that have city-specific estimates to display
  filter(X %in% c(1, 2, 3))

temp <- df_estimates1 %>%
  filter(X %in% c(1, 2, 3))

df_estimates2 <- left_join(df_estimates2, temp, by = c("X", "city_name")) %>%
  mutate(Y = Y.x,
         lower_95 = lower_95.x,
         upper_95 = upper_95.x,
         lower_50 = lower_50.x,
         upper_50 =  upper_50.x) %>%
  select(X, city_name, Y, lower_95, upper_95, lower_50, upper_50)

df_estimates <- rbind(df_estimates1, df_estimates2)

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
q <- ggplot(df_estimates) +
  theme_bw() +
  scale_x_discrete(name="", breaks = seq(1:params),
                   labels=c(bquote(gamma["intercept"]),
                            bquote(gamma["park size"]),
                            bquote(gamma["isolation"]),
                            bquote(gamma["wingspan"]),
                            bquote(gamma["migratory"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-10, 4), breaks = c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8)) +
  guides(color = guide_legend(title = "city")) +
  scale_color_manual(values=my_palette) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  ggtitle("Colonization") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip() 

# add estimates
q <- q +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 

# plot the plot
q


#-------------------------------------------------------------------------------
# persistence (phi)

for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  # cbind all cities, so we can look at city specific estimates
  # plus region name, so we can look at mean regional estimates
  cities <- c(cities, region_name)
  # length should equal the number of the cities in the region plus 1
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params_re <- 3 # intercept, park size, isolation, 
  params_fixed <- 2 #  wingspan, migratory
  params <- params_re + params_fixed
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  # create a vector of city names. 
  # We will want a df with city name and then param of interest, and then estimates of the param
  # repeated for each param of interest. and repeated for all cities in region
  city_name <- rep(cities, each=(params)) 
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (phi_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  phi_0 <- which( rownames(estimates)=="phi0" )
  phi_parksize <- which( rownames(estimates)=="mu_phi_park_size" )
  phi_isolation <- which( rownames(estimates)=="mu_phi_isolation" )
  phi_wingspan <- which( rownames(estimates)=="phi_wingspan" )
  phi_migratory <- which( rownames(estimates)=="phi_migratory" )
  first_phi_city <- which( rownames(estimates)=="phi_city[1]" )
  first_phi_parksize <- which( rownames(estimates)=="phi_park_size[1]" )
  first_phi_isolation <- which( rownames(estimates)=="phi_isolation[1]" )
  
  # now pull out param estimates for each city within the region
  # random effects
  for(j in 1:(n_cities-1)){
    
    city <- cities[j] # city name
    
    # index tells the loop where to store the param values
    # the index references by which city we are looping across
    # and which param we are interested in.
    index_lower <- 1 + ((j-1) * params_re) + ((j-1) * params_fixed)
    index_upper <- 1 + ((j-1) * params_re) + (params_re - 1) + ((j-1) * params_fixed)
    
    # Y is the mean response, we also plot out the 50 and 95% BCIs
    # using 50 and 95% because these are default values in stanfit summary
    # If you wanted to customize BCIs, transform the fit into a df and calculate desired quantiles.
    Y[index_lower:index_upper] <- c(
      estimates[first_phi_city+(j-1),2], # phi - intercept
      estimates[first_phi_parksize+(j-1),2],  # phi - park size
      estimates[first_phi_isolation+(j-1),2]  # phi - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[first_phi_city+(j-1),6], # phi - intercept
      estimates[first_phi_parksize+(j-1),6],  # phi - park size
      estimates[first_phi_isolation+(j-1),6]  # phi - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[first_phi_city+(j-1),7], # phi - intercept
      estimates[first_phi_parksize+(j-1),7],  # phi - park size
      estimates[first_phi_isolation+(j-1),7]  # phi - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[first_phi_city+(j-1),8], # phi - intercept
      estimates[first_phi_parksize+(j-1),8],  # phi - park size
      estimates[first_phi_isolation+(j-1),8]  # phi - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[first_phi_city+(j-1),9], # phi - intercept
      estimates[first_phi_parksize+(j-1),9],  # phi - park size
      estimates[first_phi_isolation+(j-1),9]  # phi - isolation
    )
    
  }
  
  # fixed effects and means of random effects
  # j (last city in n_cities) actually represents the "mean" city or regional average
  # for this regional average, just don't add any city specific random effects
  # that would tell us how much a city deviates from the regional mean
  for(j in n_cities){ # will only index one unit (the regional average)
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params) 
    index_upper <- length(Y)
    
    Y[index_lower:index_upper] <- c(
      estimates[phi_0,2], 
      estimates[phi_parksize,2],
      estimates[phi_isolation,2],
      estimates[phi_wingspan,2],
      estimates[phi_migratory,2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[phi_0,6],
      estimates[phi_parksize,6], 
      estimates[phi_isolation,6], 
      estimates[phi_wingspan,6], 
      estimates[phi_migratory,6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[phi_0,7],
      estimates[phi_parksize,7],
      estimates[phi_isolation,7],
      estimates[phi_wingspan,7],
      estimates[phi_migratory,7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[phi_0,8],
      estimates[phi_parksize,8],
      estimates[phi_isolation,8],
      estimates[phi_wingspan,8],
      estimates[phi_migratory,8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[phi_0,9],
      estimates[phi_parksize,9],
      estimates[phi_isolation,9],
      estimates[phi_wingspan,9],
      estimates[phi_migratory,9]
    )
    
  }
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$city_name <- as.factor(df_estimates$city_name)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  
  
}

# plot cities in alphabetical order
df_estimates$city_name <- fct_relevel(df_estimates$city_name, region_name)

df_estimates1 <- df_estimates %>%
  filter(city_name == region_name)

df_estimates2 <- df_estimates %>%
  filter(city_name != region_name) %>%
  # filter to the params that have city-specific estimates to display
  filter(X %in% c(1, 2, 3))

temp <- df_estimates1 %>%
  filter(X %in% c(1, 2, 3))

df_estimates2 <- left_join(df_estimates2, temp, by = c("X", "city_name")) %>%
  mutate(Y = Y.x,
         lower_95 = lower_95.x,
         upper_95 = upper_95.x,
         lower_50 = lower_50.x,
         upper_50 =  upper_50.x) %>%
  select(X, city_name, Y, lower_95, upper_95, lower_50, upper_50)

df_estimates <- rbind(df_estimates1, df_estimates2)

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
r <- ggplot(df_estimates) +
  theme_bw() +
  scale_x_discrete(name="", breaks = seq(1:params),
                   labels=c(bquote(phi["intercept"]),
                            bquote(phi["park size"]),
                            bquote(phi["isolation"]),
                            bquote(phi["wingspan"]),
                            bquote(phi["migratory"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-4, 10), breaks = c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10)) +
  guides(color = guide_legend(title = "city")) +
  scale_color_manual(values=my_palette) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  ggtitle("Persistence") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.text=element_text(size=10),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip() 

# add estimates
r <- r +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 

# plot the plot
r

#-------------------------------------------------------------------------------
# detection (p)


# try a way that pulls out the region means simultaneously
for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign object from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  cities <- c(cities, region_name)
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params_re <- 3
  params_fixed <- 4
  params <- params_re + params_fixed
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params)
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  city_name <- rep(cities, each=(params)) 
  
  # get indices for species random effects distributions for particular city
  p0 <- which( rownames(estimates)=="p0" )
  p_wingspan <- which( rownames(estimates)=="p_wingspan" )
  p_migratory <- which( rownames(estimates)=="p_migratory" )
  p_feature_diversity <- which( rownames(estimates)=="p_feature_diversity" )
  p_ease_of_id <- which( rownames(estimates)=="p_ease_of_id" )
  delta0 <- which( rownames(estimates)=="delta0" )
  epsilon0 <- which( rownames(estimates)=="epsilon0" )
  first_delta_regional_cluster <- which( rownames(estimates)=="delta_regional_cluster[1]" )
  first_epsilon_regional_cluster <- which( rownames(estimates)=="epsilon_regional_cluster[1]" )
  first_p_city <- which( rownames(estimates)=="p_city[1]" )
  
  # regions
  # california = 1
  # midwest = 2
  # northeast = 3
  # san fransisco = 4
  # southeast = 5
  # texas = 6
  
  #"Atlanta", 5
  #"Boston", 3
  #"Charlotte", 5
  #"Chicago", 2
  #"Dallas", 6
  #"DC", 3
  #"Denton", 6
  #"Houston", 6
  #"LA", 1
  #"Minneapolis", 2
  #"NYC", 3    
  #"Philadelphia", 3
  #"Raleigh", 5
  #"SD", 1
  #"SF" 4
  regional_cluster <- c(5,3,5,2,6,3,6,6,1,2,3,3,5,1,4)
  
  for(j in 1:(n_cities-1)){
    
    city <- cities[j]
    
    regionl_cluster_of_city <- regional_cluster[j]
    
    index_lower <- 1 + ((j-1) * params_re) + ((j-1) * params_fixed)
    index_upper <- 1 + ((j-1) * params_re) + (params_re - 1) + ((j-1) * params_fixed)
    
    Y[index_lower:index_upper] <- c(
      estimates[first_p_city+(j-1),2], # p - intercept
      estimates[first_delta_regional_cluster+(regionl_cluster_of_city-1),2],
      estimates[first_epsilon_regional_cluster+(regionl_cluster_of_city-1),2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[first_p_city+(j-1),6], # p - intercept
      estimates[first_delta_regional_cluster+(regionl_cluster_of_city-1),6],
      estimates[first_epsilon_regional_cluster+(regionl_cluster_of_city-1),6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[first_p_city+(j-1),7], # p - intercept
      estimates[first_delta_regional_cluster+(regionl_cluster_of_city-1),7],
      estimates[first_epsilon_regional_cluster+(regionl_cluster_of_city-1),7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[first_p_city+(j-1),8], # p - intercept
      estimates[first_delta_regional_cluster+(regionl_cluster_of_city-1),8],
      estimates[first_epsilon_regional_cluster+(regionl_cluster_of_city-1),8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[first_p_city+(j-1),9], # p - intercept
      estimates[first_delta_regional_cluster+(regionl_cluster_of_city-1),9],
      estimates[first_epsilon_regional_cluster+(regionl_cluster_of_city-1),9]
    )
    
  }
  
  for(j in n_cities){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params) 
    index_upper <- length(Y)
    
    Y[index_lower:index_upper] <- c(
      estimates[p0,2], 
      estimates[delta0,2],
      estimates[epsilon0,2],
      estimates[p_wingspan,2],
      estimates[p_migratory,2],
      estimates[p_feature_diversity,2],
      estimates[p_ease_of_id,2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[p0,6],
      estimates[delta0,6],
      estimates[epsilon0,6], 
      estimates[p_wingspan,6],
      estimates[p_migratory,6],
      estimates[p_feature_diversity,6],
      estimates[p_ease_of_id,6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[p0,7],
      estimates[delta0,7],
      estimates[epsilon0,7], 
      estimates[p_wingspan,7],
      estimates[p_migratory,7],
      estimates[p_feature_diversity,7],
      estimates[p_ease_of_id,7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[p0,8],
      estimates[delta0,8],
      estimates[epsilon0,8], 
      estimates[p_wingspan,8],
      estimates[p_migratory,8],
      estimates[p_feature_diversity,8],
      estimates[p_ease_of_id,8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[p0,9],
      estimates[delta0,9],
      estimates[epsilon0,9], 
      estimates[p_wingspan,9],
      estimates[p_migratory,9],
      estimates[p_feature_diversity,9],
      estimates[p_ease_of_id,9]
    )
    
  }
  
  df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$city_name <- as.factor(df_estimates$city_name)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  
  
}


df_estimates$city_name <- fct_relevel(df_estimates$city_name, region_name)

df_estimates1 <- df_estimates %>%
  filter(city_name == region_name)
  
df_estimates2 <- df_estimates %>%
  filter(city_name != region_name) %>%
  # filter to the params that have city-specific estimates to display
  filter(X %in% c(1, 2, 3))

temp <- df_estimates1 %>%
  filter(X %in% c(1, 2, 3))

df_estimates2 <- left_join(df_estimates2, temp, by = "X") %>%
  mutate(Y = ifelse(X %in% c(2, 3), (Y.x + Y.y), Y.x),
         lower_95 = ifelse(X %in% c(2, 3), (lower_95.x + lower_95.y), lower_95.x),
         upper_95 = ifelse(X %in% c(2, 3), (upper_95.x + upper_95.y), upper_95.x),
         lower_50 = ifelse(X %in% c(2, 3), (lower_50.x + lower_50.y), lower_50.x),
         upper_50 = ifelse(X %in% c(2, 3), (upper_50.x + upper_50.y), upper_50.x)) %>%
  rename("city_name" = "city_name.x") %>%
  select(X, city_name, Y, lower_95, upper_95, lower_50, upper_50)

df_estimates <- rbind(df_estimates1, df_estimates2)

## --------------------------------------------------
## Draw caterpillar plot

s <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(p["intercept"]),
                             bquote(p["peak phenology"]),
                             bquote(p["phenology decay"]),
                             bquote(p["wingspan"]),
                             bquote(p["migratory"]),
                             bquote(p["ft. diversity"]),
                             bquote(p["ease of ID"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-4, 4), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=my_palette) + 
                                        geom_hline(yintercept = 0, lty = "dashed") +
   ggtitle("Detection") +
   theme(plot.title = element_text(size = 18, face = "bold"),
         legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() 

s <- s +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
s

#-------------------------------------------------------------------------------
# plot the 4 panels on a 2x2 grid

cowplot::plot_grid(p, q, r, s, ncol = 2)


ggplot(site_data) +
  ggridges::geom_density_ridges(aes(x=log_total_green_space_area, y=city, fill = city), alpha = 0.3)

ggplot(site_data) +
  ggridges::geom_density_ridges(aes(x=log_isolation_scaled_across_all_cities, y=city, fill = city), alpha = 0.3)
