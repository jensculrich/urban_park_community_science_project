# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

# enter the region/regions you want to plot
# currently I think this will only work if you enter one single region,
# but I think eventually we want to plot multiple regionss simultaneously


# select a region
regions <- c(
  "midwest",
  "northeast",
  "southeast",
  "southeast_atlantic",
  "southeast_texas",
  "southwest"
)

region <- regions[4]

# list of city names

# midwest
if(region == regions[1]){
  city_names <- c(
    "Chicago",
    "Denver",
    "Des_Moines",
    "Detroit", 
    "Minneapolis",
    "St_Louis"
  )
}

# northeast
if(region == regions[2]){
  city_names <- c(
    "Boston", 
    "DC",
    "NYC", 
    "Philadelphia"
  )
}

# southeast
if(region == regions[3]){
  city_names <- c(
    "Atlanta",
    "Charlotte",
    "Dallas",
    "Denton",
    "Houston",
    "Raleigh"
  )
}

# southeast_atlantic
if(region == regions[4]){
  city_names <- c(
    "Atlanta",
    "Charlotte",
    "Raleigh"
  )
}

# southeast_texas
if(region == regions[5]){
  city_names <- c(
    "Dallas",
    "Denton",
    "Houston"
  )
}

# southwest
if(region == regions[6]){
  city_names <- c(
    "LA",
    "Phoenix",
    "Riverside",
    "SD",
    "SF"
  )
} 

n_cities <- length(city_names)

# handy for viewing column numbers
# this line of code won't work until you've actually read in a stan fit object
View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

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
  params <- 4  # 4 for intercept, park size, isolation, wingspan
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
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", region, ".rds"))
  fit_summary <- rstan::summary(stan_out)
  estimates <- as.data.frame(fit_summary)
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (psi1_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  psi1_0 <- which( rownames(estimates)=="psi1_0" )
  psi1_wingspan <- which( rownames(estimates)=="mu_psi1_wingspan" )
  psi1_parksize <- which( rownames(estimates)=="mu_psi1_park_size" )
  psi1_isolation <- which( rownames(estimates)=="mu_psi1_isolation" )
  first_psi1_city <- which( rownames(estimates)=="psi1_city[1]" )
  first_psi1_wingspan <- which( rownames(estimates)=="psi1_wingspan[1]" )
  first_psi1_parksize <- which( rownames(estimates)=="psi1_park_size[1]" )
  first_psi1_isolation <- which( rownames(estimates)=="psi1_isolation[1]" )
  
  # now pull out param estimates for each city within the region
  for(j in 1:(n_cities-1)){
    
    city <- cities[j] # city name
    
    # index tells the loop where to store the param values
    # the index references by which city we are looping across
    # and which param we are interested in.
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    # Y is the mean response, we also plot out the 50 and 95% BCIs
    # using 50 and 95% because these are default values in stanfit summary
    # If you wanted to customize BCIs, transform the fit into a df and calculate desired quantiles.
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[first_psi1_city+(j-1),1], # psi1 - intercept
      fit_summary$summary[first_psi1_wingspan+(j-1),1],  # psi1 - wingspan
      fit_summary$summary[first_psi1_parksize+(j-1),1],  # psi1 - park size
      fit_summary$summary[first_psi1_isolation+(j-1),1]  # psi1 - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_psi1_city+(j-1),4], # psi1 - intercept
      fit_summary$summary[first_psi1_wingspan+(j-1),4],  # psi1 - wingspan
      fit_summary$summary[first_psi1_parksize+(j-1),4],  # psi1 - park size
      fit_summary$summary[first_psi1_isolation+(j-1),4]  # psi1 - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_psi1_city+(j-1),8], # psi1 - intercept
      fit_summary$summary[first_psi1_wingspan+(j-1),8],  # psi1 - wingspan
      fit_summary$summary[first_psi1_parksize+(j-1),8],  # psi1 - park size
      fit_summary$summary[first_psi1_isolation+(j-1),8]  # psi1 - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_psi1_city+(j-1),5], # psi1 - intercept
      fit_summary$summary[first_psi1_wingspan+(j-1),5],  # psi1 - wingspan
      fit_summary$summary[first_psi1_parksize+(j-1),5],  # psi1 - park size
      fit_summary$summary[first_psi1_isolation+(j-1),5]  # psi1 - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_psi1_city+(j-1),7], # psi1 - intercept
      fit_summary$summary[first_psi1_wingspan+(j-1),7],  # psi1 - wingspan
      fit_summary$summary[first_psi1_parksize+(j-1),7],  # psi1 - park size
      fit_summary$summary[first_psi1_isolation+(j-1),7]  # psi1 - isolation
    )
    
  }
  
  # j (last city in n_cities) actually represents the "mean" city or regional average
  # for this regional average, just don't add any city specific random effects
  # that would tell us how much a city deviates from the regional mean
  for(j in n_cities){ # will only index one unit (the regional average)
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[psi1_0,1], 
      fit_summary$summary[psi1_wingspan,1],
      fit_summary$summary[psi1_parksize,1],
      fit_summary$summary[psi1_isolation,1]
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[psi1_0,4], 
      fit_summary$summary[psi1_wingspan,4], 
      fit_summary$summary[psi1_parksize,4], 
      fit_summary$summary[psi1_isolation,4]
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[psi1_0,8],
      fit_summary$summary[psi1_wingspan,8],
      fit_summary$summary[psi1_parksize,8],
      fit_summary$summary[psi1_isolation,8]
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[psi1_0,5],
      fit_summary$summary[psi1_wingspan,5],
      fit_summary$summary[psi1_parksize,5],
      fit_summary$summary[psi1_isolation,5]
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[psi1_0,7],
      fit_summary$summary[psi1_wingspan,7],
      fit_summary$summary[psi1_parksize,7],
      fit_summary$summary[psi1_isolation,7]
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
                                                     

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
p <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(psi["intercept"]),
                             bquote(psi["wingspan"]),
                             bquote(psi["park size"]),
                             bquote(psi["isolation"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3, 3), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=c("black", "#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                      "#1a5acd", "#E69F90", "#FFFF00")) + 
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

# try a way that pulls out the region means simultaneously
for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  cities <- c(cities, region_name)
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params <- 4
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params)
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  city_name <- rep(cities, each=(params)) 
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", region, ".rds"))
  fit_summary <- rstan::summary(stan_out)
  estimates <- as.data.frame(fit_summary)
  
  # get indices for species random effects distributions for particular city
  gamma0 <- which( rownames(estimates)=="gamma0" )
  gamma_wingspan <- which( rownames(estimates)=="mu_gamma_wingspan" )
  gamma_parksize <- which( rownames(estimates)=="mu_gamma_park_size" )
  gamma_isolation <- which( rownames(estimates)=="mu_gamma_isolation" )
  first_gamma_city <- which( rownames(estimates)=="gamma_city[1]" )
  first_gamma_wingspan <- which( rownames(estimates)=="gamma_wingspan[1]" )
  first_gamma_parksize <- which( rownames(estimates)=="gamma_park_size[1]" )
  first_gamma_isolation <- which( rownames(estimates)=="gamma_isolation[1]" )
  
  for(j in 1:(n_cities-1)){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[first_gamma_city+(j-1),1], # gamma - intercept
      fit_summary$summary[first_gamma_wingspan+(j-1),1],  # gamma - wingspan
      fit_summary$summary[first_gamma_parksize+(j-1),1],  # gamma - park size
      fit_summary$summary[first_gamma_isolation+(j-1),1]  # gamma - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_gamma_city+(j-1),4], # gamma - intercept
      fit_summary$summary[first_gamma_wingspan+(j-1),4],  # gamma - wingspan
      fit_summary$summary[first_gamma_parksize+(j-1),4],  # gamma - park size
      fit_summary$summary[first_gamma_isolation+(j-1),4]  # gamma - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_gamma_city+(j-1),8], # gamma - intercept
      fit_summary$summary[first_gamma_wingspan+(j-1),8],  # gamma - wingspan
      fit_summary$summary[first_gamma_parksize+(j-1),8],  # gamma - park size
      fit_summary$summary[first_gamma_isolation+(j-1),8]  # gamma - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_gamma_city+(j-1),5], # gamma - intercept
      fit_summary$summary[first_gamma_wingspan+(j-1),5],  # gamma - wingspan
      fit_summary$summary[first_gamma_parksize+(j-1),5],  # gamma - park size
      fit_summary$summary[first_gamma_isolation+(j-1),5]  # gamma - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_gamma_city+(j-1),7], # gamma - intercept
      fit_summary$summary[first_gamma_wingspan+(j-1),7],  # gamma - wingspan
      fit_summary$summary[first_gamma_parksize+(j-1),7],  # gamma - park size
      fit_summary$summary[first_gamma_isolation+(j-1),7]  # gamma - isolation
    )
    
  }
  
  for(j in n_cities){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[gamma0,1], 
      fit_summary$summary[gamma_wingspan,1],
      fit_summary$summary[gamma_parksize,1],
      fit_summary$summary[gamma_isolation,1]
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[gamma0,4], 
      fit_summary$summary[gamma_wingspan,4], 
      fit_summary$summary[gamma_parksize,4], 
      fit_summary$summary[gamma_isolation,4]
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[gamma0,8],
      fit_summary$summary[gamma_wingspan,8],
      fit_summary$summary[gamma_parksize,8],
      fit_summary$summary[gamma_isolation,8]
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[gamma0,5],
      fit_summary$summary[gamma_wingspan,5],
      fit_summary$summary[gamma_parksize,5],
      fit_summary$summary[gamma_isolation,5]
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[gamma0,7],
      fit_summary$summary[gamma_wingspan,7],
      fit_summary$summary[gamma_parksize,7],
      fit_summary$summary[gamma_isolation,7]
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


## --------------------------------------------------
## Draw caterpillar plot

q <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(gamma["intercept"]),
                             bquote(gamma["wingspan"]),
                             bquote(gamma["park size"]),
                             bquote(gamma["isolation"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-8, 5), breaks = c(-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=c("black", "#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                      "#1a5acd", "#E69F90", "#FFFF00")) + 
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

q <- q +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
q


#-------------------------------------------------------------------------------
# persistence (phi)

# try a way that pulls out the region means simultaneously
for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  cities <- c(cities, region_name)
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params <- 4
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params)
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  city_name <- rep(cities, each=(params)) 
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", region, ".rds"))
  fit_summary <- rstan::summary(stan_out)
  estimates <- as.data.frame(fit_summary)
  
  # get indices for species random effects distributions for particular city
  phi0 <- which( rownames(estimates)=="phi0" )
  phi_wingspan <- which( rownames(estimates)=="mu_phi_wingspan" )
  phi_parksize <- which( rownames(estimates)=="mu_phi_park_size" )
  phi_isolation <- which( rownames(estimates)=="mu_phi_isolation" )
  first_phi_city <- which( rownames(estimates)=="phi_city[1]" )
  first_phi_wingspan <- which( rownames(estimates)=="phi_wingspan[1]" )
  first_phi_parksize <- which( rownames(estimates)=="phi_park_size[1]" )
  first_phi_isolation <- which( rownames(estimates)=="phi_isolation[1]" )
  
  for(j in 1:(n_cities-1)){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[first_phi_city+(j-1),1], # phi - intercept
      fit_summary$summary[first_phi_wingspan+(j-1),1],  # phi - wingspan
      fit_summary$summary[first_phi_parksize+(j-1),1],  # phi - park size
      fit_summary$summary[first_phi_isolation+(j-1),1]  # phi - isolation
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_phi_city+(j-1),4], # phi - intercept
      fit_summary$summary[first_phi_wingspan+(j-1),4],  # phi - wingspan
      fit_summary$summary[first_phi_parksize+(j-1),4],  # phi - park size
      fit_summary$summary[first_phi_isolation+(j-1),4]  # phi - isolation
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_phi_city+(j-1),8], # phi - intercept
      fit_summary$summary[first_phi_wingspan+(j-1),8],  # phi - wingspan
      fit_summary$summary[first_phi_parksize+(j-1),8],  # phi - park size
      fit_summary$summary[first_phi_isolation+(j-1),8]  # phi - isolation
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_phi_city+(j-1),5], # phi - intercept
      fit_summary$summary[first_phi_wingspan+(j-1),5],  # phi - wingspan
      fit_summary$summary[first_phi_parksize+(j-1),5],  # phi - park size
      fit_summary$summary[first_phi_isolation+(j-1),5]  # phi - isolation
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_phi_city+(j-1),7], # phi - intercept
      fit_summary$summary[first_phi_wingspan+(j-1),7],  # phi - wingspan
      fit_summary$summary[first_phi_parksize+(j-1),7],  # phi - park size
      fit_summary$summary[first_phi_isolation+(j-1),7]  # phi - isolation
    )
    
  }
  
  for(j in n_cities){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[phi0,1], 
      fit_summary$summary[phi_wingspan,1],
      fit_summary$summary[phi_parksize,1],
      fit_summary$summary[phi_isolation,1]
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[phi0,4], 
      fit_summary$summary[phi_wingspan,4], 
      fit_summary$summary[phi_parksize,4], 
      fit_summary$summary[phi_isolation,4]
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[phi0,8],
      fit_summary$summary[phi_wingspan,8],
      fit_summary$summary[phi_parksize,8],
      fit_summary$summary[phi_isolation,8]
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[phi0,5],
      fit_summary$summary[phi_wingspan,5],
      fit_summary$summary[phi_parksize,5],
      fit_summary$summary[phi_isolation,5]
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[phi0,7],
      fit_summary$summary[phi_wingspan,7],
      fit_summary$summary[phi_parksize,7],
      fit_summary$summary[phi_isolation,7]
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


## --------------------------------------------------
## Draw caterpillar plot

r <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(phi["intercept"]),
                             bquote(phi["wingspan"]),
                             bquote(phi["park size"]),
                             bquote(phi["isolation"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-4, 7), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=c("black", "#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                      "#1a5acd", "#E69F90", "#FFFF00")) + 
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

r <- r +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
r

#-------------------------------------------------------------------------------
# detection (p)


# try a way that pulls out the region means simultaneously
for(i in 1:n_regions){
  
  # get region
  #region_name <- region_names[i]
  region_name <- region
  
  # list of city names # assign objectw from strings
  #cities <- eval(parse(text=paste0("cities_", region_name)))
  cities <- city_names
  
  # get number of cities from the region
  cities <- c(cities, region_name)
  n_cities <- length(cities)
  
  # ecological params
  # number of params to plot
  params <- 6
  X <- rep(seq(1:params), times=n_cities) #  ecological params of interest
  Y <- vector(length = n_cities*params)
  lower_95 <- vector(length = n_cities*params)
  upper_95 <- vector(length = n_cities*params)
  lower_50 <- vector(length = n_cities*params)
  upper_50 <- vector(length = n_cities*params)
  
  city_name <- rep(cities, each=(params)) 
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", region, ".rds"))
  fit_summary <- rstan::summary(stan_out)
  estimates <- as.data.frame(fit_summary)
  
  # get indices for species random effects distributions for particular city
  p0 <- which( rownames(estimates)=="p0" )
  p_wingspan <- which( rownames(estimates)=="p_wingspan" )
  p_feature_diversity <- which( rownames(estimates)=="p_feature_diversity" )
  p_ease_of_id <- which( rownames(estimates)=="p_ease_of_id" )
  mu_p_species_date <- which( rownames(estimates)=="mu_p_species_date" )
  mu_p_species_date_sq <- which( rownames(estimates)=="mu_p_species_date_sq" )
  first_p_city <- which( rownames(estimates)=="p_city[1]" )
  
  for(j in 1:(n_cities-1)){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[first_p_city+(j-1),1] # p - intercept
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_p_city+(j-1),4] # p - intercept
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[first_p_city+(j-1),8] # p - intercept
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_p_city+(j-1),5] # p - intercept
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[first_p_city+(j-1),7] # p - intercept
    )
    
  }
  
  for(j in n_cities){
    
    city <- cities[j]
    
    index_lower <- 1 + ((j-1) * params)
    index_upper <- 1 + ((j-1) * params) + (params - 1)
    
    Y[index_lower:index_upper] <- c(
      fit_summary$summary[p0,1], 
      fit_summary$summary[p_wingspan,1],
      fit_summary$summary[p_feature_diversity,1],
      fit_summary$summary[p_ease_of_id,1],
      fit_summary$summary[mu_p_species_date,1],
      fit_summary$summary[mu_p_species_date_sq,1]
    )
    
    lower_95[index_lower:index_upper] <- c(
      fit_summary$summary[p0,4], 
      fit_summary$summary[p_wingspan,4],
      fit_summary$summary[p_feature_diversity,4],
      fit_summary$summary[p_ease_of_id,4],
      fit_summary$summary[mu_p_species_date,4],
      fit_summary$summary[mu_p_species_date_sq,4]
    )
    
    upper_95[index_lower:index_upper] <- c(
      fit_summary$summary[p0,8], 
      fit_summary$summary[p_wingspan,8],
      fit_summary$summary[p_feature_diversity,8],
      fit_summary$summary[p_ease_of_id,8],
      fit_summary$summary[mu_p_species_date,8],
      fit_summary$summary[mu_p_species_date_sq,8]
    )
    
    lower_50[index_lower:index_upper] <- c(
      fit_summary$summary[p0,5], 
      fit_summary$summary[p_wingspan,5],
      fit_summary$summary[p_feature_diversity,5],
      fit_summary$summary[p_ease_of_id,5],
      fit_summary$summary[mu_p_species_date,5],
      fit_summary$summary[mu_p_species_date_sq,5]
    )
    
    upper_50[index_lower:index_upper] <- c(
      fit_summary$summary[p0,7], 
      fit_summary$summary[p_wingspan,7],
      fit_summary$summary[p_feature_diversity,7],
      fit_summary$summary[p_ease_of_id,7],
      fit_summary$summary[mu_p_species_date,7],
      fit_summary$summary[mu_p_species_date_sq,7]
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
  filter(X == "1")

df_estimates <- rbind(df_estimates1, df_estimates2)

## --------------------------------------------------
## Draw caterpillar plot

s <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(p["intercept"]),
                             bquote(p["wingspan"]),
                             bquote(p["ft. diversity"]),
                             bquote(p["ease of ID"]),
                             bquote(p["peak phenology"]),
                             bquote(p["phenology decay"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-4, 3), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12)) +
   guides(color = guide_legend(title = "city")) +
   scale_color_manual(values=c("black", "#E69F00", "#D12F00", "#56B4E9", "#99A4E9", 
                                      "#1a5acd", "#E69F90", "#FFFF00")) + 
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

