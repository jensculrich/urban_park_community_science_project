# plot effects of predictors on occurrence and detection

library(tidyverse)
library(cmdstanr)

ilogit <- function(x) exp(x)/(1+exp(x))

## get param estimates from the region
stan_out <- readRDS(
  "./analyses/part1_urban_butterfly_community_dynamics/model_outputs/stan_out_jun4.rds")

# summarise all variables with default and additional summary measures
estimates <- as.data.frame(stan_out$summary(
  variables = c(
    "psi1_0", 
    "mu_psi1_park_size",
    "mu_psi1_isolation",
    "psi1_wingspan",
    "psi1_migratory",
    
    "gamma0", 
    "mu_gamma_park_size",
    "mu_gamma_isolation",
    "gamma_wingspan",
    "gamma_migratory",
    
    "phi0", 
    "mu_phi_park_size",
    "mu_phi_isolation",
    "phi_wingspan",
    "phi_migratory",
    
    "p0", 
    "p_city_detections",
    "p_wingspan",
    "p_migratory",
    "p_feature_diversity",
    "p_ease_of_id",
    "delta0",
    "epsilon0"
),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates) <- estimates[, 1]

# handy for viewing column numbers
# this line of code won't work until you've actually read in a stan fit object
#View(cbind(1:nrow(estimates), estimates)) # View to see which row corresponds to the parameter of interest

rm(stan_out)
gc() 

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

  # ecological params
  # number of params to plot
  # plot the mean etsimates across all cites, 
  # which I'll just call params_fixed based on my earlier plotting conventions
  params_fixed <- 5 #  # intercept, park size, isolation,  wingspan, migratory
  params_re <- 0  # no city effects to display in this summary plot
  params <- params_re + params_fixed
  X <- rep(seq(1:params)) #  ecological params of interest
  Y <- vector(length = params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = params)
  upper_95 <- vector(length = params)
  lower_50 <- vector(length = params)
  upper_50 <- vector(length = params)
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (psi1_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  psi1_0 <- which( rownames(estimates)=="psi1_0" )
  psi1_parksize <- which( rownames(estimates)=="mu_psi1_park_size" )
  psi1_isolation <- which( rownames(estimates)=="mu_psi1_isolation" )
  psi1_wingspan <- which( rownames(estimates)=="psi1_wingspan" )
  psi1_migratory <- which( rownames(estimates)=="psi1_migratory" )
  
    # now fill in the param values from the model estimates summary table
    
    index_lower <- 1  
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
    
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

  df_estimates <- df_estimates[1:3,]

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
p <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(psi[1]["intercept"]),
                             bquote(psi[1]["park size"]),
                             bquote(psi[1]["connectivity"])#,
                             #bquote(psi["wingspan"]),
                             #bquote(psi["migratory"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-2, 2), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y),
             size = 5, alpha = 0.8) 

# plot the plot
p

#-------------------------------------------------------------------------------
# everything below replicates the above but for the other processes (colonization, persistence, detection)
# I may not have commented everything as well yet.

#-------------------------------------------------------------------------------
# colonization (gamma)

# ecological params
# number of params to plot
# plot the mean etsimates across all cites, 
# which I'll just call params_fixed based on my earlier plotting conventions
params_fixed <- 5 #  # intercept, park size, isolation,  wingspan, migratory
params_re <- 0  # no city effects to display in this summary plot
params <- params_re + params_fixed
X <- rep(seq(1:params)) #  ecological params of interest
Y <- vector(length = params) # Y = mean estimate for a param of interest
lower_95 <- vector(length = params)
upper_95 <- vector(length = params)
lower_50 <- vector(length = params)
upper_50 <- vector(length = params)

    # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (gamma_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  gamma_0 <- which( rownames(estimates)=="gamma0" )
  gamma_parksize <- which( rownames(estimates)=="mu_gamma_park_size" )
  gamma_isolation <- which( rownames(estimates)=="mu_gamma_isolation" )
  gamma_wingspan <- which( rownames(estimates)=="gamma_wingspan" )
  gamma_migratory <- which( rownames(estimates)=="gamma_migratory" )
  
  
    # now fill in the param values from the model estimates summary table
  
    index_lower <- 1  
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
    
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

  df_estimates <- df_estimates[1:3,]

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
q <- ggplot(df_estimates) +
  theme_bw() +
  scale_x_discrete(name="", breaks = seq(1:params),
                   labels=c(bquote(gamma["intercept"]),
                            bquote(gamma["park size"]),
                            bquote(gamma["connectivity"])#,
                            #bquote(gamma["wingspan"]),
                            #bquote(gamma["migratory"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-6.5, 3), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y),
             size = 5, alpha = 0.8) 

# plot the plot
q


#-------------------------------------------------------------------------------
# persistence (phi)

#-------------------------------------------------------------------------------
# colonization (gamma)

# ecological params
# number of params to plot
# plot the mean etsimates across all cites, 
# which I'll just call params_fixed based on my earlier plotting conventions
params_fixed <- 5 #  # intercept, park size, isolation,  wingspan, migratory
params_re <- 0  # no city effects to display in this summary plot
params <- params_re + params_fixed
X <- rep(seq(1:params)) #  ecological params of interest
Y <- vector(length = params) # Y = mean estimate for a param of interest
lower_95 <- vector(length = params)
upper_95 <- vector(length = params)
lower_50 <- vector(length = params)
upper_50 <- vector(length = params)

  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (phi_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  phi_0 <- which( rownames(estimates)=="phi0" )
  phi_parksize <- which( rownames(estimates)=="mu_phi_park_size" )
  phi_isolation <- which( rownames(estimates)=="mu_phi_isolation" )
  phi_wingspan <- which( rownames(estimates)=="phi_wingspan" )
  phi_migratory <- which( rownames(estimates)=="phi_migratory" )
  
  
    index_lower <- 1
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
    
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

df_estimates <- df_estimates[1:3,]

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
r <- ggplot(df_estimates) +
  theme_bw() +
  scale_x_discrete(name="", breaks = seq(1:params),
                   labels=c(bquote(phi["intercept"]),
                            bquote(phi["park size"]),
                            bquote(phi["connectivity"])#,
                            #bquote(phi["wingspan"]),
                            #bquote(phi["migratory"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-3, 6.5), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y),
             size = 5, alpha = 0.8) 

# plot the plot
r

#-------------------------------------------------------------------------------
# detection (p)


# ecological params
# number of params to plot
# plot the mean etsimates across all cites, 
# which I'll just call params_fixed based on my earlier plotting conventions
params_fixed <- 7 #  # intercept, park size, isolation,  wingspan, migratory
params_re <- 0  # no city effects to display in this summary plot
params <- params_re + params_fixed
X <- rep(seq(1:params)) #  ecological params of interest
Y <- vector(length = params) # Y = mean estimate for a param of interest
lower_95 <- vector(length = params)
upper_95 <- vector(length = params)
lower_50 <- vector(length = params)
upper_50 <- vector(length = params)

  # get indices for species random effects distributions for particular city
  p0 <- which( rownames(estimates)=="p0" )
  p_wingspan <- which( rownames(estimates)=="p_wingspan" )
  p_migratory <- which( rownames(estimates)=="p_migratory" )
  p_feature_diversity <- which( rownames(estimates)=="p_feature_diversity" )
  p_ease_of_id <- which( rownames(estimates)=="p_ease_of_id" )
  delta0 <- which( rownames(estimates)=="delta0" )
  epsilon0 <- which( rownames(estimates)=="epsilon0" )
  
    
    index_lower <- 1 
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
    

  df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)


df_estimates

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
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
                     limits = c(-4, 2), breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +
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
s <- s +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y),
             size = 5, alpha = 0.8) 

# plot the plot
s

#-------------------------------------------------------------------------------
# plot the 4 panels on a 2x2 grid

#cowplot::plot_grid(p, q, r, s, ncol = 2)

cowplot::plot_grid(p, q, r, labels = c("a)", "b)", "c)"),
                   ncol = 3, label_size = 18)

#ggplot(site_data) +
#  ggridges::geom_density_ridges(aes(x=log_total_green_space_area, y=city, fill = city), alpha = 0.3)

#ggplot(site_data) +
#  ggridges::geom_density_ridges(aes(x=log_isolation_scaled_across_all_cities, y=city, fill = city), alpha = 0.3)
