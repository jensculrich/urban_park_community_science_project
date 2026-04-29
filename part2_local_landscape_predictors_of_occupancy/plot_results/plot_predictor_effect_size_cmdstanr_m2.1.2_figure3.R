# plot effects of predictors on occurrence and detection

library(tidyverse)
library(cmdstanr)

## get param estimates from the region
stan_out <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_apr9.rds")

# summarise all variables with default and additional summary measures
estimates <- as.data.frame(stan_out$summary(
  variables = c(
    "psi_0", 
    "sigma_psi_species",
    "sigma_psi_city",
    "psi_wingspan",
    "psi_migratory",
    "mu_psi_park_size",
    "sigma_psi_park_size",
    "mu_psi_tree_cover",
    "sigma_psi_tree_cover",
    "mu_psi_plant_diversity",
    "sigma_psi_plant_diversity",
    "mu_psi_landscape_isolation",
    "sigma_psi_landscape_isolation",
    "mu_psi_landscape_grassherb",
    "sigma_psi_landscape_grassherb",
    "mu_psi_landscape_woody",
    "sigma_psi_landscape_woody",
    
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
    "sigma_p_species_date_sq"),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates) <- estimates[, 1]

#-------------------------------------------------------------------------------
# occupancy (psi)

  # ecological params
  # number of params to plot
  params_re <- 9 # intercept, local, landscape, species traits
  params_fixed <- 0 
  params <- params_re + params_fixed
  X <- rep(seq(1:params)) #  ecological params of interest
  Y <- vector(length = params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = params)
  upper_95 <- vector(length = params)
  lower_50 <- vector(length = params)
  upper_50 <- vector(length = params)
  
  # get indices for species random effects distributions for particular city
  # by indexing the row with the mean city estimate (usually I call this param mu_...)
  # and the first city random effect (psi_..._[1])
  # index of other cities will be the row of the first random effect plus some integer.
  psi_0 <- which( rownames(estimates)=="psi_0" )
  psi_parksize <- which( rownames(estimates)=="mu_psi_park_size" )
  psi_tree_cover <-  which( rownames(estimates)=="mu_psi_tree_cover" )
  psi_plant_diversity <- which( rownames(estimates)=="mu_psi_plant_diversity" )
  psi_isolation <- which( rownames(estimates)=="mu_psi_landscape_isolation" )
  psi_landscape_grassherb <- which( rownames(estimates)=="mu_psi_landscape_grassherb" )
  psi_landscape_woody <-  which( rownames(estimates)=="mu_psi_landscape_woody" )
  psi_wingspan <- which( rownames(estimates)=="psi_wingspan" )
  psi_migratory <- which( rownames(estimates)=="psi_migratory" )
  
    index_lower <- 1 
    index_upper <- params
    
    Y[index_lower:index_upper] <- c(
      estimates[psi_0,2],
      estimates[psi_parksize,2],
      estimates[psi_tree_cover,2],
      estimates[psi_plant_diversity,2],
      estimates[psi_isolation,2],
      estimates[psi_landscape_grassherb,2],
      estimates[psi_landscape_woody,2], 
      estimates[psi_wingspan,2],
      estimates[psi_migratory,2]
    )
    
    lower_95[index_lower:index_upper] <- c(
      estimates[psi_0,6], 
      estimates[psi_parksize,6],
      estimates[psi_tree_cover,6],
      estimates[psi_plant_diversity,6],
      estimates[psi_isolation,6],
      estimates[psi_landscape_grassherb,6],
      estimates[psi_landscape_woody,6], 
      estimates[psi_wingspan,6],
      estimates[psi_migratory,6]
    )
    
    upper_95[index_lower:index_upper] <- c(
      estimates[psi_0,7], 
      estimates[psi_parksize,7],
      estimates[psi_tree_cover,7],
      estimates[psi_plant_diversity,7],
      estimates[psi_isolation,7],
      estimates[psi_landscape_grassherb,7],
      estimates[psi_landscape_woody,7], 
      estimates[psi_wingspan,7],
      estimates[psi_migratory,7]
    )
    
    lower_50[index_lower:index_upper] <- c(
      estimates[psi_0,8], 
      estimates[psi_parksize,8],
      estimates[psi_tree_cover,8],
      estimates[psi_plant_diversity,8],
      estimates[psi_isolation,8],
      estimates[psi_landscape_grassherb,8],
      estimates[psi_landscape_woody,8], 
      estimates[psi_wingspan,8],
      estimates[psi_migratory,8]
    )
    
    upper_50[index_lower:index_upper] <- c(
      estimates[psi_0,9], 
      estimates[psi_parksize,9],
      estimates[psi_tree_cover,9],
      estimates[psi_plant_diversity,9],
      estimates[psi_isolation,9],
      estimates[psi_landscape_grassherb,9],
      estimates[psi_landscape_woody,9], 
      estimates[psi_wingspan,9],
      estimates[psi_migratory,9]
    )
    
  # now bind all of the param names, city names, and quantiles into a df for plotting
  category <- c("Intercept", "Local", "Local", "Local", "Landscape", "Landscape", "Landscape", "Species Trait", "Species Trait")
  df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50, category))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  df_estimates$category <- as.factor(df_estimates$category)
  
  df_estimates1 <- df_estimates[1,]
  df_estimates2 <- df_estimates[2:params,]
  
  df_estimates2$category <- fct_relevel(df_estimates2$category, "Species Trait", "Landscape", "Local")
  
  
## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
p2 <- ggplot(df_estimates2) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(from=2,to=(params)),
                    labels=c(bquote(psi["park size"]),
                             bquote(psi["tree cover"]),
                             bquote(psi["plant diversity"]),
                             bquote(psi["isolation"]),
                             bquote(psi["landsc. herb."]),
                             bquote(psi["landsc. woody"]),
                             bquote(psi["wingspan"]),
                             bquote(psi["migratory"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-0.5, 1.75), 
                      breaks = scales::pretty_breaks()) +
                      #breaks = c(-0.5, 0, 1, 2, 4, 6, 8)) +
   scale_color_manual(name = "", values=c("goldenrod2", "orchid3", "dodgerblue3")) + 
   geom_hline(yintercept = 0, lty = "dashed") +
   ggtitle("") +
   theme(legend.position = c(0.975, 0.025), # x=1 (right), y=0 (bottom)
         legend.justification = c(1, 0), # Justify the bottom-right corner of the legend box to these coordinates
         legend.text = element_text(size=16),
         plot.title = element_text(size = 18, face = "bold"),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip() 

# add estimates
p2 <- p2 +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=category, colour=category),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=category, colour=category),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=category, colour=category), 
             size = 5, alpha = 0.8) 

# plot the plot
p2


# layout the plot
p1 <- ggplot(df_estimates1) +
  theme_bw() +
  scale_x_discrete(name="", breaks = seq(from=1,to=1),
                   labels=c(bquote(psi["intercept"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-2, 0), 
                     breaks = scales::pretty_breaks()) +
  #breaks = c(-0.5, 0, 1, 2, 4, 6, 8)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  #scale_colour_manual(colour="black") +
  ggtitle("") +
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip() 

# add estimates
p1 <- p1 +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y), 
             size = 5, alpha = 0.8) 

# plot the plot
p1

panelled_plot <- cowplot::plot_grid(p2, p1, ncol=1, rel_heights = c(3, 1), align = "v", axis = "lr")
#saveRDS(panelled_plot, "./part2_local_landscape_predictors_of_occupancy/plot_results/figure3_param_estimates.rds")

