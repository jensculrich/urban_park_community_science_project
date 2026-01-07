# plot effect sizes with and without including plant diversity 
# will determine how much if any of the effects work indirectly 
# on butterflies through increasing plant diversity

library(tidyverse)
library(cmdstanr)

# for simplicity, we will only compare the mean effects (not city specific) 

## get param estimates from m2.1
stan_out_m2.1 <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_jan2.rds")

# summarise all variables with default and additional summary measures
estimates1 <- as.data.frame(stan_out_m2.1$summary(
  variables = c(
    "mu_psi_park_size",
    "mu_psi_tree_cover",
    "mu_psi_plant_diversity"
    ),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates1) <- estimates1[, 1]

# clear space
rm(stan_out_m2.1)
gc()

## get param estimates from m2.2
stan_out_m2.2 <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.2_jan6.rds")

# summarise all variables with default and additional summary measures
estimates2 <- as.data.frame(stan_out_m2.2$summary(
  variables = c(
    "mu_psi_park_size",
    "mu_psi_tree_cover"
  ),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates2) <- estimates2[, 1]

# clear space
rm(stan_out_m2.2)
gc()

# handy for viewing column numbers
# this line of code won't work until you've actually read in a stan fit object
#View(cbind(1:nrow(estimates), estimates)) # View to see which row corresponds to the parameter of interest

#-------------------------------------------------------------------------------
# occupancy (psi)
  
  # ecological params
  # number of params to plot
  params <- 5  # park size, park size, tree cover, tree cover, plant diversity
  X <- c("park size", "park size", "tree cover", "tree cover", "plant diversity")
  Model <- c("m2.1", "m2.2", "m2.1", "m2.2", "m2.1")
  Y <- vector(length = params) # Y = mean estimate for a param of interest
  lower_95 <- vector(length = params)
  upper_95 <- vector(length = params)
  lower_50 <- vector(length = params)
  upper_50 <- vector(length = params)
  
  # fill in the means and BCIs
    Y[1:params] <- c(
      estimates1[1,2], # park size m2.1
      estimates2[1,2], # park size m2.2
      estimates1[2,2], # tree cover m2.1
      estimates2[2,2], # tree cover m2.2
      estimates1[3,2] # plant diversity m2.1
    )
    
    lower_95[1:params] <- c(
      estimates1[1,6], # park size m2.1
      estimates2[1,6], # park size m2.2
      estimates1[2,6], # tree cover m2.1
      estimates2[2,6], # tree cover m2.2
      estimates1[3,6] # plant diversity m2.1
    )
    
    upper_95[1:params] <- c(
      estimates1[1,7], # park size m2.1
      estimates2[1,7], # park size m2.2
      estimates1[2,7], # tree cover m2.1
      estimates2[2,7], # tree cover m2.2
      estimates1[3,7] # plant diversity m2.1
    )
    
    lower_50[1:params] <- c(
      estimates1[1,8], # park size m2.1
      estimates2[1,8], # park size m2.2
      estimates1[2,8], # tree cover m2.1
      estimates2[2,8], # tree cover m2.2
      estimates1[3,8] # plant diversity m2.1
    )
    
    upper_50[1:params] <- c(
      estimates1[1,9], # park size m2.1
      estimates2[1,9], # park size m2.2
      estimates1[2,9], # tree cover m2.1
      estimates2[2,9], # tree cover m2.2
      estimates1[3,9] # plant diversity m2.1
    )
  
  # now bind all of the param names, city names, and quantiles into a df for plotting
  df_estimates <- as.data.frame(cbind(X, Model, Y, lower_95, upper_95, lower_50, upper_50))
  
  df_estimates$X <- as.factor(df_estimates$X)
  df_estimates$Model <- as.factor(df_estimates$Model)
  df_estimates$Y <- as.numeric(df_estimates$Y)
  df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
  df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
  df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
  df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)
  
  
# plot params in alphabetical order
df_estimates$X <- factor(df_estimates$X, 
                         levels = c("park size", "tree cover", "plant diversity"))                                                     
# plot models in alphabetical order
df_estimates$Model <- factor(df_estimates$Model, 
                         levels = c("m2.1", "m2.2"))  

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
p <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = c("park size", "tree cover", "plant diversity"),
                    labels=c(bquote(psi["park size"]),
                             bquote(psi["tree cover"]),
                             bquote(psi["plant diversity"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1, 1.5)) +
   geom_hline(yintercept = 0, lty = "dashed") +
   ggtitle("Occupancy") +
   guides(color = guide_legend(title = "Model")) +
   scale_color_manual(values=c("gray30", "gray80")) + 
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=Model, colour=Model),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=Model, colour=Model),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=Model, colour=Model), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 

# plot the plot
p
