library(tidyverse)

#-------------------------------------------------------------------------------
# now compare m2.1 and m2.2 means after removing plant diversity (following the DAG logic)

# for simplicity, we will only compare the mean effects (not city specific) 

## get param estimates from the region
stan_out <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_may25.rds")

# summarise all variables with default and additional summary measures
estimates1 <- as.data.frame(stan_out$summary(
  variables = c(
    "mu_psi_park_size",
    "mu_psi_tree_cover",
    "mu_psi_plant_diversity",
    "mu_psi_landscape_isolation",
    "mu_psi_landscape_grassherb",
    "mu_psi_landscape_woody"
  ),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

rownames(estimates1) <- estimates1[, 1]

# clear space
rm(stan_out)
gc()

## get param estimates from m2.2
stan_out_m2.2 <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.2_may27.rds")

# summarise all variables with default and additional summary measures
estimates2 <- as.data.frame(stan_out_m2.2$summary(
  variables = c(
    "mu_psi_park_size",
    "mu_psi_tree_cover",
    "mu_psi_landscape_isolation"
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
params <- 9
# park size, park size, tree cover, tree cover, plant diversity, 
# connectivity, connectivity, land grassherb, land woody, 
X <- c("park size", "park size", "tree cover", "tree cover", "plant diversity",
       "connectivity", "connectivity", "landscape herb.", "landscape woody")
Model <- c("m2.1", "m2.2", "m2.1", "m2.2", "m2.1",
           "m2.1", "m2.2", "m2.1", "m2.1")
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
  estimates1[3,2], # plant diversity m2.1
  estimates1[4,2], # connectivity m2.1
  estimates2[3,2], # connectivity m2.2
  estimates1[5,2], # land herb m2.1
  estimates1[6,2] # land woody m2.1
)

lower_95[1:params] <- c(
  estimates1[1,6], # park size m2.1
  estimates2[1,6], # park size m2.2
  estimates1[2,6], # tree cover m2.1
  estimates2[2,6], # tree cover m2.2
  estimates1[3,6], # plant diversity m2.1
  estimates1[4,6], # connectivity m2.1
  estimates2[3,6], # connectivity m2.2
  estimates1[5,6], # land herb m2.1
  estimates1[6,6] # land woody m2.1
)

upper_95[1:params] <- c(
  estimates1[1,7], # park size m2.1
  estimates2[1,7], # park size m2.2
  estimates1[2,7], # tree cover m2.1
  estimates2[2,7], # tree cover m2.2
  estimates1[3,7], # plant diversity m2.1
  estimates1[4,7], # connectivity m2.1
  estimates2[3,7], # connectivity m2.2
  estimates1[5,7], # land herb m2.1
  estimates1[6,7] # land woody m2.1
)

lower_50[1:params] <- c(
  estimates1[1,8], # park size m2.1
  estimates2[1,8], # park size m2.2
  estimates1[2,8], # tree cover m2.1
  estimates2[2,8], # tree cover m2.2
  estimates1[3,8], # plant diversity m2.1
  estimates1[4,8], # connectivity m2.1
  estimates2[3,8], # connectivity m2.2
  estimates1[5,8], # land herb m2.1
  estimates1[6,8] # land woody m2.1
)

upper_50[1:params] <- c(
  estimates1[1,9], # park size m2.1
  estimates2[1,9], # park size m2.2
  estimates1[2,9], # tree cover m2.1
  estimates2[2,9], # tree cover m2.2
  estimates1[3,9], # plant diversity m2.1
  estimates1[4,9], # connectivity m2.1
  estimates2[3,9], # connectivity m2.2
  estimates1[5,9], # land herb m2.1
  estimates1[6,9] # land woody m2.1
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
                         levels = c("park size", "tree cover", "plant diversity",
                                    "connectivity", "landscape herb.", "landscape woody"))                                                     
# plot models in alphabetical order
df_estimates$Model <- factor(df_estimates$Model, 
                             levels = c("m2.1", "m2.2"))  

## --------------------------------------------------
## Draw caterpillar plot

# layout the plot
r <- ggplot(df_estimates) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c("park size", "tree cover", "plant diversity",
                                       "connectivity", "landscape herb.", "landscape woody"),
                   labels=c(bquote(psi["park size"]),
                            bquote(psi["tree cover"]),
                            bquote(psi["plant diversity"]),
                            bquote(psi["connectivity"]),
                            bquote(psi["landscape herb."]),
                            bquote(psi["landscape woody"])
                   )) +
  scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                     limits = c(-0.5, 1), breaks = c(-0.5, 0, 0.5, 1, 1.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  ggtitle("") +
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
r <- r +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=Model, colour=Model),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=Model, colour=Model),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=Model, colour=Model), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 

# plot the plot
r

#-------------------------------------------------------------------------------
# plot the 3 panels on a 2x2 grid

cowplot::plot_grid(p, r, ncol = 2, rel_widths = c(1, 0.67),
                   labels = c('a)', 'b)'),
                   label_size = 20)
