# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

# params to plot
param_names <- c("rho",
                 "sigma_species[1]",
                 "sigma_species[2]",
                 
                 "psi1_0", 
                 #"sigma_psi1_species",
                 "psi1_wingspan",
                 "psi1_park_size",
                 "psi1_connectivity",
                 "psi1_plant_genera",
                 "psi1_tree_cover",
                 
                 "gamma0", 
                 "sigma_gamma_species",
                 "gamma_wingspan",
                 "gamma_park_size",
                 "gamma_connectivity",
                 "gamma_plant_genera",
                 "gamma_tree_cover",
                 #"gamma_wingspan_connectivity",
                 
                 "phi0", 
                 "sigma_phi_species",
                 "phi_wingspan",
                 "phi_park_size",
                 "phi_connectivity",
                 "phi_plant_genera",
                 "phi_tree_cover",
                 
                 "p0", 
                 #"sigma_p_species",
                 "p_wingspan",
                 "p_feature_diversity",
                 "p_ease_of_id",
                 "mu_p_species_date",
                 "sigma_p_species_date",
                 "mu_p_species_date_sq",
                 "sigma_p_species_date_sq")

# create empty list of length n_cities
# each element of the list holds the posterior distributions for 
# all parameters for each individual city 
df <- data.frame()

# read and summarize data across loop
for(city_number in 1:n_cities){
  
  city <- city_names[city_number]
  
  temp <- cbind(as.data.frame(
    readRDS(paste0(
      "./model_outputs/stan_out_", city, "_2km_connectivity_family_100buffers.rds"))
  ), city)
  
  df <- rbind(df, as.data.frame(temp[,1:length(param_names)]))
  
}

df <- as.matrix(df)

mean = vector(length=ncol(df))
lower_50 = vector(length=ncol(df))
upper_50 = vector(length=ncol(df))
lower_95 = vector(length=ncol(df))
upper_95 = vector(length=ncol(df))

for(i in 1:ncol(df)){
    
    quants = as.vector(quantile(df[,i], probs = c(0.05, 0.25, 0.50, 0.75, 0.95)))
    
    mean[i] = quants[3]
    lower_50[i] = quants[2]
    upper_50[i] = quants[4]
    lower_95[i] = quants[1]
    upper_95[i] = quants[5]

}

X <- seq(1:length(param_names)) #  ecological params of interest
df_estimates <- as.data.frame(cbind(X, param_names, mean, lower_50, upper_50, lower_95, upper_95))

df_estimates$X <- as.factor(df_estimates$X)
df_estimates$Y <- as.numeric(df_estimates$mean)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)


## --------------------------------------------------
## Draw caterpillar plot

psi1_estimates <- df_estimates[4:9,]

(p <- ggplot(psi1_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(from=4,to=9),
                    labels=c(bquote(psi["intercept"]),
                             bquote(psi["wingspan"]),
                             bquote(psi["park size"]),
                             bquote(psi["isolation"]),
                             bquote(psi["plant richness"]),
                             bquote(psi["tree cover"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3, 4)) +
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
)

p <- p +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p



#-------------------------------------------------------------------------------
# initial occurrence (psi1)

# ecological params
# number of params to plot
params <- 12
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
       fit_summary$summary[4,1], # psi1 - intercept
       fit_summary$summary[5,1], # psi1 - wingspan
       fit_summary$summary[6,1], # psi1 - park size
       fit_summary$summary[7,1], # psi1 - connectivity
       fit_summary$summary[8,1], # psi1 - plant richness
       fit_summary$summary[9,1], # psi1 - tree cover
       # NYC level inference
       fit_summary2$summary[4,1], # psi1 - intercept
       fit_summary2$summary[5,1], # psi1 - wingspan
       fit_summary2$summary[6,1], # psi1 - park size
       fit_summary2$summary[7,1], # psi1 - connectivity
       fit_summary2$summary[8,1], # psi1 - plant richness
       fit_summary2$summary[9,1] # psi1 - tree cover
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[4,4], # psi1 - intercept
  fit_summary$summary[5,4], # psi1 - wingspan
  fit_summary$summary[6,4], # psi1 - park size
  fit_summary$summary[7,4], # psi1 - connectivity
  fit_summary$summary[8,4], # psi1 - plant richness
  fit_summary$summary[9,4], # psi1 - tree cover
  # NYC level inference
  fit_summary2$summary[4,4], # psi1 - intercept
  fit_summary2$summary[5,4], # psi1 - wingspan
  fit_summary2$summary[6,4], # psi1 - park size
  fit_summary2$summary[7,4], # psi1 - connectivity
  fit_summary2$summary[8,4], # psi1 - plant richness
  fit_summary2$summary[9,4] # psi1 - tree cover
)

upper_95 <- c(# LA level inference
  fit_summary$summary[4,8], # psi1 - intercept
  fit_summary$summary[5,8], # psi1 - wingspan
  fit_summary$summary[6,8], # psi1 - park size
  fit_summary$summary[7,8], # psi1 - connectivity
  fit_summary$summary[8,8], # psi1 - plant richness
  fit_summary$summary[9,8], # psi1 - tree cover
  # NYC level inference
  fit_summary2$summary[4,8], # psi1 - intercept
  fit_summary2$summary[5,8], # psi1 - wingspan
  fit_summary2$summary[6,8], # psi1 - park size
  fit_summary2$summary[7,8], # psi1 - connectivity
  fit_summary2$summary[8,8], # psi1 - plant richness
  fit_summary2$summary[9,8] # psi1 - tree cover
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[4,5], # psi1 - intercept
  fit_summary$summary[5,5], # psi1 - wingspan
  fit_summary$summary[6,5], # psi1 - park size
  fit_summary$summary[7,5], # psi1 - connectivity
  fit_summary$summary[8,5], # psi1 - plant richness
  fit_summary$summary[9,5], # psi1 - tree cover
  # NYC level inference
  fit_summary2$summary[4,5], # psi1 - intercept
  fit_summary2$summary[5,5], # psi1 - wingspan
  fit_summary2$summary[6,5], # psi1 - park size
  fit_summary2$summary[7,5], # psi1 - connectivity
  fit_summary2$summary[8,5], # psi1 - plant richness
  fit_summary2$summary[9,5] # psi1 - tree cover
)

upper_50 <- c(# LA level inference
  fit_summary$summary[4,7], # psi1 - intercept
  fit_summary$summary[5,7], # psi1 - wingspan
  fit_summary$summary[6,7], # psi1 - park size
  fit_summary$summary[7,7], # psi1 - connectivity
  fit_summary$summary[8,7], # psi1 - plant richness
  fit_summary$summary[9,7], # psi1 - tree cover
  # NYC level inference
  fit_summary2$summary[4,7], # psi1 - intercept
  fit_summary2$summary[5,7], # psi1 - wingspan
  fit_summary2$summary[6,7], # psi1 - park size
  fit_summary2$summary[7,7], # psi1 - connectivity
  fit_summary2$summary[8,7], # psi1 - plant richness
  fit_summary2$summary[9,7] # psi1 - tree cover
)

df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))

df_estimates$X <- as.factor(df_estimates$X)
df_estimates$city_name <- as.factor(df_estimates$city_name)
df_estimates$Y <- as.numeric(df_estimates$Y)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

## --------------------------------------------------
## Draw caterpillar plot

(p <- ggplot(df_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = seq(1:(0.5*params)),
                     labels=c(bquote(psi["intercept"]),
                              bquote(psi["wingspan"]),
                              bquote(psi["park size"]),
                              bquote(psi["isolation"]),
                              bquote(psi["plant richness"]),
                              bquote(psi["tree cover"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-3, 4)) +
    guides(color = guide_legend(title = "city")) +
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
)

p <- p +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, group=city_name, colour=city_name), 
                 position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
p


#-------------------------------------------------------------------------------
# colonization (gamma)

# ecological params
# number of params to plot
params <- 12
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[10,1], # gamma - intercept
  fit_summary$summary[12,1], # gamma - wingspan
  fit_summary$summary[13,1], # gamma - park size
  fit_summary$summary[14,1], # gamma - connectivity
  fit_summary$summary[15,1], # gamma - plant richness
  fit_summary$summary[16,1], # gamma - tree cover
  # NYC level inference
  fit_summary2$summary[10,1], # gamma - intercept
  fit_summary2$summary[12,1], # gamma - wingspan
  fit_summary2$summary[13,1], # gamma - park size
  fit_summary2$summary[14,1], # gamma - connectivity
  fit_summary2$summary[15,1], # gamma - plant richness
  fit_summary2$summary[16,1] # gamma - tree cover
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[10,4], # gamma - intercept
  fit_summary$summary[12,4], # gamma - wingspan
  fit_summary$summary[13,4], # gamma - park size
  fit_summary$summary[14,4], # gamma - connectivity
  fit_summary$summary[15,4], # gamma - plant richness
  fit_summary$summary[16,4], # gamma - tree cover
  # NYC level inference
  fit_summary2$summary[10,4], # gamma - intercept
  fit_summary2$summary[12,4], # gamma - wingspan
  fit_summary2$summary[13,4], # gamma - park size
  fit_summary2$summary[14,4], # gamma - connectivity
  fit_summary2$summary[15,4], # gamma - plant richness
  fit_summary2$summary[16,4] # gamma - tree cover
)

upper_95 <- c(# LA level inference
  fit_summary$summary[10,8], # gamma - intercept
  fit_summary$summary[12,8], # gamma - wingspan
  fit_summary$summary[13,8], # gamma - park size
  fit_summary$summary[14,8], # gamma - connectivity
  fit_summary$summary[15,8], # gamma - plant richness
  fit_summary$summary[16,8], # gamma - tree cover
  # NYC level inference
  fit_summary2$summary[10,8], # gamma - intercept
  fit_summary2$summary[12,8], # gamma - wingspan
  fit_summary2$summary[13,8], # gamma - park size
  fit_summary2$summary[14,8], # gamma - connectivity
  fit_summary2$summary[15,8], # gamma - plant richness
  fit_summary2$summary[16,8] # gamma - tree cover
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[10,5], # gamma - intercept
  fit_summary$summary[12,5], # gamma - wingspan
  fit_summary$summary[13,5], # gamma - park size
  fit_summary$summary[14,5], # gamma - connectivity
  fit_summary$summary[15,5], # gamma - plant richness
  fit_summary$summary[16,5], # gamma - tree cover
  # NYC level inference
  fit_summary2$summary[10,5], # gamma - intercept
  fit_summary2$summary[12,5], # gamma - wingspan
  fit_summary2$summary[13,5], # gamma - park size
  fit_summary2$summary[14,5], # gamma - connectivity
  fit_summary2$summary[15,5], # gamma - plant richness
  fit_summary2$summary[16,5] # gamma - tree cover
)

upper_50 <- c(# LA level inference
  fit_summary$summary[10,7], # gamma - intercept
  fit_summary$summary[12,7], # gamma - wingspan
  fit_summary$summary[13,7], # gamma - park size
  fit_summary$summary[14,7], # gamma - connectivity
  fit_summary$summary[15,7], # gamma - plant richness
  fit_summary$summary[16,7], # gamma - tree cover
  # NYC level inference
  fit_summary2$summary[10,7], # gamma - intercept
  fit_summary2$summary[12,7], # gamma - wingspan
  fit_summary2$summary[13,7], # gamma - park size
  fit_summary2$summary[14,7], # gamma - connectivity
  fit_summary2$summary[15,7], # gamma - plant richness
  fit_summary2$summary[16,7] # gamma - tree cover
)

df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))

df_estimates$X <- as.factor(df_estimates$X)
df_estimates$city_name <- as.factor(df_estimates$city_name)
df_estimates$Y <- as.numeric(df_estimates$Y)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

## --------------------------------------------------
## Draw caterpillar plot

(q <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:(0.5*params)),
                    labels=c(bquote(gamma["intercept"]),
                             bquote(gamma["wingspan"]),
                             bquote(gamma["park size"]),
                             bquote(gamma["isolation"]),
                             bquote(gamma["plant richness"]),
                             bquote(gamma["tree cover"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-5, 4), breaks = c(-4, -2, 0, 2)) +
   guides(color = guide_legend(title = "city")) +
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
)

q <- q +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
q


#-------------------------------------------------------------------------------
# persistence (phi)

# ecological params
# number of params to plot
params <- 12
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[17,1], # phi - intercept
  fit_summary$summary[19,1], # phi - wingspan
  fit_summary$summary[20,1], # phi - park size
  fit_summary$summary[21,1], # phi - connectivity
  fit_summary$summary[22,1], # phi - plant richness
  fit_summary$summary[23,1], # phi - tree cover
  # NYC level inference
  fit_summary2$summary[17,1], # phi - intercept
  fit_summary2$summary[19,1], # phi - wingspan
  fit_summary2$summary[20,1], # phi - park size
  fit_summary2$summary[21,1], # phi - connectivity
  fit_summary2$summary[22,1], # phi - plant richness
  fit_summary2$summary[23,1] # phi - tree cover
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[17,4], # phi - intercept
  fit_summary$summary[19,4], # phi - wingspan
  fit_summary$summary[20,4], # phi - park size
  fit_summary$summary[21,4], # phi - connectivity
  fit_summary$summary[22,4], # phi - plant richness
  fit_summary$summary[23,4], # phi - tree cover
  # NYC level inference
  fit_summary2$summary[17,4], # phi - intercept
  fit_summary2$summary[19,4], # phi - wingspan
  fit_summary2$summary[20,4], # phi - park size
  fit_summary2$summary[21,4], # phi - connectivity
  fit_summary2$summary[22,4], # phi - plant richness
  fit_summary2$summary[23,4] # phi - tree cover
)

upper_95 <- c(# LA level inference
  fit_summary$summary[17,8], # phi - intercept
  fit_summary$summary[19,8], # phi - wingspan
  fit_summary$summary[20,8], # phi - park size
  fit_summary$summary[21,8], # phi - connectivity
  fit_summary$summary[22,8], # phi - plant richness
  fit_summary$summary[23,8], # phi - tree cover
  # NYC level inference
  fit_summary2$summary[17,8], # phi - intercept
  fit_summary2$summary[19,8], # phi - wingspan
  fit_summary2$summary[20,8], # phi - park size
  fit_summary2$summary[21,8], # phi - connectivity
  fit_summary2$summary[22,8], # phi - plant richness
  fit_summary2$summary[23,8] # phi - tree cover
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[17,5], # phi - intercept
  fit_summary$summary[19,5], # phi - wingspan
  fit_summary$summary[20,5], # phi - park size
  fit_summary$summary[21,5], # phi - connectivity
  fit_summary$summary[22,5], # phi - plant richness
  fit_summary$summary[23,5], # phi - tree cover
  # NYC level inference
  fit_summary2$summary[17,5], # phi - intercept
  fit_summary2$summary[19,5], # phi - wingspan
  fit_summary2$summary[20,5], # phi - park size
  fit_summary2$summary[21,5], # phi - connectivity
  fit_summary2$summary[22,5], # phi - plant richness
  fit_summary2$summary[23,5] # phi - tree cover
)

upper_50 <- c(# LA level inference
  fit_summary$summary[17,7], # phi - intercept
  fit_summary$summary[19,7], # phi - wingspan
  fit_summary$summary[20,7], # phi - park size
  fit_summary$summary[21,7], # phi - connectivity
  fit_summary$summary[22,7], # phi - plant richness
  fit_summary$summary[23,7], # phi - tree cover
  # NYC level inference
  fit_summary2$summary[17,7], # phi - intercept
  fit_summary2$summary[19,7], # phi - wingspan
  fit_summary2$summary[20,7], # phi - park size
  fit_summary2$summary[21,7], # phi - connectivity
  fit_summary2$summary[22,7], # phi - plant richness
  fit_summary2$summary[23,7] # phi - tree cover
)

df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))

df_estimates$X <- as.factor(df_estimates$X)
df_estimates$city_name <- as.factor(df_estimates$city_name)
df_estimates$Y <- as.numeric(df_estimates$Y)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

## --------------------------------------------------
## Draw caterpillar plot

(r <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:(0.5*params)),
                    labels=c(bquote(phi["intercept"]),
                             bquote(phi["wingspan"]),
                             bquote(phi["park size"]),
                             bquote(phi["isolation"]),
                             bquote(phi["plant richness"]),
                             bquote(phi["tree cover"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3, 4)) +
   guides(color = guide_legend(title = "city")) +
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
)

r <- r +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
r

#-------------------------------------------------------------------------------
# detection (p)

# ecological params
# number of params to plot
params <- 14
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[24,1], # p - intercept
  fit_summary$summary[25,1], # p - wingspan
  fit_summary$summary[26,1], # p - feature diversity
  fit_summary$summary[27,1], # p - ease of id
  fit_summary$summary[28,1], # p - peak phenology
  fit_summary$summary[30,1], # p - phenology decline,
  fit_summary$summary[1,1], # p - detection-occurrence correlation
  
  # NYC level inference
  fit_summary2$summary[24,1], # p - intercept
  fit_summary2$summary[25,1], # p - wingspan
  fit_summary2$summary[26,1], # p - feature diversity
  fit_summary2$summary[27,1], # p - ease of id
  fit_summary2$summary[28,1], # p - peak phenology
  fit_summary2$summary[30,1], # p - phenology decline
  fit_summary2$summary[1,1] # p - detection-occurrence correlation
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[24,4], # p - intercept
  fit_summary$summary[25,4], # p - wingspan
  fit_summary$summary[26,4], # p - feature diversity
  fit_summary$summary[27,4], # p - ease of id
  fit_summary$summary[28,4], # p - peak phenology
  fit_summary$summary[30,4], # p - phenology decline
  fit_summary$summary[1,4], # p - detection-occurrence correlation
  # NYC level inference
  fit_summary2$summary[24,4], # p - intercept
  fit_summary2$summary[25,4], # p - wingspan
  fit_summary2$summary[26,4], # p - feature diversity
  fit_summary2$summary[27,4], # p - ease of id
  fit_summary2$summary[28,4], # p - peak phenology
  fit_summary2$summary[30,4], # p - phenology decline
  fit_summary2$summary[1,4] # p - detection-occurrence correlation
)

upper_95 <- c(# LA level inference
  fit_summary$summary[24,8], # p - intercept
  fit_summary$summary[25,8], # p - wingspan
  fit_summary$summary[26,8], # p - feature diversity
  fit_summary$summary[27,8], # p - ease of id
  fit_summary$summary[28,8], # p - peak phenology
  fit_summary$summary[30,8], # p - phenology decline
  fit_summary$summary[1,8], # p - detection-occurrence correlation
  # NYC level inference
  fit_summary2$summary[24,8], # p - intercept
  fit_summary2$summary[25,8], # p - wingspan
  fit_summary2$summary[26,8], # p - feature diversity
  fit_summary2$summary[27,8], # p - ease of id
  fit_summary2$summary[28,8], # p - peak phenology
  fit_summary2$summary[30,8], # p - phenology decline
  fit_summary2$summary[1,8] # p - detection-occurrence correlation
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[24,5], # p - intercept
  fit_summary$summary[25,5], # p - wingspan
  fit_summary$summary[26,5], # p - feature diversity
  fit_summary$summary[27,5], # p - ease of id
  fit_summary$summary[28,5], # p - peak phenology
  fit_summary$summary[30,5], # p - phenology decline
  fit_summary$summary[1,5], # p - detection-occurrence correlation
  # NYC level inference
  fit_summary2$summary[24,5], # p - intercept
  fit_summary2$summary[25,5], # p - wingspan
  fit_summary2$summary[26,5], # p - feature diversity
  fit_summary2$summary[27,5], # p - ease of id
  fit_summary2$summary[28,5], # p - peak phenology
  fit_summary2$summary[30,5], # p - phenology decline
  fit_summary2$summary[1,5] # p - detection-occurrence correlation
)

upper_50 <- c(# LA level inference
  fit_summary$summary[24,7], # p - intercept
  fit_summary$summary[25,7], # p - wingspan
  fit_summary$summary[26,7], # p - feature diversity
  fit_summary$summary[27,7], # p - ease of id
  fit_summary$summary[28,7], # p - peak phenology
  fit_summary$summary[30,7], # p - phenology decline
  fit_summary$summary[1,7], # p - detection-occurrence correlation
  # NYC level inference
  fit_summary2$summary[24,7], # p - intercept
  fit_summary2$summary[25,7], # p - wingspan
  fit_summary2$summary[26,7], # p - feature diversity
  fit_summary2$summary[27,7], # p - ease of id
  fit_summary2$summary[28,7], # p - peak phenology
  fit_summary2$summary[30,7], # p - phenology decline
  fit_summary2$summary[1,7] # p - detection-occurrence correlation
)

df_estimates <- as.data.frame(cbind(X, city_name, Y, lower_95, upper_95, lower_50, upper_50))

df_estimates$X <- as.factor(df_estimates$X)
df_estimates$city_name <- as.factor(df_estimates$city_name)
df_estimates$Y <- as.numeric(df_estimates$Y)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)

## --------------------------------------------------
## Draw caterpillar plot

(s <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:(0.5*params)),
                    labels=c(bquote(p["intercept"]),
                             bquote(p["wingspan"]),
                             bquote(p["ft. diversity"]),
                             bquote(p["ease of ID"]),
                             bquote(p["peak phenology"]),
                             bquote(p["phenology decay"]),
                             bquote(rho["occurrence-detection"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3.5, 3.5), breaks = c(-2, 0, 2)) +
   guides(color = guide_legend(title = "city")) +
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
)

s <- s +
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50, group=city_name, colour=city_name),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y, colour=city_name), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
s

cowplot::plot_grid(p, q, r, s, ncol = 2)
