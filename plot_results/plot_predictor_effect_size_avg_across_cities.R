# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

# list of city names
city_names <- c(
  "Boston", # 1
  "Dallas", # 2
  "Houston", # 3
  "LA", # 4
  "NYC", # 5
  "Riverside", # 6
  "SF" # 7
)

n_cities <- length(city_names)

# params to plot
param_names <- c("psi1_0", 
                 "sigma_psi1_species",
                 "psi1_wingspan",
                 "psi1_park_size",
                 "psi1_connectivity",
                 #"psi1_plant_genera",
                 "psi1_tree_cover",
                 
                 "gamma0", 
                 "sigma_gamma_species",
                 "gamma_wingspan",
                 "gamma_park_size",
                 "gamma_connectivity",
                 #"gamma_plant_genera",
                 "gamma_tree_cover",
                 #"gamma_wingspan_connectivity",
                 
                 "phi0", 
                 "sigma_phi_species",
                 "phi_wingspan",
                 "phi_park_size",
                 "phi_connectivity",
                 #"phi_plant_genera",
                 "phi_tree_cover",
                 
                 "p0", 
                 "sigma_p_species",
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

df_estimates$X <- as.factor(as.numeric(df_estimates$X))
df_estimates$Y <- as.numeric(df_estimates$mean)
df_estimates$lower_95 <- as.numeric(df_estimates$lower_95)
df_estimates$upper_95 <- as.numeric(df_estimates$upper_95)
df_estimates$lower_50 <- as.numeric(df_estimates$lower_50)
df_estimates$upper_50 <- as.numeric(df_estimates$upper_50)


## --------------------------------------------------
## Draw caterpillar plot for psi1

psi1_estimates <- df_estimates[c(1, 3, 4, 5, 6),]

(p <- ggplot(psi1_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = c(1, 3, 4, 5, 6),
                    labels=c(bquote(psi["intercept"]),
                             bquote(psi["wingspan"]),
                             bquote(psi["park size"]),
                             bquote(psi["isolation"]),
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


## --------------------------------------------------
## Draw caterpillar plot for gamma

gamma_estimates <- df_estimates[c(7, 9, 10, 11, 12),]

(q <- ggplot(gamma_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(7, 9, 10, 11, 12),
                     labels=c(bquote(gamma["intercept"]),
                              bquote(gamma["wingspan"]),
                              bquote(gamma["park size"]),
                              bquote(gamma["isolation"]),
                              bquote(gamma["tree cover"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-5, 2.5), breaks = c(-4, -2, 0, 2)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
q


## --------------------------------------------------
## Draw caterpillar plot for phi

phi_estimates <- df_estimates[c(13, 15, 16, 17, 18),]

(r <- ggplot(phi_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(13, 15, 16, 17, 18),
                     labels=c(bquote(phi["intercept"]),
                              bquote(phi["wingspan"]),
                              bquote(phi["park size"]),
                              bquote(phi["isolation"]),
                              bquote(phi["tree cover"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-2.5, 4), breaks = c(-2, 0, 2, 4)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
r


## --------------------------------------------------
## Draw caterpillar plot for p

p_estimates <- df_estimates[c(19, 21, 22, 23, 24, 26),]

(s <- ggplot(p_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = c(19, 21, 22, 23, 24, 26),
                     labels=c(bquote(p["intercept"]),
                              bquote(p["wingspan"]),
                              bquote(p["ft. diversity"]),
                              bquote(p["ease of ID"]),
                              bquote(p["peak phenology"]),
                              bquote(p["phenology decay"]),
                              bquote(rho["occurrence-detection"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-3, 3), breaks = c(-2, 0, 2)) +
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
  geom_errorbar(aes(x=X, ymin=lower_95, ymax=upper_95),
                position=position_dodge(width=0.5),width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                position=position_dodge(width=0.5),width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y), 
             position=position_dodge(width=0.5),
             size = 5, alpha = 0.8) 
s

## --------------------------------------------------
## Draw caterpillar plot for all

cowplot::plot_grid(p, q, r, s, ncol = 2)
