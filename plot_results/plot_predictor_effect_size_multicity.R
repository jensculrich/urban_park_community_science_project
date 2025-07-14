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

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

# ecological params
# number of params to plot
params <- 5
X <- rep(seq(1:params), 2) #  ecological params of interest
Y <- vector(length = n_cities*params)
lower_95 <- vector(length = n_cities*params)
upper_95 <- vector(length = n_cities*params)
lower_50 <- vector(length = n_cities*params)
upper_50 <- vector(length = n_cities*params)

city_name <- rep(city_names, each=(params)) 

for(i in 1:n_cities){
  city <- city_names[i]
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", 
    city, "_2km_connectivity_family_100buffers.rds"))
  fit_summary <- rstan::summary(stan_out)
  
  index_lower <- 1 + ((i-1) * params)
  index_upper <- 1 + ((i-1) * params) + (params - 1)
  
  Y[index_lower:index_upper] <- c(
    fit_summary$summary[1,1], # psi1 - intercept
    fit_summary$summary[3,1], # psi1 - wingspan
    fit_summary$summary[4,1], # psi1 - park size
    fit_summary$summary[5,1], # psi1 - connectivity
    #fit_summary$summary[8,1], # psi1 - plant richness
    fit_summary$summary[6,1] # psi1 - tree cover
  )
  
  lower_95[index_lower:index_upper] <- c(
    fit_summary$summary[1,4], # psi1 - intercept
    fit_summary$summary[3,4], # psi1 - wingspan
    fit_summary$summary[4,4], # psi1 - park size
    fit_summary$summary[5,4], # psi1 - connectivity
    #fit_summary$summary[8,4], # psi1 - plant richness
    fit_summary$summary[6,4] # psi1 - tree cover
  )
  
  upper_95[index_lower:index_upper] <- c(
    fit_summary$summary[1,8], # psi1 - intercept
    fit_summary$summary[3,8], # psi1 - wingspan
    fit_summary$summary[4,8], # psi1 - park size
    fit_summary$summary[5,8], # psi1 - connectivity
    #fit_summary$summary[8,8], # psi1 - plant richness
    fit_summary$summary[6,8] # psi1 - tree cover
  )
  
  lower_50[index_lower:index_upper] <- c(
    fit_summary$summary[1,5], # psi1 - intercept
    fit_summary$summary[3,5], # psi1 - wingspan
    fit_summary$summary[4,5], # psi1 - park size
    fit_summary$summary[5,5], # psi1 - connectivity
    #fit_summary$summary[8,5, # psi1 - plant richness
    fit_summary$summary[6,5] # psi1 - tree cover
  )
  
  upper_50[index_lower:index_upper] <- c(
    fit_summary$summary[1,7], # psi1 - intercept
    fit_summary$summary[3,7], # psi1 - wingspan
    fit_summary$summary[4,7], # psi1 - park size
    fit_summary$summary[5,7], # psi1 - connectivity
    #fit_summary$summary[8,7], # psi1 - plant richness
    fit_summary$summary[6,7] # psi1 - tree cover
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

## --------------------------------------------------
## Draw caterpillar plot

(p <- ggplot(df_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = seq(1:params),
                     labels=c(bquote(psi["intercept"]),
                              bquote(psi["wingspan"]),
                              bquote(psi["park size"]),
                              bquote(psi["isolation"]),
                              #bquote(psi["plant richness"]),
                              bquote(psi["tree cover"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-5, 5), breaks = c(-4, -2, 0, 2, 4)) +
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
params <- 5
X <- rep(seq(1:params), 2) #  ecological params of interest
Y <- vector(length = n_cities*params)
lower_95 <- vector(length = n_cities*params)
upper_95 <- vector(length = n_cities*params)
lower_50 <- vector(length = n_cities*params)
upper_50 <- vector(length = n_cities*params)

city_name <- rep(city_names, each=(params)) 

for(i in 1:n_cities){
  city <- city_names[i]
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", 
    city, "_2km_connectivity_family_100buffers.rds"))
  fit_summary <- rstan::summary(stan_out)
  
  index_lower <- 1 + ((i-1) * params)
  index_upper <- 1 + ((i-1) * params) + (params - 1)
  
  Y[index_lower:index_upper] <- c(
    fit_summary$summary[7,1], # psi1 - intercept
    fit_summary$summary[9,1], # psi1 - wingspan
    fit_summary$summary[10,1], # psi1 - park size
    fit_summary$summary[11,1], # psi1 - connectivity
    #fit_summary$summary[8,1], # psi1 - plant richness
    fit_summary$summary[12,1] # psi1 - tree cover
  )
  
  lower_95[index_lower:index_upper] <- c(
    fit_summary$summary[7,4], # psi1 - intercept
    fit_summary$summary[9,4], # psi1 - wingspan
    fit_summary$summary[10,4], # psi1 - park size
    fit_summary$summary[11,4], # psi1 - connectivity
    #fit_summary$summary[8,4], # psi1 - plant richness
    fit_summary$summary[12,4] # psi1 - tree cover
  )
  
  upper_95[index_lower:index_upper] <- c(
    fit_summary$summary[7,8], # psi1 - intercept
    fit_summary$summary[9,8], # psi1 - wingspan
    fit_summary$summary[10,8], # psi1 - park size
    fit_summary$summary[11,8], # psi1 - connectivity
    #fit_summary$summary[8,8], # psi1 - plant richness
    fit_summary$summary[12,8] # psi1 - tree cover
  )
  
  lower_50[index_lower:index_upper] <- c(
    fit_summary$summary[7,5], # psi1 - intercept
    fit_summary$summary[9,5], # psi1 - wingspan
    fit_summary$summary[10,5], # psi1 - park size
    fit_summary$summary[11,5], # psi1 - connectivity
    #fit_summary$summary[8,5, # psi1 - plant richness
    fit_summary$summary[12,5] # psi1 - tree cover
  )
  
  upper_50[index_lower:index_upper] <- c(
    fit_summary$summary[7,7], # psi1 - intercept
    fit_summary$summary[9,7], # psi1 - wingspan
    fit_summary$summary[10,7], # psi1 - park size
    fit_summary$summary[11,7], # psi1 - connectivity
    #fit_summary$summary[8,7], # psi1 - plant richness
    fit_summary$summary[12,7] # psi1 - tree cover
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

## --------------------------------------------------
## Draw caterpillar plot

(q <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(gamma["intercept"]),
                             bquote(gamma["wingspan"]),
                             bquote(gamma["park size"]),
                             bquote(gamma["isolation"]),
                             #bquote(gamma["plant richness"]),
                             bquote(gamma["tree cover"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-6, 4.5), breaks = c(-6, -4, -2, 0, 2, 4)) +
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
params <- 5
X <- rep(seq(1:params), 2) #  ecological params of interest
Y <- vector(length = n_cities*params)
lower_95 <- vector(length = n_cities*params)
upper_95 <- vector(length = n_cities*params)
lower_50 <- vector(length = n_cities*params)
upper_50 <- vector(length = n_cities*params)

city_name <- rep(city_names, each=(params)) 

for(i in 1:n_cities){
  city <- city_names[i]
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", 
    city, "_2km_connectivity_family_100buffers.rds"))
  fit_summary <- rstan::summary(stan_out)
  
  index_lower <- 1 + ((i-1) * params)
  index_upper <- 1 + ((i-1) * params) + (params - 1)
  
  Y[index_lower:index_upper] <- c(
    fit_summary$summary[13,1], # psi1 - intercept
    fit_summary$summary[15,1], # psi1 - wingspan
    fit_summary$summary[16,1], # psi1 - park size
    fit_summary$summary[17,1], # psi1 - connectivity
    #fit_summary$summary[8,1], # psi1 - plant richness
    fit_summary$summary[18,1] # psi1 - tree cover
  )
  
  lower_95[index_lower:index_upper] <- c(
    fit_summary$summary[13,4], # psi1 - intercept
    fit_summary$summary[15,4], # psi1 - wingspan
    fit_summary$summary[16,4], # psi1 - park size
    fit_summary$summary[17,4], # psi1 - connectivity
    #fit_summary$summary[8,4], # psi1 - plant richness
    fit_summary$summary[18,4] # psi1 - tree cover
  )
  
  upper_95[index_lower:index_upper] <- c(
    fit_summary$summary[13,8], # psi1 - intercept
    fit_summary$summary[15,8], # psi1 - wingspan
    fit_summary$summary[16,8], # psi1 - park size
    fit_summary$summary[17,8], # psi1 - connectivity
    #fit_summary$summary[8,8], # psi1 - plant richness
    fit_summary$summary[18,8] # psi1 - tree cover
  )
  
  lower_50[index_lower:index_upper] <- c(
    fit_summary$summary[13,5], # psi1 - intercept
    fit_summary$summary[15,5], # psi1 - wingspan
    fit_summary$summary[16,5], # psi1 - park size
    fit_summary$summary[17,5], # psi1 - connectivity
    #fit_summary$summary[8,5, # psi1 - plant richness
    fit_summary$summary[18,5] # psi1 - tree cover
  )
  
  upper_50[index_lower:index_upper] <- c(
    fit_summary$summary[13,7], # psi1 - intercept
    fit_summary$summary[15,7], # psi1 - wingspan
    fit_summary$summary[16,7], # psi1 - park size
    fit_summary$summary[17,7], # psi1 - connectivity
    #fit_summary$summary[8,7], # psi1 - plant richness
    fit_summary$summary[18,7] # psi1 - tree cover
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

## --------------------------------------------------
## Draw caterpillar plot

(r <- ggplot(df_estimates) +
   theme_bw() +
   scale_x_discrete(name="", breaks = seq(1:params),
                    labels=c(bquote(phi["intercept"]),
                             bquote(phi["wingspan"]),
                             bquote(phi["park size"]),
                             bquote(phi["isolation"]),
                             #bquote(phi["plant richness"]),
                             bquote(phi["tree cover"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-5, 5), breaks=c(-4, -2, 0, 2, 4)) +
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

# detection params
# number of params to plot
params <- 6
X <- rep(seq(1:params), 2) #  ecological params of interest
Y <- vector(length = n_cities*params)
lower_95 <- vector(length = n_cities*params)
upper_95 <- vector(length = n_cities*params)
lower_50 <- vector(length = n_cities*params)
upper_50 <- vector(length = n_cities*params)

city_name <- rep(city_names, each=(params)) 

for(i in 1:n_cities){
  city <- city_names[i]
  
  stan_out <- readRDS(paste0(
    "./model_outputs/stan_out_", 
    city, "_2km_connectivity_family_100buffers.rds"))
  fit_summary <- rstan::summary(stan_out)
  
  index_lower <- 1 + ((i-1) * params)
  index_upper <- 1 + ((i-1) * params) + (params - 1)
  
  Y[index_lower:index_upper] <- c(
    fit_summary$summary[19,1], # p - intercept
    fit_summary$summary[21,1], # p - wingspan
    fit_summary$summary[22,1], # p - feature diversity
    fit_summary$summary[23,1], # p - ease of id
    fit_summary$summary[24,1], # p - date
    fit_summary$summary[26,1] # p - date ^2
  )
  
  lower_95[index_lower:index_upper] <- c(
    fit_summary$summary[19,4], # p - intercept
    fit_summary$summary[21,4], # p - wingspan
    fit_summary$summary[22,4], # p - feature diversity
    fit_summary$summary[23,4], # p - ease of id
    fit_summary$summary[24,4], # p - date
    fit_summary$summary[26,4] # p - date ^2
  )
  
  upper_95[index_lower:index_upper] <- c(
    fit_summary$summary[19,8], # p - intercept
    fit_summary$summary[21,8], # p - wingspan
    fit_summary$summary[22,8], # p - feature diversity
    fit_summary$summary[23,8], # p - ease of id
    fit_summary$summary[24,8], # p - date
    fit_summary$summary[26,8] # p - date ^2
  )
  
  lower_50[index_lower:index_upper] <- c(
    fit_summary$summary[19,5], # p - intercept
    fit_summary$summary[21,5], # p - wingspan
    fit_summary$summary[22,5], # p - feature diversity
    fit_summary$summary[23,5], # p - ease of id
    fit_summary$summary[24,5], # p - date
    fit_summary$summary[26,5] # p - date ^2
  )
  
  upper_50[index_lower:index_upper] <- c(
    fit_summary$summary[19,7], # p - intercept
    fit_summary$summary[21,7], # p - wingspan
    fit_summary$summary[22,7], # p - feature diversity
    fit_summary$summary[23,7], # p - ease of id
    fit_summary$summary[24,7], # p - date
    fit_summary$summary[26,7] # p - date ^2
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

## --------------------------------------------------
## Draw caterpillar plot

(s <- ggplot(df_estimates) +
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
                      limits = c(-3.5, 3.5), breaks = c(-4, -2, 0, 2)) +
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
