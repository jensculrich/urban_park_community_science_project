# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

stan_out <- readRDS("./model_outputs/stan_out_LA_2km_connectivity_family.rds")
stan_out2 <- readRDS("./model_outputs/stan_out_NYC_2km_connectivity_family.rds")
fit_summary <- rstan::summary(stan_out)
fit_summary2 <- rstan::summary(stan_out2)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

#-------------------------------------------------------------------------------
# initial occurrence (psi1)

# ecological params
# number of params to plot
params <- 10
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
       fit_summary$summary[1,1], # psi1 - intercept
       fit_summary$summary[3,1], # psi1 - wingspan
       fit_summary$summary[4,1], # psi1 - park size
       fit_summary$summary[5,1], # psi1 - connectivity
       fit_summary$summary[6,1], # psi1 - plant richness
       # NYC level inference
       fit_summary2$summary[1,1], # psi1 - intercept
       fit_summary2$summary[3,1], # psi1 - wingspan
       fit_summary2$summary[4,1], # psi1 - park size
       fit_summary2$summary[5,1], # psi1 - connectivity
       fit_summary2$summary[6,1] # psi1 - plant richness
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[1,4], # psi1 - intercept
  fit_summary$summary[3,4], # psi1 - wingspan
  fit_summary$summary[4,4], # psi1 - park size
  fit_summary$summary[5,4], # psi1 - connectivity
  fit_summary$summary[6,4], # psi1 - plant richness
  # NYC level inference
  fit_summary2$summary[1,4], # psi1 - intercept
  fit_summary2$summary[3,4], # psi1 - wingspan
  fit_summary2$summary[4,4], # psi1 - park size
  fit_summary2$summary[5,4], # psi1 - connectivity
  fit_summary2$summary[6,4] # psi1 - plant richness
)

upper_95 <- c(# LA level inference
  fit_summary$summary[1,8], # psi1 - intercept
  fit_summary$summary[3,8], # psi1 - wingspan
  fit_summary$summary[4,8], # psi1 - park size
  fit_summary$summary[5,8], # psi1 - connectivity
  fit_summary$summary[6,8], # psi1 - plant richness
  # NYC level inference
  fit_summary2$summary[1,8], # psi1 - intercept
  fit_summary2$summary[3,8], # psi1 - wingspan
  fit_summary2$summary[4,8], # psi1 - park size
  fit_summary2$summary[5,8], # psi1 - connectivity
  fit_summary2$summary[6,8] # psi1 - plant richness
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[1,5], # psi1 - intercept
  fit_summary$summary[3,5], # psi1 - wingspan
  fit_summary$summary[4,5], # psi1 - park size
  fit_summary$summary[5,5], # psi1 - connectivity
  fit_summary$summary[6,5], # psi1 - plant richness
  # NYC level inference
  fit_summary2$summary[1,5], # psi1 - intercept
  fit_summary2$summary[3,5], # psi1 - wingspan
  fit_summary2$summary[4,5], # psi1 - park size
  fit_summary2$summary[5,5], # psi1 - connectivity
  fit_summary2$summary[6,5] # psi1 - plant richness
)

upper_50 <- c(# LA level inference
  fit_summary$summary[1,7], # psi1 - intercept
  fit_summary$summary[3,7], # psi1 - wingspan
  fit_summary$summary[4,7], # psi1 - park size
  fit_summary$summary[5,7], # psi1 - connectivity
  fit_summary$summary[6,7], # psi1 - plant richness
  # NYC level inference
  fit_summary2$summary[1,7], # psi1 - intercept
  fit_summary2$summary[3,7], # psi1 - wingspan
  fit_summary2$summary[4,7], # psi1 - park size
  fit_summary2$summary[5,7], # psi1 - connectivity
  fit_summary2$summary[6,7] # psi1 - plant richness
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
                              bquote(psi["connectivity"]),
                              bquote(psi["plant richness"])
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
params <- 10
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[7,1], # gamma - intercept
  fit_summary$summary[9,1], # gamma - wingspan
  fit_summary$summary[10,1], # gamma - park size
  fit_summary$summary[11,1], # gamma - connectivity
  fit_summary$summary[12,1], # gamma - plant richness
  # NYC level inference
  fit_summary2$summary[7,1], # gamma - intercept
  fit_summary2$summary[9,1], # gamma - wingspan
  fit_summary2$summary[10,1], # gamma - park size
  fit_summary2$summary[11,1], # gamma - connectivity
  fit_summary2$summary[12,1] # gamma - plant richness
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[7,4], # gamma - intercept
  fit_summary$summary[9,4], # gamma - wingspan
  fit_summary$summary[10,4], # gamma - park size
  fit_summary$summary[11,4], # gamma - connectivity
  fit_summary$summary[12,4], # gamma - plant richness
  # NYC level inference
  fit_summary2$summary[7,4], # gamma - intercept
  fit_summary2$summary[9,4], # gamma - wingspan
  fit_summary2$summary[10,4], # gamma - park size
  fit_summary2$summary[11,4], # gamma - connectivity
  fit_summary2$summary[12,4] # gamma - plant richness
)

upper_95 <- c(# LA level inference
  fit_summary$summary[7,8], # gamma - intercept
  fit_summary$summary[9,8], # gamma - wingspan
  fit_summary$summary[10,8], # gamma - park size
  fit_summary$summary[11,8], # gamma - connectivity
  fit_summary$summary[12,8], # gamma - plant richness
  # NYC level inference
  fit_summary2$summary[7,8], # gamma - intercept
  fit_summary2$summary[9,8], # gamma - wingspan
  fit_summary2$summary[10,8], # gamma - park size
  fit_summary2$summary[11,8], # gamma - connectivity
  fit_summary2$summary[12,8] # gamma - plant richness
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[7,5], # gamma - intercept
  fit_summary$summary[9,5], # gamma - wingspan
  fit_summary$summary[10,5], # gamma - park size
  fit_summary$summary[11,5], # gamma - connectivity
  fit_summary$summary[12,5], # gamma - plant richness
  # NYC level inference
  fit_summary2$summary[7,5], # gamma - intercept
  fit_summary2$summary[9,5], # gamma - wingspan
  fit_summary2$summary[10,5], # gamma - park size
  fit_summary2$summary[11,5], # gamma - connectivity
  fit_summary2$summary[12,5] # gamma - plant richness
)

upper_50 <- c(# LA level inference
  fit_summary$summary[7,7], # gamma - intercept
  fit_summary$summary[9,7], # gamma - wingspan
  fit_summary$summary[10,7], # gamma - park size
  fit_summary$summary[11,7], # gamma - connectivity
  fit_summary$summary[12,7], # gamma - plant richness
  # NYC level inference
  fit_summary2$summary[7,7], # gamma - intercept
  fit_summary2$summary[9,7], # gamma - wingspan
  fit_summary2$summary[10,7], # gamma - park size
  fit_summary2$summary[11,7], # gamma - connectivity
  fit_summary2$summary[12,7] # gamma - plant richness
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
                             bquote(gamma["connectivity"]),
                             bquote(gamma["plant richness"])
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
params <- 10
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[13,1], # phi - intercept
  fit_summary$summary[15,1], # phi - wingspan
  fit_summary$summary[16,1], # phi - park size
  fit_summary$summary[17,1], # phi - connectivity
  fit_summary$summary[18,1], # phi - plant richness
  # NYC level inference
  fit_summary2$summary[13,1], # phi - intercept
  fit_summary2$summary[15,1], # phi - wingspan
  fit_summary2$summary[16,1], # phi - park size
  fit_summary2$summary[17,1], # phi - connectivity
  fit_summary2$summary[18,1] # phi - plant richness
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[13,4], # phi - intercept
  fit_summary$summary[15,4], # phi - wingspan
  fit_summary$summary[16,4], # phi - park size
  fit_summary$summary[17,4], # phi - connectivity
  fit_summary$summary[18,4], # phi - plant richness
  # NYC level inference
  fit_summary2$summary[13,4], # phi - intercept
  fit_summary2$summary[15,4], # phi - wingspan
  fit_summary2$summary[16,4], # phi - park size
  fit_summary2$summary[17,4], # phi - connectivity
  fit_summary2$summary[18,4] # phi - plant richness
)

upper_95 <- c(# LA level inference
  fit_summary$summary[13,8], # phi - intercept
  fit_summary$summary[15,8], # phi - wingspan
  fit_summary$summary[16,8], # phi - park size
  fit_summary$summary[17,8], # phi - connectivity
  fit_summary$summary[18,8], # phi - plant richness
  # NYC level inference
  fit_summary2$summary[13,8], # phi - intercept
  fit_summary2$summary[15,8], # phi - wingspan
  fit_summary2$summary[16,8], # phi - park size
  fit_summary2$summary[17,8], # phi - connectivity
  fit_summary2$summary[18,8] # phi - plant richness
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[13,5], # phi - intercept
  fit_summary$summary[15,5], # phi - wingspan
  fit_summary$summary[16,5], # phi - park size
  fit_summary$summary[17,5], # phi - connectivity
  fit_summary$summary[18,5], # phi - plant richness
  # NYC level inference
  fit_summary2$summary[13,5], # phi - intercept
  fit_summary2$summary[15,5], # phi - wingspan
  fit_summary2$summary[16,5], # phi - park size
  fit_summary2$summary[17,5], # phi - connectivity
  fit_summary2$summary[18,5] # phi - plant richness
)

upper_50 <- c(# LA level inference
  fit_summary$summary[13,7], # phi - intercept
  fit_summary$summary[15,7], # phi - wingspan
  fit_summary$summary[16,7], # phi - park size
  fit_summary$summary[17,7], # phi - connectivity
  fit_summary$summary[18,7], # phi - plant richness
  # NYC level inference
  fit_summary2$summary[13,7], # phi - intercept
  fit_summary2$summary[15,7], # phi - wingspan
  fit_summary2$summary[16,7], # phi - park size
  fit_summary2$summary[17,7], # phi - connectivity
  fit_summary2$summary[18,7] # phi - plant richness
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
                             bquote(phi["connectivity"]),
                             bquote(phi["plant richness"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-2, 7)) +
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
params <- 12
X <- rep(seq(1:(0.5*params)), 2) #  ecological params of interest
city_name <- rep(c("LA", "NYC"), each=(0.5*params)) 

# mean of eco params
Y <- c(# LA level inference
  fit_summary$summary[19,1], # p - intercept
  fit_summary$summary[21,1], # p - wingspan
  fit_summary$summary[22,1], # p - feature diversity
  fit_summary$summary[23,1], # p - ease of id
  fit_summary$summary[24,1], # p - peak phenology
  fit_summary$summary[26,1], # p - phenology decline
  # NYC level inference
  fit_summary2$summary[19,1], # p - intercept
  fit_summary2$summary[21,1], # p - wingspan
  fit_summary2$summary[22,1], # p - feature diversity
  fit_summary2$summary[23,1], # p - ease of id
  fit_summary2$summary[24,1], # p - peak phenology
  fit_summary2$summary[26,1] # p - phenology decline
)

# confidence intervals
lower_95 <- c(# LA level inference
  fit_summary$summary[19,4], # p - intercept
  fit_summary$summary[21,4], # p - wingspan
  fit_summary$summary[22,4], # p - feature diversity
  fit_summary$summary[23,4], # p - ease of id
  fit_summary$summary[24,4], # p - peak phenology
  fit_summary$summary[26,4], # p - phenology decline
  # NYC level inference
  fit_summary2$summary[19,4], # p - intercept
  fit_summary2$summary[21,4], # p - wingspan
  fit_summary2$summary[22,4], # p - feature diversity
  fit_summary2$summary[23,4], # p - ease of id
  fit_summary2$summary[24,4], # p - peak phenology
  fit_summary2$summary[26,4] # p - phenology decline
)

upper_95 <- c(# LA level inference
  fit_summary$summary[19,8], # p - intercept
  fit_summary$summary[21,8], # p - wingspan
  fit_summary$summary[22,8], # p - feature diversity
  fit_summary$summary[23,8], # p - ease of id
  fit_summary$summary[24,8], # p - peak phenology
  fit_summary$summary[26,8], # p - phenology decline
  # NYC level inference
  fit_summary2$summary[19,8], # p - intercept
  fit_summary2$summary[21,8], # p - wingspan
  fit_summary2$summary[22,8], # p - feature diversity
  fit_summary2$summary[23,8], # p - ease of id
  fit_summary2$summary[24,8], # p - peak phenology
  fit_summary2$summary[26,8] # p - phenology decline
)

# confidence intervals
lower_50 <- c(# LA level inference
  fit_summary$summary[19,5], # p - intercept
  fit_summary$summary[21,5], # p - wingspan
  fit_summary$summary[22,5], # p - feature diversity
  fit_summary$summary[23,5], # p - ease of id
  fit_summary$summary[24,5], # p - peak phenology
  fit_summary$summary[26,5], # p - phenology decline
  # NYC level inference
  fit_summary2$summary[19,5], # p - intercept
  fit_summary2$summary[21,5], # p - wingspan
  fit_summary2$summary[22,5], # p - feature diversity
  fit_summary2$summary[23,5], # p - ease of id
  fit_summary2$summary[24,5], # p - peak phenology
  fit_summary2$summary[26,5] # p - phenology decline
)

upper_50 <- c(# LA level inference
  fit_summary$summary[19,7], # p - intercept
  fit_summary$summary[21,7], # p - wingspan
  fit_summary$summary[22,7], # p - feature diversity
  fit_summary$summary[23,7], # p - ease of id
  fit_summary$summary[24,7], # p - peak phenology
  fit_summary$summary[26,7], # p - phenology decline
  # NYC level inference
  fit_summary2$summary[19,7], # p - intercept
  fit_summary2$summary[21,7], # p - wingspan
  fit_summary2$summary[22,7], # p - feature diversity
  fit_summary2$summary[23,7], # p - ease of id
  fit_summary2$summary[24,7], # p - peak phenology
  fit_summary2$summary[26,7] # p - phenology decline
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
                             bquote(p["phenology decay"])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3, 2), breaks = c(-2, 0, 2)) +
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
