# plot effects of predictors on occurrence and detection

library(tidyverse)
library(rstan)

stan_out <- readRDS("./model_outputs/stan_out.rds")
fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

# parameter means
params <- 6

# mu_alpha0, mu_alpha1, mu_beta0, mu_beta1
X <- seq(1:params) # 10 ecological params of interest

# mean of eco params
Y <- c(fit_summary$summary[8,1], # psi - wingspan
       fit_summary$summary[9,1], # psi - park size
       fit_summary$summary[17,1], # p - wingspan
       fit_summary$summary[18,1], # p - feature diversity
       fit_summary$summary[19,1], # p - ease of id
       fit_summary$summary[20,1] # p - park size
)

# confidence intervals
lower_95 <- c(fit_summary$summary[8,4], # psi - wingspan
              fit_summary$summary[9,4], # psi - park size
              fit_summary$summary[17,4], # p - wingspan
              fit_summary$summary[18,4], # p - feature diversity
              fit_summary$summary[19,4], # p - ease of id
              fit_summary$summary[20,4] # p - park size
)

upper_95 <- c(fit_summary$summary[8,8], # psi - wingspan
              fit_summary$summary[9,8], # psi - park size
              fit_summary$summary[17,8], # p - wingspan
              fit_summary$summary[18,8], # p - feature diversity
              fit_summary$summary[19,8], # p - ease of id
              fit_summary$summary[20,8] # p - park size
)

# confidence intervals
lower_50 <- c(fit_summary$summary[8,5], # psi - wingspan
              fit_summary$summary[9,5], # psi - park size
              fit_summary$summary[17,5], # p - wingspan
              fit_summary$summary[18,5], # p - feature diversity
              fit_summary$summary[19,5], # p - ease of id
              fit_summary$summary[20,5] # p - park size
)

upper_50 <- c(fit_summary$summary[8,7], # psi - wingspan
              fit_summary$summary[9,7], # psi - park size
              fit_summary$summary[17,7], # p - wingspan
              fit_summary$summary[18,7], # p - feature diversity
              fit_summary$summary[19,7], # p - ease of id
              fit_summary$summary[20,7] # p - park size
)

df_estimates <- as.data.frame(cbind(X, Y, lower_95, upper_95, lower_50, upper_50))

df_estimates$X <- as.factor(df_estimates$X)

## --------------------------------------------------
## Draw caterpillar plot

(p <- ggplot(df_estimates) +
    theme_bw() +
    scale_x_discrete(name="", breaks = seq(1:params),
                     labels=c(bquote(psi[wingspan]),
                              bquote(psi["park size"]),
                              bquote(p[wingspan]),
                              bquote(p["ft. diversity"]),
                              bquote(p["ease of id"]),
                              bquote(p["park size"])
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-1, 2.5)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(plot.title = element_text(size = 32, face = "bold"),
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
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X, ymin=lower_50, ymax=upper_50),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X, y=Y),
             size = 5, alpha = 0.8) 
p
