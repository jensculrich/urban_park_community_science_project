# https://avehtari.github.io/modelselection/collinear.html
# running through above example^

library(rstanarm)
library(tidyverse)
library(bayesplot)
library(projpred)

SEED <- 19

# all this data generation is from Cade 2015
# doesn't matter what this is -- if you use a different number your results will be different from mine.
set.seed(SEED)
df <- tibble(
  pos.tot = runif(200,min=0.8,max=1.0),
  urban.tot = pmin(runif(200,min=0.0,max=0.02),1.0 - pos.tot),
  neg.tot = (1.0 - pmin(pos.tot + urban.tot,1)),
  x1= pmax(pos.tot - runif(200,min=0.05,max=0.30),0),
  x3= pmax(neg.tot - runif(200,min=0.0,max=0.10),0),
  x2= pmax(pos.tot - x1 - x3/2,0),
  x4= pmax(1 - x1 - x2 - x3 - urban.tot,0))
# true model and 200 Poisson observations
mean.y <- exp(-5.8 + 6.3*df$x1 + 15.2*df$x2)
df$y <- rpois(200,mean.y)

pairs(df,diag=list(continuous="barDiag"))

fitg <- stan_glm(y ~ x1 + x2 + x3 + x4, data = df, na.action = na.fail, family=poisson(), seed=SEED)
summary(fitg)

fitg <- stan_glm(y ~ x1 + x2 + x3 + x4, data = df, na.action = na.fail, family=poisson(), QR=TRUE, seed=SEED)
summary(fitg)

mcmc_areas(as.matrix(fitg),prob_outer = .99)
mcmc_pairs(as.matrix(fitg),pars = c("x1","x2","x3","x4"))

# In case of even more variables with some being relevant and some irrelevant, 
# it will be difficult to analyse joint posterior to see which variables are 
# jointly relevant. We can easily test whether any of the covariates are 
# useful by using cross-validation to compare to a null model,
fitg0 <- stan_glm(y ~ 1, data = df, na.action = na.fail, family=poisson(), seed=SEED)
(loog <- loo(fitg))
(loog0 <- loo(fitg0))
loo_compare(loog0, loog)
# Based on cross-validation covariates together contain significant information to improve predictions.

fitg_cv <- cv_varsel(fitg, method='forward', cv_method='LOO')
plot(fitg_cv, stats = c('elpd', 'rmse'))

# And we get a LOO based recommendation for the model size to choose
(nsel <- suggest_size(fitg_cv, alpha=0.1))
(vsel <- solution_terms(fitg_cv)[1:nsel])

# Next we form the projected posterior for the chosen model.
projg <- project(fitg_cv, nv = nsel, ns = 4000)
projdraws <- as.matrix(projg)
round(colMeans(projdraws),1)
round(posterior_interval(projdraws),1)
mcmc_areas(projdraws, pars=c("(Intercept)",vsel))
