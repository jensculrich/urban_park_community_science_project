library(rstan)
library(tidyverse)
library(gridExtra)

#stan_out <- readRDS("./model_outputs/stan_out4.rds")
fit_summary <- rstan::summary(stan_out)

View(cbind(1:nrow(fit_summary$summary), fit_summary$summary)) # View to see which row corresponds to the parameter of interest

n_species <- nrow(species_names)
species_names <- cbind(species_names, as.data.frame(family_lookup))
#View(species_names)

# parameter means
params = 4 # init, colonization, persistence, detection

x <- (rep(1:params, each=n_species)) # parameter reference
y = (rep(1:n_species, times=params)) # species reference

#x <- 1:n_species

estimate <- matrix(nrow = n_species, ncol = params)
lower <- matrix(nrow = n_species, ncol = params)
upper <- matrix(nrow = n_species, ncol = params)

start <- 165 # row of first species param

for(i in 1:n_species){
  
  estimate[i,1] <- c(
    # param 1 - psi1
    rev(fit_summary$summary[165+(i -1),1])
  )
  estimate[i,2] <- c(
    # param 2 - colonization
    rev(fit_summary$summary[316+(i -1),1])
  )
  estimate[i,3] <- c(
    # param 3 - persistence
    rev(fit_summary$summary[467+(i -1),1])
  )
  estimate[i,4] <- c(
    # param 4 - detection
    rev(fit_summary$summary[618+(i -1),1] + 
          # add family effect
          fit_summary$summary[364+(species_names$family_lookup[i] - 1),1]) 
  )
  
  lower[i,1] <- c(
    # param 1 - psi1
    rev(fit_summary$summary[165+(i -1),4])
  )
  lower[i,2] <- c(
    # param 2 - colonization
    rev(fit_summary$summary[316+(i -1),4])
  )
  lower[i,3] <- c(
    # param 3 - persistence
    rev(fit_summary$summary[467+(i -1),4])
  )
  lower[i,4] <- c(
    # param 4 - detection
    rev(fit_summary$summary[618+(i -1),4] + 
          # add family effect
          fit_summary$summary[769+(species_names$family_lookup[i] - 1),4]) 
  )
  
  upper[i,1] <- c(
    # param 1 - psi1
    rev(fit_summary$summary[165+(i -1),8])
  )
  upper[i,2] <- c(
    # param 2 - colonization
    rev(fit_summary$summary[316+(i -1),8])
  )
  upper[i,3] <- c(
    # param 3 - persistence
    rev(fit_summary$summary[467+(i -1),8])
  )
  upper[i,4] <- c(
    # param 4 - detection
    rev(fit_summary$summary[618+(i -1),8] + 
          # add family effect
          fit_summary$summary[769+(species_names$family_lookup[i] - 1),8]) 
  )

}

View(as.data.frame(estimate))

estimate <- as.numeric(estimate)
lower <- as.numeric(lower)
upper <- as.numeric(upper)

df = as.data.frame(cbind(species_names$species, species_names$family, 
                         x, y, 
                         estimate, lower, upper)) %>%
  rename("species" = "V1",
         "family" = "V2") %>%
  mutate(species = as.factor(species),
         family = as.factor(family),
         x = as.factor(x),
         estimate = as.numeric(estimate),
         lower = as.numeric(lower),
         upper = as.numeric(upper))
  
# flip species names
species_names_label <- species_names$species %>%
  as.data.frame(.) %>%
  map_df(., rev) %>%
  pull(.)

num_per_page = 35
i = 2

df_filtered <- df %>%
  mutate(y_num = as.integer(y)) %>%
  filter(y_num >= num_per_page*i - num_per_page) %>%
  filter(y_num < num_per_page*i)

start <- 100
end <- 140
df_filtered <- df %>%
  mutate(y_num = as.integer(y)) %>%
  filter(y_num >= start) %>%
  filter(y_num < end)

p <- ggplot(df_filtered, aes(x, as.factor(as.numeric(y)), width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(1, 2, 3, 4),
                   labels=c("init. occ.", "colonization", 
                            "persistence", "detection")
                   ) +
  scale_y_discrete(name="", breaks = rep(1:n_species),
                   labels=species_names_label) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df_filtered, 
  #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = df_filtered, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
p


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## make a table with top and bottom (most versus least likely _)
## let's do this for detectability

df_filtered <- filter(df, x == 4) # x == 4 is detection
top <- top_n(df_filtered, 16, estimate) %>%
  mutate(y = as.factor(as.numeric(y)))
bottom <- top_n(df_filtered, -16, estimate) %>%
  mutate(y = as.factor(as.numeric(y)))

top <- top[order(top$estimate,decreasing=TRUE),]
bottom <- bottom[order(bottom$estimate,decreasing=TRUE),]

# join persistence intercepts (x == 3)
df_persistence <- filter(df, x == 3) %>%
  mutate(y = as.factor(as.numeric(y)))
# grab the persistence for the top detec.
test <- filter(df_persistence, y %in% top$y) 
top <- rbind(top, test)
# grab the persistence for the bottom detec.
test <- filter(df_persistence, y %in% bottom$y) 
bottom <- rbind(bottom, test)

top_and_bottom <- rbind(top, bottom) %>% map_df(., rev)

# join species names
species_names_df <- species_names$species %>%
  #rev(.) %>%
  as.data.frame(.) %>%
  mutate(row_id=row_number()) %>%
  mutate(row_id = as.factor(row_id)) %>%
  rename("y" = "row_id",
         "species_name" = ".")

top_and_bottom <- left_join(top_and_bottom, species_names_df, by = "y")

# make a single plot
p1 <- ggplot(top_and_bottom, aes(x, y, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(3, 4),
                   labels=c(
                     "persistence", "detection"
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(species_names_df)),
                   labels=species_names_df$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3")) +
  #geom_text(data = df_filtered, 
  #        aes(x = x, y = y, label = signif(estimate, 2)), size = 3.5) +
  
  geom_text(data = top_and_bottom, 
            aes(x = x, y = y, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", signif(lower,2), ", ", signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

# make a column for the detection rate
temp <- top_and_bottom %>% filter(x == 4) %>%
  mutate(row_id=row_number()) %>%
  mutate(row_id = as.factor(row_id))

# grab row id's for later to match rows
rows <- temp %>%
  select(species_name, row_id)

p1.2 <- ggplot(temp, aes(x, row_id, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(4),
                   labels=c("detection \nrate"
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(temp)),
                   labels=temp$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                       midpoint = median(temp$estimate)) +

  geom_text(data = temp, 
            aes(x = x, y = row_id, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", sprintf("%.1f",lower), ", ",
              #"[", signif(lower,2), ", ", 
              sprintf("%.1f",upper), "]")),
            #signif(upper,2), "]")),
            size = 3.5) +
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

# row order needs to match row order of slope column
temp2 <- top_and_bottom %>% filter(x == 3)
temp2 <- left_join(temp2, rows)

p1.1 <- ggplot(temp2, aes(x, row_id, width=1, height=1)) +
  geom_tile(aes(fill = estimate)) +
  theme_bw() +
  scale_x_discrete(name="", breaks = c(3),
                   labels=c("persiste.. \nrate"
                   )) +
  scale_y_discrete(name="", breaks = rep(1:nrow(temp)),
                   labels=temp$species_name) +
  scale_fill_gradient2(low = ("firebrick3"), high = ("dodgerblue3"),
                       midpoint = median(temp2$estimate)) +
  geom_text(data = temp2, 
            aes(x = x, y = row_id, label = paste0(
              #signif(estimate, 2),"\n(", 
              "[", sprintf("%.1f",lower), ", ",
              #"[", signif(lower,2), ", ", 
              sprintf("%.1f",upper), "]")),
            #signif(upper,2), "]")),
            size = 3.5) +
  
  theme(legend.position = "none",
        #legend.text=element_text(size=14),
        #legend.title=element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 45, hjust=1),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())

grid.arrange(p1.2, p1.1, ncol=2)

