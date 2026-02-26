m3.1 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.1_plot.rds")
m3.2 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.2_plot.rds")
m3.3 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.3_plot.rds")
m3.4 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.4_plot.rds")
legend <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/m3_legend.rds")

cowplot::plot_grid(m3.1, m3.2, m3.3, m3.4, legend, ncol = 1, 
                   labels = c("a)","b)","c)","d)", ""), label_size = 12,
                   rel_heights = c(1,1,1,1,0.2))

# save as 20 X 30 portrait?

cowplot::plot_grid(m3.1, m3.2, ncol = 1, 
                   labels = c("a)","b)", ""), label_size = 20,
                   rel_heights = c(1,1,0.1))

cowplot::plot_grid(m3.3, m3.4, ncol = 1, 
                   labels = c("a)","b)", ""), label_size = 20,
                   rel_heights = c(1,1,0.1))
