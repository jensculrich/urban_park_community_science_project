
f5.1 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/figure5.1.rds")
f5.2 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/figure5.2.rds")
f5.3 <- readRDS("./part3_citywide_drivers_of_diversity/figures/m3_plots/figure5.3.rds")

cowplot::plot_grid(f5.1, f5.2, f5.3, ncol = 1
)


