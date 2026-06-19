library(ggdag)
library(tidyverse)

# ----------------------------------------------------------------------------
# now recreate for how local and landscape urban park factors may impact butterfly occurrence

# create a dag outline causal hypotheses about how local and landscape park characteristics impact butterflies
local <- c("parkSize", "plantDiversity", "percTree")
landscape <- c("landscapeConnectivity", "landscapeHerbaceous", "landscapeWoody")

occurrence_dag <- dagify(occurrence ~ plantDiversity,
                               occurrence ~ parkSize,
                               occurrence ~ percTree,
                               plantDiversity ~ parkSize,
                               plantDiversity ~ percTree,
                               occurrence ~ landscapeConnectivity,
                               landscapeWoody ~ landscapeConnectivity,
                               landscapeHerbaceous ~ landscapeConnectivity,
                               occurrence ~ landscapeHerbaceous,
                               occurrence ~ landscapeWoody,
                               labels = c(
                                 "occurrence" = "Butterfly\n Occupancy",
                                 "plantDiversity" = "Plant\n Diversity",
                                 "parkSize" = "Park Size",
                                 "percTree" = "% Tree\n Cover",
                                 "landscapeConnectivity" = "Landscape\n Connectivity",
                                 "landscapeHerbaceous" = "Landscape\n Herb. Vegetation\n Area",
                                 "landscapeWoody" = "Landscape\n Woody Vegetation\n Area"
                               ),
                               exposure = "parkSize",
                               outcome = "occurrence"
) %>%
  tidy_dagitty() %>%
  mutate(colour = ifelse(name %in% local, "Local Predictor", "Outcome"),
         colour = ifelse(name %in% landscape, "Landscape Predictor", colour))

occurrence_dag %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(colour = colour), size = 30) +
  geom_dag_edges(edge_alpha = 1, edge_colour = 'grey') +
  geom_dag_text(aes(label = label), colour = "black", size = 6) +
  theme_dag() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=18))

#ggdag(local_occurrence_dag, text = FALSE, use_labels = "label")

ggdag_paths(occurrence_dag, text = FALSE, use_labels = "label", shadow = TRUE)
# if we're interest in effect of plant species density on butterfly occurrence
# there is 1 backdoor path through habitat diversity and another through park size

# identify the minimum adjustment set
ggdag_adjustment_set(occurrence_dag, text = FALSE, use_labels = "label", shadow = TRUE)
ggdag::ggdag_dconnected(occurrence_dag)
