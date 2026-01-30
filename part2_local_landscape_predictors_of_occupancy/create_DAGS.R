library(ggdag)
library(tidyverse)

----------------------------------------------------------------------------
# now recreate for how local and landscape urban park factors may impact butterfly occurrence

# create a dag outline causal hypotheses about how smoking impacts cardiac arrest
local <- c("parkSize", "plantDensity", "habitatDiversity", "percTree")
landscape <- c("landscapeIsolation", "landscapeGreen", "landscapeDiversity")

occurrence_dag <- dagify(occurrence ~ plantDensity,
                               occurrence ~ parkSize,
                               occurrence ~ percTree,
                               occurrence ~ habitatDiversity,
                               habitatDiversity ~ parkSize,
                               plantDensity ~ habitatDiversity,
                               plantDensity ~ percTree,
                               occurrence ~ landscapeIsolation,
                               occurrence ~ landscapeGreen,
                               occurrence ~ landscapeDiversity,
                               labels = c(
                                 "occurrence" = "Butterfly\n Occurrence",
                                 "plantDensity" = "Plant Species\n Density",
                                 "habitatDiversity" = "Habitat\n Diversity\n (latent)",
                                 "parkSize" = "Park Size",
                                 "percTree" = "% Tree\n Cover",
                                 "landscapeIsolation" = "Landscape\n Isolation",
                                 "landscapeGreen" = "Landscape\n Vegetation",
                                 "landscapeDiversity" = "Landscape\n Diversity"
                               ),
                               latent = "habitatDiversity",
                               exposure = "parkSize",
                               outcome = "occurrence"
) %>%
  tidy_dagitty() %>%
  mutate(colour = ifelse(name %in% local, "Local Predictor", "Outcome"),
         colour = ifelse(name %in% landscape, "Landscape Predictor", colour))

occurrence_dag %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(colour = colour)) +
  geom_dag_edges(edge_alpha = 1, edge_colour = 'grey') +
  geom_dag_text(aes(label = label), colour = "black", size = 4) +
  theme_dag() +
  theme(legend.title = element_blank())

#ggdag(local_occurrence_dag, text = FALSE, use_labels = "label")

ggdag_paths(occurrence_dag, text = FALSE, use_labels = "label", shadow = TRUE)
# if we're interest in effect of plant species density on butterfly occurrence
# there is 1 backdoor path through habitat diversity and another through park size

# identify the minimum adjustment set
ggdag_adjustment_set(occurrence_dag, text = FALSE, use_labels = "label", shadow = TRUE)
ggdag::ggdag_dconnected(occurrence_dag)
