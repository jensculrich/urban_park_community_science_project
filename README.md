# urban_park_community_science_project

right now I've started working towards a framework for identifying the effects of urban park habitat on insect biodiversity.

I've started by using LA county as a case study (based on large size, long running community science monitoring programs and some familiarity with the group in terms of knowledge of local species and geographical/socieoeconomic context).

I've tried utilizing only data from the city nature challenge program (annual 4 day survey event that takes place each April in cities all across the country) as well as utilizing all of the iNaturalist data available. I've tried defining spatial patches using a simple grid and alternatively using a layer with urban park polygons + buffer area.

I've applied a dynamic occupancy model to the data that allows for estimating species colonization and persistence rates as well as species-specific detectability. Ecological covariates such as park size, plant diversity, vegetation height, landscape context, etc., could be added to such a model. Additional detection covariates could also be added.

So far I've treated detection as a continuous process, i.e., assuming that all species are searched for at all sites at all times. It might be more appropriate to identify community sampling events to then infer non-detections.