# urban_park_community_science_project

Data and code to accompany "Community science synthesis across 22 major US cities reveals 
consistent positive effects of urban park size and tree cover on butterfly diversity".

With this project we wanted to examine how within and between city differences impact
butterfly diversity in urban parks. We used community science data from iNaturalist
to tackle this objective.

The community science data used in the analyses are available directly from GBIF at the cited DOIs. 
Land cover data were obtained from previously published publicly available sources cited in the text:
UrbanWatch data (fine scale land cover data within city boundaries only), 
NLCD data (coarse scale land cover data within and outside of city boundaries), and
ParkServe data (shapefile that outlines delineated park spaces across the U.S.).  
We used some previously published data on butterfly traits to support our analyses, 
however, we also collected some new trait data to supplement the existing sources (see below). 

The analyses were broken into 3 main parts:

## Part 1: Effects of Park Size and Connectivity on Butterfly Community Dynamics

Using dynamic multi-species occupancy models, we examined effects of park size and connectivity on 
initial occupancy, colonization, and persistence of butterflies within parks. We quantified 
park size as the total green area within a park boundary and connectivity as the area-weighted 
distance from other parks within a 2 km buffer of the parks.

Use the "/run_model/run_model_m1_cmdstanr.R" file to conduct the analysis (dynamic occupancy model). The file will first call additional
R files in the same directory to gather and prepare the data (or alternatively load previously prepared data).
The file will then call cmdstanr program to fit the dynamic occupancy model.

Use the "/make_figures/" subdirectory to process and plot results from the model fit.

## Part 2: Influence of Additional Local and Landscape Park Properties on Butterfly Species Occupancy

Because our dynamic models estimated that overall probability of persistence was high and 
overall colonization probability was low, we opted to use simplified static occupancy models 
to test how additional local and landscape park properties contribute to patterns of butterfly 
occupancy. Using the simpler models helped us to incorporate the additional random effects predictors. 
We used iNaturalist detection data to quantify the diversity of flowering plant genera within parks, 
UrbanWatch's 1 meter resolution remote sensing data to calculate the proportion of tree cover within 
parks, and NLCD 30 meter resolution remote sensing data to calculate the total area of 
herbaceous and woody vegetation within a 2 km buffer of the parks. 

Use the "/run_model/run_model_2.1.R" file to conduct the primary analysis (static occupancy model with all 
covariates to estimate direct effects). The file will first call additional
R files in the same directory to gather and prepare the data (or alternatively load previously prepared data).
The file will then call cmdstanr program to fit the dynamic occupancy model.
Instead use the "/run_model/run_model_2.2.R" file to conduct fit a reduced static occupancy model 
that allowed us to estimate the total effects of some predictors that were hypothesized to
be mediated by other predictors.

Use the "/make_figures/" subdirectory to process and plot results from the model fit.

## Part 3: Impacts of Between-City Differences on Multi-scale Patterns of Butterfly Diversity

Cities showed differences in park and landscape features, such as the distribution of park sizes
or the total area of vegetation spread across the city and in how these features impact butterfly 
species occupancy. To understand how these between-city differences impact patterns of alpha, 
beta, and gamma diversity, we used our fitted occupancy model from Part 2 to 
reconstruct butterfly species occupancy in all parks in each of the 22 cities. We then calculated 
local to landscape scale metrics of city-wide diversity and quantified relationships between 
diversity and city-wide features including median park size, city-wide park connectivity, 
and city-wide area of tree or herbaceous vegetation cover. We accounted for latitude as a 
predictor of alpha and gamma diversity measures, and total 
city area as a predictor of gamma diversity. Given the limited sample size of 22 
independent cities, our final models included only the predictors that showed at least
marginal evidence of an association with the diversity outcome (i.e., at least the 50% BCI departing from 0).

Use the "predict_diversity.R" file to simulate occupancy in urban parks, and then derive 
various diversity metrics from the presence/absence occupancy simulations. The predict diversity
function propagates uncertainty from the occupancy model output (m2.1)

Use the "run_m3.[...].R" files to assess the relationships between city-wide variables 
and the key diversity metrics that we monitored:
- m3.1: species richness of the median park in a city 
- m3.2: average richness of disturbance or edge avoidant species in a community in a city
- m3.3: jaccard index of community dissimilarity among communities in a city
- m3.3 submodels: contributions of species nestedness and turnover to the jaccard index
- m3.4: total number of species occurring across all parks in a city.

## Data Preparation

We use pre-existing datasets from [iNaturalist Research-grade Observations](https://doi.org/10.15468/ab3s5x) 
(accessed via GBIF) for occurrence data, [Urbanwatch](https://urbanwatch.charlotte.edu/) for land cover 
data, and [ParkServe](https://www.tpl.org/parkserve) for park shapefiles. Scripts to download and wrangle 
these into intermediate products for downstream analysis can be found [here](https://github.com/jensculrich/urban_park_community_science_project/tree/main/script/01_data_wrangling).

The major intermediate products generated per city are:

1. A dataframe of iNaturalist Lepidoptera occurrence records joined to park boundaries, with park-level attributes (park size and land cover type areas) appended per record.
2. A dataframe of iNaturalist Lepidoptera occurrence records within a defined regional boundary.
3. A dataframe of iNaturalist flowering plant occurrence records joined to park boundaries, with park-level attributes (park size and land cover type areas) appended per record.
4. A dataframe of park-level isolation values.

The derived park site covariate data for all parks in all cities are located in a table here: 
"./data/derived_park_site_covariate_data/derived_park_site_covariate_data.csv" 

Across all cities, we also generate:

5. A dataframe of city-level IIC connectivity values.
6. A dataframe of city-level mean isolation values.
7. A dataframe of city-level land cover composition, including total area per cover type and land cover diversity (number of distinct cover types).
8. A dataframe of city-level mean and median park sizes.

The derived city covariate data are located in a table here: 
"./data/city_wide_data/derived_city_wide_data.csv"

## Lepidoptera Trait Data
Our analyses incorporated some butterfly species traits as predictors of occupancy and detection.
The data are contained in the "./data/lepidoptera_trait_data/" folder.

Ease of ID ("./ease_of_id/identifiability_by_genus.csv") was the only data collected from scratch for this project.
For these data we searched the number of research grade versus total detections (within the U.S.)
for each butterfly genus included in the study.

The Species list for traits file ("./lepidoptera_trait_data/SpeciesListForTraits.csv") was obtained from
Goldstein et al., 2024. We added any species missing for the trait data using the same methods provided. 
This is the file used to access species trait data for all traits other than ease of ID for the analyses.

To add data for species/traits missing from Goldstein et al., 2024, we used the following resources:
The migratory status was obtained from Chowhury et al., 2021.
Wingspan information was taken from the Leptraits V1.0 database.
Feature Diversity was manually parsed from ButterfliesandMothsofNorthAmerica.com (BAMONA website)