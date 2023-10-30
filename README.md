Code and data used to generate results and figures of "Substantial urbanization-driven declines of larval and adult moths in a subtropical environment" manuscript.
If you download or clone this repository, you can fully replicate results and figures. We also release the raw datasets for light and frass sampling. Scored trait data can be downloaded in the data/traits subdirectory. 

## A brief overview of what is found in each directory and subdirectory

## data
Contains data that was used to estimate phenology metrics and data used in modeling framework.

### data_products
Data used in analyses that were the product of data cleaning or simple analyses.  
- lightData.csv is the average amount of artificial light in lux within a 1-km radius of each site. 
- lunarIllumination_lunarRPackage.csv is the lunar illumination for each date in 2019 and 2020.
- richness.csv provides a list with unique moth species found at each site.
- tempgradient.csv is how much cooler a site is relative to the most urban site (BACA). 
- traits.csv is the associated life history traits for each species in our analysis.
- urbanization_gradient is the proportion of land cover within 1-km of each site that is developed. 

### phylogeny
This subdirectory contains a single .tre file that is a subtree representing the identified macro-moths in our study. 

### rawMeasurements
Unmanipulated data that represents raw field collected samples. 
- adultDataSet_validNames.csv is the associated photo number and lowest taxonomic identification of each collected macro-moth.
- frass_biomass.csv is the amount of frass in mg captured per trap per sampling event.
- surveyDateSheet.csv is the total number of macro-moths and micro-moths captered per light trap sampling event.

### weatherData
Daily weather data for study region for 2019 and 2020. 

## scripts
Contains R scripts used to fit models and generate figures. 

### analyses
