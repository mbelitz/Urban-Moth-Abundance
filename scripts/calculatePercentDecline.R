library(tidyverse)

# calculate total % declines

### rich
rich <- read.csv("data/data_products/richness.csv") %>% 
  rename(validName = richness) %>% 
  group_by(Site) %>% 
  summarise(richness = length(unique(validName)))

min(rich$richness) / max(rich$richness) 

(max(rich$richness) - min(rich$richness)) / max(rich$richness) * 100


######## adult
##macro
adult_df <- read.csv("data/rawMeasurements/surveyDateSheet.csv") %>% 
  group_by(location) %>% 
  summarise(totalMacro = sum(macroMoths),
            totalMicro = sum(microMoths))

(max(adult_df$totalMacro) - min(adult_df$totalMacro)) / max(adult_df$totalMacro) * 100

##micro
(6344 - 1277) / 6344 * 100


## frass
frass_df <- read.csv("data/data_products/frass_biomass_reweigh_doyDiff.csv") %>% 
  group_by(Site) %>% 
  reframe(totalFrass = sum(Mass, na.rm = T))

(11.9230 - 7.0990) / 7.0990 * 100
