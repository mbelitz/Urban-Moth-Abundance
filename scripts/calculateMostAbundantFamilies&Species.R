library(tidyverse)

# read in site by species matrix in long form
mdf <- read.csv("data/data_products/siteXspeciesMatrix.csv")


gb <- mdf %>% 
  group_by(validName) %>% 
  summarise(totalAbundance = sum(abundance))


filter(ungroup(gb), totalAbundance < 5)$validName %>% length()

# how many families?
moth_df <- read.csv("data/data_products/adultDataSet_validNames.csv") %>% 
  distinct(id, .keep_all = T) %>% 
  mutate(year = year(mdy(eventDate))) %>% 
  mutate(doy = if_else(
    condition = year == 2019,
    true = yday(mdy(eventDate)),
    false = 365 + yday(mdy(eventDate))
  ))

moth_df <- moth_df %>% 
  filter(Family != "cantID",
         Family != "retake_photo",
         Family != "rephotograph",
         Family != "retake_photos")

head(moth_df)

# families.
unique(moth_fam_df$Family)
c("Cossidae", "Crambidae", "Erebidae", "Geometridae", "Lasiocampidae",
  "Limocodidae", "Megalopygidae", "Mimallonidae", "Noctuidae", "Notodontidae", "Pyralidae",
  "Saturniidae", "Sphingidae", "Depressariidae", "Tineidae") # 15 total families

gb <- moth_fam_df %>% 
  group_by(Family) %>% 
  summarise(totalAbundance = sum(count))