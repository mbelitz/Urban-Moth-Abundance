library(tidyverse)

# how many total morpho species?
# read in mannually controlled richness file
rich <- read.csv("data/data_products/richness.csv")
unique(rich$richness) %>% length()

# now examine most abundant species
# read in species-specific data
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

abund_df <- moth_df %>% 
  group_by(Site, validName) %>% 
  summarise(abundance = n())

abund_df <- ungroup(abund_df) %>% 
  mutate(genus = word(validName, 1 ,1),
         species = word(validName,2,2))

abund_df <- abund_df %>% 
  filter(species != "NA",
         species != "",
         genus != "Datana")

abund_df <- abund_df %>% 
  complete(validName,Site, fill = list(abundance = 0))

# rename as model data frame
mdf <- abund_df
#read in traits
traits <- read.csv("data/data_products//traits.csv") %>% 
  select(-species, -max_lat, -min_lat, -med_lat, -geoPF, -notes.1, -X)

mdf <- left_join(mdf,traits)

mdf <- mdf %>% 
  mutate(hsp = case_when(
    dietBreadth == "multiFamily" | dietBreadth == "detritus" | 
      dietBreadth == "Bracket Fungi" ~ 1,
    dietBreadth == "Family" ~ 2,
    dietBreadth == "Genus" | dietBreadth == "genus" | 
      dietBreadth == "Species" | dietBreadth == "species" ~ 3,
    dietBreadth == "unk" ~ 4,
    is.na(dietBreadth) ~ 5
  )) %>% 
  mutate(
    totalLength = if_else(notes == "WS",
                          true = maxWingspan/2,
                          false = as.double(maxWingspan))) %>% 
  filter(totalLength >= 10)

mdf <- mdf %>% 
  mutate(hsp = if_else(condition = hsp == 4, 
                       true = 5, false = hsp)) %>% 
  mutate(hsp = na_if(hsp, 5))

mdf <- filter(mdf, !is.na(hsp))

urb <- read.csv("data/data_products/urbanization_gradient.csv")

mdf <- mdf %>% left_join(urb)  

meanDev <- mean(mdf$Dev_1)
sdDev <- sd(mdf$Dev_1)

mdf_scaled_hsp <- mdf %>% 
  mutate(Dev_1 = scale(Dev_1),
         Dev_1 = scale(Dev_1),
         bio1_sd = scale(bio1_sd),
         bio1_mean = scale(bio1_mean),
         totalLength = scale(totalLength))

# opening results paragraph
# abundance for identified species
gb <- mdf_scaled_hsp %>% 
  group_by(validName) %>% 
  summarise(totalAbundance = sum(abundance))

filter(ungroup(gb), totalAbundance < 2)$validName %>% length()

write.csv(x = gb, file = "tabOutputs/identifiedMacroMothsAbundance.csv", row.names = F)

# how many families?
moth_df <- moth_df %>% 
  filter(Family != "cantID",
         Family != "retake_photo",
         Family != "rephotograph",
         Family != "retake_photos")

head(moth_df)

# families.
unique(moth_df$Family)
fams <- c("Cossidae", "Crambidae", "Erebidae", "Geometridae", "Lasiocampidae",
  "Limacodidae", "Megalopygidae", "Mimallonidae", "Noctuidae", "Notodontidae", "Pyralidae",
  "Saturniidae", "Sphingidae", "Depressariidae", "Tineidae") # 15 total families

moth_fam_df <- filter(moth_df, Family %in% fams)

gb2 <- moth_fam_df %>% 
  group_by(Family) %>% 
  summarise(totalAbundance = n()) %>% 
  ungroup()

write.csv(x = gb2, file="tabOutputs/abundanceByFamily.csv", row.names = F)
