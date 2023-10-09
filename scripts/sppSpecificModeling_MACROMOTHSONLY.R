library(tidyverse)
library(brms)
library(ape)
# read in site by species matrix in long form
mdf <- read.csv("data/data_products//siteXspeciesMatrix.csv")
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

mdf_scaled <- mdf %>% 
  mutate(Dev_1 = scale(Dev_1),
         Dev_1 = scale(Dev_1),
         bio1_sd = scale(bio1_sd),
         bio1_mean = scale(bio1_mean),
         meanLight = scale(log(meanLight + 0.01)),
         mean_temp = scale((mean_temp - 1.02) * -1),
         totalLength = scale(totalLength))

## cool let's also look at the subset of data for which we have hsp data
mdf_scaled_hsp <- mdf_scaled %>% 
  filter(!is.na(hsp))

mdf_scaled_hsp$hsp <- as.factor(mdf_scaled_hsp$hsp)

# function to capitalize spp names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tt <- read.tree("data/phylogeny/insect_tree_wBranches_allSpecies.tre")
tt$tip.label <- stringr::str_replace(tt$tip.label, pattern = "_", " ")
tt$tip.label <- firstup(tt$tip.label)


mdf_phylo <- left_join(mdf_scaled_hsp, data.frame(validName = tt$tip.label, Phylo = "Yes"))

sppNotInAnlysis <- data.frame(validName = tt$tip.label) %>% 
  filter(!validName %in% mdf_phylo$validName)

mdf_phylo <- mdf_phylo %>% 
  filter(!is.na(Phylo))

tt <- ape::drop.tip(tt, tip = sppNotInAnlysis$validName)
A <- ape::vcv.phylo(tt)

# what families are in analysis?
moth_df <- read.csv("data/data_products/adultDataSet_validNames.csv") %>% 
  distinct(id, .keep_all = T) %>% 
  mutate(year = year(mdy(eventDate))) %>% 
  mutate(doy = if_else(
    condition = year == 2019,
    true = yday(mdy(eventDate)),
    false = 365 + yday(mdy(eventDate))
  )) %>% 
  distinct(validName, Family)

# what moth validNames have multipe Families?
familyCount <- moth_df %>% 
  group_by(validName) %>% 
  summarise(count = n())

moth_df <- moth_df %>% 
  distinct(validName, .keep_all = TRUE)

mdf_phylo_fam <- left_join(mdf_phylo, moth_df, by = c("validName" = "validName"))

unique(mdf_phylo_fam$Family)

# what should be kept, based on 
# https://royalsocietypublishing.org/doi/full/10.1098/rspb.2019.3014

#' Cossidae, Erebidae, Geometridae, Lasiocampidae, Limacodidae, Noctuidae, 
#' Mimallonidae, Notodontidae, Saturniidae, Sphingidae
#' 
#' Removed, Crambidae, Megalopygidae, Tineidae

# both filterings work
#t1 <- mdf_phylo_fam %>% 
#  filter(Family %in% c("Cossidae", "Erebidae",    
#                       "Geometridae",   "Lasiocampidae", "Limacodidae",   
#                       "Mimallonidae",  "Noctuidae",     
#                       "Notodontidae",  "Saturniidae",   "Sphingidae"))
#t2 <- mdf_phylo_fam %>% 
#  filter(!Family %in% c("Crambidae", "Megalopygidae", "Tineidae"))

mdf_phylo <- mdf_phylo_fam %>% 
  filter(Family %in% c("Cossidae", "Erebidae",    
                       "Geometridae",   "Lasiocampidae", "Limacodidae",   
                       "Mimallonidae",  "Noctuidae",     
                       "Notodontidae",  "Saturniidae",   "Sphingidae"))

# fit new spp-specific model with only macro moth data
library(cmdstanr)
set_cmdstan_path()

sppSpecificDev <- brm(formula = bf(abundance ~ Dev_1  + mean_temp + 
                                     totalLength + bio1_mean + hsp + 
                                     Dev_1:totalLength +
                                     mean_temp:bio1_mean + 
                                     Dev_1:hsp +
                                     (1|validName),
                                   zi ~ Dev_1 + hsp + totalLength),
                      data = mdf_phylo,
                      family = zero_inflated_negbinomial(),
                      chains = 4, iter = 2400, warmup = 1000,
                      control = list(adapt_delta = 0.95),
                      cores = 4, seed = 1234, 
                      threads = threading(2),
                      backend = "cmdstanr", 
)

# examine model assumptions
plot(sppSpecificDev)
pp_check(sppSpecificDev)
pp_check(sppSpecificDev, type = "stat", stat = "mean")

summary(sppSpecificDev, prob = 0.89)

plot(conditional_effects(sppSpecificDev, terms = "mean_temp:bio1_mean"))


mdf_phylo <- mdf_phylo %>% 
  mutate(phyloSpp = validName)


# Add phylogenetic term
sppSpecificDev_phylo <- brm(formula = bf(abundance ~ Dev_1  + mean_temp + 
                                           totalLength + bio1_mean + hsp + 
                                           Dev_1:totalLength +
                                           mean_temp:bio1_mean + 
                                           Dev_1:hsp +
                                           (1|validName) + (1|gr(phyloSpp, cov = A)),
                                         zi ~ Dev_1 + hsp + totalLength),
                            data = mdf_phylo,
                            data2 = list(A = A),
                            family = zero_inflated_negbinomial(),
                            chains = 4, iter = 2400, warmup = 1000,
                            control = list(adapt_delta = 0.95),
                            cores = 4, seed = 1234, 
                            threads = threading(2),
                            backend = "cmdstanr", 
)

# examine model assumptions
plot(sppSpecificDev_phylo)
pp_check(sppSpecificDev_phylo)
pp_check(sppSpecificDev_phylo, type = "stat", stat = "mean")

summary(sppSpecificDev_phylo)

# plot intercations
ce_tl <- conditional_effects(x = sppSpecificDev_phylo, effects = "Dev_1:totalLength", prob = 0.89, method = "fitted")
ce_tl_df <- ce_tl$`Dev_1:totalLength`

a <- ggplot() +
  geom_line(ce_tl_df, mapping = aes(x = Dev_1, y = estimate__, color = effect2__), size = 1.2) +
  geom_ribbon(ce_tl_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__, fill = effect2__), alpha = 0.3) +
  geom_path(ce_tl_df, mapping = aes(x = Dev_1, y =  lower__, color = effect2__), size = 0.25, linetype = 2) +
  geom_path(ce_tl_df, mapping = aes(x = Dev_1, y = upper__, color = effect2__), size = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Urban development", y = "Abundance",
       color = "Body size", fill = "Body size") +
  scale_color_manual(values = c("#9e2a2b","#e09f3e","#223d44"),
                     labels = c("1.08 (Large)", "0.04 (Average)", "-1 (Small)")) +
  scale_fill_manual(values = c("#9e2a2b","#e09f3e","#223d44"),
                    labels = c("1.08 (Large)", "0.04 (Average)", "-1 (Small)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(legend.position = "bottom")


# hsp interaction plot
ce_hsp <- conditional_effects(x = sppSpecificDev_phylo, effects = "Dev_1:hsp", prob = 0.89)
ce_hsp_df <- ce_hsp$`Dev_1:hsp`

b <- ggplot(ce_hsp_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), size = 1.05) +
  geom_path(mapping = aes(y =  lower__, color = effect2__), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = upper__, color = effect2__), size = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#D89A9E","#2E4057","#519872"))  +
  scale_fill_manual(values = c("#D89A9E", "#2E4057", "#519872")) +
  labs(y = "Abundance", x = "Urban development", 
       color = "Host plant \n specialization", fill = "Host plant \n specialization") +
  theme_classic() +
  theme(legend.position = "bottom")

# tempNiche interaction plot
ce_tn <- conditional_effects(x = sppSpecificDev_phylo, effects = "mean_temp:bio1_mean", prob = 0.89)
ce_tn_df <- ce_tn$`mean_temp:bio1_mean`

c <- ggplot(ce_tn_df, mapping = aes(x = mean_temp, y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), size = 1.05) +
  geom_path(mapping = aes(y =  lower__, color = effect2__), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = upper__, color = effect2__), size = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c( "#E76F51", "#D5A220","#2A9D8F"),
                    labels = c("0.99 (Warm-adapted)", "-0.02 (Average)", "-1.04 (Cold-adapted)")) +
  scale_color_manual(values = c( "#E76F51", "#D5A220","#2A9D8F"),
                     labels = c("0.99 (Warm-adapted)", "-0.02 (Average)", "-1.04 (Cold-adapted)"))  +
  labs(y = "Abundance", x = "Relative temperature of site", 
       color = "Temperature niche", fill = "Temperature niche") +
  theme_classic() +
  theme(legend.position = "bottom")

library(ggpubr)
cpp <- ggarrange(
  a + theme(legend.position = c(.83,.85), legend.background = element_blank()),
  b + theme(legend.position = c(.83,.85)),
  c + theme(legend.position = c(.83,.9)),
  ncol = 1, labels = LETTERS
)
ggsave(cpp, filename = "figOutputs/sppSpecificAbundance_MACROMOTHSONLY.png", height = 10, width = 6)
