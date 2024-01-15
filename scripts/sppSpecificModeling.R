library(tidyverse)
library(brms)
library(ape)
library(cmdstanr)
set_cmdstan_path()

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

## cool let's also look at the subset of data for which we have hsp data
mdf_scaled_hsp$hsp <- as.factor(mdf_scaled_hsp$hsp)

# run non phylogenetic model
# run model
slopePrior <- prior(normal(0,10), class = b)

set.seed(124)
sppSpecificDev <- brm(formula = bf(abundance ~ Dev_1 + 
                                     totalLength + bio1_mean + hsp + 
                                     Dev_1:totalLength +
                                     Dev_1:bio1_mean + 
                                     Dev_1:hsp +
                                     (1|validName),
                                   zi ~ Dev_1 + hsp + totalLength),
                      data = mdf_scaled_hsp,
                      family = zero_inflated_negbinomial(),
                      chains = 4, iter = 2400, warmup = 1000,
                      control = list(adapt_delta = 0.99),
                      cores = 4, seed = 1234, 
                      threads = threading(2),
                      backend = "cmdstanr", 
                      prior = slopePrior
)


# examine model assumptions
plot(sppSpecificDev)
pp_check(sppSpecificDev, ndraws = 100)
pp_check(sppSpecificDev, type = "stat", stat = "mean")

summary(sppSpecificDev, prob = 0.89)

plot(conditional_effects(sppSpecificDev, terms = "mean_temp:bio1_mean"))

m_sum_noPhy <- summary(sppSpecificDev, prob = 0.89) 
fixed_noPhy <- m_sum_noPhy$fixed %>% 
  tibble::rownames_to_column()
random_species_noPhy <- m_sum_noPhy$random$validName %>% 
  tibble::rownames_to_column()

m_sum_allSites_noPhy <- bind_rows(fixed_noPhy, random_species_noPhy)

# tabOutputs of non phylogenetic model
write.csv(m_sum_allSites_noPhy, "tabOutputs/speciesSpecificResults_noPhylogeney.csv", row.names = F)

##### Phylogenetic model

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

# Add phylogenetic term
mdf_phylo <- mdf_phylo %>% 
  mutate(phyloSpp = validName)

set.seed(1234)
sppSpecificDev_phylo <- brm(formula = bf(abundance ~ Dev_1   + 
                                           totalLength + bio1_mean + hsp + 
                                           Dev_1:totalLength +
                                           Dev_1:bio1_mean + 
                                           Dev_1:hsp +
                                           (1|validName) + (1|gr(phyloSpp, cov = A)),
                                         zi ~ Dev_1 + hsp + totalLength),
                            data = mdf_phylo,
                            data2 = list(A = A),
                            family = zero_inflated_negbinomial(),
                            chains = 4, iter = 2400, warmup = 1000,
                            control = list(adapt_delta = 0.99),
                            cores = 4, seed = 1234, 
                            threads = threading(2),
                            backend = "cmdstanr", 
                            prior = slopePrior
)

# examine model assumptions
plot(sppSpecificDev_phylo)
pp_check(sppSpecificDev_phylo, ndraws = 100)
pp_check(sppSpecificDev_phylo, type = "stat", stat = "mean")

summary(sppSpecificDev_phylo, prob = 0.89)
m_sum <- summary(sppSpecificDev_phylo, prob = 0.89) 
fixed <- m_sum$fixed %>% 
  tibble::rownames_to_column()
random_phy <- m_sum$random$phyloSpp %>% 
  tibble::rownames_to_column()
random_species <- m_sum$random$validName %>% 
  tibble::rownames_to_column()

m_sum_allSites <- bind_rows(fixed, random_phy, random_species)

# plot intercations
ce_tl <- conditional_effects(x = sppSpecificDev_phylo, 
                             effects = "Dev_1:totalLength", 
                             prob = 0.89, method = "fitted", int_conditions = list(bodyLength = c(-1,0,1)))
ce_tl_df <- ce_tl$`Dev_1:totalLength`
head(ce_tl_df)

# pred_df <- expand.grid(totalLength = c(-1,0,1), Dev_1 = seq(-1,2,by= 0.01), phyloSpp = NA, validName = NA,
#                        mean_temp = 0, bio1_mean = 0, hsp = 1)
# out <- predict(sppSpecificDev_phylo, pred_df, prob = c(.055, .945))
# head(out)
# out2 <- data.frame(out, pred_df)
# 
# ggplot(data = out2) +
#   aes(x = Dev_1, color = factor(totalLength), fill = factor(totalLength)) +
#   geom_ribbon(aes(ymin = Q5.5, ymax = Q94.5), alpha = 0.2) +
#   geom_line(aes(y=Estimate))

a <- ggplot() +
  geom_line(ce_tl_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__, color = effect2__), linewidth = 1.2) +
  geom_ribbon(ce_tl_df, mapping = aes(x = meanDev + (sdDev*Dev_1), ymax = upper__, ymin = lower__, fill = effect2__), alpha = 0.3) +
  geom_path(ce_tl_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y =  lower__, color = effect2__), linewidth = 0.25, linetype = 2) +
  geom_path(ce_tl_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = upper__, color = effect2__), linewidth = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Urban development", y = "Abundance",
       color = "Body size", fill = "Body size") +
  scale_color_manual(values = c("#9e2a2b","#e09f3e","#223d44"),
                     labels = c("1.04 (Large)", "0.01 (Average)", "-1.03 (Small)")) +
  scale_fill_manual(values = c("#9e2a2b","#e09f3e","#223d44"),
                    labels = c("1.04 (Large)", "0.01 (Average)", "-1.03 (Small)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  theme(legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title = element_text(size = 14))


# hsp interaction plot
ce_hsp <- conditional_effects(x = sppSpecificDev_phylo, effects = "Dev_1:hsp", prob = 0.89)
ce_hsp_df <- ce_hsp$`Dev_1:hsp`

b <- ggplot(ce_hsp_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), linewidth = 1.05) +
  geom_path(mapping = aes(y =  lower__, color = effect2__), linewidth = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = upper__, color = effect2__), linewidth = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c("#D89A9E","#2E4057","#519872"), 
                     labels = c("Multi-family", "Family", "Genus/Species"))  +
  scale_fill_manual(values = c("#D89A9E", "#2E4057", "#519872"), 
                    labels = c("Multi-family", "Family", "Genus/Species")) +
  labs(y = "Abundance", x = "Urban development", 
       color = "Host plant \n specialization", fill = "Host plant \n specialization") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title = element_text(size = 14))

# tempNiche interaction plot
ce_tn <- conditional_effects(x = sppSpecificDev_phylo, effects = "Dev_1:bio1_mean", prob = 0.89)
ce_tn_df <- ce_tn$`Dev_1:bio1_mean`

c <- ggplot(ce_tn_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), size = 1.05) +
  geom_path(mapping = aes(y =  lower__, color = effect2__), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = upper__, color = effect2__), size = 0.25, linetype = 2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c( "#E76F51", "#D5A220","#2A9D8F"),
                    labels = c("0.99 (Warm-adapted)", "-0.01 (Average)", "-1 (Cold-adapted)")) +
  scale_color_manual(values = c( "#E76F51", "#D5A220","#2A9D8F"),
                     labels = c("0.99 (Warm-adapted)", "-0.01 (Average)", "-1 (Cold-adapted)"))  +
  labs(y = "Abundance", x = "Urban development", 
       color = "Temperature niche", fill = "Temperature niche") +
  theme_classic() +
  theme(legend.position = "bottom", 
        axis.text=element_text(size=12),
        axis.title = element_text(size = 14))

library(ggpubr)
cpp <- ggarrange(
  a + theme(legend.position = c(.83,.85), legend.background = element_blank()),
  b + theme(legend.position = c(.83,.85)), 
  c + theme(legend.position = c(.83,.85), legend.background = element_blank())
  ,
  ncol = 1, labels = LETTERS
)
ggsave(cpp, filename = "figOutputs/sppSpecificAbundance2_tempNicheIncluded.png",
       height = 11, width = 6)


# make results table for no BACA
mdf_phylo_noBaca <- mdf_phylo %>% 
  filter(Site != "Baca")


set.seed(50)
sppSpecificDev_phylo_noBaca <- brm(formula = bf(abundance ~ Dev_1   + 
                                                  totalLength + bio1_mean + hsp + 
                                                  Dev_1:totalLength +
                                                  Dev_1:bio1_mean + 
                                                  Dev_1:hsp +
                                                  (1|validName) + (1|gr(phyloSpp, cov = A)),
                                                zi ~ Dev_1 + hsp + totalLength),
                                   data = mdf_phylo_noBaca,
                                   data2 = list(A = A),
                                   family = zero_inflated_negbinomial(),
                                   chains = 4, iter = 2400, warmup = 1000,
                                   control = list(adapt_delta = 0.99),
                                   cores = 4, seed = 1234, 
                                   threads = threading(2),
                                   backend = "cmdstanr", 
                                   prior = slopePrior
)

# examine model assumptions
plot(sppSpecificDev_phylo_noBaca)
pp_check(sppSpecificDev_phylo_noBaca)
pp_check(sppSpecificDev_phylo_noBaca, type = "stat", stat = "mean")

summary(sppSpecificDev_phylo_noBaca, prob = 0.89)
m_sum_noBaca <- summary(sppSpecificDev_phylo_noBaca, prob = 0.89) 
fixed_noBaca <- m_sum_noBaca$fixed %>% 
  tibble::rownames_to_column()
random_phy_noBaca <- m_sum_noBaca$random$phyloSpp %>% 
  tibble::rownames_to_column()
random_species_noBaca <- m_sum_noBaca$random$validName %>% 
  tibble::rownames_to_column()

m_sum_allSites_noBaca <- bind_rows(fixed_noBaca, random_phy_noBaca, random_species_noBaca)

# tabOutputs
write.csv(m_sum_allSites, "tabOutputs/speciesSpecificResults.csv", row.names = F)
write.csv(m_sum_allSites_noBaca, "tabOutputs/speciesSpecificResults_noBaca.csv", row.names = F)
