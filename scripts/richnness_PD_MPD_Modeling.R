library(tidyverse)
library(brms)
library(ape)
library(picante)
library(cmdstanr)
set_cmdstan_path()
# read in mannually controlled richness file
rich <- read.csv("data/data_products/richness.csv") %>% 
  rename(validName = richness) %>% 
  group_by(Site) %>% 
  summarise(richness = length(unique(validName)))

# combine dataset with urbanization values
urb <- read.csv("data/data_products/urbanization_gradient.csv")

## get richness values
rich <- rich %>% 
  left_join(urb)

# make my one, well thought-out model
set.seed(seed = 5)
rich_dev <- brm(formula = bf(richness ~ Dev_1),
                 data = rich,
                 family = gaussian(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.93),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
)



# examine model assumptions
plot(rich_dev)
pp_check(rich_dev, ndraws = 100)
pp_check(rich_dev, type = "stat", stat = "mean")

summary(rich_dev, prob = 0.89)
rich_dev_sum <- summary(rich_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce <- conditional_effects(x = rich_dev, effects = "Dev_1", prob = 0.89)
ce_df <- ce$Dev_1

rich_plot <- ggplot() +
  geom_jitter(rich, mapping = aes(x = Dev_1, y = richness), alpha = 0.5) +
  geom_line(ce_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Richness") +
  scale_y_continuous(breaks = c(0,25,50,75,100,125,150,175))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

rich_plot  

# do richness, no baca
rich_noBaca <- rich %>% 
  filter(Site != "Baca")

set.seed(6)
rich_dev_noBaca <- brm(formula = bf(richness ~ Dev_1),
                data = rich_noBaca,
                family = gaussian(),
                chains = 4, iter = 2400, warmup = 1000,
                control = list(adapt_delta = 0.93),
                cores = 4, seed = 1234, 
                threads = threading(2),
                backend = "cmdstanr", 
)



# examine model assumptions
plot(rich_dev_noBaca)
pp_check(rich_dev_noBaca)
pp_check(rich_dev_noBaca, type = "stat", stat = "mean")

summary(rich_dev_noBaca, prob = 0.89)
rich_dev_sum_noBaca <- summary(rich_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_noBaca <- conditional_effects(x = rich_dev_noBaca, effects = "Dev_1", prob = 0.89)
ce_noBaca_df <- ce_noBaca$Dev_1

rich_noBaca_plot <- ggplot() +
  geom_jitter(rich_noBaca, mapping = aes(x = Dev_1, y = richness), alpha = 0.5) +
  geom_line(ce_noBaca_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_noBaca_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Richness") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

rich_noBaca_plot  

## PD Modeling
# function to capitalize spp names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# read in site by species matrix in long form
mdf <- read.csv("data/data_products/siteXspeciesMatrix.csv")

mdf <- mdf %>% 
  select(Site, validName, abundance)

# pd calculations
p <- pivot_wider(mdf, names_from = validName, values_from = abundance) 
p2 <- p %>% tibble::column_to_rownames(var = "Site")
p2

tree <- read.tree("data/phylogeny/insect_tree_wBranches_allSpecies.tre")
tree$tip.label
tree$tip.label <- stringr::str_replace(tree$tip.label, pattern = "_", " ")
tree$tip.label <- firstup(tree$tip.label)

tree <- prune.sample(p2, tree)

comm <-p[,tree$tip.label]
comm

plot(tree)
for (i in rownames(comm)) {
  plot(tree, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(comm[i, ] > 0), pch = 19, cex = 2, col ="red")
  legend("topleft", i, bty = "n")
  
}

# get vectors to join to data.frames
Site <- c(rownames(p2), "totalPhy")

#PD
moth.pd <- pd(comm, tree, include.root = TRUE)
moth.pd

#make a row in the matrix with each species so we can caclulate proportional pd
r <- rep(1,104)
comm_pdp <- rbind(comm,r)
moth.pdp <- pd(comm_pdp, tree)
moth.pdp
moth.pdp <- cbind(moth.pdp, Site)
moth.pdp <- moth.pdp %>% 
  mutate(pd_p = PD/5441.055)

urb <- read.csv("data/data_products/urbanization_gradient.csv")

moth.pdp <- moth.pdp %>% 
  left_join(urb) %>% 
  na.omit()

# make my one, well thought-out model
set.seed(7)
pdp_dev <- brm(formula = bf(PD ~ Dev_1),
                data = moth.pdp,
                family = gaussian(),
                chains = 4, iter = 2400, warmup = 1000,
                control = list(adapt_delta = 0.93),
                cores = 4, seed = 1234, 
                threads = threading(2),
                backend = "cmdstanr", 
)

# examine model assumptions
plot(pdp_dev)
pp_check(pdp_dev)
pp_check(pdp_dev, type = "stat", stat = "mean")

summary(pdp_dev, prob = 0.89)

pdp_dev_sum <- summary(pdp_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_pdp <- conditional_effects(x = pdp_dev, effects = "Dev_1", prob = 0.89)
ce_pdp_df <- ce_pdp$Dev_1

pdp_plot <- ggplot() +
  geom_jitter(moth.pdp, mapping = aes(x = Dev_1, y = PD), alpha = 0.5) +
  geom_line(ce_pdp_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_pdp_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Phylogenetic \n diversity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

pdp_plot  

# do pdp, no baca
pdp_noBaca <- moth.pdp %>% 
  filter(Site != "Baca")

set.seed(8)
pdp_dev_noBaca <- brm(formula = bf(PD ~ Dev_1),
                       data = pdp_noBaca,
                       family = gaussian(),
                       chains = 4, iter = 2400, warmup = 1000,
                       control = list(adapt_delta = 0.93),
                       cores = 4, seed = 1234, 
                       threads = threading(2),
                       backend = "cmdstanr", 
)

# examine model assumptions
plot(pdp_dev_noBaca)
pp_check(pdp_dev_noBaca)
pp_check(pdp_dev_noBaca, type = "stat", stat = "mean")

summary(pdp_dev_noBaca, prob = 0.89)
pdp_dev_sum_noBaca <- summary(pdp_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_pdp_noBaca <- conditional_effects(x = pdp_dev_noBaca, effects = "Dev_1", prob = 0.89)
ce_pdp_noBaca_df <- ce_pdp_noBaca$Dev_1

pdp_noBaca_plot <- ggplot() +
  geom_jitter(pdp_noBaca, mapping = aes(x = Dev_1, y = PD), alpha = 0.5) +
  geom_line(ce_pdp_noBaca_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_pdp_noBaca_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "PD") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

pdp_noBaca_plot  

###############################################################################
###############################################################################
###############################MPD#############################################
###############################################################################
###############################################################################
# MPD modeling
Site <- c(rownames(p2))
phy.dist<-cophenetic(tree)
phydist <- as.matrix(phy.dist)
comm.m <- as.matrix(comm)
moth.mpd.aw <-mpd(comm.m, phy.dist, abundance.weighted = TRUE)
moth.mpd.aw

moth.mpd <-mpd(comm.m, phy.dist, abundance.weighted = FALSE)
moth.mpd

mpd <- data.frame(Site = Site,
                  mpd.raw = moth.mpd, 
                  abundance.weighted = moth.mpd.aw)

mpd <- mpd %>% 
  left_join(urb)

mpd

# make model
set.seed(9)
mpd_dev <- brm(formula = bf(abundance.weighted ~ Dev_1),
               data = mpd,
               family = gaussian(),
               chains = 4, iter = 2400, warmup = 1000,
               control = list(adapt_delta = 0.93),
               cores = 4, seed = 1234, 
               threads = threading(2),
               backend = "cmdstanr", 
)

# examine model assumptions
plot(mpd_dev)
pp_check(mpd_dev)
pp_check(mpd_dev, type = "stat", stat = "mean")

summary(mpd_dev, prob = 0.89)
mpd_dev_sum <- summary(mpd_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_mpd <- conditional_effects(x = mpd_dev, effects = "Dev_1", prob = 0.89)
ce_mpd_df <- ce_mpd$Dev_1

mpd_plot <- ggplot() +
  geom_jitter(mpd, mapping = aes(x = Dev_1, y = abundance.weighted), alpha = 0.5) +
  geom_line(ce_mpd_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_mpd_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Mean pairwise \n distance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

mpd_plot  

# remove baca
mpd_noBaca <- mpd %>% 
  filter(Site != "Baca")

# make model
set.seed(10)
mpd_dev_noBaca <- brm(formula = bf(abundance.weighted ~ Dev_1),
               data = mpd_noBaca,
               family = gaussian(),
               chains = 4, iter = 2400, warmup = 1000,
               control = list(adapt_delta = 0.93),
               cores = 4, seed = 1234, 
               threads = threading(2),
               backend = "cmdstanr", 
)

# examine model assumptions
plot(mpd_dev_noBaca)
pp_check(mpd_dev_noBaca)
pp_check(mpd_dev_noBaca, type = "stat", stat = "mean")

summary(mpd_dev_noBaca, prob = 0.89)
mpd_dev_sum_noBaca <- summary(mpd_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_mpd_noBaca <- conditional_effects(x = mpd_dev_noBaca, effects = "Dev_1", prob = 0.89)
ce_mpd_noBaca_df <- ce_mpd$Dev_1

mpd_plot_noBaca <- ggplot() +
  geom_jitter(mpd_noBaca, mapping = aes(x = Dev_1, y = abundance.weighted), alpha = 0.5) +
  geom_line(ce_mpd_noBaca_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_mpd_noBaca_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Mean pairwise \n distance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

mpd_plot_noBaca  


# plot modeling outputs
cp <- cowplot::plot_grid(rich_plot, pdp_plot, mpd_plot, 
                         labels = c("A", "B", "C"),
                         nrow = 3, ncol = 1)

if(dir.exists("figOutputs")){
  ggsave(plot = cp, filename = "figOutputs/rich_PD_MDP.png", dpi = 500,
         width = 4, height = 6)
} else{
  dir.create("figOutputs")
  ggsave(plot = cp, filename = "figOutputs/rich_PD_MPD.png", dpi = 500,
         width = 4, height = 6)
}

# Table outputs
rich_dev_sum_noBaca <- rich_dev_sum_noBaca %>% 
  mutate(model = "Richness No Baca")
rich_dev_sum <- rich_dev_sum %>% 
  mutate(model = "Richness")

pdp_dev_sum <- pdp_dev_sum %>% 
  mutate(model = "PDP")
pdp_dev_sum_noBaca <-  pdp_dev_sum_noBaca %>% 
  mutate(model = "PDP No Baca")

mpd_dev_sum <- mpd_dev_sum %>% 
  mutate(model = "MPD")
mpd_dev_sum_noBaca <- mpd_dev_sum_noBaca %>% 
  mutate(model = "MPD No Baca")


allModels <- bind_rows(rich_dev_sum, rich_dev_sum_noBaca,
                       pdp_dev_sum, pdp_dev_sum_noBaca,
                       mpd_dev_sum, mpd_dev_sum_noBaca)

write.csv(allModels, "tabOutputs/Richness_PD_MPD_ModelResults.csv", row.names = F)
