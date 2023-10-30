library(tidyverse)
library(brms)
library(cmdstanr)
set_cmdstan_path()

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

moth_fam_df <- moth_df %>% 
  group_by(Site, Family) %>% 
  summarise(count = n())

head(moth_fam_df)

# how many sites are these families observed at?
sites <- moth_fam_df %>% 
  group_by(Family) %>% 
  summarise(nSites = n())

# remove families observed at less than 5 sites
removeSites <- sites %>% 
  filter(nSites < 6)

moth_fam_df <- moth_fam_df %>% 
  filter(!Family %in% removeSites$Family)

# now impute sites with zeros
mdf <- ungroup(moth_fam_df) %>% 
  tidyr::complete(Site,Family,
                  fill = list(count = 0))


# combine dataset with urbanization values
urb <- read.csv("data/data_products/urbanization_gradient.csv")

mdf <- mdf %>% 
  left_join(urb) 

fam_dev <- brm(formula = bf(count ~ Dev_1 +
                              (1 + Dev_1|Family)),
               data = mdf,
               family = negbinomial(),
               chains = 4, iter = 2400, warmup = 1000,
               control = list(adapt_delta = 0.93),
               cores = 4, seed = 1234, 
               threads = threading(2),
               backend = "cmdstanr", 
)

# examine model assumptions
plot(fam_dev)
pp_check(fam_dev)
pp_check(fam_dev, type = "stat", stat = "mean")

library(tidybayes)

fam_dev %>%
  spread_draws(b_Dev_1, r_Family[Family,Intercept]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_Dev_1 + r_Family) %>%
  ungroup() %>%
  mutate(Family = str_replace_all(Family, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(Family, mu))) +
  geom_vline(xintercept = fixef(fam_dev)[2, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(fam_dev)[2, 3:4], color = "#839496", linetype = 2) +
  stat_halfeye(.width = .5, size = 2/3, fill = "#859900") +
  labs(x = expression("Sensitivity to urban development"),
       y = "") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 


fam_dev <- brm(formula = bf(count ~ Dev_1 * Family),
                 data = mdf,
                 family = negbinomial(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.93),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
)

# examine model assumptions
plot(fam_dev)
pp_check(fam_dev)
pp_check(fam_dev, type = "stat", stat = "mean")

summary(fam_dev, prob = 0.89)
fam_dev_sum <- summary(fam_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

plot(conditional_effects(fam_dev, effects = "Dev_1:Family", prob = 0.89)) +
  theme_classic()

# make plot
ce_fam <- conditional_effects(x = fam_dev, effects = "Dev_1:Family", prob = 0.89)
ce_fam_df <- ce_fam$`Dev_1:Family`

fam_plot <- ggplot(ce_fam_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), size = 1.05, alpha = 0.5) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Abundance", x = "Urban development", 
       color = "Family", fill = "Family") +
  scale_color_discrete() +
  scale_fill_discrete() +
  theme_classic() +
  theme(legend.position = "right")

fam_plot

fam_plot_log <- ggplot(ce_fam_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(mapping = aes(ymin = lower__, ymax = upper__, fill = effect2__), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = effect2__), size = 1.05, alpha = 0.5) +
  scale_y_log10(expand = c(0,0),
                     breaks = c(0.1,1,10,100),
                     labels = c(0.1,1,10,100)) +
  labs(y = "Abundance", x = "Urban development", 
       color = "Family", fill = "Family") +
  scale_color_discrete() +
  scale_fill_discrete() +
  theme_classic() +
  theme(legend.position = "right")

fam_plot_log

library(ggpubr)

gp <- ggarrange(fam_plot, fam_plot_log, nrow = 2, ncol = 1, common.legend = TRUE, legend = "right")

ggsave(plot = gp, filename = "figOutputs/FamilyResponses.png", dpi = 400, width = 6, height = 6)

# run model without baca
mdf_noBaca <- filter(mdf, Site != "Baca")

fam_dev_noBaca <- brm(formula = bf(count ~ Dev_1 * Family),
               data = mdf_noBaca,
               family = negbinomial(),
               chains = 4, iter = 2400, warmup = 1000,
               control = list(adapt_delta = 0.93),
               cores = 4, seed = 1234, 
               threads = threading(2),
               backend = "cmdstanr", 
)

summary(fam_dev_noBaca, prob = 0.89)
fam_dev_sum_noBaca <- summary(fam_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

## tabOutputs
fam_dev_sum <- fam_dev_sum %>% 
  mutate(data = "Full")

fam_dev_sum_noBaca <- fam_dev_sum_noBaca %>% 
  mutate(data = "No Baca")

write.csv(fam_dev_sum, file = "tabOutputs/FamilyResults.csv", row.names = F)
write.csv(fam_dev_sum_noBaca, file = "tabOutputs/FamilyResults_NoBaca.csv", row.names = F)
