library(dplyr)
library(ggplot2)
library(lubridate)
library(brms)

# read in adult data
adult_df <- read.csv("data/rawMeasurements/surveyDateSheet.csv") %>% 
  mutate(doy = yday(mdy(eventDate)),
         year = year(mdy(eventDate))) %>% 
  mutate(doy2 = if_else(condition = year == 2019,
                        true = doy,
                        false = doy + 365)) %>% 
  rename(Site = location) %>% 
  mutate(SiteDate = paste(Site, eventDate, sep = "_")) %>% 
  mutate(Date = as_date(mdy(eventDate)))

## create models for each cumulative meaasure
# combine dataset with urbanization values
urb <- read.csv("data/data_products/urbanization_gradient_SITE.csv")
light <- read.csv("data/data_products/lightData_SITE.csv")
temp <- read.csv("data/data_products/temp_gradient_SITE.csv")
weather <- read.csv("data/weatherData/daymetDat.csv")
weather <- weather %>% 
  rename(doy = data.yday,
         year = data.year,
         prcp = data.prcp..mm.day.,
         tmin = data.tmin..deg.c.) %>% 
  select(doy, year, prcp, tmin)
lunarIllumination <- read.csv("data/data_products/lunarIllumination_lunarRPackage.csv")
lunarIllumination$Date <- as.Date(lunarIllumination$Date)
# join urb values with moth_df
adult_df <- adult_df%>% 
  left_join(urb) %>% 
  left_join(light) %>% 
  left_join(temp) %>% 
  left_join(weather) %>% 
  left_join(lunarIllumination, by = c("Date" = "Date")) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1) %>% 
  mutate(meanLight = log(meanLight + 0.01))

adult_df <- adult_df %>% 
  mutate(rel_temp = scale(rel_temp),
         meanLight = scale(meanLight),
         tmin = scale(tmin),
         prcp = scale(prcp),
         Dev_1 = scale(Dev_1),
         lunarIllumination = scale(lunarIllumination))

# make my one, well thought-out model
macro_dev <- brm(formula = bf( macroMoths ~ Dev_1  + rel_temp + lunarIllumination + prcp + tmin +
                                    (1|Site),
                                  zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                            prior = c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                                      set_prior("normal(0, 1)", class = "b")),
                            data = adult_df,
                            family = zero_inflated_negbinomial(),
                            chains = 4, iter = 2400, warmup = 1000,
                            control = list(adapt_delta = 0.93),
                            cores = 4, seed = 1234, 
                            threads = threading(2),
                            backend = "cmdstanr", 
  )



# examine model assumptions
plot(macro_dev)
pp_check(macro_dev)
pp_check(macro_dev, type = "stat", stat = "mean")


summary(macro_dev, prob = 0.89)

ce <- conditional_effects(x = macro_dev, effects = "Dev_1", prob = 0.89)
ce_df <- ce$Dev_1

macro_plot <- ggplot() +
  geom_jitter(adult_df, mapping = aes(x = Dev_1, y = macroMoths), alpha = 0.5) +
  geom_line(ce_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Total abundance") +
  ggtitle("Macro moths") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

macro_plot


adult_df_noBaca <- adult_df %>% 
  filter(Site != "BACA")
# make my one, well thought-out model
macro_dev_noBaca <- brm(formula = bf( macroMoths ~ Dev_1  + rel_temp + lunarIllumination + prcp + tmin +
                                 (1|Site),
                               zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                 prior = c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                           set_prior("normal(0, 1)", class = "b")),
                 data = adult_df_noBaca,
                 family = zero_inflated_negbinomial(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.95),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
)



# examine model assumptions
pp_check(macro_dev_noBaca)
pp_check(macro_dev_noBaca, type = "stat", stat = "mean")


summary(macro_dev_noBaca, prob = 0.89)


################################################################################
#######################MICROMOTHS###############################################
################################################################################
# make my one, well thought-out model
micro_dev <- brm(formula = bf( microMoths ~ Dev_1  + rel_temp + lunarIllumination + prcp + tmin +
                                 (1|Site),
                               zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                 prior = c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                           set_prior("normal(0, 1)", class = "b")),
                 data = adult_df,
                 family = zero_inflated_negbinomial(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.95),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
)

# examine model assumptions
plot(micro_dev)
pp_check(micro_dev)
pp_check(micro_dev, type = "stat", stat = "mean")


summary(micro_dev, prob = 0.89)

ce_micro <- conditional_effects(x = micro_dev, effects = "Dev_1", prob = 0.89)
ce_micro_df <- ce_micro$Dev_1

micro_plot <- ggplot() +
  geom_jitter(adult_df, mapping = aes(x = Dev_1, y = microMoths), alpha = 0.5) +
  geom_line(ce_micro_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_micro_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Total abundance") +
  ggtitle("Micro moths") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

micro_plot


###############################################################################
################################FRASS##########################################
###############################################################################
## total abundance -- frass
# read in frass data
frass_df <- read.csv("data/data_products/frass_biomass_reweigh_doyDiff.csv")
frass_df <- frass_df %>% 
  mutate(MassPerDay = log((Mass/diff) + 0.0001))
## total abundance -- frass
frass_weather <- left_join(weather, frass_df)

l <- unique(filter(frass_weather, !is.na(Date))$Date)

meanWeatherFun <- function(x){
  
  x <- x
  doyd <- filter(frass_df, Date == x)$doy[1]
  yeard <- filter(frass_df, Date == x)$year[1]
  d <- filter(frass_df, Date == x)$diff[1]
  
  meanWeekWeather <- weather %>% 
    filter(year == yeard & doy > (doyd - d) & doy <= doyd)
  
  out_df <- data.frame(Date = x, 
                       doy = doyd,
                       year = yeard,
                       mean_prcp = mean(meanWeekWeather$prcp),
                       mean_tmin = mean(meanWeekWeather$tmin))
  
}

ldf <- lapply(l, meanWeatherFun)

meanWeather_df <- bind_rows(ldf)
# do this for lunar date 
d_2019 <- weather %>% 
  filter(year == 2019) %>% 
  mutate(Date = as.Date(doy, origin = "2019-01-01") - 1)

d_2020 <- weather %>% 
  filter(year == 2020) %>% 
  mutate(Date = as.Date(doy, origin = "2020-01-01") - 1)

d_df <- rbind(d_2019, d_2020)

unique_dates <- unique(d_df$Date)

lunarIllumination <- read.csv('data/data_products/lunarIllumination_lunarRPackage.csv')

lunarDateFun <- function(x){
  
  x <- x
  doyd <- filter(frass_df, Date == x)$doy[1]
  yeard <- filter(frass_df, Date == x)$year[1]
  d <- filter(frass_df, Date == x)$diff[1]
  
  meanWeekLunar <- lunarIllumination %>% 
    mutate(year = year(Date), doy = yday(Date)) %>% 
    filter(year == yeard & doy > (doyd - d) & doy <= doyd)
  
  out_df <- data.frame(Date = x, 
                       doy = doyd,
                       year = yeard,
                       meanLunarIllumination = mean(meanWeekLunar$lunarIllumination))
  
  
}

surveyDates <- unique(frass_df$Date)

lp_ldf <- lapply(surveyDates, lunarDateFun)

lunarDates_df <- bind_rows(lp_ldf)

frass_df <- frass_df %>% 
  left_join(urb) %>% 
  left_join(light) %>% 
  left_join(temp) %>% 
  left_join(meanWeather_df) %>% 
  left_join(lunarDates_df) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1) %>% 
  mutate(meanLight = log(meanLight + 0.01))

frass_df <- frass_df %>% 
  mutate(rel_temp = scale(rel_temp),
         meanLight = scale(meanLight),
         mean_tmin = scale(mean_tmin),
         mean_prcp = scale(mean_prcp),
         Dev_1 = scale(Dev_1),
         meanLunarIllumination = scale(meanLunarIllumination)
  )

# per site?
frass_df_perSite <- frass_df %>% 
  group_by(Site, Date) %>% 
  summarise(MeanMassPerDay = mean(MassPerDay, na.rm = T))

frass_df_perSite <- ungroup(frass_df_perSite) %>% 
  left_join(urb) %>% 
  left_join(light) %>% 
  left_join(temp) %>% 
  left_join(meanWeather_df) %>% 
  left_join(lunarDates_df) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1) %>% 
  mutate(meanLight = log(meanLight + 0.01))

frass_df_perSite <- frass_df_perSite %>% 
  mutate(rel_temp = scale(rel_temp),
         meanLight = scale(meanLight),
         mean_tmin = scale(mean_tmin),
         mean_prcp = scale(mean_prcp),
         Dev_1 = scale(Dev_1),
         meanLunarIllumination = scale(meanLunarIllumination)
  )


frass_dev <- brm(formula = bf(MeanMassPerDay ~ Dev_1  + rel_temp + meanLunarIllumination + mean_tmin + mean_prcp +
                                (1|Site)),
                 prior = c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                           set_prior("normal(0, 1)", class = "b")),
                 data = frass_df_perSite,
                 family = gaussian(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.95),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
)

# examine model assumptions
plot(frass_dev)
pp_check(frass_dev)
pp_check(frass_dev, type = "stat", stat = "mean")

summary(frass_dev, prob = 0.89)

ce_frass <- conditional_effects(x = frass_dev, effects = "Dev_1", prob = 0.89)
ce_frass_df <- ce_frass$Dev_1

frass_plot <- ggplot() +
  geom_jitter(frass_df_perSite, mapping = aes(x = Dev_1, y = MeanMassPerDay), alpha = 0.5) +
  geom_line(ce_frass_df, mapping = aes(x = Dev_1, y = estimate__)) +
  geom_ribbon(ce_frass_df, mapping = aes(x = Dev_1, ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Mean frass mass per day \n log(x + 0.001)") +
  ggtitle("Caterpillars") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

frass_plot

# plot modeling outputs
cp <- cowplot::plot_grid(macro_plot, micro_plot, frass_plot, 
                         labels = c("A", "B", "C"),
                         nrow = 3, ncol = 1)

if(dir.exists("figOutputs")){
  ggsave(plot = cp, filename = "figOutputs/pooledAbundance.png", dpi = 500,
         width = 6, height = 10)
} else{
  dir.create("figOutputs")
  ggsave(plot = cp, filename = "figOutputs/pooledAbundance.png", dpi = 500,
         width = 6, height = 10)
}


