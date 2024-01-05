library(dplyr)
library(ggplot2)
library(lubridate)
library(brms)
library(equatiomatic)
library(lme4)
library(cmdstanr)
set_cmdstan_path()

# read in adult data
adult_df <- read.csv("data/rawMeasurements/surveyDateSheet.csv") %>% 
  mutate(doy = yday(mdy(eventDate)),
         year = year(mdy(eventDate))) %>% 
  mutate(doy2 = if_else(condition = year == 2019,
                        true = doy,
                        false = doy + 365)) %>% 
  rename(Site = location) %>% 
  mutate(SiteDate = paste(Site, eventDate, sep = "_")) %>% 
  mutate(Date = as_date(mdy(eventDate))) %>% 
  mutate(Site = stringr::str_to_title(Site))

abundanceDF <- adult_df %>% 
  group_by(Date, Site) %>% 
  summarise(macroMoths = sum(macroMoths),
            microMoths = sum(microMoths))

sampleDates <- dplyr::distinct(adult_df, Date, doy2, year, doy)

adult_df <- left_join(abundanceDF, sampleDates)

## create models for each cumulative measure
# combine dataset with urbanization values
urb <- read.csv("data/data_products/urbanization_gradient.csv")
temp <- read.csv("data/data_products/temp_gradient.csv")
weather <- read.csv("data/weatherData/daymetData.csv")
weather <- weather %>% 
  rename(doy = yday,
         year = year,
         prcp = prcp..mm.day.,
         tmin = tmin..deg.c.) %>% 
  select(doy, year, prcp, tmin)
lunarIllumination <- read.csv("data/data_products/lunarIllumination_lunarRPackage.csv")
lunarIllumination$Date <- as.Date(lunarIllumination$Date)
# join urb values with moth_df
adult_df <- adult_df%>% 
  left_join(urb) %>% 
  left_join(temp) %>% 
  left_join(weather) %>% 
  left_join(lunarIllumination, by = c("Date" = "Date")) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1) 

meanDev <- mean(adult_df$Dev_1)
sdDev <- sd(adult_df$Dev_1)

adult_df <- ungroup(adult_df) %>% 
  mutate(rel_temp = scale(rel_temp),
         tmin = scale(tmin),
         prcp = scale(prcp),
         Dev_1 = scale(Dev_1),
         lunarIllumination = scale(lunarIllumination))

# get priors
get_prior(macroMoths ~ Dev_1  + lunarIllumination + prcp + tmin +
            (1|Site), data = adult_df, family = zero_inflated_negbinomial())

slopePrior <- prior(normal(0,10), class = b)
# make my one, well thought-out model
set.seed(2963)
macro_dev <- brm(formula = bf(macroMoths ~ Dev_1  + lunarIllumination + prcp + tmin +
                                    (1|Site),
                                  zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                            data = adult_df,
                            family = zero_inflated_negbinomial(),
                            chains = 4, iter = 2400, warmup = 1000,
                            control = list(adapt_delta = 0.99),
                            cores = 4, seed = 1234, 
                            threads = threading(2),
                            backend = "cmdstanr", 
                 prior = slopePrior)



# examine model assumptions
plot(macro_dev)
pp_check(macro_dev, ndraws = 100)
pp_check(macro_dev, type = "stat", stat = "mean")

# check for spatial autocorrelation
resids <- residuals(macro_dev)

resids_df <- adult_df %>% 
  mutate(resids = resids[,1]) %>% 
  group_by(Site) %>% 
  summarise(meanResid = mean(resids)) %>% 
  ungroup()

library(sf)
coordsDF <- read.csv("data/data_products/siteCoordinates.csv")
coordsDF_m <- st_as_sf(coordsDF, coords = c("long", "lat"), crs = "WGS84") %>% 
  st_transform(crs = "+proj=laea +lon_0=-78.7495368 +lat_0=32.7548488 +datum=WGS84 +units=m +no_defs")
x <- st_coordinates(coordsDF_m)[,1]
y <- st_coordinates(coordsDF_m)[,2]
  
resids_df <- left_join(resids_df, coordsDF)
library(ncf)
spAuCoTest <- ncf::correlog(x = x, y, resids_df$meanResid, latlon = F, increment = 100) 
plot(spAuCoTest) # no signal of spatial autocorrelation in residuals

summary(macro_dev, prob = 0.89)
conditional_effects(x = macro_dev, effects = "prcp", prob = 0.89, dpar = "zi")

macro_sum <- summary(macro_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce <- conditional_effects(x = macro_dev, effects = "Dev_1", prob = 0.89)
ce_df <- ce$Dev_1

# plot on original urb development scale 
macro_plot <- ggplot() +
  geom_jitter(adult_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = macroMoths), alpha = 0.5) +
  geom_line(ce_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__)) +
  geom_ribbon(ce_df, mapping = aes(x = meanDev + (sdDev*Dev_1), ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Total abundance") +
  ggtitle("Macro moths") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

macro_plot

# rerun model without baca
adult_df_noBaca <- adult_df %>% 
  filter(Site != "Baca")
# make my one, well thought-out model
macro_dev_noBaca <- brm(formula = bf( macroMoths ~ Dev_1  + 
                                        lunarIllumination + prcp + tmin +
                                 (1|Site),
                               zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                 data = adult_df_noBaca,
                 family = zero_inflated_negbinomial(),
                 chains = 4, iter = 2800, warmup = 1000,
                 control = list(adapt_delta = 0.999),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
                 prior = slopePrior
)

# examine model assumptions
pp_check(macro_dev_noBaca)
pp_check(macro_dev_noBaca, type = "stat", stat = "mean")

summary(macro_dev_noBaca, prob = 0.89)

macro_sum_noBaca <- summary(macro_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

################################################################################
#######################MICROMOTHS###############################################
################################################################################
# make my one, well thought-out model
set.seed(42)
micro_dev <- brm(formula = bf( microMoths ~ Dev_1  +
                                 lunarIllumination + prcp + tmin +
                                 (1|Site),
                               zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                 data = adult_df,
                 family = zero_inflated_negbinomial(),
                 chains = 4, iter = 2700, warmup = 1000,
                 control = list(adapt_delta = 0.99),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
                 prior = slopePrior
)

# examine model assumptions
plot(micro_dev)
pp_check(micro_dev,ndraws = 100)
pp_check(micro_dev, type = "stat", stat = "mean")

# check for spatial autocorrelation
resids <- residuals(micro_dev)

resids_df <- adult_df %>% 
  mutate(resids = resids[,1]) %>% 
  group_by(Site) %>% 
  summarise(meanResid = mean(resids)) %>% 
  ungroup()

resids_df <- left_join(resids_df, coordsDF)

spAuCoTest <- ncf::correlog(x = x, y, resids_df$meanResid, latlon = F, increment = 100) 
plot(spAuCoTest) # no signal of spatial autocorrelation in residuals

ggplot(resids_df, mapping = aes(x = long, y = lat, color = meanResid)) +
  geom_point(size = 2) +
  scale_color_gradient2() 

# results
summary(micro_dev, prob = 0.89)
micro_sum <- summary(micro_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_micro <- conditional_effects(x = micro_dev, effects = "Dev_1", prob = 0.89)
ce_micro_df <- ce_micro$Dev_1

micro_plot <- ggplot() +
  geom_jitter(adult_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = microMoths), alpha = 0.5) +
  geom_line(ce_micro_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__)) +
  geom_ribbon(ce_micro_df, mapping = aes(x = meanDev + (sdDev*Dev_1), ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Total abundance") +
  ggtitle("Micro moths") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

micro_plot

# micro model no baca
set.seed(21)
micro_dev_noBaca <- brm(formula = bf(microMoths ~ Dev_1  +  
                                       lunarIllumination + prcp + tmin +
                                 (1|Site),
                               zi ~ Dev_1 + lunarIllumination + prcp + tmin + (1|Site)),
                 data = adult_df_noBaca,
                 family = zero_inflated_negbinomial(),
                 chains = 4, iter = 3200, warmup = 1000,
                 control = list(adapt_delta = 0.99999999),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
                 prior = slopePrior
)

# examine model assumptions
plot(micro_dev_noBaca)
pp_check(micro_dev_noBaca)
pp_check(micro_dev_noBaca, type = "stat", stat = "mean")

summary(micro_dev_noBaca, prob = 0.89)
micro_sum_noBaca <- summary(micro_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

###############################################################################
################################FRASS##########################################
###############################################################################
## total abundance -- frass
# read in frass data
frass_df <- read.csv("data/rawMeasurements/frass_biomass.csv")
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

# per site
frass_df_perSite <- frass_df %>% 
  mutate(Site = stringr::str_to_title(Site)) %>% 
  group_by(Site, Date) %>% 
  summarise(MeanMassPerDay = mean(MassPerDay, na.rm = T))

frass_df_perSite <- ungroup(frass_df_perSite) %>% 
  left_join(urb) %>% 
  left_join(temp) %>% 
  left_join(meanWeather_df) %>% 
  left_join(lunarDates_df) %>% 
  mutate(rel_temp = (mean_temp - 1.02)*-1)

meanDev <- mean(frass_df_perSite$Dev_1)
sdDev <- sd(frass_df_perSite$Dev_1)

frass_df_perSite <- frass_df_perSite %>% 
  mutate(rel_temp = scale(rel_temp),
         mean_tmin = scale(mean_tmin),
         mean_prcp = scale(mean_prcp),
         Dev_1 = scale(Dev_1),
         meanLunarIllumination = scale(meanLunarIllumination)
  )


set.seed(789)
frass_dev <- brm(formula = bf(MeanMassPerDay ~ Dev_1  +
                                meanLunarIllumination + mean_tmin + mean_prcp +
                                (1|Site)),
                 data = frass_df_perSite,
                 family = gaussian(),
                 chains = 4, iter = 2400, warmup = 1000,
                 control = list(adapt_delta = 0.999),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
                 prior = slopePrior
)

# examine model assumptions
plot(frass_dev)
pp_check(frass_dev, ndraws = 100)
pp_check(frass_dev, type = "stat", stat = "mean")

summary(frass_dev, prob = 0.89)

frass_sum <- summary(frass_dev, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

ce_frass <- conditional_effects(x = frass_dev, effects = "Dev_1", prob = 0.89)
ce_frass_df <- ce_frass$Dev_1

frass_plot <- ggplot() +
  geom_jitter(frass_df_perSite, mapping = aes(x = meanDev + (sdDev*Dev_1), y = MeanMassPerDay), alpha = 0.5) +
  geom_line(ce_frass_df, mapping = aes(x = meanDev + (sdDev*Dev_1), y = estimate__)) +
  geom_ribbon(ce_frass_df, mapping = aes(x = meanDev + (sdDev*Dev_1), ymax = upper__, ymin = lower__), alpha = 0.3) +
  labs(x = "Urban development", y = "Mean frass mass per day \n log(x + 0.001)") +
  ggtitle("Caterpillars") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

frass_plot

# frass model no baca
frass_df_perSite_noBaca <- filter(frass_df_perSite, Site != "Baca")
set.seed(987)
frass_dev_noBaca <- brm(formula = bf(MeanMassPerDay ~ Dev_1  +
                                       meanLunarIllumination + mean_tmin + mean_prcp +
                                (1|Site)),
                 data = frass_df_perSite_noBaca,
                 family = gaussian(),
                 chains = 4, iter = 2800, warmup = 1000,
                 control = list(adapt_delta = 0.9999),
                 cores = 4, seed = 1234, 
                 threads = threading(2),
                 backend = "cmdstanr", 
                 prior = slopePrior
)

# examine model assumptions
plot(frass_dev_noBaca)
pp_check(frass_dev_noBaca)
pp_check(frass_dev_noBaca, type = "stat", stat = "mean")

summary(frass_dev_noBaca, prob = 0.89)

frass_sum_noBaca <- summary(frass_dev_noBaca, prob = 0.89)$fixed %>% 
  tibble::rownames_to_column()

# plot modeling outputs
cp <- cowplot::plot_grid(macro_plot, micro_plot, frass_plot, 
                         labels = c("A", "B", "C"),
                         nrow = 3, ncol = 1)

if(dir.exists("figOutputs")){
  ggsave(plot = cp, filename = "figOutputs/pooledAbundance.png", dpi = 500,
         width = 4, height = 7)
} else{
  dir.create("figOutputs")
  ggsave(plot = cp, filename = "figOutputs/pooledAbundance.png", dpi = 500,
         width = 4, height = 7)
}

## Table outputs
if(dir.exists("tabOutputs")){
  print("Directory exists, wahoo")
} else{
  dir.create("tabOutputs")
}

write.csv(macro_sum, file = "tabOutputs/pooledMacroMothResults.csv", row.names = F)
write.csv(macro_sum_noBaca, file = "tabOutputs/pooledMacroMothResults_noBaca.csv", row.names = F)

write.csv(micro_sum, file = "tabOutputs/pooledMicroMothResults.csv", row.names = F)
write.csv(micro_sum_noBaca, file = "tabOutputs/pooledMicroMothResults_noBaca.csv", row.names = F)

write.csv(frass_sum, file = "tabOutputs/pooledCaterpillarsResults.csv", row.names = F)
write.csv(frass_sum_noBaca, file = "tabOutputs/pooledCaterpillarsResults_noBaca.csv", row.names = F)
