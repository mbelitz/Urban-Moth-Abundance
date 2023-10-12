library(tidyverse)
library(plotly)

## make 3d plot of data
dev <- read.csv("data/data_products/urbanization_gradient.csv")
light <- read.csv("data/data_products/lightData.csv")
temp <- read.csv("data/data_products/temp_gradient.csv")

tdf <- dev %>% 
  left_join(light) %>% 
  left_join(temp) %>% 
  mutate(Class = case_when(
    Site %in% c("Baca", "Joma", "Cofr") ~ "Urban",
    Site %in% c("Demi", "Biva", "Bowa") ~ "Suburban",
    Site %in% c("Auca", "Rist", "Prcr") ~ "Rural"
  ))

# make sub categories categorical
plot_ly(x=tdf$Dev_1, y=((tdf$mean_temp - 1.02) * 2), z=log(tdf$meanLight + 0.01), 
        type="scatter3d", mode="markers", color=tdf$Class, ) 

