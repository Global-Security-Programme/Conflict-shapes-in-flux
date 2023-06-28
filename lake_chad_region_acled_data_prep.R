### Preamble ###############################################################
# Armed conflict in the Lake Chad region - ACLED: Data preparation
# Supplementary material for the article Conflict shapes in flux by Idler and Tkacova (2023)
# Date created: 27 March 2023

rm(list= ls())

library(tidyverse)

data <- read.csv("data/red-nigeria-acled-4-no-outliers-after-friction-added.csv")

data <- data %>% filter(event_relevant != 0 &
                          event_rel_out == 0 &
                          geo_precision <= 2)

data <- data %>% dplyr::select(data_id, year, actor1, actor2, 
                               longitude, latitude, geo_precision,
                               fatalities)

saveRDS(data, "data/lake-chad-region-acled.rds")
