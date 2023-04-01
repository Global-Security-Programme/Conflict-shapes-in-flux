### Preamble ###############################################################
# Calculations for Table 7 in the Extended appendix
# Supplementary material for the article Conflict shapes in flux by Idler and Tkacova (2023)
# Date created: 27 March 2023

rm(list= ls())

library(devtools)
library(purrr)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(tidyverse)


# Nigeria -----------------------------------------------------------------
rm(list= ls())

data <- readRDS("data/nigeria.rds")

# Number of engagements 
data.net <- data %>% 
  dplyr::select(side_a, side_b, year, id, dyad_name) %>% 
  mutate(From = side_a, To = side_b, Time = year, ID.match = id)

calculate_engag <- function(g) {
  data.frame(N_engagements = igraph::degree(g), total_N_engagements = gsize(g)) }

Dsplit <- split(data.net, data.net$Time)
lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))

res_engag <- lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))
engagements <- map_df(res_engag, ~.x, .id="Time")
engagements <- tibble::rownames_to_column(engagements, "Side")
engagements$Side <- str_replace(engagements$Side,"\\..*","")
engagements <- engagements %>% mutate(engagements_prop = round(N_engagements/total_N_engagements, 1))
engagements <- engagements %>% filter(Side == "Government of Nigeria", 
                                      Time >=2011 & Time <= 2016)

write_csv(engagements, "out/nigeria-ext-app-tab7-engagements.csv")

# Number of actors 
data.net.uniq <- data.net %>% distinct(dyad_name, Time, .keep_all = T)

calculate_actors <- function(g) {
  data.frame(N_actors = igraph::degree(g), total_N_actors = gorder(g)) }

Dsplit <- split(data.net.uniq, data.net.uniq$Time)
lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))

res_act <- lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))
actors <- map_df(res_act, ~.x, .id="Time")
actors <- tibble::rownames_to_column(actors, "Side")
actors$Side <- str_replace(actors$Side,"\\..*","")
actors <- actors %>% mutate(actors_prop = round(N_actors/total_N_actors, 1))
actors <- actors %>% filter(Side == "Government of Nigeria", 
                                      Time >=2011 & Time <= 2016)

write_csv(actors, "out/nigeria-ext-app-tab7-actors.csv")


# Colombia ----------------------------------------------------------

# Data
rm(list= ls())

data <- readRDS("data/colombia.rds")

# Number of engagements 
data.net <- data %>% 
  dplyr::select(side_a, side_b, year, id, dyad_name) %>% 
  mutate(From = side_a, To = side_b, Time = year, ID.match = id)

calculate_engag <- function(g) {
  data.frame(N_engagements = igraph::degree(g), total_N_engagements = gsize(g)) }

Dsplit <- split(data.net, data.net$Time)
lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))

res_engag <- lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))
engagements <- map_df(res_engag, ~.x, .id="Time")
engagements <- tibble::rownames_to_column(engagements, "Side")
engagements$Side <- str_replace(engagements$Side,"\\..*","")
engagements <- engagements %>% mutate(engagements_prop = round(N_engagements/total_N_engagements, 1))

engagements <- engagements %>% filter(Side == "Government of Colombia", 
                            Time >=2006 & Time <= 2016)

write_csv(engagements, "out/colombia-ext-app-tab7-engagements.csv")

# Number of actors 
data.net.uniq <- data.net %>% distinct(dyad_name, Time, .keep_all = T)

calculate_actors <- function(g) {
  data.frame(N_actors = igraph::degree(g), total_N_actors = gorder(g)) }

Dsplit <- split(data.net.uniq, data.net.uniq$Time)
lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))

res_act <- lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))
actors <- map_df(res_act, ~.x, .id="Time")
actors <- tibble::rownames_to_column(actors, "Side")
actors$Side <- str_replace(actors$Side,"\\..*","")
actors <- actors %>% mutate(actors_prop = round(N_actors/total_N_actors, 1))

actors <- actors %>% filter(Side == "Government of Colombia", 
                                      Time >=2006 & Time <= 2016)

write_csv(actors, "out/colombia-ext-app-tab7-actors.csv")


# Afghanistan-Pakistan ----------------------------------------------------

rm(list= ls())

data <- readRDS("data/afghanistan-pakistan.rds")

# Number of engagements 
data.net <- data %>% 
  dplyr::select(side_a, side_b, year, id, dyad_name) %>% 
  mutate(From = side_a, To = side_b, Time = year, ID.match = id)

calculate_engag <- function(g) {
  data.frame(N_engagements = igraph::degree(g), total_N_engagements = gsize(g)) }

Dsplit <- split(data.net, data.net$Time)
lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))

res_engag <- lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))
engagements <- map_df(res_engag, ~.x, .id="Time")
engagements <- tibble::rownames_to_column(engagements, "Side")
engagements$Side <- str_replace(engagements$Side,"\\..*","")
engagements <- engagements %>% mutate(engagements_prop = round(N_engagements/total_N_engagements, 1))

engagements <- engagements %>% filter(Side == "Government of Afghanistan") %>% 
  filter(Time >=2006 & Time <= 2008)
                            
write_csv(engagements, "out/afghanistan-pakistan-ext-app-tab7-engagements.csv")

# Number of actors 
data.net.uniq <- data.net %>% distinct(dyad_name, Time, .keep_all = T)

calculate_actors <- function(g) {
  data.frame(N_actors = igraph::degree(g), total_N_actors = gorder(g)) }

Dsplit <- split(data.net.uniq, data.net.uniq$Time)
lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))

res_act <- lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))
actors <- map_df(res_act, ~.x, .id="Time")
actors <- tibble::rownames_to_column(actors, "Side")
actors$Side <- str_replace(actors$Side,"\\..*","")
actors <- actors %>% mutate(actors_prop = round(N_actors/total_N_actors, 1))

actors <- actors %>% filter(Side == "Government of Afghanistan") %>% 
  filter(Time >=2006 & Time <= 2008)

write_csv(actors, "out/afghanistan-pakistan-ext-app-tab7-actors.csv")


# Syria-Iraq --------------------------------------------------------------

rm(list= ls())

data <- readRDS("data/syria-iraq.rds")

# Number of engagements 
data.net <- data %>% 
  dplyr::select(side_a, side_b, year, id, dyad_name) %>% 
  mutate(From = side_a, To = side_b, Time = year, ID.match = id)

calculate_engag <- function(g) {
  data.frame(N_engagements = igraph::degree(g), total_N_engagements = gsize(g)) }

Dsplit <- split(data.net, data.net$Time)
lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))

res_engag <- lapply(Dsplit, function(x) calculate_engag(graph_from_data_frame(x)))
engagements <- map_df(res_engag, ~.x, .id="Time")
engagements <- tibble::rownames_to_column(engagements, "Side")
engagements$Side <- str_replace(engagements$Side,"\\..*","")
engagements <- engagements %>% mutate(engagements_prop = round(N_engagements/total_N_engagements, 1))

engagements <- engagements %>% filter(Side == "Government of Iraq" | Side == "Government of Syria") %>% 
  filter(Time >=2010 & Time <= 2016)

write_csv(engagements, "out/syria-iraq-ext-app-tab7-engagements.csv")

# Number of actors 
data.net.uniq <- data.net %>% distinct(dyad_name, Time, .keep_all = T)

calculate_actors <- function(g) {
  data.frame(N_actors = igraph::degree(g), total_N_actors = gorder(g)) }

Dsplit <- split(data.net.uniq, data.net.uniq$Time)
lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))

res_act <- lapply(Dsplit, function(x) calculate_actors(graph_from_data_frame(x)))
actors <- map_df(res_act, ~.x, .id="Time")
actors <- tibble::rownames_to_column(actors, "Side")
actors$Side <- str_replace(actors$Side,"\\..*","")
actors <- actors %>% mutate(actors_prop = round(N_actors/total_N_actors, 1))

actors <- actors %>% filter(Side == "Government of Iraq" | Side == "Government of Syria") %>% 
  filter(Time >=2010 & Time <= 2016)

write_csv(actors, "out/syria-iraq-ext-app-tab7-actors.csv")
