### Preamble ###############################################################
# Figure 5 in the Appendix
# Supplementary material for the article Conflict shapes in flux by Idler and Tkacova (2023)
# Date created: 27 March 2023

rm(list= ls())

library(ggplot2)
library(concaveman)
library(GGally)
library(igraph)
library(sna)
library(centiserve)
library(stats)
library(raster)
library(MASS)
library(rgeos)
library(rnaturalearth)
library(sf)
sf::sf_use_s2(FALSE)
library(spdep)
library(units)
library(tidyverse)


# Figure 5.1 --------------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

data <- readRDS("data/nigeria.rds")

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326

points.sf.2011 <- points.sf %>% filter(year == 2011)

# wzoneData
wzone.2011 <- readOGR(paste0('data/replication-kikuta/2011/2011_12_31.shp'), stringsAsFactors = F)
wzone.2011.sel <- wzone.2011[wzone.2011$dyad_id == 640 |        
                               wzone.2011$dyad_id == 5581 |      
                               wzone.2011$dyad_id == 5558,]      

wzone.2011.sel.un <- unionSpatialPolygons(wzone.2011.sel, rep(1, length(wzone.2011.sel)))
wzone.2011.sel.un.sf <- st_as_sf(wzone.2011.sel.un, crs = st_crs(4326))

# map
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2011,
          aes(color = dyad_name, show.legend = "point"), alpha = 0.6, shape = 1) +
  geom_sf(data = wzone.2011.sel.un.sf, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Islamist insurgency, 2011: Wzone and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3) +
  theme(legend.position="right") +
  labs(color = "Conflict dyad")
ggsave("figs/nigeria-2011-wzones-events.png", width = 8, height = 4)



# Figure 5.2 --------------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

# wzoneData
wzone.2016 <- readOGR(paste0('data/replication-kikuta/2016/2016_12_31.shp'), stringsAsFactors = F)
wzone.2016.sel <- wzone.2016[wzone.2016$dyad_id == 640 |        
                               wzone.2016$dyad_id == 12422 |      
                               wzone.2016$dyad_id == 14666 |      
                               wzone.2016$dyad_id == 14671 |      
                               wzone.2016$dyad_id == 14668 |      
                               wzone.2016$dyad_id == 14669 |      
                               wzone.2016$dyad_id == 14667 ,]     

wzone.2016.sel.un <- unionSpatialPolygons(wzone.2016.sel, rep(1, length(wzone.2016.sel)))
wzone.2016.sel.un.sf <- st_as_sf(wzone.2016.sel.un, crs = st_crs(4326))

# map
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2016,
          aes(color = dyad_name, show.legend = "point"), alpha = 0.6, shape = 1) +
  geom_sf(data = wzone.2016.sel.un.sf, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Islamist insurgency, 2016: Wzone and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3) +
  theme(legend.position="right") +
  labs(color = "Conflict dyad")
ggsave("figs/nigeria-2016-wzones-events.png", width = 8, height = 4)


# Figure 5.3 --------------------------------------------------------------

rm(list= ls())

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

data <- readRDS("data/syria-iraq.rds")

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326
points.sf.2013 <- points.sf %>% filter(year == 2013)

# wzoneData
wzone.2013 <- readOGR(paste0('data/replication-kikuta/2013/2013_12_31.shp'), stringsAsFactors = F)
wzone.2013.sel <- wzone.2013[wzone.2013$dyad_id == 754 |        
                               wzone.2013$dyad_id == 524 |      
                               wzone.2013$dyad_id == 781 |      
                               wzone.2013$dyad_id == 12337 |      
                               wzone.2013$dyad_id == 11973 |      
                               wzone.2013$dyad_id == 13767 |      
                               wzone.2013$dyad_id == 13793 |      
                               wzone.2013$dyad_id == 13813 |      
                               wzone.2013$dyad_id == 13821 |      
                               wzone.2013$dyad_id == 14341 |      
                               wzone.2013$dyad_id == 14620 ,]      


wzone.2013.sel.un <- unionSpatialPolygons(wzone.2013.sel, rep(1, length(wzone.2013.sel)))
wzone.2013.sel.un.sf <- st_as_sf(wzone.2013.sel.un, crs = st_crs(4326))

points.sf.2013$without_wzones <- "Yes"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Iran - PJAK"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Iraq - IS"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Turkey - PKK"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Hezbollah - Jabhat Fateh al-Sham"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Syria - Syrian insurgents"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "IS, Jabhat Fateh al-Sham - PYD"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "al-Tawhid Brigade - PYD"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "IS - Jabhat Fateh al-Sham"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Syria - PYD"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "IS - Islamic Front"] <- "No"
points.sf.2013$without_wzones[points.sf.2013$dyad_name == "Government of Syria - IS"] <- "No"

table(points.sf.2013$without_wzones)

# map
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2013,
          aes(color = without_wzones, show.legend = "point"), alpha = 0.6, shape = 1) +
  geom_sf(data = wzone.2013.sel.un.sf, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq, 2013: Wzone and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(31, 51), ylim = c(29, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3) +
  theme(legend.position="right") +
  labs(color = "Without wzone")
ggsave("figs/syria-iraq-2013-wzones-events.png", width = 8, height = 4)
