### Preamble ###############################################################
# Armed conflict in the Afghan-Pakistani borderlands - Wzones: Spatial change analysis
# Supplementary material for the article Conflict shapes in flux by Idler and Tkacova (2023)
# Date created: 27 March 2023

rm(list= ls())

library(ggplot2)
library(concaveman)
library(GGally)
library(igraph)
library(sna)
library(stats)
library(raster)
library(MASS)
library(rgeos)
library(rnaturalearth)
library(sf)
sf::sf_use_s2(FALSE)
library(spdep)
library(rgeos)
library(maptools)
library(units)
library(tidyverse)


# Data and map preparation ------------------------------------------------

data <- readRDS("data/afghanistan-pakistan.rds")

data <- data %>% mutate(deaths_total = deaths_battle + deaths_civilians + 
                          deaths_unknown)

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))


# 2006: Conflict shape and hotspots --------------------------------------------------------

points.sf.2006 <- points.sf %>% filter(year == 2006)

# wzoneData
wzone.2006 <- readOGR(paste0('data/replication-kikuta/2006/2006_12_31.shp'), stringsAsFactors = F)
wzone.2006.sel <- wzone.2006[wzone.2006$dyad_id == 735 |        
                               wzone.2006$dyad_id == 726,]     

wzone.2006.sel.un <- unionSpatialPolygons(wzone.2006.sel, rep(1, length(wzone.2006.sel)))
wzone.2006.sel.un.sf <- st_as_sf(wzone.2006.sel.un, crs = st_crs(4326))

grid.clip.2006 <- st_make_grid(wzone.2006.sel.un.sf, 
                               cellsize = c(0.5,0.5),  
                               what="polygons") %>%
  st_intersection(wzone.2006.sel.un.sf)

points.count.2006 <- st_intersects(grid.clip.2006, points.sf.2006)
grid.count.2006 <- st_sf(pt_count = lengths(points.count.2006), 
                         geometry = st_cast(grid.clip.2006, "MULTIPOLYGON"))

neighbors.2006 <- poly2nb(grid.count.2006, queen = T)
weighted.neighbors.2006 <- nb2listw(include.self(neighbors.2006), 
                                    style = "B", 
                                    zero.policy = T)
grid.count.2006$HOTSPOT <- as.vector(localG_perm(grid.count.2006$pt_count, 
                                                 weighted.neighbors.2006))

grid.count.2006$z_cat <- "Not significant or cold spot"
grid.count.2006$z_cat[grid.count.2006$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2006$z_cat[grid.count.2006$HOTSPOT < 3.091 & grid.count.2006$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2006$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2006, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands - Wzones, 2006: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-wzones2006-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2006.sel <- grid.count.2006[grid.count.2006$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2006$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006 <- st_union(grid.count.2006.sel)
st_write(hotspot.2006, "out/afghanistan-pakistan-wzones2006-hotspot.shp")


# 2008: Conflict shape and hotspots --------------------------------------------------------

points.sf.2008 <- points.sf %>% filter(year == 2008)

# wzoneData
wzone.2008 <- readOGR(paste0('data/replication-kikuta/2008/2008_12_31.shp'), stringsAsFactors = F)
wzone.2008.sel <- wzone.2008[wzone.2008$dyad_id == 735 |        
                               wzone.2008$dyad_id == 726 |     
                               wzone.2008$dyad_id == 11336 |   
                               wzone.2008$dyad_id == 857 |      
                               wzone.2008$dyad_id == 5198 |     
                               wzone.2008$dyad_id == 5269|      
                               wzone.2008$dyad_id == 5277 |    
                               wzone.2008$dyad_id == 5189 |     
                               wzone.2008$dyad_id == 5550 |     
                               wzone.2008$dyad_id == 5523 |     
                               wzone.2008$dyad_id == 5257,]     

wzone.2008.sel.un <- unionSpatialPolygons(wzone.2008.sel, rep(1, length(wzone.2008.sel)))
wzone.2008.sel.un.sf <- st_as_sf(wzone.2008.sel.un, crs = st_crs(4326))

grid.clip.2008 <- st_make_grid(wzone.2008.sel.un.sf, 
                               cellsize = c(0.5,0.5), #cell size in degrees
                               what="polygons") %>%
  st_intersection(wzone.2008.sel.un.sf)

points.count.2008 <- st_intersects(grid.clip.2008, points.sf.2008)
grid.count.2008 <- st_sf(pt_count = lengths(points.count.2008), 
                         geometry = st_cast(grid.clip.2008, "MULTIPOLYGON"))

neighbors.2008 <- poly2nb(grid.count.2008, queen = T)
weighted.neighbors.2008 <- nb2listw(include.self(neighbors.2008), 
                                    style = "B", 
                                    zero.policy = T)

grid.count.2008$HOTSPOT <- as.vector(localG_perm(grid.count.2008$pt_count, 
                                                 weighted.neighbors.2008))
grid.count.2008$z_cat <- "Not significant or cold spot"
grid.count.2008$z_cat[grid.count.2008$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2008$z_cat[grid.count.2008$HOTSPOT < 3.091 & grid.count.2008$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2008$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2008, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands - Wzones, 2008: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-wzones2008-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2008.sel <- grid.count.2008[grid.count.2008$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2008$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008 <- st_union(grid.count.2008.sel)
st_write(hotspot.2008, "out/afghanistan-pakistan-wzones2008-hotspot.shp")


# 2006-2008: Spatial change ---------------------------------------------------------------------

label.2006 <- data.frame(
  x = 64, 
  y = 28,
  text = "2006")
label.2008 <- data.frame(
  x = 73, 
  y = 36.5, 
  text = "2008")

# map
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2006, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2006.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  
  geom_sf(data = hotspot.2008, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2008.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands - Wzones, 2006-2008: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2006, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2008, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/afghanistan-pakistan-wzones2006-2008-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2006.2008 <- st_area(wzone.2006.sel.un.sf)/100
new.perc.2006.2008 <- st_area(wzone.2008.sel.un.sf)/one.perc.2006.2008
new.perc.2006.2008 <- drop_units(new.perc.2006.2008)
perc.change.2006.2008 <- new.perc.2006.2008 - 100
perc.change.2006.2008 # 83.212 (expansion)

# overlap: hotspots
hotspot.2006.sp <- as(st_geometry(hotspot.2006), "Spatial")
hotspot.2008.sp <- as(st_geometry(hotspot.2008), "Spatial")
intersection.h.2006.2008 <- gIntersection(hotspot.2006.sp, hotspot.2008.sp)
intersection.h.2006.2008.sf <- st_as_sf(intersection.h.2006.2008)
overl_hotspot.2006.2008 <- st_area(intersection.h.2006.2008.sf)/(st_area(hotspot.2006)/100)  
overl_hotspot.2006.2008 # 42.96087

# overlap: conflict shape
buffer.2006.sp <- as(st_geometry(wzone.2006.sel.un.sf), "Spatial")
buffer.2008.sp <- as(st_geometry(wzone.2008.sel.un.sf), "Spatial")
intersection.2006.2008 <- gIntersection(buffer.2006.sp, buffer.2008.sp)
intersection.2006.2008.sf <- st_as_sf(intersection.2006.2008)
overl_shape.2006.2008 <- st_area(intersection.2006.2008.sf)/(st_area(wzone.2006.sel.un.sf)/100) 
overl_shape.2006.2008 # 100

# spatial change
spatial_ch.2006.2008 <- (overl_hotspot.2006.2008 + overl_shape.2006.2008)/2
spatial_ch.2006.2008 # 71.48044 (no shift)
