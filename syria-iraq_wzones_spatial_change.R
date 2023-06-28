### Preamble ###############################################################
# Armed conflict in Syria and Iraq - Wzones - Wzones: Spatial change analysis
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
library(gridExtra)
library(ggpolypath)

# Data and map preparation ------------------------------------------------
data <- readRDS("data/syria-iraq.rds")

data <- data %>% mutate(deaths_total = deaths_battle + deaths_civilians + 
                          deaths_unknown)

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))


# 2010: Wzones and hotspots ----------------------------------------------------
sf::sf_use_s2(TRUE)

points.sf.2010 <- points.sf %>% filter(year == 2010)

# wzoneData
wzone.2010 <- readOGR(paste0('data/replication-kikuta/2010/2010_12_31.shp'), stringsAsFactors = F)
wzone.2010.sel <- wzone.2010[wzone.2010$dyad_id == 754 |       
                               wzone.2010$dyad_id == 781,]     


wzone.2010.sel.un <- unionSpatialPolygons(wzone.2010.sel, rep(1, length(wzone.2010.sel)))

wzone.2010.sel.un.sf <- st_as_sf(wzone.2010.sel.un, crs = st_crs(4326))

grid.clip.2010 <- st_make_grid(wzone.2010.sel.un.sf, 
                               cellsize = c(0.5,0.5), #cell size in degrees
                               what="polygons") %>%
  st_intersection(wzone.2010.sel.un.sf)

points.count.2010 <- st_intersects(grid.clip.2010, points.sf.2010)

grid.count.2010 <- st_sf(pt_count = lengths(points.count.2010), 
                         geometry = st_cast(grid.clip.2010, "MULTIPOLYGON"))

neighbors.2010 <- poly2nb(grid.count.2010, queen = T)
weighted.neighbors.2010 <- nb2listw(include.self(neighbors.2010), 
                                    style = "B", 
                                    zero.policy = T)
grid.count.2010$HOTSPOT <- as.vector(localG_perm(grid.count.2010$pt_count, 
                                                 weighted.neighbors.2010))

grid.count.2010$z_cat <- "Not significant or cold spot"
grid.count.2010$z_cat[grid.count.2010$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010$z_cat[grid.count.2010$HOTSPOT < 3.091 & grid.count.2010$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2010, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in Syria and Iraq - Wzones, 2010: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-wzones-2010-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2010.sel <- grid.count.2010[grid.count.2010$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2010$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010 <- st_union(grid.count.2010.sel)
st_write(hotspot.2010, "out/syria-iraq-wzones-2010-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2010,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = wzone.2010.sel.un.sf, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq - Wzones, 2010: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-wzones-2010-events.png", width = 8, height = 4)



# 2013: Wzones and hotspots --------------------------------------------------------

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

grid.clip.2013 <- st_make_grid(wzone.2013.sel.un.sf, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(wzone.2013.sel.un.sf)

points.count.2013 <- st_intersects(grid.clip.2013, points.sf.2013)
grid.count.2013 <- st_sf(pt_count = lengths(points.count.2013), 
                         geometry = st_cast(grid.clip.2013, "MULTIPOLYGON"))

neighbors.2013 <- poly2nb(grid.count.2013, queen = T)
weighted.neighbors.2013 <- nb2listw(include.self(neighbors.2013), 
                                    style = "B", 
                                    zero.policy = T)

grid.count.2013$HOTSPOT <- as.vector(localG_perm(grid.count.2013$pt_count, 
                                                 weighted.neighbors.2013))
grid.count.2013$z_cat <- "Not significant or cold spot"
grid.count.2013$z_cat[grid.count.2013$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013$z_cat[grid.count.2013$HOTSPOT < 3.091 & grid.count.2013$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013$z_cat)

# wzones and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2013, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in Syria and Iraq - Wzones, 2013: Conflict shape and Hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-wzones-2013-conflict-shape-hotspot.png", width = 8, height = 4)

# hostpots
grid.count.2013.sel <- grid.count.2013[grid.count.2013$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2013$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013 <- st_union(grid.count.2013.sel)
st_write(hotspot.2013, "out/syria-iraq-wzones-2013-hotspot.shp")



# 2010-2013: Spatial change -----------------------------------------------

# map
label.2010 <- data.frame(
  x = 37, 
  y = 40,
  text = "2010")
label.2013 <- data.frame(
  x = 33, 
  y = 33, 
  text = "2013")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2010, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2010.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2013, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2013.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq - Wzones, 2010-2013: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(31, 51), ylim = c(29, 42), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2010, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2013, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/syria-iraq-wzones-2010-2013-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2010.2013 <- st_area(wzone.2010.sel.un.sf)/100
new.perc.2010.2013 <- st_area(wzone.2013.sel.un.sf)/one.perc.2010.2013
new.perc.2010.2013 <- drop_units(new.perc.2010.2013)
perc.change.2010.2013 <- new.perc.2010.2013 - 100
perc.change.2010.2013 #  29.86264 (expansion)

# overlap: hotspots
hotspot.2010.sp <- as(st_geometry(hotspot.2010), "Spatial")
hotspot.2013.sp <- as(st_geometry(hotspot.2013), "Spatial")
intersection.h.2010.2013 <- gIntersection(hotspot.2010.sp, hotspot.2013.sp)
intersection.h.2010.2013.sf <- st_as_sf(intersection.h.2010.2013) # no overlap
overl_hotspot.2011.2013 <- 0
overl_hotspot.2011.2013

# overlap shape
buffer.2010.sp <- as(st_geometry(wzone.2010.sel.un.sf), "Spatial")
buffer.2013.sp <- as(st_geometry(wzone.2013.sel.un.sf), "Spatial")
intersection.2010.2013 <- gIntersection(buffer.2010.sp, buffer.2013.sp)
intersection.2010.2013.sf <- st_as_sf(intersection.2010.2013)
overl_shape.2011.2013 <- st_area(intersection.2010.2013.sf)/(st_area(wzone.2010.sel.un.sf)/100)
overl_shape.2011.2013 # 75.08904

# spatial change
spatial_ch.2011.2013 <- (overl_hotspot.2011.2013 + drop_units(overl_shape.2011.2013))/2
spatial_ch.2011.2013 # 37.54452 (shift)


# 2014: Wzones and hotspots --------------------------------------------------------

points.sf.2014 <- points.sf %>% filter(year == 2014)

# wzoneData
wzone.2014 <- readOGR(paste0('data/replication-kikuta/2014/2014_12_31.shp'), stringsAsFactors = F)

wzone.2014.sel <- wzone.2014[wzone.2014$dyad_id == 524 |       
                               wzone.2014$dyad_id == 781, ]    

wzone.2014.sel.un <- unionSpatialPolygons(wzone.2014.sel, rep(1, length(wzone.2014.sel)))

wzone.2014.sel.un.sf <- st_as_sf(wzone.2014.sel.un, crs = st_crs(4326))

grid.clip.2014 <- st_make_grid(wzone.2014.sel.un.sf, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(wzone.2014.sel.un.sf)

points.count.2014 <- st_intersects(grid.clip.2014, points.sf.2014)
grid.count.2014 <- st_sf(pt_count = lengths(points.count.2014), 
                         geometry = st_cast(grid.clip.2014, "MULTIPOLYGON"))

neighbors.2014 <- poly2nb(grid.count.2014, queen = T)
weighted.neighbors.2014 <- nb2listw(include.self(neighbors.2014), 
                                    style = "B", 
                                    zero.policy = T)

grid.count.2014$HOTSPOT <- as.vector(localG_perm(grid.count.2014$pt_count, 
                                                 weighted.neighbors.2014))
grid.count.2014$z_cat <- "Not significant or cold spot"
grid.count.2014$z_cat[grid.count.2014$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2014$z_cat[grid.count.2014$HOTSPOT < 3.091 & grid.count.2014$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2014$z_cat)

# wzones and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2014, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed Conflict in Syria and Iraq - Wzones, 2014: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-wzones-2014-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2014.sel <- grid.count.2014[grid.count.2014$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2014$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014 <- st_union(grid.count.2014.sel)
st_write(hotspot.2014, "out/syria-iraq-wzones-2014-hotspot.shp")


# 2016: Wzones and hotspots --------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

wzone.2016 <- readOGR(paste0('data/replication-kikuta/2016/2016_12_31.shp'), stringsAsFactors = F)

wzone.2016.sel <- wzone.2016[wzone.2016$dyad_id == 406 |        
                               wzone.2016$dyad_id == 524 |      
                               wzone.2016$dyad_id == 754 |      
                               wzone.2016$dyad_id == 781 |     
                               wzone.2016$dyad_id == 891 |      
                               wzone.2016$dyad_id == 11973 |      
                               wzone.2016$dyad_id == 14620 |      
                               wzone.2016$dyad_id == 14637 |      
                               wzone.2016$dyad_id == 14701,]      

wzone.2016.sel.un <- unionSpatialPolygons(wzone.2016.sel, rep(1, length(wzone.2016.sel)))
wzone.2016.sel.un.sf <- st_as_sf(wzone.2016.sel.un, crs = st_crs(4326))

grid.clip.2016 <- st_make_grid(wzone.2016.sel.un.sf, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(wzone.2016.sel.un.sf)


points.count.2016 <- st_intersects(grid.clip.2016, points.sf.2016)
grid.count.2016 <- st_sf(pt_count = lengths(points.count.2016), 
                         geometry = st_cast(grid.clip.2016, "MULTIPOLYGON"))

neighbors.2016 <- poly2nb(grid.count.2016, queen = T)
weighted.neighbors.2016 <- nb2listw(include.self(neighbors.2016), 
                                    style = "B", 
                                    zero.policy = T)

grid.count.2016$HOTSPOT <- as.vector(localG_perm(grid.count.2016$pt_count, 
                                                 weighted.neighbors.2016))

grid.count.2016$z_cat <- "Not significant or cold spot"
grid.count.2016$z_cat[grid.count.2016$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2016$z_cat[grid.count.2016$HOTSPOT < 3.091 & grid.count.2016$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2016$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2016, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed Conflict in Syria and Iraq - Wzones, 2016: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-wzones-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/syria-iraq-wzones-2016-hotspot.shp")


# 2014-2016: Spatial change -----------------------------------------------
label.2014 <- data.frame(
  x = 42, 
  y = 40,
  text = "2014")
label.2016 <- data.frame(
  x = 33, 
  y = 34, 
  text = "2016")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2014, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2014.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2016, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2016.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed Conflict in Syria and Iraq - Wzones, 2014-2016: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(31, 51), ylim = c(29, 42), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2014, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2016, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/syria-iraq-wzones-2014-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2014.2016 <- st_area(wzone.2014.sel.un.sf)/100
new.perc.2014.2016 <- st_area(wzone.2016.sel.un.sf)/one.perc.2014.2016
new.perc.2014.2016 <- drop_units(new.perc.2014.2016)
perc.change.2014.2016 <- new.perc.2014.2016 - 100
perc.change.2014.2016 # 50.17533 (expansion)

# overlap: hotspots
hotspot.2014.sp <- as(st_geometry(hotspot.2014), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2014.2016 <- gIntersection(hotspot.2014.sp, hotspot.2016.sp)
intersection.h.2014.2016.sf <- st_as_sf(intersection.h.2014.2016)
overl_hotspot.2014.2016 <- st_area(intersection.h.2014.2016.sf)/(st_area(hotspot.2014)/100) 
overl_hotspot.2014.2016  # 20.45896 

# overlap: conflict shape
buffer.2014.sp <- as(st_geometry(wzone.2014.sel.un.sf), "Spatial")
buffer.2016.sp <- as(st_geometry(wzone.2016.sel.un.sf), "Spatial")
intersection.2014.2016 <- gIntersection(buffer.2014.sp, buffer.2016.sp)
intersection.2014.2016.sf <- st_as_sf(intersection.2014.2016)
overl_shape.2014.2016 <- st_area(intersection.2014.2016.sf)/(st_area(wzone.2014.sel.un.sf)/100) 
overl_shape.2014.2016  # 87.61377 

# spatial change
spatial_ch.2014.2016 <- (overl_hotspot.2014.2016 + overl_shape.2014.2016)/2
spatial_ch.2014.2016  # 54.0106 (no shift)

