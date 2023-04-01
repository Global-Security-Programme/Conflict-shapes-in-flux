### Preamble ###############################################################
# Armed conflict in Colombia - WZONES: Spatial change analysis
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
library(rgdal)
library(maptools)
library(units)
library(tidyverse)

# Data and map preparation ------------------------------------------------
data <- readRDS("data/colombia.rds")

data <- data %>% mutate(deaths_total = deaths_battle + deaths_civilians + 
                          deaths_unknown)

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))


# 2006: Wzones and hotspots ----------------------------------------------------

points.sf.2006 <- points.sf %>% filter(year == 2006)

# wzoneData
wzone.2006 <- readOGR(paste0('data/replication-kikuta/2006/2006_12_31.shp'), stringsAsFactors = F)
wzone.2006.sel <- wzone.2006[wzone.2006$dyad_id == 623 |        
                               wzone.2006$dyad_id == 624 |     
                               wzone.2006$dyad_id == 15558 |     
                               wzone.2006$dyad_id == 5196,]      

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
  ggtitle("Armed conflict in Colombia - Wzones, 2006: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/colombia-wzones-2006-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2006.sel <- grid.count.2006[grid.count.2006$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2006$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006 <- st_union(grid.count.2006.sel)
st_write(hotspot.2006, "out/colombia-wzones-2006-hotspot.shp")


# 2011: Conflict shape and hotspots --------------------------------------------------------

points.sf.2011 <- points.sf %>% filter(year == 2011)

# wzoneData
wzone.2011 <- readOGR(paste0('data/replication-kikuta/2011/2011_12_31.shp'), stringsAsFactors = F)
wzone.2011.sel <- wzone.2011[wzone.2011$dyad_id == 623,]      

wzone.2011.sel.un <- unionSpatialPolygons(wzone.2011.sel, rep(1, length(wzone.2011.sel)))
wzone.2011.sel.un.sf <- st_as_sf(wzone.2011.sel.un, crs = st_crs(4326))

grid.clip.2011 <- st_make_grid(wzone.2011.sel.un.sf, 
                               cellsize = c(0.5,0.5),
                               what="polygons") %>%
  st_intersection(wzone.2011.sel.un.sf)

points.count.2011 <- st_intersects(grid.clip.2011, points.sf.2011)
grid.count.2011 <- st_sf(pt_count = lengths(points.count.2011), 
                         geometry = st_cast(grid.clip.2011, "MULTIPOLYGON"))

neighbors.2011 <- poly2nb(grid.count.2011, queen = T)
weighted.neighbors.2011 <- nb2listw(include.self(neighbors.2011), 
                                    style = "B", 
                                    zero.policy = T)

grid.count.2011$HOTSPOT <- as.vector(localG_perm(grid.count.2011$pt_count, 
                                                 weighted.neighbors.2011))
grid.count.2011$z_cat <- "Not significant or cold spot"
grid.count.2011$z_cat[grid.count.2011$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2011$z_cat[grid.count.2011$HOTSPOT < 3.091 & grid.count.2011$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2011$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2011, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in Colombia - Wzones, 2011: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/colombia-wzones-2011-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2011.sel <- grid.count.2011[grid.count.2011$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2011$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011 <- st_union(grid.count.2011.sel)
st_write(hotspot.2011, "out/colombia-wzones-2011-hotspot.shp")


# 2006-2011: Spatial change -----------------------------------------------

# map
label.2006 <- data.frame(
  x = -75, 
  y = 12,
  text = "2006")
label.2011 <- data.frame(
  x = -80.5, 
  y = 2, 
  text = "2011")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2006, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2006.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2011, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2011.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in Colombia - Wzones, 2006-2011: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2006, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2011, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/colombia-wzones-2006-2011-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2006.2011 <- st_area(wzone.2006.sel.un.sf)/100
new.perc.2006.2011 <- st_area(wzone.2011.sel.un.sf)/one.perc.2006.2011
new.perc.2006.2011 <- drop_units(new.perc.2006.2011)
perc.change.2006.2011 <- new.perc.2006.2011 - 100
perc.change.2006.2011 # -71.2088(contraction)

# overlap: hotspots
hotspot.2006.sp <- as(st_geometry(hotspot.2006), "Spatial")
hotspot.2011.sp <- as(st_geometry(hotspot.2011), "Spatial")
intersection.h.2006.2011 <- gIntersection(hotspot.2006.sp, hotspot.2011.sp)
intersection.h.2006.2011.sf <- st_as_sf(intersection.h.2006.2011)
overl_hotspot.2006.2011 <- st_area(intersection.h.2006.2011.sf)/(st_area(hotspot.2011)/100) 
overl_hotspot.2006.2011 # 90.04717

# overlap: conflict shape
buffer.2006.sp <- as(st_geometry(wzone.2006.sel.un.sf), "Spatial")
buffer.2011.sp <- as(st_geometry(wzone.2011.sel.un.sf), "Spatial")
intersection.2006.2011 <- gIntersection(buffer.2006.sp, buffer.2011.sp)
intersection.2006.2011.sf <- st_as_sf(intersection.2006.2011)
overl_shape.2006.2011 <- st_area(intersection.2006.2011.sf)/(st_area(wzone.2011.sel.un.sf)/100) # 88.80934 [1]
overl_shape.2006.2011 # 100

# spatial change
spatial_ch.2006.2011 <- (overl_hotspot.2006.2011 + overl_shape.2006.2011)/2
spatial_ch.2006.2011 # 95.02358 (no shift)



# 2012: Wzones and hotspots  --------------------------------------

points.sf.2012 <- points.sf %>% filter(year == 2012)

# wzoneData
wzone.2012 <- readOGR(paste0('data/replication-kikuta/2012/2012_12_31.shp'), stringsAsFactors = F)
wzone.2012.sel <- wzone.2012[wzone.2012$dyad_id == 623,]      

wzone.2012.sel.un <- unionSpatialPolygons(wzone.2012.sel, rep(1, length(wzone.2012.sel)))
wzone.2012.sel.un.sf <- st_as_sf(wzone.2012.sel.un, crs = st_crs(4326))

grid.clip.2012 <- st_make_grid(wzone.2012.sel.un.sf, 
                               cellsize = c(0.5,0.5), #cell size in degrees
                               what="polygons") %>%
  st_intersection(wzone.2012.sel.un.sf)

points.count.2012 <- st_intersects(grid.clip.2012, points.sf.2012)
grid.count.2012 <- st_sf(pt_count = lengths(points.count.2012), 
                         geometry = st_cast(grid.clip.2012, "MULTIPOLYGON"))

neighbors.2012 <- poly2nb(grid.count.2012, queen = T)
weighted.neighbors.2012 <- nb2listw(include.self(neighbors.2012), 
                                    style = "B", 
                                    zero.policy = T)
grid.count.2012$HOTSPOT <- as.vector(localG_perm(grid.count.2012$pt_count, 
                                                 weighted.neighbors.2012))

grid.count.2012$z_cat <- "Not significant or cold spot"
grid.count.2012$z_cat[grid.count.2012$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2012$z_cat[grid.count.2012$HOTSPOT < 3.091 & grid.count.2012$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2012$z_cat)

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2012, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in Colombia - Wzones, 2012: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/colombia-wzones-2012-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2012.sel <- grid.count.2012[grid.count.2012$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2012$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2012 <- st_union(grid.count.2012.sel)
st_write(hotspot.2012, "out/colombia-wzones-2012-hotspot.shp")


# 2016: Wzones and hotspots --------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

# wzoneData
wzone.2016 <- readOGR(paste0('data/replication-kikuta/2016/2016_12_31.shp'), stringsAsFactors = F)
wzone.2016.sel <- wzone.2016[wzone.2016$dyad_id == 623 |        
                               wzone.2016$dyad_id == 624 |      
                               wzone.2016$dyad_id == 15558,]    

wzone.2016.sel.un <- unionSpatialPolygons(wzone.2016.sel, rep(1, length(wzone.2016.sel)))
wzone.2016.sel.un.sf <- st_as_sf(wzone.2016.sel.un, crs = st_crs(4326))

grid.clip.2016 <- st_make_grid(wzone.2016.sel.un.sf, 
                               cellsize = c(0.5,0.5), #cell size in degrees
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
  scale_fill_manual(values = c("orange", "grey90")) +
  ggtitle("Armed conflict in Colombia - Wzones, 2016: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/colombia-wzones-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/colombia-wzones-2016-hotspot.shp")


# 2012-2016: Spatial change -----------------------------------------------

# map
label.2012 <- data.frame(
  x = -75, 
  y = 12,
  text = "2012")
label.2016 <- data.frame(
  x = -80.5, 
  y = 2, 
  text = "2016")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2012, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2012.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2016, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2016.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in Colombia - Wzones, 2012-2016: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(-82, -65), ylim = c(-3, 15), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2012, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2016, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/colombia-wzones-2012-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2012.2016 <- st_area(wzone.2012.sel.un.sf)/100
new.perc.2012.2016 <- st_area(wzone.2016.sel.un.sf)/one.perc.2012.2016
new.perc.2012.2016 <- drop_units(new.perc.2012.2016)
perc.change.2012.2016 <- new.perc.2012.2016 - 100
perc.change.2012.2016 # 12.58297 (expansion)

# overlap: hotspots
hotspot.2012.sp <- as(st_geometry(hotspot.2012), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2012.2016 <- gIntersection(hotspot.2012.sp, hotspot.2016.sp)
intersection.h.2012.2016.sf <- st_as_sf(intersection.h.2012.2016) # no overlap
overl_hotspot.2012.2016 <- 0
overl_hotspot.2012.2016 # 0

# overlap: conflict shape
buffer.2012.sp <- as(st_geometry(wzone.2012.sel.un.sf), "Spatial")
buffer.2016.sp <- as(st_geometry(wzone.2016.sel.un.sf), "Spatial")
intersection.2012.2016 <- gIntersection(buffer.2012.sp, buffer.2016.sp)
intersection.2012.2016.sf <- st_as_sf(intersection.2012.2016)
overl_shape.2012.2016 <- st_area(intersection.2012.2016.sf)/(st_area(wzone.2012.sel.un.sf)/100) 
overl_shape.2012.2016 # 65.10216

# spatial change
spatial_ch.2012.2016 <- (overl_hotspot.2012.2016 + drop_units(overl_shape.2012.2016))/2
spatial_ch.2012.2016 # 32.55108 (shift)

