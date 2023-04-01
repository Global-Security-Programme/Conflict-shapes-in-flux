### Preamble ###############################################################
# Islamist insurgency in Nigeria - Wzones: Spatial change analysis
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
data <- readRDS("data/nigeria.rds")

data <- data %>% mutate(deaths_total = deaths_battle + deaths_civilians + 
                          deaths_unknown)

data$side_a[data$side_a_new_id == 1051] <- "Boko Haram"
data$side_b[data$side_b_new_id == 1051] <- "Boko Haram"

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))



# 2011: Wzones and hotspots ----------------------------------------------------

points.sf.2011 <- points.sf %>% filter(year == 2011)

# wzoneData
wzone.2011 <- readOGR(paste0('data/replication-kikuta/2011/2011_12_31.shp'), stringsAsFactors = F)
wzone.2011.sel <- wzone.2011[wzone.2011$dyad_id == 640 |        
                               wzone.2011$dyad_id == 5581 |      
                               wzone.2011$dyad_id == 5558,]      

wzone.2011.sel.un <- unionSpatialPolygons(wzone.2011.sel, rep(1, length(wzone.2011.sel)))
wzone.2011.sel.un.sf <- st_as_sf(wzone.2011.sel.un, crs = st_crs(4326))

grid.clip.2011 <- st_make_grid(wzone.2011.sel.un.sf, 
                               cellsize = c(0.5,0.5), #cell size in degrees
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
  ggtitle("Islamist insurgency, 2011: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/nigeria-wzones-2011-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2011.sel <- grid.count.2011[grid.count.2011$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2011$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011 <- st_union(grid.count.2011.sel)
st_write(hotspot.2011, "out/nigeria-wzones-2011-hotspot.shp")


# 2016: Wzones and hotspots --------------------------------------------------------

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
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Islamist insurgency, 2016: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/nigeria-wzones-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/nigeria-wzones-2016-hotspot.shp")


# 2011-2016: Spatial change -----------------------------------------------

# map
label.2011 <- data.frame(
  x = 7, 
  y = 8,
  text = "2011")
label.2016 <- data.frame(
  x = 14.5, 
  y = 15, 
  text = "2016")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2011, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = wzone.2011.sel.un.sf, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2016, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = wzone.2016.sel.un.sf, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Islamist insurgency, 2011-2016: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2011, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2016, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/nigeria-wzones-2011-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2011.2016 <- st_area(wzone.2011.sel.un.sf)/100
new.perc.2011.2016 <- st_area(wzone.2016.sel.un.sf)/one.perc.2011.2016
new.perc.2011.2016 <- drop_units(new.perc.2011.2016)
perc.change.2011.2016 <- new.perc.2011.2016 - 100
perc.change.2011.2016 # 60.45466 (expansion)

# overlap: hotspots
hotspot.2011.sp <- as(st_geometry(hotspot.2011), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2011.2016 <- gIntersection(hotspot.2011.sp, hotspot.2016.sp)
intersection.h.2011.2016.sf <- st_as_sf(intersection.h.2011.2016)
overl_hotspot.2011.2016 <- st_area(intersection.h.2011.2016.sf)/(st_area(hotspot.2016)/100) 
overl_hotspot.2011.2016 # 38.04669

# overlap: conflict shape
buffer.2011.sp <- as(st_geometry(wzone.2011.sel.un.sf), "Spatial")
buffer.2016.sp <- as(st_geometry(wzone.2016.sel.un.sf), "Spatial")
intersection.2011.2016 <- gIntersection(buffer.2011.sp, buffer.2016.sp)
intersection.2011.2016.sf <- st_as_sf(intersection.2011.2016)
overl_shape.2011.2016 <- st_area(intersection.2011.2016.sf)/(st_area(wzone.2011.sel.un.sf)/100) 
overl_shape.2011.2016 # 65.29603

# spatial change
spatial_ch.2011.2016 <- (overl_hotspot.2011.2016 + overl_shape.2011.2016)/2
spatial_ch.2011.2016 # 51.67136 (no shift)

