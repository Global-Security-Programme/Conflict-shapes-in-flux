### Preamble ###############################################################
# Armed conflict in the Afghan-Pakistani borderlands: Spatial change analysis
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


# 2006: Dominant actors -----------------------------------------------------------

data.2006 <- data %>% filter(year == 2006)

data.2006 <- data.2006 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2006 <- data.2006 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2006 <- edges.2006 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2006 <- edges.2006 %>% distinct(from, to, weight)

side.a.2006 <- data.2006 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2006 <- data.2006 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2006 <- rbind(side.a.2006, side.b.2006)
nodes.2006 <- nodes.2006 %>% distinct(id, label)

net.2006 <- graph_from_data_frame(edges.2006, 
                                  vertices = nodes.2006,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2006$nodes_n <- igraph::degree(net.2006)
nodes.2006$degree <- strength(net.2006)
nodes.2006$eigenvec <- round(igraph::evcent(net.2006)$vector,2)
nodes.2006$katz <- round(katzcent(net.2006),2)
write_csv(nodes.2006, "out/afghanistan-pakistan-2006-centrality-measures.csv")

net.2006 <- graph_from_data_frame(edges.2006, 
                                  vertices = nodes.2006,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2006,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 10, 
                         area = 100),
       layout.exp = 0.5,
       label = paste(V(net.2006)$label,":", V(net.2006)$degree),
       label.size = 4,
       color = "blue",
       alpha = 0.3,
       size = 5, 
       legend.position = "none",
       edge.size = 0.5,
       edge.color = "grey",
       edge.label = "weight",
       edge.label.size = 4,
       edge.label.color = "blue",
       edge.label.alpha = 0.6) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2006")

ggsave("figs/afghanistan-pakistan-2006-network.png", width = 10, height = 3)

# 2006: Conflict shape and hotspots --------------------------------------------------------

points.sf.2006 <- points.sf %>% filter(year == 2006)

# conflict shape
concave.hull.2006 <- concaveman(points.sf.2006, concavity = 2, length_threshold = 0)
buffer.2006 <- st_buffer(concave.hull.2006, dist = 0.5)
st_write(buffer.2006, "out/afghanistan-pakistan-2006-conflict-shape.shp")

grid.clip.2006 <- st_make_grid(buffer.2006, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(buffer.2006)

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
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2006: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-2006-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2006.sel <- grid.count.2006[grid.count.2006$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2006$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006 <- st_union(grid.count.2006.sel)
st_write(hotspot.2006, "out/afghanistan-pakistan-2006-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2006,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2006, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2006: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-2006-events.png", width = 8, height = 4)


# 2008: Dominant actors -----------------------------------------------------------

data.2008 <- data %>% filter(year == 2008)

data.2008 <- data.2008 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2008 <- data.2008 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2008 <- edges.2008 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2008 <- edges.2008 %>% distinct(from, to, weight)

side.a.2008 <- data.2008 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2008 <- data.2008 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2008 <- rbind(side.a.2008, side.b.2008)
nodes.2008 <- nodes.2008 %>% distinct(id, label)

net.2008 <- graph_from_data_frame(edges.2008, 
                                  vertices = nodes.2008,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2008$nodes_n <- igraph::degree(net.2008)
nodes.2008$degree <- strength(net.2008)
nodes.2008$eigenvec <- round(igraph::evcent(net.2008)$vector,2)
nodes.2008$katz <- round(katzcent(net.2008),2)
write_csv(nodes.2008, "out/afghanistan-pakistan-2008-centrality-measures.csv")

net.2008 <- graph_from_data_frame(edges.2008, 
                                  vertices = nodes.2008,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2008,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 10, 
                         area = 100),
       layout.exp = 0.5,
       label = paste(V(net.2008)$label,":", V(net.2008)$degree),
       label.size = 4,
       color = "blue",
       alpha = 0.3,
       size = 5, 
       legend.position = "none",
       edge.size = 0.5,
       edge.color = "grey",
       edge.label = "weight",
       edge.label.size = 4,
       edge.label.color = "blue",
       edge.label.alpha = 0.6) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands,2008")
ggsave("figs/afghanistan-pakistan-2008-network.png", width = 10, height = 4)


# 2008: Conflict shape and hotspots --------------------------------------------------------

points.sf.2008 <- points.sf %>% filter(year == 2008)

# conflict shape
concave.hull.2008 <- concaveman(points.sf.2008, concavity = 2, length_threshold = 0)
buffer.2008 <- st_buffer(concave.hull.2008, dist = 0.5)
st_write(buffer.2008, "out/afghanistan-pakistan-2008-conflict-shape.shp")

grid.clip.2008 <- st_make_grid(buffer.2008, 
                               cellsize = c(0.5,0.5), #
                               what="polygons") %>%
  st_intersection(buffer.2008)

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
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2008: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-2008-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2008.sel <- grid.count.2008[grid.count.2008$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2008$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008 <- st_union(grid.count.2008.sel)
st_write(hotspot.2008, "out/afghanistan-pakistan-2008-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2008,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2008, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2008: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(58, 76), ylim = c(24, 42), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/afghanistan-pakistan-2008-events.png", width = 8, height = 4)


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
  geom_sf(data = buffer.2006, alpha = 0.3, fill = "blue", color = NA) +
  
  geom_sf(data = hotspot.2008, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = buffer.2008, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in the Afghan-Pakistani borderlands, 2006-2008: Spatial change") +
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
ggsave("figs/afghanistan-pakistan-2006-2008-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2006.2008 <- st_area(buffer.2006)/100
new.perc.2006.2008 <- st_area(buffer.2008)/one.perc.2006.2008
new.perc.2006.2008 <- drop_units(new.perc.2006.2008)
perc.change.2006.2008 <- new.perc.2006.2008 - 100
perc.change.2006.2008 # 87.50662 (expansion)

# overlap: hotspots
hotspot.2006.sp <- as(st_geometry(hotspot.2006), "Spatial")
hotspot.2008.sp <- as(st_geometry(hotspot.2008), "Spatial")
intersection.h.2006.2008 <- gIntersection(hotspot.2006.sp, hotspot.2008.sp)
intersection.h.2006.2008.sf <- st_as_sf(intersection.h.2006.2008)
overl_hotspot.2006.2008 <- st_area(intersection.h.2006.2008.sf)/(st_area(hotspot.2006)/100)  
overl_hotspot.2006.2008 # 56.09073

# overlap: conflict shape
buffer.2006.sp <- as(st_geometry(buffer.2006), "Spatial")
buffer.2008.sp <- as(st_geometry(buffer.2008), "Spatial")
intersection.2006.2008 <- gIntersection(buffer.2006.sp, buffer.2008.sp)
intersection.2006.2008.sf <- st_as_sf(intersection.2006.2008)
overl_shape.2006.2008 <- st_area(intersection.2006.2008.sf)/(st_area(buffer.2006)/100) 
overl_shape.2006.2008 # 96.35396

# spatial change
spatial_ch.2006.2008 <- (overl_hotspot.2006.2008 + overl_shape.2006.2008)/2
spatial_ch.2006.2008 # 76.22234 (no shift)


# 2006.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2006.GRID75 <- points.sf %>% filter(year == 2006)

# conflict shape
concave.hull.2006.GRID75 <- concaveman(points.sf.2006.GRID75, concavity = 2, length_threshold = 0)
buffer.2006.GRID75 <- st_buffer(concave.hull.2006.GRID75, dist = 0.5)
st_write(buffer.2006.GRID75, "out/afghanistan-pakistan-conflict-2006_GRID75-shape.shp")

grid.clip.2006.GRID75 <- st_make_grid(buffer.2006.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2006.GRID75)

points.count.2006.GRID75 <- st_intersects(grid.clip.2006.GRID75, points.sf.2006.GRID75)
grid.count.2006.GRID75 <- st_sf(pt_count = lengths(points.count.2006.GRID75), 
                                geometry = st_cast(grid.clip.2006.GRID75, "MULTIPOLYGON"))

neighbors.2006.GRID75 <- poly2nb(grid.count.2006.GRID75, queen = T)
weighted.neighbors.2006.GRID75 <- nb2listw(include.self(neighbors.2006.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2006.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2006.GRID75$pt_count, 
                                                        weighted.neighbors.2006.GRID75))

grid.count.2006.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2006.GRID75$z_cat[grid.count.2006.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2006.GRID75$z_cat[grid.count.2006.GRID75$HOTSPOT < 3.091 & grid.count.2006.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2006.GRID75$z_cat)

# hotspots
grid.count.2006.GRID75.sel <- grid.count.2006.GRID75[grid.count.2006.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2006.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006.GRID75 <- st_union(grid.count.2006.GRID75.sel)
st_write(hotspot.2006.GRID75, "out/afghanistan-pakistan-2006_GRID75-hotspot.shp")


# 2008.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2008.GRID75 <- points.sf %>% filter(year == 2008)

# conflict shape
concave.hull.2008.GRID75 <- concaveman(points.sf.2008.GRID75, concavity = 2, length_threshold = 0)
buffer.2008.GRID75 <- st_buffer(concave.hull.2008.GRID75, dist = 0.5)
st_write(buffer.2008.GRID75, "out/afghanistan-pakistan-conflict-2008_GRID75-shape.shp")

grid.clip.2008.GRID75 <- st_make_grid(buffer.2008.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2008.GRID75)

points.count.2008.GRID75 <- st_intersects(grid.clip.2008.GRID75, points.sf.2008.GRID75)
grid.count.2008.GRID75 <- st_sf(pt_count = lengths(points.count.2008.GRID75), 
                                geometry = st_cast(grid.clip.2008.GRID75, "MULTIPOLYGON"))

neighbors.2008.GRID75 <- poly2nb(grid.count.2008.GRID75, queen = T)
weighted.neighbors.2008.GRID75 <- nb2listw(include.self(neighbors.2008.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2008.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2008.GRID75$pt_count, 
                                                        weighted.neighbors.2008.GRID75))

grid.count.2008.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2008.GRID75$z_cat[grid.count.2008.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2008.GRID75$z_cat[grid.count.2008.GRID75$HOTSPOT < 3.091 & grid.count.2008.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2008.GRID75$z_cat)

# hotspots
grid.count.2008.GRID75.sel <- grid.count.2008.GRID75[grid.count.2008.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2008.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008.GRID75 <- st_union(grid.count.2008.GRID75.sel)
st_write(hotspot.2008.GRID75, "out/afghanistan-pakistan-2008_GRID75-hotspot.shp")


# 2006.GRID75-2008.GRID75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2006.GRID75.2008.GRID75 <- st_area(buffer.2006.GRID75)/100
new.perc.2006.GRID75.2008.GRID75 <- st_area(buffer.2008.GRID75)/one.perc.2006.GRID75.2008.GRID75
new.perc.2006.GRID75.2008.GRID75 <- drop_units(new.perc.2006.GRID75.2008.GRID75)
perc.change.2006.GRID75.2008.GRID75 <- new.perc.2006.GRID75.2008.GRID75 - 100
perc.change.2006.GRID75.2008.GRID75 #  87.50662 (expansion)

# overlap: hotspots
hotspot.2006.GRID75.sp <- as(st_geometry(hotspot.2006.GRID75), "Spatial")
hotspot.2008.GRID75.sp <- as(st_geometry(hotspot.2008.GRID75), "Spatial")
intersection.h.2006.GRID75.2008.GRID75 <- gIntersection(hotspot.2006.GRID75.sp, hotspot.2008.GRID75.sp)
intersection.h.2006.GRID75.2008.GRID75.sf <- st_as_sf(intersection.h.2006.GRID75.2008.GRID75)
overl_hotspot.2006.GRID75.2008.GRID75 <- st_area(intersection.h.2006.GRID75.2008.GRID75.sf)/(st_area(hotspot.2006.GRID75)/100) 
overl_hotspot.2006.GRID75.2008.GRID75 # 51.90085

# overlap: conflict shape
buffer.2006.GRID75.sp <- as(st_geometry(buffer.2006.GRID75), "Spatial")
buffer.2008.GRID75.sp <- as(st_geometry(buffer.2008.GRID75), "Spatial")
intersection.2006.GRID75.2008.GRID75 <- gIntersection(buffer.2006.GRID75.sp, buffer.2008.GRID75.sp)
intersection.2006.GRID75.2008.GRID75.sf <- st_as_sf(intersection.2006.GRID75.2008.GRID75)
overl_shape.2006.GRID75.2008.GRID75 <- st_area(intersection.2006.GRID75.2008.GRID75.sf)/(st_area(buffer.2006.GRID75)/100) # 88.80934 [1]
overl_shape.2006.GRID75.2008.GRID75 # 96.35396 

# spatial change
spatial_ch.2006.GRID75.2008.GRID75 <- (overl_hotspot.2006.GRID75.2008.GRID75 + overl_shape.2006.GRID75.2008.GRID75)/2
spatial_ch.2006.GRID75.2008.GRID75 # 74.1274 (no shift)


# 2006.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2006.GRID100 <- points.sf %>% filter(year == 2006)

# conflict shape
concave.hull.2006.GRID100 <- concaveman(points.sf.2006.GRID100, concavity = 2, length_threshold = 0)
buffer.2006.GRID100 <- st_buffer(concave.hull.2006.GRID100, dist = 0.5)
st_write(buffer.2006.GRID100, "out/afghanistan-pakistan-conflict-2006_GRID100-shape.shp")

grid.clip.2006.GRID100 <- st_make_grid(buffer.2006.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2006.GRID100)

points.count.2006.GRID100 <- st_intersects(grid.clip.2006.GRID100, points.sf.2006.GRID100)
grid.count.2006.GRID100 <- st_sf(pt_count = lengths(points.count.2006.GRID100), 
                                 geometry = st_cast(grid.clip.2006.GRID100, "MULTIPOLYGON"))

neighbors.2006.GRID100 <- poly2nb(grid.count.2006.GRID100, queen = T)
weighted.neighbors.2006.GRID100 <- nb2listw(include.self(neighbors.2006.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2006.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2006.GRID100$pt_count, 
                                                         weighted.neighbors.2006.GRID100))

grid.count.2006.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2006.GRID100$z_cat[grid.count.2006.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2006.GRID100$z_cat[grid.count.2006.GRID100$HOTSPOT < 3.091 & grid.count.2006.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2006.GRID100$z_cat)

# hotspots
grid.count.2006.GRID100.sel <- grid.count.2006.GRID100[grid.count.2006.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2006.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006.GRID100 <- st_union(grid.count.2006.GRID100.sel)
st_write(hotspot.2006.GRID100, "out/afghanistan-pakistan-2006_GRID100-hotspot.shp")


# 2008.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2008.GRID100 <- points.sf %>% filter(year == 2008)

# conflict shape
concave.hull.2008.GRID100 <- concaveman(points.sf.2008.GRID100, concavity = 2, length_threshold = 0)
buffer.2008.GRID100 <- st_buffer(concave.hull.2008.GRID100, dist = 0.5)
st_write(buffer.2008.GRID100, "out/afghanistan-pakistan-conflict-2008_GRID100-shape.shp")

grid.clip.2008.GRID100 <- st_make_grid(buffer.2008.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2008.GRID100)

points.count.2008.GRID100 <- st_intersects(grid.clip.2008.GRID100, points.sf.2008.GRID100)
grid.count.2008.GRID100 <- st_sf(pt_count = lengths(points.count.2008.GRID100), 
                                 geometry = st_cast(grid.clip.2008.GRID100, "MULTIPOLYGON"))

neighbors.2008.GRID100 <- poly2nb(grid.count.2008.GRID100, queen = T)
weighted.neighbors.2008.GRID100 <- nb2listw(include.self(neighbors.2008.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2008.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2008.GRID100$pt_count, 
                                                         weighted.neighbors.2008.GRID100))

grid.count.2008.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2008.GRID100$z_cat[grid.count.2008.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2008.GRID100$z_cat[grid.count.2008.GRID100$HOTSPOT < 3.091 & grid.count.2008.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2008.GRID100$z_cat)

# hotspots
grid.count.2008.GRID100.sel <- grid.count.2008.GRID100[grid.count.2008.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2008.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008.GRID100 <- st_union(grid.count.2008.GRID100.sel)
st_write(hotspot.2008.GRID100, "out/afghanistan-pakistan-2008_GRID100-hotspot.shp")


# 2006.GRID100-2008.GRID100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2006.GRID100.2008.GRID100 <- st_area(buffer.2006.GRID100)/100
new.perc.2006.GRID100.2008.GRID100 <- st_area(buffer.2008.GRID100)/one.perc.2006.GRID100.2008.GRID100
new.perc.2006.GRID100.2008.GRID100 <- drop_units(new.perc.2006.GRID100.2008.GRID100)
perc.change.2006.GRID100.2008.GRID100 <- new.perc.2006.GRID100.2008.GRID100 - 100
perc.change.2006.GRID100.2008.GRID100 #  87.50662 (expansion)

# overlap: hotspots
hotspot.2006.GRID100.sp <- as(st_geometry(hotspot.2006.GRID100), "Spatial")
hotspot.2008.GRID100.sp <- as(st_geometry(hotspot.2008.GRID100), "Spatial")
intersection.h.2006.GRID100.2008.GRID100 <- gIntersection(hotspot.2006.GRID100.sp, hotspot.2008.GRID100.sp)
intersection.h.2006.GRID100.2008.GRID100.sf <- st_as_sf(intersection.h.2006.GRID100.2008.GRID100)
overl_hotspot.2006.GRID100.2008.GRID100 <- st_area(intersection.h.2006.GRID100.2008.GRID100.sf)/(st_area(hotspot.2006.GRID100)/100) 
overl_hotspot.2006.GRID100.2008.GRID100 # 28.57645 

# overlap: conflict shape
buffer.2006.GRID100.sp <- as(st_geometry(buffer.2006.GRID100), "Spatial")
buffer.2008.GRID100.sp <- as(st_geometry(buffer.2008.GRID100), "Spatial")
intersection.2006.GRID100.2008.GRID100 <- gIntersection(buffer.2006.GRID100.sp, buffer.2008.GRID100.sp)
intersection.2006.GRID100.2008.GRID100.sf <- st_as_sf(intersection.2006.GRID100.2008.GRID100)
overl_shape.2006.GRID100.2008.GRID100 <- st_area(intersection.2006.GRID100.2008.GRID100.sf)/(st_area(buffer.2006.GRID100)/100) # 88.80934 [1]
overl_shape.2006.GRID100.2008.GRID100 # 96.35396

# spatial change
spatial_ch.2006.GRID100.2008.GRID100 <- (overl_hotspot.2006.GRID100.2008.GRID100 + overl_shape.2006.GRID100.2008.GRID100)/2
spatial_ch.2006.GRID100.2008.GRID100 # 62.46521 (no shift)


# 2006.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2006.BUFF75 <- points.sf %>% filter(year == 2006)

# conflict shape
concave.hull.2006.BUFF75 <- concaveman(points.sf.2006.BUFF75, concavity = 2, length_threshold = 0)
buffer.2006.BUFF75 <- st_buffer(concave.hull.2006.BUFF75, dist = 0.75)
st_write(buffer.2006.BUFF75, "out/afghanistan-pakistan-conflict-2006_BUFF75-shape.shp")

grid.clip.2006.BUFF75 <- st_make_grid(buffer.2006.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2006.BUFF75)

points.count.2006.BUFF75 <- st_intersects(grid.clip.2006.BUFF75, points.sf.2006.BUFF75)
grid.count.2006.BUFF75 <- st_sf(pt_count = lengths(points.count.2006.BUFF75), 
                                geometry = st_cast(grid.clip.2006.BUFF75, "MULTIPOLYGON"))

neighbors.2006.BUFF75 <- poly2nb(grid.count.2006.BUFF75, queen = T)
weighted.neighbors.2006.BUFF75 <- nb2listw(include.self(neighbors.2006.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2006.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2006.BUFF75$pt_count, 
                                                        weighted.neighbors.2006.BUFF75))

grid.count.2006.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2006.BUFF75$z_cat[grid.count.2006.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2006.BUFF75$z_cat[grid.count.2006.BUFF75$HOTSPOT < 3.091 & grid.count.2006.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2006.BUFF75$z_cat)

# hotspots
grid.count.2006.BUFF75.sel <- grid.count.2006.BUFF75[grid.count.2006.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2006.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006.BUFF75 <- st_union(grid.count.2006.BUFF75.sel)
st_write(hotspot.2006.BUFF75, "out/afghanistan-pakistan-2006_BUFF75-hotspot.shp")


# 2008.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2008.BUFF75 <- points.sf %>% filter(year == 2008)

# conflict shape
concave.hull.2008.BUFF75 <- concaveman(points.sf.2008.BUFF75, concavity = 2, length_threshold = 0)
buffer.2008.BUFF75 <- st_buffer(concave.hull.2008.BUFF75, dist = 0.75)
st_write(buffer.2008.BUFF75, "out/afghanistan-pakistan-conflict-2008_BUFF75-shape.shp")

grid.clip.2008.BUFF75 <- st_make_grid(buffer.2008.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2008.BUFF75)

points.count.2008.BUFF75 <- st_intersects(grid.clip.2008.BUFF75, points.sf.2008.BUFF75)
grid.count.2008.BUFF75 <- st_sf(pt_count = lengths(points.count.2008.BUFF75), 
                                geometry = st_cast(grid.clip.2008.BUFF75, "MULTIPOLYGON"))

neighbors.2008.BUFF75 <- poly2nb(grid.count.2008.BUFF75, queen = T)
weighted.neighbors.2008.BUFF75 <- nb2listw(include.self(neighbors.2008.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2008.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2008.BUFF75$pt_count, 
                                                        weighted.neighbors.2008.BUFF75))

grid.count.2008.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2008.BUFF75$z_cat[grid.count.2008.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2008.BUFF75$z_cat[grid.count.2008.BUFF75$HOTSPOT < 3.091 & grid.count.2008.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2008.BUFF75$z_cat)

# hotspots
grid.count.2008.BUFF75.sel <- grid.count.2008.BUFF75[grid.count.2008.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2008.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008.BUFF75 <- st_union(grid.count.2008.BUFF75.sel)
st_write(hotspot.2008.BUFF75, "out/afghanistan-pakistan-2008_BUFF75-hotspot.shp")


# 2006.BUFF75-2008.BUFF75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2006.BUFF75.2008.BUFF75 <- st_area(buffer.2006.BUFF75)/100
new.perc.2006.BUFF75.2008.BUFF75 <- st_area(buffer.2008.BUFF75)/one.perc.2006.BUFF75.2008.BUFF75
new.perc.2006.BUFF75.2008.BUFF75 <- drop_units(new.perc.2006.BUFF75.2008.BUFF75)
perc.change.2006.BUFF75.2008.BUFF75 <- new.perc.2006.BUFF75.2008.BUFF75 - 100
perc.change.2006.BUFF75.2008.BUFF75 #  73.93985 (expansion)

# overlap: hotspots
hotspot.2006.BUFF75.sp <- as(st_geometry(hotspot.2006.BUFF75), "Spatial")
hotspot.2008.BUFF75.sp <- as(st_geometry(hotspot.2008.BUFF75), "Spatial")
intersection.h.2006.BUFF75.2008.BUFF75 <- gIntersection(hotspot.2006.BUFF75.sp, hotspot.2008.BUFF75.sp)
intersection.h.2006.BUFF75.2008.BUFF75.sf <- st_as_sf(intersection.h.2006.BUFF75.2008.BUFF75)
overl_hotspot.2006.BUFF75.2008.BUFF75 <- st_area(intersection.h.2006.BUFF75.2008.BUFF75.sf)/(st_area(hotspot.2006.BUFF75)/100) 
overl_hotspot.2006.BUFF75.2008.BUFF75 # 75.05605

# overlap: conflict shape
buffer.2006.BUFF75.sp <- as(st_geometry(buffer.2006.BUFF75), "Spatial")
buffer.2008.BUFF75.sp <- as(st_geometry(buffer.2008.BUFF75), "Spatial")
intersection.2006.BUFF75.2008.BUFF75 <- gIntersection(buffer.2006.BUFF75.sp, buffer.2008.BUFF75.sp)
intersection.2006.BUFF75.2008.BUFF75.sf <- st_as_sf(intersection.2006.BUFF75.2008.BUFF75)
overl_shape.2006.BUFF75.2008.BUFF75 <- st_area(intersection.2006.BUFF75.2008.BUFF75.sf)/(st_area(buffer.2006.BUFF75)/100) # 88.80934 [1]
overl_shape.2006.BUFF75.2008.BUFF75 # 97.02872

# spatial change
spatial_ch.2006.BUFF75.2008.BUFF75 <- (overl_hotspot.2006.BUFF75.2008.BUFF75 + overl_shape.2006.BUFF75.2008.BUFF75)/2
spatial_ch.2006.BUFF75.2008.BUFF75 # 86.04238 (no shift)


# 2006.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2006.BUFF100 <- points.sf %>% filter(year == 2006)

# conflict shape
concave.hull.2006.BUFF100 <- concaveman(points.sf.2006.BUFF100, concavity = 2, length_threshold = 0)
buffer.2006.BUFF100 <- st_buffer(concave.hull.2006.BUFF100, dist = 1)
st_write(buffer.2006.BUFF100, "out/afghanistan-pakistan-conflict-2006_BUFF100-shape.shp")

grid.clip.2006.BUFF100 <- st_make_grid(buffer.2006.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2006.BUFF100)

points.count.2006.BUFF100 <- st_intersects(grid.clip.2006.BUFF100, points.sf.2006.BUFF100)
grid.count.2006.BUFF100 <- st_sf(pt_count = lengths(points.count.2006.BUFF100), 
                                 geometry = st_cast(grid.clip.2006.BUFF100, "MULTIPOLYGON"))

neighbors.2006.BUFF100 <- poly2nb(grid.count.2006.BUFF100, queen = T)
weighted.neighbors.2006.BUFF100 <- nb2listw(include.self(neighbors.2006.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2006.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2006.BUFF100$pt_count, 
                                                         weighted.neighbors.2006.BUFF100))

grid.count.2006.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2006.BUFF100$z_cat[grid.count.2006.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2006.BUFF100$z_cat[grid.count.2006.BUFF100$HOTSPOT < 3.091 & grid.count.2006.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2006.BUFF100$z_cat)

# hotspots
grid.count.2006.BUFF100.sel <- grid.count.2006.BUFF100[grid.count.2006.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2006.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2006.BUFF100 <- st_union(grid.count.2006.BUFF100.sel)
st_write(hotspot.2006.BUFF100, "out/afghanistan-pakistan-2006_BUFF100-hotspot.shp")


# 2008.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2008.BUFF100 <- points.sf %>% filter(year == 2008)

# conflict shape
concave.hull.2008.BUFF100 <- concaveman(points.sf.2008.BUFF100, concavity = 2, length_threshold = 0)
buffer.2008.BUFF100 <- st_buffer(concave.hull.2008.BUFF100, dist = 1)
st_write(buffer.2008.BUFF100, "out/afghanistan-pakistan-conflict-2008_BUFF100-shape.shp")

grid.clip.2008.BUFF100 <- st_make_grid(buffer.2008.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2008.BUFF100)

points.count.2008.BUFF100 <- st_intersects(grid.clip.2008.BUFF100, points.sf.2008.BUFF100)
grid.count.2008.BUFF100 <- st_sf(pt_count = lengths(points.count.2008.BUFF100), 
                                 geometry = st_cast(grid.clip.2008.BUFF100, "MULTIPOLYGON"))

neighbors.2008.BUFF100 <- poly2nb(grid.count.2008.BUFF100, queen = T)
weighted.neighbors.2008.BUFF100 <- nb2listw(include.self(neighbors.2008.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2008.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2008.BUFF100$pt_count, 
                                                         weighted.neighbors.2008.BUFF100))

grid.count.2008.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2008.BUFF100$z_cat[grid.count.2008.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2008.BUFF100$z_cat[grid.count.2008.BUFF100$HOTSPOT < 3.091 & grid.count.2008.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2008.BUFF100$z_cat)

# hotspots
grid.count.2008.BUFF100.sel <- grid.count.2008.BUFF100[grid.count.2008.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2008.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2008.BUFF100 <- st_union(grid.count.2008.BUFF100.sel)
st_write(hotspot.2008.BUFF100, "out/afghanistan-pakistan-2008_BUFF100-hotspot.shp")


# 2006.BUFF100-2008.BUFF100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2006.BUFF100.2008.BUFF100 <- st_area(buffer.2006.BUFF100)/100
new.perc.2006.BUFF100.2008.BUFF100 <- st_area(buffer.2008.BUFF100)/one.perc.2006.BUFF100.2008.BUFF100
new.perc.2006.BUFF100.2008.BUFF100 <- drop_units(new.perc.2006.BUFF100.2008.BUFF100)
perc.change.2006.BUFF100.2008.BUFF100 <- new.perc.2006.BUFF100.2008.BUFF100 - 100
perc.change.2006.BUFF100.2008.BUFF100 #  64.49447 (expansion)

# overlap: hotspots
hotspot.2006.BUFF100.sp <- as(st_geometry(hotspot.2006.BUFF100), "Spatial")
hotspot.2008.BUFF100.sp <- as(st_geometry(hotspot.2008.BUFF100), "Spatial")
intersection.h.2006.BUFF100.2008.BUFF100 <- gIntersection(hotspot.2006.BUFF100.sp, hotspot.2008.BUFF100.sp)
intersection.h.2006.BUFF100.2008.BUFF100.sf <- st_as_sf(intersection.h.2006.BUFF100.2008.BUFF100)
overl_hotspot.2006.BUFF100.2008.BUFF100 <- st_area(intersection.h.2006.BUFF100.2008.BUFF100.sf)/(st_area(hotspot.2006.BUFF100)/100) 
overl_hotspot.2006.BUFF100.2008.BUFF100 # 58.26419

# overlap: conflict shape
buffer.2006.BUFF100.sp <- as(st_geometry(buffer.2006.BUFF100), "Spatial")
buffer.2008.BUFF100.sp <- as(st_geometry(buffer.2008.BUFF100), "Spatial")
intersection.2006.BUFF100.2008.BUFF100 <- gIntersection(buffer.2006.BUFF100.sp, buffer.2008.BUFF100.sp)
intersection.2006.BUFF100.2008.BUFF100.sf <- st_as_sf(intersection.2006.BUFF100.2008.BUFF100)
overl_shape.2006.BUFF100.2008.BUFF100 <- st_area(intersection.2006.BUFF100.2008.BUFF100.sf)/(st_area(buffer.2006.BUFF100)/100) # 88.80934 [1]
overl_shape.2006.BUFF100.2008.BUFF100 # 97.57801

# spatial change
spatial_ch.2006.BUFF100.2008.BUFF100 <- (overl_hotspot.2006.BUFF100.2008.BUFF100 + overl_shape.2006.BUFF100.2008.BUFF100)/2
spatial_ch.2006.BUFF100.2008.BUFF100 # 77.9211 (no shift)



# End of script -----------------------------------------------------------
