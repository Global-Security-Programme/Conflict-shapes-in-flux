### Preamble ###############################################################
# Islamist insurgency in Nigeria: Spatial change analysis
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


# 2011: Dominant actors -----------------------------------------------------------

data.2011 <- data %>% filter(year == 2011)

data.2011 <- data.2011 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2011 <- data.2011 %>% dplyr::select(side_a_new_id, side_b_new_id, weight)
edges.2011 <- edges.2011 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2011 <- edges.2011 %>% distinct(from, to, weight)

side.a.2011 <- data.2011 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2011 <- data.2011 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2011 <- rbind(side.a.2011, side.b.2011)
nodes.2011 <- nodes.2011 %>% distinct(id, label)

net.2011 <- graph_from_data_frame(edges.2011, 
                                  vertices = nodes.2011,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2011$nodes_n <- igraph::degree(net.2011)
nodes.2011$degree <- strength(net.2011)
nodes.2011$eigenvec <- round(igraph::evcent(net.2011)$vector,2)
nodes.2011$katz <- round(katzcent(net.2011),2)
write_csv(nodes.2011, "out/nigeria-2011-centrality-measures.csv")

net.2011 <- graph_from_data_frame(edges.2011, 
                                  vertices = nodes.2011,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2011,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 0.1, 
                         area = 5),
       label = paste(V(net.2011)$label,":", V(net.2011)$degree),
       label.size = 4,
       color = "blue",
       alpha = 0.3,
       size = "degree", 
       legend.position = "none",
       edge.size = 0.5,
       edge.color = "grey",
       edge.label = "weight",
       edge.label.size = 4,
       edge.label.color = "blue",
       edge.label.alpha = 0.6) +
  ggtitle("Islamist insurgency in Nigeria, 2011")
ggsave("figs/nigeria-2011-network.png", width = 8, height = 4)


# 2011: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011 <- concaveman(points.sf.2011, concavity = 2, length_threshold = 0)
buffer.2011 <- st_buffer(concave.hull.2011, dist = 0.5)
st_write(buffer.2011, "out/nigeria-2011-conflict-shape.shp")

grid.clip.2011 <- st_make_grid(buffer.2011, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(buffer.2011)

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
  scale_fill_manual(values = c("red", "grey90")) +
  ggtitle("Islamist insurgency, 2011: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/nigeria-2011-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2011.sel <- grid.count.2011[grid.count.2011$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2011$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011 <- st_union(grid.count.2011.sel)
st_write(hotspot.2011, "out/nigeria-2011-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2011,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2011, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Islamist insurgency, 2011: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/nigeria-2011-events.png", width = 8, height = 4)


# 2016: Dominant actors ------------------------------------------------------------

data.2016 <- data %>% filter(year == 2016)

data.2016 <- data.2016 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2016 <- data.2016 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2016 <- edges.2016 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2016 <- edges.2016 %>% distinct(from, to, weight)

side.a.2016 <- data.2016 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2016 <- data.2016 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)

nodes.2016 <- rbind(side.a.2016, side.b.2016)
nodes.2016 <- nodes.2016 %>% distinct(id, label)

net.2016 <- graph_from_data_frame(edges.2016, 
                                  vertices = nodes.2016,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2016$nodes_n <- igraph::degree(net.2016)
nodes.2016$degree <- strength(net.2016)
nodes.2016$eigenvec <- round(igraph::evcent(net.2016)$vector,2)
nodes.2016$katz <- round(katzcent(net.2016),2)
write_csv(nodes.2016, "out/nigeria-2016-centrality-measures.csv")

net.2016 <- graph_from_data_frame(edges.2016, 
                                  vertices = nodes.2016,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2016,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 0.1, 
                         area = 5),
       label = paste(V(net.2016)$label,":", V(net.2016)$degree),
       label.size = 4,
       legend.position = "none",
       color = "red",
       alpha = 0.3,
       size = "degree",
       edge.size = 0.5,
       edge.color = "grey",
       edge.label = "weight",
       edge.label.size = 4,
       edge.label.color = "red",
       edge.label.alpha = 0.6) +
  ggtitle("Islamist insurgency in Nigeria, 2016")
ggsave("figs/nigeria-2016-network.png", width = 8, height = 4)


# 2016: Conflict shape and hotspots --------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016 <- concaveman(points.sf.2016, concavity = 2, length_threshold = 0)
buffer.2016 <- st_buffer(concave.hull.2016, dist = 0.5)
st_write(buffer.2016, "out/nigeria-2016-conflict-shape.shp")

grid.clip.2016 <- st_make_grid(buffer.2016, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(buffer.2016)

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
ggsave("figs/nigeria-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/nigeria-2016-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2016,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2016, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Islamist insurgency, 2016: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/nigeria-2016-events.png", width = 8, height = 4)


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
  geom_sf(data = buffer.2011, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2016, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = buffer.2016, alpha = 0.3, fill = "red", color = NA) +
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
ggsave("figs/nigeria-2011-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2011.2016 <- st_area(buffer.2011)/100
new.perc.2011.2016 <- st_area(buffer.2016)/one.perc.2011.2016
new.perc.2011.2016 <- drop_units(new.perc.2011.2016)
perc.change.2011.2016 <- new.perc.2011.2016 - 100
perc.change.2011.2016 # -48.79866 (contraction)

# overlap: hotspots
hotspot.2011.sp <- as(st_geometry(hotspot.2011), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2011.2016 <- gIntersection(hotspot.2011.sp, hotspot.2016.sp)
intersection.h.2011.2016.sf <- st_as_sf(intersection.h.2011.2016)
overl_hotspot.2011.2016 <- st_area(intersection.h.2011.2016.sf)/(st_area(hotspot.2016)/100) 
overl_hotspot.2011.2016 # 42.72037 

# overlap: conflict shape
buffer.2011.sp <- as(st_geometry(buffer.2011), "Spatial")
buffer.2016.sp <- as(st_geometry(buffer.2016), "Spatial")
intersection.2011.2016 <- gIntersection(buffer.2011.sp, buffer.2016.sp)
intersection.2011.2016.sf <- st_as_sf(intersection.2011.2016)
overl_shape.2011.2016 <- st_area(intersection.2011.2016.sf)/(st_area(buffer.2016)/100) 
overl_shape.2011.2016 # 48.33015

# spatial change
spatial_ch.2011.2016 <- (overl_hotspot.2011.2016 + overl_shape.2011.2016)/2
spatial_ch.2011.2016 # 45.52526 (shift)


# 2011.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011.GRID75 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011.GRID75 <- concaveman(points.sf.2011.GRID75, concavity = 2, length_threshold = 0)
buffer.2011.GRID75 <- st_buffer(concave.hull.2011.GRID75, dist = 0.5)
st_write(buffer.2011.GRID75, "out/nigeria-conflict-2011_GRID75-shape.shp")

grid.clip.2011.GRID75 <- st_make_grid(buffer.2011.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2011.GRID75)

points.count.2011.GRID75 <- st_intersects(grid.clip.2011.GRID75, points.sf.2011.GRID75)
grid.count.2011.GRID75 <- st_sf(pt_count = lengths(points.count.2011.GRID75), 
                                geometry = st_cast(grid.clip.2011.GRID75, "MULTIPOLYGON"))

neighbors.2011.GRID75 <- poly2nb(grid.count.2011.GRID75, queen = T)
weighted.neighbors.2011.GRID75 <- nb2listw(include.self(neighbors.2011.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2011.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2011.GRID75$pt_count, 
                                                        weighted.neighbors.2011.GRID75))

grid.count.2011.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2011.GRID75$z_cat[grid.count.2011.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2011.GRID75$z_cat[grid.count.2011.GRID75$HOTSPOT < 3.091 & grid.count.2011.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2011.GRID75$z_cat)

# hotspots
grid.count.2011.GRID75.sel <- grid.count.2011.GRID75[grid.count.2011.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2011.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011.GRID75 <- st_union(grid.count.2011.GRID75.sel)
st_write(hotspot.2011.GRID75, "out/nigeria-2011_GRID75-hotspot.shp")


# 2016.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.GRID75 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.GRID75 <- concaveman(points.sf.2016.GRID75, concavity = 2, length_threshold = 0)
buffer.2016.GRID75 <- st_buffer(concave.hull.2016.GRID75, dist = 0.5)
st_write(buffer.2016.GRID75, "out/nigeria-conflict-2016_GRID75-shape.shp")

grid.clip.2016.GRID75 <- st_make_grid(buffer.2016.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2016.GRID75)

points.count.2016.GRID75 <- st_intersects(grid.clip.2016.GRID75, points.sf.2016.GRID75)
grid.count.2016.GRID75 <- st_sf(pt_count = lengths(points.count.2016.GRID75), 
                                geometry = st_cast(grid.clip.2016.GRID75, "MULTIPOLYGON"))

neighbors.2016.GRID75 <- poly2nb(grid.count.2016.GRID75, queen = T)
weighted.neighbors.2016.GRID75 <- nb2listw(include.self(neighbors.2016.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2016.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2016.GRID75$pt_count, 
                                                        weighted.neighbors.2016.GRID75))

grid.count.2016.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2016.GRID75$z_cat[grid.count.2016.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2016.GRID75$z_cat[grid.count.2016.GRID75$HOTSPOT < 3.091 & grid.count.2016.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2016.GRID75$z_cat)

# hotspots
grid.count.2016.GRID75.sel <- grid.count.2016.GRID75[grid.count.2016.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2016.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016.GRID75 <- st_union(grid.count.2016.GRID75.sel)
st_write(hotspot.2016.GRID75, "out/nigeria-2016_GRID75-hotspot.shp")


# 2011.GRID75-2016.GRID75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2011.GRID75.2016.GRID75 <- st_area(buffer.2011.GRID75)/100
new.perc.2011.GRID75.2016.GRID75 <- st_area(buffer.2016.GRID75)/one.perc.2011.GRID75.2016.GRID75
new.perc.2011.GRID75.2016.GRID75 <- drop_units(new.perc.2011.GRID75.2016.GRID75)
perc.change.2011.GRID75.2016.GRID75 <- new.perc.2011.GRID75.2016.GRID75 - 100
perc.change.2011.GRID75.2016.GRID75 #  -48.79866 (contraction)

# overlap: hotspots
hotspot.2011.GRID75.sp <- as(st_geometry(hotspot.2011.GRID75), "Spatial")
hotspot.2016.GRID75.sp <- as(st_geometry(hotspot.2016.GRID75), "Spatial")
intersection.h.2011.GRID75.2016.GRID75 <- gIntersection(hotspot.2011.GRID75.sp, hotspot.2016.GRID75.sp)
intersection.h.2011.GRID75.2016.GRID75.sf <- st_as_sf(intersection.h.2011.GRID75.2016.GRID75)
overl_hotspot.2011.GRID75.2016.GRID75 <- st_area(intersection.h.2011.GRID75.2016.GRID75.sf)/(st_area(hotspot.2016.GRID75)/100) 
overl_hotspot.2011.GRID75.2016.GRID75 # 51.8529

# overlap: conflict shape
buffer.2011.GRID75.sp <- as(st_geometry(buffer.2011.GRID75), "Spatial")
buffer.2016.GRID75.sp <- as(st_geometry(buffer.2016.GRID75), "Spatial")
intersection.2011.GRID75.2016.GRID75 <- gIntersection(buffer.2011.GRID75.sp, buffer.2016.GRID75.sp)
intersection.2011.GRID75.2016.GRID75.sf <- st_as_sf(intersection.2011.GRID75.2016.GRID75)
overl_shape.2011.GRID75.2016.GRID75 <- st_area(intersection.2011.GRID75.2016.GRID75.sf)/(st_area(buffer.2016.GRID75)/100) # 88.80934 [1]
overl_shape.2011.GRID75.2016.GRID75 # 48.33015

# spatial change
spatial_ch.2011.GRID75.2016.GRID75 <- (overl_hotspot.2011.GRID75.2016.GRID75 + overl_shape.2011.GRID75.2016.GRID75)/2
spatial_ch.2011.GRID75.2016.GRID75 # 50.09153 (no shift)


# 2011.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011.GRID100 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011.GRID100 <- concaveman(points.sf.2011.GRID100, concavity = 2, length_threshold = 0)
buffer.2011.GRID100 <- st_buffer(concave.hull.2011.GRID100, dist = 0.5)
st_write(buffer.2011.GRID100, "out/nigeria-conflict-2011_GRID100-shape.shp")

grid.clip.2011.GRID100 <- st_make_grid(buffer.2011.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2011.GRID100)

points.count.2011.GRID100 <- st_intersects(grid.clip.2011.GRID100, points.sf.2011.GRID100)
grid.count.2011.GRID100 <- st_sf(pt_count = lengths(points.count.2011.GRID100), 
                                 geometry = st_cast(grid.clip.2011.GRID100, "MULTIPOLYGON"))

neighbors.2011.GRID100 <- poly2nb(grid.count.2011.GRID100, queen = T)
weighted.neighbors.2011.GRID100 <- nb2listw(include.self(neighbors.2011.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2011.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2011.GRID100$pt_count, 
                                                         weighted.neighbors.2011.GRID100))

grid.count.2011.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2011.GRID100$z_cat[grid.count.2011.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2011.GRID100$z_cat[grid.count.2011.GRID100$HOTSPOT < 3.091 & grid.count.2011.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2011.GRID100$z_cat)

# hotspots
grid.count.2011.GRID100.sel <- grid.count.2011.GRID100[grid.count.2011.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2011.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011.GRID100 <- st_union(grid.count.2011.GRID100.sel)
st_write(hotspot.2011.GRID100, "out/nigeria-2011_GRID100-hotspot.shp")


# 2016.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.GRID100 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.GRID100 <- concaveman(points.sf.2016.GRID100, concavity = 2, length_threshold = 0)
buffer.2016.GRID100 <- st_buffer(concave.hull.2016.GRID100, dist = 0.5)
st_write(buffer.2016.GRID100, "out/nigeria-conflict-2016_GRID100-shape.shp")

grid.clip.2016.GRID100 <- st_make_grid(buffer.2016.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2016.GRID100)

points.count.2016.GRID100 <- st_intersects(grid.clip.2016.GRID100, points.sf.2016.GRID100)
grid.count.2016.GRID100 <- st_sf(pt_count = lengths(points.count.2016.GRID100), 
                                 geometry = st_cast(grid.clip.2016.GRID100, "MULTIPOLYGON"))

neighbors.2016.GRID100 <- poly2nb(grid.count.2016.GRID100, queen = T)
weighted.neighbors.2016.GRID100 <- nb2listw(include.self(neighbors.2016.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2016.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2016.GRID100$pt_count, 
                                                         weighted.neighbors.2016.GRID100))

grid.count.2016.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2016.GRID100$z_cat[grid.count.2016.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2016.GRID100$z_cat[grid.count.2016.GRID100$HOTSPOT < 3.091 & grid.count.2016.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2016.GRID100$z_cat)

# hotspots
grid.count.2016.GRID100.sel <- grid.count.2016.GRID100[grid.count.2016.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2016.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016.GRID100 <- st_union(grid.count.2016.GRID100.sel)
st_write(hotspot.2016.GRID100, "out/nigeria-2016_GRID100-hotspot.shp")


# 2011.GRID100-2016.GRID100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2011.GRID100.2016.GRID100 <- st_area(buffer.2011.GRID100)/100
new.perc.2011.GRID100.2016.GRID100 <- st_area(buffer.2016.GRID100)/one.perc.2011.GRID100.2016.GRID100
new.perc.2011.GRID100.2016.GRID100 <- drop_units(new.perc.2011.GRID100.2016.GRID100)
perc.change.2011.GRID100.2016.GRID100 <- new.perc.2011.GRID100.2016.GRID100 - 100
perc.change.2011.GRID100.2016.GRID100 #  -48.79866 (contraction)

# overlap: hotspots
hotspot.2011.GRID100.sp <- as(st_geometry(hotspot.2011.GRID100), "Spatial")
hotspot.2016.GRID100.sp <- as(st_geometry(hotspot.2016.GRID100), "Spatial")
intersection.h.2011.GRID100.2016.GRID100 <- gIntersection(hotspot.2011.GRID100.sp, hotspot.2016.GRID100.sp)
intersection.h.2011.GRID100.2016.GRID100.sf <- st_as_sf(intersection.h.2011.GRID100.2016.GRID100)
overl_hotspot.2011.GRID100.2016.GRID100 <- st_area(intersection.h.2011.GRID100.2016.GRID100.sf)/(st_area(hotspot.2016.GRID100)/100) 
overl_hotspot.2011.GRID100.2016.GRID100 # 23.58659

# overlap: conflict shape
buffer.2011.GRID100.sp <- as(st_geometry(buffer.2011.GRID100), "Spatial")
buffer.2016.GRID100.sp <- as(st_geometry(buffer.2016.GRID100), "Spatial")
intersection.2011.GRID100.2016.GRID100 <- gIntersection(buffer.2011.GRID100.sp, buffer.2016.GRID100.sp)
intersection.2011.GRID100.2016.GRID100.sf <- st_as_sf(intersection.2011.GRID100.2016.GRID100)
overl_shape.2011.GRID100.2016.GRID100 <- st_area(intersection.2011.GRID100.2016.GRID100.sf)/(st_area(buffer.2016.GRID100)/100) # 88.80934 [1]
overl_shape.2011.GRID100.2016.GRID100 # 48.33015

# spatial change
spatial_ch.2011.GRID100.2016.GRID100 <- (overl_hotspot.2011.GRID100.2016.GRID100 + overl_shape.2011.GRID100.2016.GRID100)/2
spatial_ch.2011.GRID100.2016.GRID100 # 35.95837 (shift)


# 2011.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011.BUFF75 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011.BUFF75 <- concaveman(points.sf.2011.BUFF75, concavity = 2, length_threshold = 0)
buffer.2011.BUFF75 <- st_buffer(concave.hull.2011.BUFF75, dist = 0.75)
st_write(buffer.2011.BUFF75, "out/nigeria-conflict-2011_BUFF75-shape.shp")

grid.clip.2011.BUFF75 <- st_make_grid(buffer.2011.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2011.BUFF75)

points.count.2011.BUFF75 <- st_intersects(grid.clip.2011.BUFF75, points.sf.2011.BUFF75)
grid.count.2011.BUFF75 <- st_sf(pt_count = lengths(points.count.2011.BUFF75), 
                                geometry = st_cast(grid.clip.2011.BUFF75, "MULTIPOLYGON"))

neighbors.2011.BUFF75 <- poly2nb(grid.count.2011.BUFF75, queen = T)
weighted.neighbors.2011.BUFF75 <- nb2listw(include.self(neighbors.2011.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2011.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2011.BUFF75$pt_count, 
                                                        weighted.neighbors.2011.BUFF75))

grid.count.2011.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2011.BUFF75$z_cat[grid.count.2011.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2011.BUFF75$z_cat[grid.count.2011.BUFF75$HOTSPOT < 3.091 & grid.count.2011.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2011.BUFF75$z_cat)

# hotspots
grid.count.2011.BUFF75.sel <- grid.count.2011.BUFF75[grid.count.2011.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2011.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011.BUFF75 <- st_union(grid.count.2011.BUFF75.sel)
st_write(hotspot.2011.BUFF75, "out/nigeria-2011_BUFF75-hotspot.shp")


# 2016.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.BUFF75 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.BUFF75 <- concaveman(points.sf.2016.BUFF75, concavity = 2, length_threshold = 0)
buffer.2016.BUFF75 <- st_buffer(concave.hull.2016.BUFF75, dist = 0.75)
st_write(buffer.2016.BUFF75, "out/nigeria-conflict-2016_BUFF75-shape.shp")

grid.clip.2016.BUFF75 <- st_make_grid(buffer.2016.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2016.BUFF75)

points.count.2016.BUFF75 <- st_intersects(grid.clip.2016.BUFF75, points.sf.2016.BUFF75)
grid.count.2016.BUFF75 <- st_sf(pt_count = lengths(points.count.2016.BUFF75), 
                                geometry = st_cast(grid.clip.2016.BUFF75, "MULTIPOLYGON"))

neighbors.2016.BUFF75 <- poly2nb(grid.count.2016.BUFF75, queen = T)
weighted.neighbors.2016.BUFF75 <- nb2listw(include.self(neighbors.2016.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2016.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2016.BUFF75$pt_count, 
                                                        weighted.neighbors.2016.BUFF75))

grid.count.2016.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2016.BUFF75$z_cat[grid.count.2016.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2016.BUFF75$z_cat[grid.count.2016.BUFF75$HOTSPOT < 3.091 & grid.count.2016.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2016.BUFF75$z_cat)

# hotspots
grid.count.2016.BUFF75.sel <- grid.count.2016.BUFF75[grid.count.2016.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2016.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016.BUFF75 <- st_union(grid.count.2016.BUFF75.sel)
st_write(hotspot.2016.BUFF75, "out/nigeria-2016_BUFF75-hotspot.shp")


# 2011.BUFF75-2016.BUFF75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2011.BUFF75.2016.BUFF75 <- st_area(buffer.2011.BUFF75)/100
new.perc.2011.BUFF75.2016.BUFF75 <- st_area(buffer.2016.BUFF75)/one.perc.2011.BUFF75.2016.BUFF75
new.perc.2011.BUFF75.2016.BUFF75 <- drop_units(new.perc.2011.BUFF75.2016.BUFF75)
perc.change.2011.BUFF75.2016.BUFF75 <- new.perc.2011.BUFF75.2016.BUFF75 - 100
perc.change.2011.BUFF75.2016.BUFF75 #  -46.9609 (contraction)

# overlap: hotspots
hotspot.2011.BUFF75.sp <- as(st_geometry(hotspot.2011.BUFF75), "Spatial")
hotspot.2016.BUFF75.sp <- as(st_geometry(hotspot.2016.BUFF75), "Spatial")
intersection.h.2011.BUFF75.2016.BUFF75 <- gIntersection(hotspot.2011.BUFF75.sp, hotspot.2016.BUFF75.sp)
intersection.h.2011.BUFF75.2016.BUFF75.sf <- st_as_sf(intersection.h.2011.BUFF75.2016.BUFF75)
overl_hotspot.2011.BUFF75.2016.BUFF75 <- st_area(intersection.h.2011.BUFF75.2016.BUFF75.sf)/(st_area(hotspot.2016.BUFF75)/100) 
overl_hotspot.2011.BUFF75.2016.BUFF75 # 32.48339

# overlap: conflict shape
buffer.2011.BUFF75.sp <- as(st_geometry(buffer.2011.BUFF75), "Spatial")
buffer.2016.BUFF75.sp <- as(st_geometry(buffer.2016.BUFF75), "Spatial")
intersection.2011.BUFF75.2016.BUFF75 <- gIntersection(buffer.2011.BUFF75.sp, buffer.2016.BUFF75.sp)
intersection.2011.BUFF75.2016.BUFF75.sf <- st_as_sf(intersection.2011.BUFF75.2016.BUFF75)
overl_shape.2011.BUFF75.2016.BUFF75 <- st_area(intersection.2011.BUFF75.2016.BUFF75.sf)/(st_area(buffer.2016.BUFF75)/100) # 88.80934 [1]
overl_shape.2011.BUFF75.2016.BUFF75 # 53.3102

# spatial change
spatial_ch.2011.BUFF75.2016.BUFF75 <- (overl_hotspot.2011.BUFF75.2016.BUFF75 + overl_shape.2011.BUFF75.2016.BUFF75)/2
spatial_ch.2011.BUFF75.2016.BUFF75 # 42.8968 (shift)


# 2011.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011.BUFF100 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011.BUFF100 <- concaveman(points.sf.2011.BUFF100, concavity = 2, length_threshold = 0)
buffer.2011.BUFF100 <- st_buffer(concave.hull.2011.BUFF100, dist = 1)
st_write(buffer.2011.BUFF100, "out/nigeria-conflict-2011_BUFF100-shape.shp")

grid.clip.2011.BUFF100 <- st_make_grid(buffer.2011.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2011.BUFF100)

points.count.2011.BUFF100 <- st_intersects(grid.clip.2011.BUFF100, points.sf.2011.BUFF100)
grid.count.2011.BUFF100 <- st_sf(pt_count = lengths(points.count.2011.BUFF100), 
                                 geometry = st_cast(grid.clip.2011.BUFF100, "MULTIPOLYGON"))

neighbors.2011.BUFF100 <- poly2nb(grid.count.2011.BUFF100, queen = T)
weighted.neighbors.2011.BUFF100 <- nb2listw(include.self(neighbors.2011.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2011.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2011.BUFF100$pt_count, 
                                                         weighted.neighbors.2011.BUFF100))

grid.count.2011.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2011.BUFF100$z_cat[grid.count.2011.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2011.BUFF100$z_cat[grid.count.2011.BUFF100$HOTSPOT < 3.091 & grid.count.2011.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2011.BUFF100$z_cat)

# hotspots
grid.count.2011.BUFF100.sel <- grid.count.2011.BUFF100[grid.count.2011.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2011.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011.BUFF100 <- st_union(grid.count.2011.BUFF100.sel)
st_write(hotspot.2011.BUFF100, "out/nigeria-2011_BUFF100-hotspot.shp")


# 2016.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.BUFF100 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.BUFF100 <- concaveman(points.sf.2016.BUFF100, concavity = 2, length_threshold = 0)
buffer.2016.BUFF100 <- st_buffer(concave.hull.2016.BUFF100, dist = 1)
st_write(buffer.2016.BUFF100, "out/nigeria-conflict-2016_BUFF100-shape.shp")

grid.clip.2016.BUFF100 <- st_make_grid(buffer.2016.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2016.BUFF100)

points.count.2016.BUFF100 <- st_intersects(grid.clip.2016.BUFF100, points.sf.2016.BUFF100)
grid.count.2016.BUFF100 <- st_sf(pt_count = lengths(points.count.2016.BUFF100), 
                                 geometry = st_cast(grid.clip.2016.BUFF100, "MULTIPOLYGON"))

neighbors.2016.BUFF100 <- poly2nb(grid.count.2016.BUFF100, queen = T)
weighted.neighbors.2016.BUFF100 <- nb2listw(include.self(neighbors.2016.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2016.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2016.BUFF100$pt_count, 
                                                         weighted.neighbors.2016.BUFF100))

grid.count.2016.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2016.BUFF100$z_cat[grid.count.2016.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2016.BUFF100$z_cat[grid.count.2016.BUFF100$HOTSPOT < 3.091 & grid.count.2016.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2016.BUFF100$z_cat)

# hotspots
grid.count.2016.BUFF100.sel <- grid.count.2016.BUFF100[grid.count.2016.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2016.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016.BUFF100 <- st_union(grid.count.2016.BUFF100.sel)
st_write(hotspot.2016.BUFF100, "out/nigeria-2016_BUFF100-hotspot.shp")


# 2011.BUFF100-2016.BUFF100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2011.BUFF100.2016.BUFF100 <- st_area(buffer.2011.BUFF100)/100
new.perc.2011.BUFF100.2016.BUFF100 <- st_area(buffer.2016.BUFF100)/one.perc.2011.BUFF100.2016.BUFF100
new.perc.2011.BUFF100.2016.BUFF100 <- drop_units(new.perc.2011.BUFF100.2016.BUFF100)
perc.change.2011.BUFF100.2016.BUFF100 <- new.perc.2011.BUFF100.2016.BUFF100 - 100
perc.change.2011.BUFF100.2016.BUFF100 #  -45.58158 (contraction)

# overlap: hotspots
hotspot.2011.BUFF100.sp <- as(st_geometry(hotspot.2011.BUFF100), "Spatial")
hotspot.2016.BUFF100.sp <- as(st_geometry(hotspot.2016.BUFF100), "Spatial")
intersection.h.2011.BUFF100.2016.BUFF100 <- gIntersection(hotspot.2011.BUFF100.sp, hotspot.2016.BUFF100.sp)
intersection.h.2011.BUFF100.2016.BUFF100.sf <- st_as_sf(intersection.h.2011.BUFF100.2016.BUFF100)
overl_hotspot.2011.BUFF100.2016.BUFF100 <- st_area(intersection.h.2011.BUFF100.2016.BUFF100.sf)/(st_area(hotspot.2016.BUFF100)/100) 
overl_hotspot.2011.BUFF100.2016.BUFF100 # 36.16869

# overlap: conflict shape
buffer.2011.BUFF100.sp <- as(st_geometry(buffer.2011.BUFF100), "Spatial")
buffer.2016.BUFF100.sp <- as(st_geometry(buffer.2016.BUFF100), "Spatial")
intersection.2011.BUFF100.2016.BUFF100 <- gIntersection(buffer.2011.BUFF100.sp, buffer.2016.BUFF100.sp)
intersection.2011.BUFF100.2016.BUFF100.sf <- st_as_sf(intersection.2011.BUFF100.2016.BUFF100)
overl_shape.2011.BUFF100.2016.BUFF100 <- st_area(intersection.2011.BUFF100.2016.BUFF100.sf)/(st_area(buffer.2016.BUFF100)/100) # 88.80934 [1]
overl_shape.2011.BUFF100.2016.BUFF100 # 57.46443

# spatial change
spatial_ch.2011.BUFF100.2016.BUFF100 <- (overl_hotspot.2011.BUFF100.2016.BUFF100 + overl_shape.2011.BUFF100.2016.BUFF100)/2
spatial_ch.2011.BUFF100.2016.BUFF100 # 46.81656 (shift)


# End of script -----------------------------------------------------------


