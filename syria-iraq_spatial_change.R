### Preamble ###############################################################
# Armed conflict in Syria and Iraq: Spatial change analysis
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


# 2010: Dominant actors -----------------------------------------------------------

data.2010 <- data %>% filter(year == 2010)

data.2010 <- data.2010 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()
edges.2010 <- data.2010 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2010 <- edges.2010 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2010 <- edges.2010 %>% distinct(from, to, weight)

side.a.2010 <- data.2010 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2010 <- data.2010 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2010 <- rbind(side.a.2010, side.b.2010)
nodes.2010 <- nodes.2010 %>% distinct(id, label)

net.2010 <- graph_from_data_frame(edges.2010, 
                                  vertices = nodes.2010,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2010$nodes_n <- igraph::degree(net.2010)
nodes.2010$degree <- strength(net.2010)
nodes.2010$eigenvec <- round(igraph::evcent(net.2010)$vector,2)
nodes.2010$katz <- round(katzcent(net.2010),2)
write_csv(nodes.2010, "out/syria-iraq-2010-centrality-measures.csv")

net.2010 <- graph_from_data_frame(edges.2010, 
                                  vertices = nodes.2010,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2010,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 10, 
                         area = 100),
       layout.exp = 0.5,
       label = paste(V(net.2010)$label,":", V(net.2010)$degree),
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
  ggtitle("Armed conflict in Syria and Iraq, 2010")
ggsave("figs/syria-iraq-2010-network.png", width = 10, height = 2)


# 2010: Conflict shape and hotspots ----------------------------------------------------

points.sf.2010 <- points.sf %>% filter(year == 2010)

# conflict shape
concave.hull.2010 <- concaveman(points.sf.2010, concavity = 2, length_threshold = 0)
buffer.2010 <- st_buffer(concave.hull.2010, dist = 0.5)
st_write(buffer.2010, "out/syria-iraq-2010-conflict-shape.shp")

grid.clip.2010 <- st_make_grid(buffer.2010, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(buffer.2010)

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
  ggtitle("Armed conflict in Syria and Iraq, 2010: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2010-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2010.sel <- grid.count.2010[grid.count.2010$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2010$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010 <- st_union(grid.count.2010.sel)
st_write(hotspot.2010, "out/syria-iraq-2010-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2010,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2010, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq, 2010: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2010-events.png", width = 8, height = 4)


# 2013: Dominant actors -----------------------------------------------------------

data.2013 <- data %>% filter(year == 2013)

data.2013 <- data.2013 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2013 <- data.2013 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2013 <- edges.2013 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2013 <- edges.2013 %>% distinct(from, to, weight)

side.a.2013 <- data.2013 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2013 <- data.2013 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2013 <- rbind(side.a.2013, side.b.2013)
nodes.2013 <- nodes.2013 %>% distinct(id, label)

net.2013 <- graph_from_data_frame(edges.2013, 
                                  vertices = nodes.2013,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2013$nodes_n <- igraph::degree(net.2013)
nodes.2013$degree <- strength(net.2013)
nodes.2013$eigenvec <- round(igraph::evcent(net.2013)$vector,2)
nodes.2013$katz <- round(katzcent(net.2013),2)
write_csv(nodes.2013, "out/syria-iraq-2013-centrality-measures.csv")

net.2013 <- graph_from_data_frame(edges.2013, 
                                  vertices = nodes.2013,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2013,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 10, 
                         area = 100),
       layout.exp = 0.5,
       label = paste(str_trunc(V(net.2013)$label, 40, "center"),":", V(net.2013)$degree),
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
  ggtitle("Armed conflict in Syria and Iraq, 2013")
ggsave("figs/syria-iraq-2013-network.png", width = 10, height = 6)


# 2013: Conflict shapes and hotspots --------------------------------------------------------

points.sf.2013 <- points.sf %>% filter(year == 2013)

# conflict shape
concave.hull.2013 <- concaveman(points.sf.2013, concavity = 2, length_threshold = 0)
buffer.2013 <- st_buffer(concave.hull.2013, dist = 0.5)
st_write(buffer.2013, "out/syria-iraq-2013-conflict-shape.shp")

grid.clip.2013 <- st_make_grid(buffer.2013, 
                               cellsize = c(0.5,0.5), #cell size in degrees
                               what="polygons") %>%
  st_intersection(buffer.2013)

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

# confict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2013, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed conflict in Syria and Iraq, 2013: Conflict shape and Hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2013-conflict-shape-hotspot.png", width = 8, height = 4)

# hostpots
grid.count.2013.sel <- grid.count.2013[grid.count.2013$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2013$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013 <- st_union(grid.count.2013.sel)
st_write(hotspot.2013, "out/syria-iraq-2013-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2013,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2013, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("AArmed conflict in Syria and Iraq, 2013: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2013-events.png", width = 8, height = 4)


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
  geom_sf(data = buffer.2010, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2013, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = buffer.2013, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq, 2010-2013: Spatial change") +
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
ggsave("figs/syria-iraq-2010-2013-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2010.2013 <- st_area(buffer.2010)/100
new.perc.2010.2013 <- st_area(buffer.2013)/one.perc.2010.2013
new.perc.2010.2013 <- drop_units(new.perc.2010.2013)
perc.change.2010.2013 <- new.perc.2010.2013 - 100
perc.change.2010.2013 #  41.00099 (expansion)

# overlap: hotspots
hotspot.2010.sp <- as(st_geometry(hotspot.2010), "Spatial")
hotspot.2013.sp <- as(st_geometry(hotspot.2013), "Spatial")
intersection.h.2010.2013 <- gIntersection(hotspot.2010.sp, hotspot.2013.sp)
intersection.h.2010.2013.sf <- st_as_sf(intersection.h.2010.2013) # no overlap
overl_hotspot.2011.2013 <- 0
overl_hotspot.2011.2013

# overlap shape
buffer.2010.sp <- as(st_geometry(buffer.2010), "Spatial")
buffer.2013.sp <- as(st_geometry(buffer.2013), "Spatial")
intersection.2010.2013 <- gIntersection(buffer.2010.sp, buffer.2013.sp)
intersection.2010.2013.sf <- st_as_sf(intersection.2010.2013)
overl_shape.2011.2013 <- st_area(intersection.2010.2013.sf)/(st_area(buffer.2010)/100)
overl_shape.2011.2013 # 58.86652

# spatial change
spatial_ch.2011.2013 <- (overl_hotspot.2011.2013 + drop_units(overl_shape.2011.2013))/2
spatial_ch.2011.2013 # 29.43326 (shift)


# 2014: Dominant actors -----------------------------------------------------------

data.2014 <- data %>% filter(year == 2014)

data.2014 <- data.2014 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2014 <- data.2014 %>% dplyr::select(side_a_new_id,side_b_new_id, weight)
edges.2014 <- edges.2014 %>% rename(from = side_a_new_id,
                                    to = side_b_new_id)
edges.2014 <- edges.2014 %>% distinct(from, to, weight)

side.a.2014 <- data.2014 %>% distinct(side_a_new_id, side_a) %>% 
  rename(id = side_a_new_id, label = side_a)
side.b.2014 <- data.2014 %>% distinct(side_b_new_id, side_b) %>% 
  rename(id = side_b_new_id, label = side_b)
nodes.2014 <- rbind(side.a.2014, side.b.2014)
nodes.2014 <- nodes.2014 %>% distinct(id, label)

net.2014 <- graph_from_data_frame(edges.2014, 
                                  vertices = nodes.2014,
                                  directed = F)

# degree centrality and other network centrality measures
nodes.2014$nodes_n <- igraph::degree(net.2014)
nodes.2014$degree <- strength(net.2014)
nodes.2014$eigenvec <- round(igraph::evcent(net.2014)$vector,2)
nodes.2014$katz <- round(katzcent(net.2014),2)
write_csv(nodes.2014, "out/syria-iraq-2014-centrality-measures.csv")

net.2014 <- graph_from_data_frame(edges.2014, 
                                  vertices = nodes.2014,
                                  directed = F)


# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(18)

ggnet2(net.2014,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 5, 
                         area = 1000),
       label = paste(V(net.2014)$label,":", V(net.2014)$degree),
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
  ggtitle("Armed conflict in Syria and Iraq, 2014")
ggsave("figs/syria-iraq-2014-network.png", width = 8, height = 4)



# 2014: Conflict shape and hotspots --------------------------------------------------------

points.sf.2014 <- points.sf %>% filter(year == 2014)

# conflict shape
concave.hull.2014 <- concaveman(points.sf.2014, concavity = 2, length_threshold = 0)
buffer.2014 <- st_buffer(concave.hull.2014, dist = 0.5)
st_write(buffer.2014, "out/syria-iraq-2014-conflict-shape.shp")

grid.clip.2014 <- st_make_grid(buffer.2014, 
                               cellsize = c(0.5,0.5), 
                               what="polygons") %>%
  st_intersection(buffer.2014)

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

# conflict shape and hotspots
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = grid.count.2014, 
          aes(fill = z_cat), 
          color = "grey90",
          alpha = 0.6) + 
  scale_fill_manual(values = c("red", "orange", "grey90")) +
  ggtitle("Armed Conflict in Syria and Iraq, 2014: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2014-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2014.sel <- grid.count.2014[grid.count.2014$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2014$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014 <- st_union(grid.count.2014.sel)
st_write(hotspot.2014, "out/syria-iraq-2014-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2014,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2014, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed Conflict in Syria and Iraq, 2014: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2014-events.png", width = 8, height = 4)


# 2016: Dominant actors -----------------------------------------------------------

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
write_csv(nodes.2016, "out/syria-iraq-2016-centrality-measures.csv")

net.2016 <- graph_from_data_frame(edges.2016, 
                                  vertices = nodes.2016,
                                  directed = F)

# network graph
theme_update(plot.title = element_text(hjust = 0.5))
set.seed(20)

ggnet2(net.2016,
       mode = "fruchtermanreingold", 
       layout.par = list(repulse.rad = 5, 
                         area = 1000),
       label = paste(V(net.2016)$label,":", V(net.2016)$degree),
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
  ggtitle("2016")
ggsave("figs/syria-iraq-2016-network.png", width = 8, height = 4)


# 2016: Conflict shape and hotspots --------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016 <- concaveman(points.sf.2016, concavity = 2, length_threshold = 0)
buffer.2016 <- st_buffer(concave.hull.2016, dist = 0.5)
st_write(buffer.2016, "out/syria-iraq-2016-conflict-shape.shp")

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
  ggtitle("Armed Conflict in Syria and Iraq, 2016: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/syria-iraq-2016-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2016,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2016, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed Conflict in Syria and Iraq, 2016: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2016-events.png", width = 8, height = 4)


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
  geom_sf(data = buffer.2014, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2016, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = buffer.2016, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed Conflict in Syria and Iraq, 2014-2016: Spatial change") +
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
ggsave("figs/syria-iraq-2014-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2014.2016 <- st_area(buffer.2014)/100
new.perc.2014.2016 <- st_area(buffer.2016)/one.perc.2014.2016
new.perc.2014.2016 <- drop_units(new.perc.2014.2016)
perc.change.2014.2016 <- new.perc.2014.2016 - 100
perc.change.2014.2016 # 34.21461 (expansion)

# overlap: hotspots
hotspot.2014.sp <- as(st_geometry(hotspot.2014), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2014.2016 <- gIntersection(hotspot.2014.sp, hotspot.2016.sp)
intersection.h.2014.2016.sf <- st_as_sf(intersection.h.2014.2016)
overl_hotspot.2014.2016 <- st_area(intersection.h.2014.2016.sf)/(st_area(hotspot.2014)/100) 
overl_hotspot.2014.2016  # 74.20413 

# overlap: conflict shape
buffer.2014.sp <- as(st_geometry(buffer.2014), "Spatial")
buffer.2016.sp <- as(st_geometry(buffer.2016), "Spatial")
intersection.2014.2016 <- gIntersection(buffer.2014.sp, buffer.2016.sp)
intersection.2014.2016.sf <- st_as_sf(intersection.2014.2016)
overl_shape.2014.2016 <- st_area(intersection.2014.2016.sf)/(st_area(buffer.2014)/100) 
overl_shape.2014.2016  # 93.0838 

# spatial change
spatial_ch.2014.2016 <- (overl_hotspot.2014.2016 + overl_shape.2014.2016)/2
spatial_ch.2014.2016  # 83.64396 (no shift)


# 2010.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2010.GRID75 <- points.sf %>% filter(year == 2010)

# conflict shape
concave.hull.2010.GRID75 <- concaveman(points.sf.2010.GRID75, concavity = 2, length_threshold = 0)
buffer.2010.GRID75 <- st_buffer(concave.hull.2010.GRID75, dist = 0.5)
st_write(buffer.2010.GRID75, "out/syria-iraq-conflict-2010_GRID75-shape.shp")

grid.clip.2010.GRID75 <- st_make_grid(buffer.2010.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2010.GRID75)

points.count.2010.GRID75 <- st_intersects(grid.clip.2010.GRID75, points.sf.2010.GRID75)
grid.count.2010.GRID75 <- st_sf(pt_count = lengths(points.count.2010.GRID75), 
                                geometry = st_cast(grid.clip.2010.GRID75, "MULTIPOLYGON"))

neighbors.2010.GRID75 <- poly2nb(grid.count.2010.GRID75, queen = T)
weighted.neighbors.2010.GRID75 <- nb2listw(include.self(neighbors.2010.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2010.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2010.GRID75$pt_count, 
                                                        weighted.neighbors.2010.GRID75))

grid.count.2010.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2010.GRID75$z_cat[grid.count.2010.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010.GRID75$z_cat[grid.count.2010.GRID75$HOTSPOT < 3.091 & grid.count.2010.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010.GRID75$z_cat)

# hotspots
grid.count.2010.GRID75.sel <- grid.count.2010.GRID75[grid.count.2010.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2010.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010.GRID75 <- st_union(grid.count.2010.GRID75.sel)
st_write(hotspot.2010.GRID75, "out/syria-iraq-2010_GRID75-hotspot.shp")


# 2013.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2013.GRID75 <- points.sf %>% filter(year == 2013)

# conflict shape
concave.hull.2013.GRID75 <- concaveman(points.sf.2013.GRID75, concavity = 2, length_threshold = 0)
buffer.2013.GRID75 <- st_buffer(concave.hull.2013.GRID75, dist = 0.5)
st_write(buffer.2013.GRID75, "out/syria-iraq-conflict-2013_GRID75-shape.shp")

grid.clip.2013.GRID75 <- st_make_grid(buffer.2013.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2013.GRID75)

points.count.2013.GRID75 <- st_intersects(grid.clip.2013.GRID75, points.sf.2013.GRID75)
grid.count.2013.GRID75 <- st_sf(pt_count = lengths(points.count.2013.GRID75), 
                                geometry = st_cast(grid.clip.2013.GRID75, "MULTIPOLYGON"))

neighbors.2013.GRID75 <- poly2nb(grid.count.2013.GRID75, queen = T)
weighted.neighbors.2013.GRID75 <- nb2listw(include.self(neighbors.2013.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2013.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2013.GRID75$pt_count, 
                                                        weighted.neighbors.2013.GRID75))

grid.count.2013.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2013.GRID75$z_cat[grid.count.2013.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013.GRID75$z_cat[grid.count.2013.GRID75$HOTSPOT < 3.091 & grid.count.2013.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013.GRID75$z_cat)

# hotspots
grid.count.2013.GRID75.sel <- grid.count.2013.GRID75[grid.count.2013.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2013.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013.GRID75 <- st_union(grid.count.2013.GRID75.sel)
st_write(hotspot.2013.GRID75, "out/syria-iraq-2013_GRID75-hotspot.shp")


# 2010.GRID75-2013.GRID75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2010.GRID75.2013.GRID75 <- st_area(buffer.2010.GRID75)/100
new.perc.2010.GRID75.2013.GRID75 <- st_area(buffer.2013.GRID75)/one.perc.2010.GRID75.2013.GRID75
new.perc.2010.GRID75.2013.GRID75 <- drop_units(new.perc.2010.GRID75.2013.GRID75)
perc.change.2010.GRID75.2013.GRID75 <- new.perc.2010.GRID75.2013.GRID75 - 100
perc.change.2010.GRID75.2013.GRID75 #  41.00099 (expansion)

# overlap: hotspots
hotspot.2010.GRID75.sp <- as(st_geometry(hotspot.2010.GRID75), "Spatial")
hotspot.2013.GRID75.sp <- as(st_geometry(hotspot.2013.GRID75), "Spatial")
intersection.h.2010.GRID75.2013.GRID75 <- gIntersection(hotspot.2010.GRID75.sp, hotspot.2013.GRID75.sp)
intersection.h.2010.GRID75.2013.GRID75.sf <- st_as_sf(intersection.h.2010.GRID75.2013.GRID75) # no overlap
overl_hotspot.2010.GRID75.2013.GRID75 <- 0
overl_hotspot.2010.GRID75.2013.GRID75

# overlap shape
buffer.2010.GRID75.sp <- as(st_geometry(buffer.2010.GRID75), "Spatial")
buffer.2013.GRID75.sp <- as(st_geometry(buffer.2013.GRID75), "Spatial")
intersection.2010.GRID75.2013.GRID75 <- gIntersection(buffer.2010.GRID75.sp, buffer.2013.GRID75.sp)
intersection.2010.GRID75.2013.GRID75.sf <- st_as_sf(intersection.2010.GRID75.2013.GRID75)
overl_shape.2010.GRID75.2013.GRID75 <- st_area(intersection.2010.GRID75.2013.GRID75.sf)/(st_area(buffer.2010.GRID75)/100)
overl_shape.2010.GRID75.2013.GRID75 # 58.86652

# spatial change
spatial_ch.2010.GRID75.2013.GRID75 <- (overl_hotspot.2010.GRID75.2013.GRID75 + drop_units(overl_shape.2010.GRID75.2013.GRID75))/2
spatial_ch.2010.GRID75.2013.GRID75 # 29.43326 (shift)


# 2014.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2014.GRID75 <- points.sf %>% filter(year == 2014)

# conflict shape
concave.hull.2014.GRID75 <- concaveman(points.sf.2014.GRID75, concavity = 2, length_threshold = 0)
buffer.2014.GRID75 <- st_buffer(concave.hull.2014.GRID75, dist = 0.5)
st_write(buffer.2014.GRID75, "out/syria-iraq-conflict-2014_GRID75-shape.shp")

grid.clip.2014.GRID75 <- st_make_grid(buffer.2014.GRID75, 
                                      cellsize = c(0.75,0.75), 
                                      what="polygons") %>%
  st_intersection(buffer.2014.GRID75)

points.count.2014.GRID75 <- st_intersects(grid.clip.2014.GRID75, points.sf.2014.GRID75)
grid.count.2014.GRID75 <- st_sf(pt_count = lengths(points.count.2014.GRID75), 
                                geometry = st_cast(grid.clip.2014.GRID75, "MULTIPOLYGON"))

neighbors.2014.GRID75 <- poly2nb(grid.count.2014.GRID75, queen = T)
weighted.neighbors.2014.GRID75 <- nb2listw(include.self(neighbors.2014.GRID75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2014.GRID75$HOTSPOT <- as.vector(localG_perm(grid.count.2014.GRID75$pt_count, 
                                                        weighted.neighbors.2014.GRID75))

grid.count.2014.GRID75$z_cat <- "Not significant or cold spot"
grid.count.2014.GRID75$z_cat[grid.count.2014.GRID75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2014.GRID75$z_cat[grid.count.2014.GRID75$HOTSPOT < 3.091 & grid.count.2014.GRID75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2014.GRID75$z_cat)

# hotspots
grid.count.2014.GRID75.sel <- grid.count.2014.GRID75[grid.count.2014.GRID75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2014.GRID75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014.GRID75 <- st_union(grid.count.2014.GRID75.sel)
st_write(hotspot.2014.GRID75, "out/syria-iraq-2014_GRID75-hotspot.shp")


# 2016.GRID75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.GRID75 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.GRID75 <- concaveman(points.sf.2016.GRID75, concavity = 2, length_threshold = 0)
buffer.2016.GRID75 <- st_buffer(concave.hull.2016.GRID75, dist = 0.5)
st_write(buffer.2016.GRID75, "out/syria-iraq-conflict-2016_GRID75-shape.shp")

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
st_write(hotspot.2016.GRID75, "out/syria-iraq-2016_GRID75-hotspot.shp")


# 2014.GRID75-2016.GRID75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2014.GRID75.2016.GRID75 <- st_area(buffer.2014.GRID75)/100
new.perc.2014.GRID75.2016.GRID75 <- st_area(buffer.2016.GRID75)/one.perc.2014.GRID75.2016.GRID75
new.perc.2014.GRID75.2016.GRID75 <- drop_units(new.perc.2014.GRID75.2016.GRID75)
perc.change.2014.GRID75.2016.GRID75 <- new.perc.2014.GRID75.2016.GRID75 - 100
perc.change.2014.GRID75.2016.GRID75 #  34.21461 (expansion)

# overlap: hotspots
hotspot.2014.GRID75.sp <- as(st_geometry(hotspot.2014.GRID75), "Spatial")
hotspot.2016.GRID75.sp <- as(st_geometry(hotspot.2016.GRID75), "Spatial")
intersection.h.2014.GRID75.2016.GRID75 <- gIntersection(hotspot.2014.GRID75.sp, hotspot.2016.GRID75.sp)
intersection.h.2014.GRID75.2016.GRID75.sf <- st_as_sf(intersection.h.2014.GRID75.2016.GRID75) # no overlap
overl_hotspot.2014.GRID75.2016.GRID75 <- st_area(intersection.h.2014.GRID75.2016.GRID75.sf)/(st_area(hotspot.2014.GRID75)/100) 
overl_hotspot.2014.GRID75.2016.GRID75 # 86.41572

# overlap shape
buffer.2014.GRID75.sp <- as(st_geometry(buffer.2014.GRID75), "Spatial")
buffer.2016.GRID75.sp <- as(st_geometry(buffer.2016.GRID75), "Spatial")
intersection.2014.GRID75.2016.GRID75 <- gIntersection(buffer.2014.GRID75.sp, buffer.2016.GRID75.sp)
intersection.2014.GRID75.2016.GRID75.sf <- st_as_sf(intersection.2014.GRID75.2016.GRID75)
overl_shape.2014.GRID75.2016.GRID75 <- st_area(intersection.2014.GRID75.2016.GRID75.sf)/(st_area(buffer.2014.GRID75)/100)
overl_shape.2014.GRID75.2016.GRID75 # 93.0838

# spatial change
spatial_ch.2014.GRID75.2016.GRID75 <- (overl_hotspot.2014.GRID75.2016.GRID75 + overl_shape.2014.GRID75.2016.GRID75)/2
spatial_ch.2014.GRID75.2016.GRID75 # 89.74976 (no shift)


# 2010.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2010.GRID100 <- points.sf %>% filter(year == 2010)

# conflict shape
concave.hull.2010.GRID100 <- concaveman(points.sf.2010.GRID100, concavity = 2, length_threshold = 0)
buffer.2010.GRID100 <- st_buffer(concave.hull.2010.GRID100, dist = 0.5)
st_write(buffer.2010.GRID100, "out/syria-iraq-conflict-2010_GRID100-shape.shp")

grid.clip.2010.GRID100 <- st_make_grid(buffer.2010.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2010.GRID100)

points.count.2010.GRID100 <- st_intersects(grid.clip.2010.GRID100, points.sf.2010.GRID100)
grid.count.2010.GRID100 <- st_sf(pt_count = lengths(points.count.2010.GRID100), 
                                 geometry = st_cast(grid.clip.2010.GRID100, "MULTIPOLYGON"))

neighbors.2010.GRID100 <- poly2nb(grid.count.2010.GRID100, queen = T)
weighted.neighbors.2010.GRID100 <- nb2listw(include.self(neighbors.2010.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2010.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2010.GRID100$pt_count, 
                                                         weighted.neighbors.2010.GRID100))

grid.count.2010.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2010.GRID100$z_cat[grid.count.2010.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010.GRID100$z_cat[grid.count.2010.GRID100$HOTSPOT < 3.091 & grid.count.2010.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010.GRID100$z_cat)

# hotspots
grid.count.2010.GRID100.sel <- grid.count.2010.GRID100[grid.count.2010.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2010.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010.GRID100 <- st_union(grid.count.2010.GRID100.sel)
st_write(hotspot.2010.GRID100, "out/syria-iraq-2010_GRID100-hotspot.shp")


# 2013.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2013.GRID100 <- points.sf %>% filter(year == 2013)

# conflict shape
concave.hull.2013.GRID100 <- concaveman(points.sf.2013.GRID100, concavity = 2, length_threshold = 0)
buffer.2013.GRID100 <- st_buffer(concave.hull.2013.GRID100, dist = 0.5)
st_write(buffer.2013.GRID100, "out/syria-iraq-conflict-2013_GRID100-shape.shp")

grid.clip.2013.GRID100 <- st_make_grid(buffer.2013.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2013.GRID100)

points.count.2013.GRID100 <- st_intersects(grid.clip.2013.GRID100, points.sf.2013.GRID100)
grid.count.2013.GRID100 <- st_sf(pt_count = lengths(points.count.2013.GRID100), 
                                 geometry = st_cast(grid.clip.2013.GRID100, "MULTIPOLYGON"))

neighbors.2013.GRID100 <- poly2nb(grid.count.2013.GRID100, queen = T)
weighted.neighbors.2013.GRID100 <- nb2listw(include.self(neighbors.2013.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2013.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2013.GRID100$pt_count, 
                                                         weighted.neighbors.2013.GRID100))

grid.count.2013.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2013.GRID100$z_cat[grid.count.2013.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013.GRID100$z_cat[grid.count.2013.GRID100$HOTSPOT < 3.091 & grid.count.2013.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013.GRID100$z_cat)

# hotspots
grid.count.2013.GRID100.sel <- grid.count.2013.GRID100[grid.count.2013.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2013.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013.GRID100 <- st_union(grid.count.2013.GRID100.sel)
st_write(hotspot.2013.GRID100, "out/syria-iraq-2013_GRID100-hotspot.shp")


# 2010.GRID100-2013.GRID100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2010.GRID100.2013.GRID100 <- st_area(buffer.2010.GRID100)/100
new.perc.2010.GRID100.2013.GRID100 <- st_area(buffer.2013.GRID100)/one.perc.2010.GRID100.2013.GRID100
new.perc.2010.GRID100.2013.GRID100 <- drop_units(new.perc.2010.GRID100.2013.GRID100)
perc.change.2010.GRID100.2013.GRID100 <- new.perc.2010.GRID100.2013.GRID100 - 100
perc.change.2010.GRID100.2013.GRID100 #  41.00099 (expansion)

# overlap: hotspots
hotspot.2010.GRID100.sp <- as(st_geometry(hotspot.2010.GRID100), "Spatial")
hotspot.2013.GRID100.sp <- as(st_geometry(hotspot.2013.GRID100), "Spatial")
intersection.h.2010.GRID100.2013.GRID100 <- gIntersection(hotspot.2010.GRID100.sp, hotspot.2013.GRID100.sp)
intersection.h.2010.GRID100.2013.GRID100.sf <- st_as_sf(intersection.h.2010.GRID100.2013.GRID100) # no overlap
overl_hotspot.2010.GRID100.2013.GRID100 <- 0
overl_hotspot.2010.GRID100.2013.GRID100

# overlap shape
buffer.2010.GRID100.sp <- as(st_geometry(buffer.2010.GRID100), "Spatial")
buffer.2013.GRID100.sp <- as(st_geometry(buffer.2013.GRID100), "Spatial")
intersection.2010.GRID100.2013.GRID100 <- gIntersection(buffer.2010.GRID100.sp, buffer.2013.GRID100.sp)
intersection.2010.GRID100.2013.GRID100.sf <- st_as_sf(intersection.2010.GRID100.2013.GRID100)
overl_shape.2010.GRID100.2013.GRID100 <- st_area(intersection.2010.GRID100.2013.GRID100.sf)/(st_area(buffer.2010.GRID100)/100)
overl_shape.2010.GRID100.2013.GRID100 # 58.86652

# spatial change
spatial_ch.2010.GRID100.2013.GRID100 <- (overl_hotspot.2010.GRID100.2013.GRID100 + drop_units(overl_shape.2010.GRID100.2013.GRID100))/2
spatial_ch.2010.GRID100.2013.GRID100 # 29.43326 (shift)


# 2014.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2014.GRID100 <- points.sf %>% filter(year == 2014)

# conflict shape
concave.hull.2014.GRID100 <- concaveman(points.sf.2014.GRID100, concavity = 2, length_threshold = 0)
buffer.2014.GRID100 <- st_buffer(concave.hull.2014.GRID100, dist = 0.5)
st_write(buffer.2014.GRID100, "out/syria-iraq-conflict-2014_GRID100-shape.shp")

grid.clip.2014.GRID100 <- st_make_grid(buffer.2014.GRID100, 
                                       cellsize = c(1,1), 
                                       what="polygons") %>%
  st_intersection(buffer.2014.GRID100)

points.count.2014.GRID100 <- st_intersects(grid.clip.2014.GRID100, points.sf.2014.GRID100)
grid.count.2014.GRID100 <- st_sf(pt_count = lengths(points.count.2014.GRID100), 
                                 geometry = st_cast(grid.clip.2014.GRID100, "MULTIPOLYGON"))

neighbors.2014.GRID100 <- poly2nb(grid.count.2014.GRID100, queen = T)
weighted.neighbors.2014.GRID100 <- nb2listw(include.self(neighbors.2014.GRID100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2014.GRID100$HOTSPOT <- as.vector(localG_perm(grid.count.2014.GRID100$pt_count, 
                                                         weighted.neighbors.2014.GRID100))

grid.count.2014.GRID100$z_cat <- "Not significant or cold spot"
grid.count.2014.GRID100$z_cat[grid.count.2014.GRID100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2014.GRID100$z_cat[grid.count.2014.GRID100$HOTSPOT < 3.091 & grid.count.2014.GRID100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2014.GRID100$z_cat)

# hotspots
grid.count.2014.GRID100.sel <- grid.count.2014.GRID100[grid.count.2014.GRID100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2014.GRID100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014.GRID100 <- st_union(grid.count.2014.GRID100.sel)
st_write(hotspot.2014.GRID100, "out/syria-iraq-2014_GRID100-hotspot.shp")


# 2016.GRID100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.GRID100 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.GRID100 <- concaveman(points.sf.2016.GRID100, concavity = 2, length_threshold = 0)
buffer.2016.GRID100 <- st_buffer(concave.hull.2016.GRID100, dist = 0.5)
st_write(buffer.2016.GRID100, "out/syria-iraq-conflict-2016_GRID100-shape.shp")

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
st_write(hotspot.2016.GRID100, "out/syria-iraq-2016_GRID100-hotspot.shp")


# 2014.GRID100-2016.GRID100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2014.GRID100.2016.GRID100 <- st_area(buffer.2014.GRID100)/100
new.perc.2014.GRID100.2016.GRID100 <- st_area(buffer.2016.GRID100)/one.perc.2014.GRID100.2016.GRID100
new.perc.2014.GRID100.2016.GRID100 <- drop_units(new.perc.2014.GRID100.2016.GRID100)
perc.change.2014.GRID100.2016.GRID100 <- new.perc.2014.GRID100.2016.GRID100 - 100
perc.change.2014.GRID100.2016.GRID100 #  34.21461 (expansion)

# overlap: hotspots
hotspot.2014.GRID100.sp <- as(st_geometry(hotspot.2014.GRID100), "Spatial")
hotspot.2016.GRID100.sp <- as(st_geometry(hotspot.2016.GRID100), "Spatial")
intersection.h.2014.GRID100.2016.GRID100 <- gIntersection(hotspot.2014.GRID100.sp, hotspot.2016.GRID100.sp)
intersection.h.2014.GRID100.2016.GRID100.sf <- st_as_sf(intersection.h.2014.GRID100.2016.GRID100) # no overlap
overl_hotspot.2014.GRID100.2016.GRID100 <- st_area(intersection.h.2014.GRID100.2016.GRID100.sf)/(st_area(hotspot.2014.GRID100)/100) 
overl_hotspot.2014.GRID100.2016.GRID100 # 67.97917 

# overlap shape
buffer.2014.GRID100.sp <- as(st_geometry(buffer.2014.GRID100), "Spatial")
buffer.2016.GRID100.sp <- as(st_geometry(buffer.2016.GRID100), "Spatial")
intersection.2014.GRID100.2016.GRID100 <- gIntersection(buffer.2014.GRID100.sp, buffer.2016.GRID100.sp)
intersection.2014.GRID100.2016.GRID100.sf <- st_as_sf(intersection.2014.GRID100.2016.GRID100)
overl_shape.2014.GRID100.2016.GRID100 <- st_area(intersection.2014.GRID100.2016.GRID100.sf)/(st_area(buffer.2014.GRID100)/100)
overl_shape.2014.GRID100.2016.GRID100 # 93.0838

# spatial change
spatial_ch.2014.GRID100.2016.GRID100 <- (overl_hotspot.2014.GRID100.2016.GRID100 + overl_shape.2014.GRID100.2016.GRID100)/2
spatial_ch.2014.GRID100.2016.GRID100 # 80.53148 (no shift)


# 2010.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2010.BUFF75 <- points.sf %>% filter(year == 2010)

# conflict shape
concave.hull.2010.BUFF75 <- concaveman(points.sf.2010.BUFF75, concavity = 2, length_threshold = 0)
buffer.2010.BUFF75 <- st_buffer(concave.hull.2010.BUFF75, dist = 0.75)
st_write(buffer.2010.BUFF75, "out/syria-iraq-conflict-2010_BUFF75-shape.shp")

grid.clip.2010.BUFF75 <- st_make_grid(buffer.2010.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2010.BUFF75)

points.count.2010.BUFF75 <- st_intersects(grid.clip.2010.BUFF75, points.sf.2010.BUFF75)
grid.count.2010.BUFF75 <- st_sf(pt_count = lengths(points.count.2010.BUFF75), 
                                geometry = st_cast(grid.clip.2010.BUFF75, "MULTIPOLYGON"))

neighbors.2010.BUFF75 <- poly2nb(grid.count.2010.BUFF75, queen = T)
weighted.neighbors.2010.BUFF75 <- nb2listw(include.self(neighbors.2010.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2010.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2010.BUFF75$pt_count, 
                                                        weighted.neighbors.2010.BUFF75))

grid.count.2010.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2010.BUFF75$z_cat[grid.count.2010.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010.BUFF75$z_cat[grid.count.2010.BUFF75$HOTSPOT < 3.091 & grid.count.2010.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010.BUFF75$z_cat)

# hotspots
grid.count.2010.BUFF75.sel <- grid.count.2010.BUFF75[grid.count.2010.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2010.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010.BUFF75 <- st_union(grid.count.2010.BUFF75.sel)
st_write(hotspot.2010.BUFF75, "out/syria-iraq-2010_BUFF75-hotspot.shp")


# 2013.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2013.BUFF75 <- points.sf %>% filter(year == 2013)

# conflict shape
concave.hull.2013.BUFF75 <- concaveman(points.sf.2013.BUFF75, concavity = 2, length_threshold = 0)
buffer.2013.BUFF75 <- st_buffer(concave.hull.2013.BUFF75, dist = 0.75)
st_write(buffer.2013.BUFF75, "out/syria-iraq-conflict-2013_BUFF75-shape.shp")

grid.clip.2013.BUFF75 <- st_make_grid(buffer.2013.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2013.BUFF75)

points.count.2013.BUFF75 <- st_intersects(grid.clip.2013.BUFF75, points.sf.2013.BUFF75)
grid.count.2013.BUFF75 <- st_sf(pt_count = lengths(points.count.2013.BUFF75), 
                                geometry = st_cast(grid.clip.2013.BUFF75, "MULTIPOLYGON"))

neighbors.2013.BUFF75 <- poly2nb(grid.count.2013.BUFF75, queen = T)
weighted.neighbors.2013.BUFF75 <- nb2listw(include.self(neighbors.2013.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2013.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2013.BUFF75$pt_count, 
                                                        weighted.neighbors.2013.BUFF75))

grid.count.2013.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2013.BUFF75$z_cat[grid.count.2013.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013.BUFF75$z_cat[grid.count.2013.BUFF75$HOTSPOT < 3.091 & grid.count.2013.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013.BUFF75$z_cat)

# hotspots
grid.count.2013.BUFF75.sel <- grid.count.2013.BUFF75[grid.count.2013.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2013.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013.BUFF75 <- st_union(grid.count.2013.BUFF75.sel)
st_write(hotspot.2013.BUFF75, "out/syria-iraq-2013_BUFF75-hotspot.shp")


# 2010.BUFF75-2013.BUFF75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2010.BUFF75.2013.BUFF75 <- st_area(buffer.2010.BUFF75)/100
new.perc.2010.BUFF75.2013.BUFF75 <- st_area(buffer.2013.BUFF75)/one.perc.2010.BUFF75.2013.BUFF75
new.perc.2010.BUFF75.2013.BUFF75 <- drop_units(new.perc.2010.BUFF75.2013.BUFF75)
perc.change.2010.BUFF75.2013.BUFF75 <- new.perc.2010.BUFF75.2013.BUFF75 - 100
perc.change.2010.BUFF75.2013.BUFF75 #  34.79026 (expansion)

# overlap: hotspots
hotspot.2010.BUFF75.sp <- as(st_geometry(hotspot.2010.BUFF75), "Spatial")
hotspot.2013.BUFF75.sp <- as(st_geometry(hotspot.2013.BUFF75), "Spatial")
intersection.h.2010.BUFF75.2013.BUFF75 <- gIntersection(hotspot.2010.BUFF75.sp, hotspot.2013.BUFF75.sp)
intersection.h.2010.BUFF75.2013.BUFF75.sf <- st_as_sf(intersection.h.2010.BUFF75.2013.BUFF75) # no overlap
overl_hotspot.2010.BUFF75.2013.BUFF75 <- 0
overl_hotspot.2010.BUFF75.2013.BUFF75

# overlap shape
buffer.2010.BUFF75.sp <- as(st_geometry(buffer.2010.BUFF75), "Spatial")
buffer.2013.BUFF75.sp <- as(st_geometry(buffer.2013.BUFF75), "Spatial")
intersection.2010.BUFF75.2013.BUFF75 <- gIntersection(buffer.2010.BUFF75.sp, buffer.2013.BUFF75.sp)
intersection.2010.BUFF75.2013.BUFF75.sf <- st_as_sf(intersection.2010.BUFF75.2013.BUFF75)
overl_shape.2010.BUFF75.2013.BUFF75 <- st_area(intersection.2010.BUFF75.2013.BUFF75.sf)/(st_area(buffer.2010.BUFF75)/100)
overl_shape.2010.BUFF75.2013.BUFF75 # 64.3497

# spatial change
spatial_ch.2010.BUFF75.2013.BUFF75 <- (overl_hotspot.2010.BUFF75.2013.BUFF75 + drop_units(overl_shape.2010.BUFF75.2013.BUFF75))/2
spatial_ch.2010.BUFF75.2013.BUFF75 # 32.17485 (shift)


# 2014.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2014.BUFF75 <- points.sf %>% filter(year == 2014)

# conflict shape
concave.hull.2014.BUFF75 <- concaveman(points.sf.2014.BUFF75, concavity = 2, length_threshold = 0)
buffer.2014.BUFF75 <- st_buffer(concave.hull.2014.BUFF75, dist = 0.75)
st_write(buffer.2014.BUFF75, "out/syria-iraq-conflict-2014_BUFF75-shape.shp")

grid.clip.2014.BUFF75 <- st_make_grid(buffer.2014.BUFF75, 
                                      cellsize = c(0.5,0.5), 
                                      what="polygons") %>%
  st_intersection(buffer.2014.BUFF75)

points.count.2014.BUFF75 <- st_intersects(grid.clip.2014.BUFF75, points.sf.2014.BUFF75)
grid.count.2014.BUFF75 <- st_sf(pt_count = lengths(points.count.2014.BUFF75), 
                                geometry = st_cast(grid.clip.2014.BUFF75, "MULTIPOLYGON"))

neighbors.2014.BUFF75 <- poly2nb(grid.count.2014.BUFF75, queen = T)
weighted.neighbors.2014.BUFF75 <- nb2listw(include.self(neighbors.2014.BUFF75), 
                                           style = "B", 
                                           zero.policy = T)
grid.count.2014.BUFF75$HOTSPOT <- as.vector(localG_perm(grid.count.2014.BUFF75$pt_count, 
                                                        weighted.neighbors.2014.BUFF75))

grid.count.2014.BUFF75$z_cat <- "Not significant or cold spot"
grid.count.2014.BUFF75$z_cat[grid.count.2014.BUFF75$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2014.BUFF75$z_cat[grid.count.2014.BUFF75$HOTSPOT < 3.091 & grid.count.2014.BUFF75$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2014.BUFF75$z_cat)

# hotspots
grid.count.2014.BUFF75.sel <- grid.count.2014.BUFF75[grid.count.2014.BUFF75$z_cat == "Hotspot (P-value: 0.001)" |
                                                       grid.count.2014.BUFF75$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014.BUFF75 <- st_union(grid.count.2014.BUFF75.sel)
st_write(hotspot.2014.BUFF75, "out/syria-iraq-2014_BUFF75-hotspot.shp")


# 2016.BUFF75: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.BUFF75 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.BUFF75 <- concaveman(points.sf.2016.BUFF75, concavity = 2, length_threshold = 0)
buffer.2016.BUFF75 <- st_buffer(concave.hull.2016.BUFF75, dist = 0.75)
st_write(buffer.2016.BUFF75, "out/syria-iraq-conflict-2016_BUFF75-shape.shp")

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
st_write(hotspot.2016.BUFF75, "out/syria-iraq-2016_BUFF75-hotspot.shp")


# 2014.BUFF75-2016.BUFF75: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2014.BUFF75.2016.BUFF75 <- st_area(buffer.2014.BUFF75)/100
new.perc.2014.BUFF75.2016.BUFF75 <- st_area(buffer.2016.BUFF75)/one.perc.2014.BUFF75.2016.BUFF75
new.perc.2014.BUFF75.2016.BUFF75 <- drop_units(new.perc.2014.BUFF75.2016.BUFF75)
perc.change.2014.BUFF75.2016.BUFF75 <- new.perc.2014.BUFF75.2016.BUFF75 - 100
perc.change.2014.BUFF75.2016.BUFF75 #  30.47522 (expansion)

# overlap: hotspots
hotspot.2014.BUFF75.sp <- as(st_geometry(hotspot.2014.BUFF75), "Spatial")
hotspot.2016.BUFF75.sp <- as(st_geometry(hotspot.2016.BUFF75), "Spatial")
intersection.h.2014.BUFF75.2016.BUFF75 <- gIntersection(hotspot.2014.BUFF75.sp, hotspot.2016.BUFF75.sp)
intersection.h.2014.BUFF75.2016.BUFF75.sf <- st_as_sf(intersection.h.2014.BUFF75.2016.BUFF75) # no overlap
overl_hotspot.2014.BUFF75.2016.BUFF75 <- st_area(intersection.h.2014.BUFF75.2016.BUFF75.sf)/(st_area(hotspot.2014.BUFF75)/100) 
overl_hotspot.2014.BUFF75.2016.BUFF75 # 69.59328

# overlap shape
buffer.2014.BUFF75.sp <- as(st_geometry(buffer.2014.BUFF75), "Spatial")
buffer.2016.BUFF75.sp <- as(st_geometry(buffer.2016.BUFF75), "Spatial")
intersection.2014.BUFF75.2016.BUFF75 <- gIntersection(buffer.2014.BUFF75.sp, buffer.2016.BUFF75.sp)
intersection.2014.BUFF75.2016.BUFF75.sf <- st_as_sf(intersection.2014.BUFF75.2016.BUFF75)
overl_shape.2014.BUFF75.2016.BUFF75 <- st_area(intersection.2014.BUFF75.2016.BUFF75.sf)/(st_area(buffer.2014.BUFF75)/100)
overl_shape.2014.BUFF75.2016.BUFF75 # 93.38716

# spatial change
spatial_ch.2014.BUFF75.2016.BUFF75 <- (overl_hotspot.2014.BUFF75.2016.BUFF75 + overl_shape.2014.BUFF75.2016.BUFF75)/2
spatial_ch.2014.BUFF75.2016.BUFF75 # 81.49022 (no shift)


# 2010.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2010.BUFF100 <- points.sf %>% filter(year == 2010)

# conflict shape
concave.hull.2010.BUFF100 <- concaveman(points.sf.2010.BUFF100, concavity = 2, length_threshold = 0)
buffer.2010.BUFF100 <- st_buffer(concave.hull.2010.BUFF100, dist = 1)
st_write(buffer.2010.BUFF100, "out/syria-iraq-conflict-2010_BUFF100-shape.shp")

grid.clip.2010.BUFF100 <- st_make_grid(buffer.2010.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2010.BUFF100)

points.count.2010.BUFF100 <- st_intersects(grid.clip.2010.BUFF100, points.sf.2010.BUFF100)
grid.count.2010.BUFF100 <- st_sf(pt_count = lengths(points.count.2010.BUFF100), 
                                 geometry = st_cast(grid.clip.2010.BUFF100, "MULTIPOLYGON"))

neighbors.2010.BUFF100 <- poly2nb(grid.count.2010.BUFF100, queen = T)
weighted.neighbors.2010.BUFF100 <- nb2listw(include.self(neighbors.2010.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2010.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2010.BUFF100$pt_count, 
                                                         weighted.neighbors.2010.BUFF100))

grid.count.2010.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2010.BUFF100$z_cat[grid.count.2010.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010.BUFF100$z_cat[grid.count.2010.BUFF100$HOTSPOT < 3.091 & grid.count.2010.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010.BUFF100$z_cat)

# hotspots
grid.count.2010.BUFF100.sel <- grid.count.2010.BUFF100[grid.count.2010.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2010.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010.BUFF100 <- st_union(grid.count.2010.BUFF100.sel)
st_write(hotspot.2010.BUFF100, "out/syria-iraq-2010_BUFF100-hotspot.shp")


# 2013.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2013.BUFF100 <- points.sf %>% filter(year == 2013)

# conflict shape
concave.hull.2013.BUFF100 <- concaveman(points.sf.2013.BUFF100, concavity = 2, length_threshold = 0)
buffer.2013.BUFF100 <- st_buffer(concave.hull.2013.BUFF100, dist = 1)
st_write(buffer.2013.BUFF100, "out/syria-iraq-conflict-2013_BUFF100-shape.shp")

grid.clip.2013.BUFF100 <- st_make_grid(buffer.2013.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2013.BUFF100)

points.count.2013.BUFF100 <- st_intersects(grid.clip.2013.BUFF100, points.sf.2013.BUFF100)
grid.count.2013.BUFF100 <- st_sf(pt_count = lengths(points.count.2013.BUFF100), 
                                 geometry = st_cast(grid.clip.2013.BUFF100, "MULTIPOLYGON"))

neighbors.2013.BUFF100 <- poly2nb(grid.count.2013.BUFF100, queen = T)
weighted.neighbors.2013.BUFF100 <- nb2listw(include.self(neighbors.2013.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2013.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2013.BUFF100$pt_count, 
                                                         weighted.neighbors.2013.BUFF100))

grid.count.2013.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2013.BUFF100$z_cat[grid.count.2013.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013.BUFF100$z_cat[grid.count.2013.BUFF100$HOTSPOT < 3.091 & grid.count.2013.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013.BUFF100$z_cat)

# hotspots
grid.count.2013.BUFF100.sel <- grid.count.2013.BUFF100[grid.count.2013.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2013.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013.BUFF100 <- st_union(grid.count.2013.BUFF100.sel)
st_write(hotspot.2013.BUFF100, "out/syria-iraq-2013_BUFF100-hotspot.shp")


# 2010.BUFF100-2013.BUFF100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2010.BUFF100.2013.BUFF100 <- st_area(buffer.2010.BUFF100)/100
new.perc.2010.BUFF100.2013.BUFF100 <- st_area(buffer.2013.BUFF100)/one.perc.2010.BUFF100.2013.BUFF100
new.perc.2010.BUFF100.2013.BUFF100 <- drop_units(new.perc.2010.BUFF100.2013.BUFF100)
perc.change.2010.BUFF100.2013.BUFF100 <- new.perc.2010.BUFF100.2013.BUFF100 - 100
perc.change.2010.BUFF100.2013.BUFF100 #  30.42688 (expansion)

# overlap: hotspots
hotspot.2010.BUFF100.sp <- as(st_geometry(hotspot.2010.BUFF100), "Spatial")
hotspot.2013.BUFF100.sp <- as(st_geometry(hotspot.2013.BUFF100), "Spatial")
intersection.h.2010.BUFF100.2013.BUFF100 <- gIntersection(hotspot.2010.BUFF100.sp, hotspot.2013.BUFF100.sp)
intersection.h.2010.BUFF100.2013.BUFF100.sf <- st_as_sf(intersection.h.2010.BUFF100.2013.BUFF100) # no overlap
overl_hotspot.2010.BUFF100.2013.BUFF100 <- 0
overl_hotspot.2010.BUFF100.2013.BUFF100

# overlap shape
buffer.2010.BUFF100.sp <- as(st_geometry(buffer.2010.BUFF100), "Spatial")
buffer.2013.BUFF100.sp <- as(st_geometry(buffer.2013.BUFF100), "Spatial")
intersection.2010.BUFF100.2013.BUFF100 <- gIntersection(buffer.2010.BUFF100.sp, buffer.2013.BUFF100.sp)
intersection.2010.BUFF100.2013.BUFF100.sf <- st_as_sf(intersection.2010.BUFF100.2013.BUFF100)
overl_shape.2010.BUFF100.2013.BUFF100 <- st_area(intersection.2010.BUFF100.2013.BUFF100.sf)/(st_area(buffer.2010.BUFF100)/100)
overl_shape.2010.BUFF100.2013.BUFF100 # 68.29977

# spatial change
spatial_ch.2010.BUFF100.2013.BUFF100 <- (overl_hotspot.2010.BUFF100.2013.BUFF100 + drop_units(overl_shape.2010.BUFF100.2013.BUFF100))/2
spatial_ch.2010.BUFF100.2013.BUFF100 # 34.14988 (shift)


# 2014.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2014.BUFF100 <- points.sf %>% filter(year == 2014)

# conflict shape
concave.hull.2014.BUFF100 <- concaveman(points.sf.2014.BUFF100, concavity = 2, length_threshold = 0)
buffer.2014.BUFF100 <- st_buffer(concave.hull.2014.BUFF100, dist = 1)
st_write(buffer.2014.BUFF100, "out/syria-iraq-conflict-2014_BUFF100-shape.shp")

grid.clip.2014.BUFF100 <- st_make_grid(buffer.2014.BUFF100, 
                                       cellsize = c(0.5,0.5), 
                                       what="polygons") %>%
  st_intersection(buffer.2014.BUFF100)

points.count.2014.BUFF100 <- st_intersects(grid.clip.2014.BUFF100, points.sf.2014.BUFF100)
grid.count.2014.BUFF100 <- st_sf(pt_count = lengths(points.count.2014.BUFF100), 
                                 geometry = st_cast(grid.clip.2014.BUFF100, "MULTIPOLYGON"))

neighbors.2014.BUFF100 <- poly2nb(grid.count.2014.BUFF100, queen = T)
weighted.neighbors.2014.BUFF100 <- nb2listw(include.self(neighbors.2014.BUFF100), 
                                            style = "B", 
                                            zero.policy = T)
grid.count.2014.BUFF100$HOTSPOT <- as.vector(localG_perm(grid.count.2014.BUFF100$pt_count, 
                                                         weighted.neighbors.2014.BUFF100))

grid.count.2014.BUFF100$z_cat <- "Not significant or cold spot"
grid.count.2014.BUFF100$z_cat[grid.count.2014.BUFF100$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2014.BUFF100$z_cat[grid.count.2014.BUFF100$HOTSPOT < 3.091 & grid.count.2014.BUFF100$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2014.BUFF100$z_cat)

# hotspots
grid.count.2014.BUFF100.sel <- grid.count.2014.BUFF100[grid.count.2014.BUFF100$z_cat == "Hotspot (P-value: 0.001)" |
                                                         grid.count.2014.BUFF100$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2014.BUFF100 <- st_union(grid.count.2014.BUFF100.sel)
st_write(hotspot.2014.BUFF100, "out/syria-iraq-2014_BUFF100-hotspot.shp")


# 2016.BUFF100: Conflict shape and hotspots ----------------------------------------------------

points.sf.2016.BUFF100 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016.BUFF100 <- concaveman(points.sf.2016.BUFF100, concavity = 2, length_threshold = 0)
buffer.2016.BUFF100 <- st_buffer(concave.hull.2016.BUFF100, dist = 1)
st_write(buffer.2016.BUFF100, "out/syria-iraq-conflict-2016_BUFF100-shape.shp")

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
st_write(hotspot.2016.BUFF100, "out/syria-iraq-2016_BUFF100-hotspot.shp")


# 2014.BUFF100-2016.BUFF100: Spatial change -----------------------------------------------

# contraction/expansion
one.perc.2014.BUFF100.2016.BUFF100 <- st_area(buffer.2014.BUFF100)/100
new.perc.2014.BUFF100.2016.BUFF100 <- st_area(buffer.2016.BUFF100)/one.perc.2014.BUFF100.2016.BUFF100
new.perc.2014.BUFF100.2016.BUFF100 <- drop_units(new.perc.2014.BUFF100.2016.BUFF100)
perc.change.2014.BUFF100.2016.BUFF100 <- new.perc.2014.BUFF100.2016.BUFF100 - 100
perc.change.2014.BUFF100.2016.BUFF100 #  427.90947(expansion)

# overlap: hotspots
hotspot.2014.BUFF100.sp <- as(st_geometry(hotspot.2014.BUFF100), "Spatial")
hotspot.2016.BUFF100.sp <- as(st_geometry(hotspot.2016.BUFF100), "Spatial")
intersection.h.2014.BUFF100.2016.BUFF100 <- gIntersection(hotspot.2014.BUFF100.sp, hotspot.2016.BUFF100.sp)
intersection.h.2014.BUFF100.2016.BUFF100.sf <- st_as_sf(intersection.h.2014.BUFF100.2016.BUFF100) # no overlap
overl_hotspot.2014.BUFF100.2016.BUFF100 <- st_area(intersection.h.2014.BUFF100.2016.BUFF100.sf)/(st_area(hotspot.2014.BUFF100)/100) 
overl_hotspot.2014.BUFF100.2016.BUFF100 # 76.39964

# overlap shape
buffer.2014.BUFF100.sp <- as(st_geometry(buffer.2014.BUFF100), "Spatial")
buffer.2016.BUFF100.sp <- as(st_geometry(buffer.2016.BUFF100), "Spatial")
intersection.2014.BUFF100.2016.BUFF100 <- gIntersection(buffer.2014.BUFF100.sp, buffer.2016.BUFF100.sp)
intersection.2014.BUFF100.2016.BUFF100.sf <- st_as_sf(intersection.2014.BUFF100.2016.BUFF100)
overl_shape.2014.BUFF100.2016.BUFF100 <- st_area(intersection.2014.BUFF100.2016.BUFF100.sf)/(st_area(buffer.2014.BUFF100)/100)
overl_shape.2014.BUFF100.2016.BUFF100 # 93.7616

# spatial change
spatial_ch.2014.BUFF100.2016.BUFF100 <- (overl_shape.2014.BUFF100.2016.BUFF100 + overl_shape.2014.BUFF100.2016.BUFF100)/2
spatial_ch.2014.BUFF100.2016.BUFF100 # 93.7616 (no shift)


# End of script -----------------------------------------------------------

