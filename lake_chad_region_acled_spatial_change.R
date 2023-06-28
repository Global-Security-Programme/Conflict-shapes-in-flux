### Preamble ###############################################################
# Armed conflict in the Lake Chad region - ACLED: Spatial change analysis
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
data <- readRDS("data/lake-chad-region-acled.rds")

data$dyad_name <- paste0(as.character(data$actor1),
                         " - ", 
                         as.character(data$actor2))

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

edges.2011 <- data.2011 %>% dplyr::select(actor1,actor2, weight)
edges.2011 <- edges.2011 %>% rename(from = actor1,
                                    to = actor2)
edges.2011 <- edges.2011 %>% distinct(from, to, weight)

side.a.2011 <- data.2011 %>% distinct(actor1) %>% 
  rename(label = actor1)
side.b.2011 <- data.2011 %>% distinct(actor2) %>% 
  rename(label = actor2)
nodes.2011 <- rbind(side.a.2011, side.b.2011)
nodes.2011 <- nodes.2011 %>% distinct(label)
nodes.2011$id <- 1:nrow(nodes.2011)  

net.2011 <- graph_from_data_frame(edges.2011, 
                                  vertices = nodes.2011,
                                  directed = F)

# degree centrality 
nodes.2011$nodes_n <- igraph::degree(net.2011)
nodes.2011$degree <- strength(net.2011)
write_csv(nodes.2011, "out/lake-chad-region-acled-2011-centrality-measures.csv")


# 2011: Conflict shape and hotspots ----------------------------------------------------

points.sf.2011 <- points.sf %>% filter(year == 2011)

# conflict shape
concave.hull.2011 <- concaveman(points.sf.2011, concavity = 2, length_threshold = 0)
buffer.2011 <- st_buffer(concave.hull.2011, dist = 0.5)
st_write(buffer.2011, "out/lake-chad-region-acled-2011-conflict-shape.shp")

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
  ggtitle("Armed conflict in the Lake Chad region, 2011: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/lake-chad-region-acled-2011-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2011.sel <- grid.count.2011[grid.count.2011$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2011$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2011 <- st_union(grid.count.2011.sel)
st_write(hotspot.2011, "out/lake-chad-region-acled-2011-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2011,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2011, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in the Lake Chad region, 2011: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/lake-chad-region-acled-2011-events.png", width = 8, height = 4)


# 2016: Dominant actors ------------------------------------------------------------

data.2016 <- data %>% filter(year == 2016)

data.2016 <- data.2016 %>% group_by(dyad_name) %>% 
  mutate(weight = n()) %>% ungroup()

edges.2016 <- data.2016 %>% dplyr::select(actor1,actor2, weight)
edges.2016 <- edges.2016 %>% rename(from = actor1,
                                    to = actor2)
edges.2016 <- edges.2016 %>% distinct(from, to, weight)

side.a.2016 <- data.2016 %>% distinct(actor1) %>% 
  rename(label = actor1)
side.b.2016 <- data.2016 %>% distinct(actor2) %>% 
  rename(label = actor2)
nodes.2016 <- rbind(side.a.2016, side.b.2016)
nodes.2016 <- nodes.2016 %>% distinct(label)
nodes.2016$id <- 1:nrow(nodes.2016)  

net.2016 <- graph_from_data_frame(edges.2016, 
                                  vertices = nodes.2016,
                                  directed = F)

# degree centrality 
nodes.2016$nodes_n <- igraph::degree(net.2016)
nodes.2016$degree <- strength(net.2016)
write_csv(nodes.2016, "out/lake-chad-region-acled-2016-centrality-measures.csv")


# 2016: Conflict shape and hotspots --------------------------------------------------------

points.sf.2016 <- points.sf %>% filter(year == 2016)

# conflict shape
concave.hull.2016 <- concaveman(points.sf.2016, concavity = 2, length_threshold = 0)
buffer.2016 <- st_buffer(concave.hull.2016, dist = 0.5)
st_write(buffer.2016, "out/lake-chad-region-acled-2016-conflict-shape.shp")

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
  ggtitle("Armed conflict in the Lake Chad region, 2016: Conflict shape and hotspots") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  labs(fill = "Gi* Statistic") +
  theme(legend.position = "bottom") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/lake-chad-region-acled-2016-conflict-shape-hotspot.png", width = 8, height = 4)

# hotspots
grid.count.2016.sel <- grid.count.2016[grid.count.2016$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2016$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2016 <- st_union(grid.count.2016.sel)
st_write(hotspot.2016, "out/lake-chad-region-acled-2016-hotspot.shp")

# conflict events
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2016,
          color = "red", alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2016, alpha = 0.3,  fill = "grey", color = NA) +
  ggtitle("Armed conflict in the Lake Chad region, 2016: Conflict shape and conflict events") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(0, 20), ylim = c(3, 19), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/lake-chad-region-acled-2016-events.png", width = 8, height = 4)


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
  ggtitle("Armed conflict in the Lake Chad region - ACLED, 2011-2016: Spatial change") +
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
ggsave("figs/lake-chad-region-acled-2011-2016-spatial-change.png", width = 8, height = 4)

# contraction/expansion
one.perc.2011.2016 <- st_area(buffer.2011)/100
new.perc.2011.2016 <- st_area(buffer.2016)/one.perc.2011.2016
new.perc.2011.2016 <- drop_units(new.perc.2011.2016)
perc.change.2011.2016 <- new.perc.2011.2016 - 100
perc.change.2011.2016 # -46.73145 (contraction)

# overlap: hotspots
hotspot.2011.sp <- as(st_geometry(hotspot.2011), "Spatial")
hotspot.2016.sp <- as(st_geometry(hotspot.2016), "Spatial")
intersection.h.2011.2016 <- gIntersection(hotspot.2011.sp, hotspot.2016.sp)
intersection.h.2011.2016.sf <- st_as_sf(intersection.h.2011.2016)
overl_hotspot.2011.2016 <- st_area(intersection.h.2011.2016.sf)/(st_area(hotspot.2016)/100) 
overl_hotspot.2011.2016 # 17.5844

# overlap: conflict shape
buffer.2011.sp <- as(st_geometry(buffer.2011), "Spatial")
buffer.2016.sp <- as(st_geometry(buffer.2016), "Spatial")
intersection.2011.2016 <- gIntersection(buffer.2011.sp, buffer.2016.sp)
intersection.2011.2016.sf <- st_as_sf(intersection.2011.2016)
overl_shape.2011.2016 <- st_area(intersection.2011.2016.sf)/(st_area(buffer.2016)/100) 
overl_shape.2011.2016 # 48.34024

# spatial change
spatial_ch.2011.2016 <- (overl_hotspot.2011.2016 + overl_shape.2011.2016)/2
spatial_ch.2011.2016 # 32.96232 (shift)

