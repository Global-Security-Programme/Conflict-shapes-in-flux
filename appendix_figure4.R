### Preamble ###############################################################
# Figure 4 in the Appendix
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

world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(world) <- 4326

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

data <- readRDS("data/syria-iraq.rds")

points.sf <- st_as_sf(data, 
                      coords = c("longitude", "latitude"))
st_crs(points.sf) <- 4326


# Figure 4.1 -----------------------------------------------

# 2010
points.sf.2010 <- points.sf %>% filter(year == 2010)

concave.hull.2010 <- concaveman(points.sf.2010, concavity = 2, length_threshold = 0)
buffer.2010 <- st_buffer(concave.hull.2010, dist = 0.5)

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

grid.count.2010.sel <- grid.count.2010[grid.count.2010$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2010$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2010 <- st_union(grid.count.2010.sel)

# 2013
points.sf.2013 <- points.sf %>% filter(year == 2013)

concave.hull.2013 <- concaveman(points.sf.2013, concavity = 2, length_threshold = 0)
buffer.2013 <- st_buffer(concave.hull.2013, dist = 0.5)

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

grid.count.2013.sel <- grid.count.2013[grid.count.2013$z_cat == "Hotspot (P-value: 0.001)" |
                                         grid.count.2013$z_cat == "Hotspot (P-value: 0.05)", ] 
hotspot.2013 <- st_union(grid.count.2013.sel)


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


# Figure 4.2 -----------------------------------------------

# 2010
points.sf.2010.dyad <- points.sf.2010 %>% 
  filter(dyad_name == "Government of Iraq - IS")

concave.hull.2010.dyad <- concaveman(points.sf.2010.dyad, concavity = 2, length_threshold = 0)
buffer.2010.dyad <- st_buffer(concave.hull.2010.dyad, dist = 0.5)

grid.clip.2010.dyad <- st_make_grid(buffer.2010.dyad, 
                                    cellsize = c(0.5,0.5), 
                                    what="polygons") %>%
  st_intersection(buffer.2010.dyad)

points.count.2010.dyad <- st_intersects(grid.clip.2010.dyad, points.sf.2010.dyad)
grid.count.2010.dyad <- st_sf(pt_count = lengths(points.count.2010.dyad), 
                              geometry = st_cast(grid.clip.2010.dyad, "MULTIPOLYGON"))


neighbors.2010.dyad <- poly2nb(grid.count.2010.dyad, queen = T)
weighted.neighbors.2010.dyad <- nb2listw(include.self(neighbors.2010.dyad), 
                                         style = "B", 
                                         zero.policy = T)

grid.count.2010.dyad$HOTSPOT <- as.vector(localG_perm(grid.count.2010.dyad$pt_count, 
                                                      weighted.neighbors.2010.dyad))

grid.count.2010.dyad$z_cat <- "Not significant or cold spot"
grid.count.2010.dyad$z_cat[grid.count.2010.dyad$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2010.dyad$z_cat[grid.count.2010.dyad$HOTSPOT < 3.091 & grid.count.2010.dyad$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2010.dyad$z_cat)


grid.count.2010.dyad.sel <- grid.count.2010.dyad[grid.count.2010.dyad$z_cat == "Hotspot (P-value: 0.001)" |
                                                   grid.count.2010.dyad$z_cat == "Hotspot (P-value: 0.05)", ] 

hotspot.2010.dyad <- st_union(grid.count.2010.dyad.sel)


# 2013
points.sf.2013.dyad <- points.sf.2013 %>% 
  filter(dyad_name == "Government of Iraq - IS")

concave.hull.2013.dyad <- concaveman(points.sf.2013.dyad, concavity = 2, length_threshold = 0)
buffer.2013.dyad <- st_buffer(concave.hull.2013.dyad, dist = 0.5)

grid.clip.2013.dyad <- st_make_grid(buffer.2013.dyad, 
                                    cellsize = c(0.5,0.5), 
                                    what="polygons") %>%
  st_intersection(buffer.2013.dyad)

points.count.2013.dyad <- st_intersects(grid.clip.2013.dyad, points.sf.2013.dyad)
grid.count.2013.dyad <- st_sf(pt_count = lengths(points.count.2013.dyad), 
                              geometry = st_cast(grid.clip.2013.dyad, "MULTIPOLYGON"))


neighbors.2013.dyad <- poly2nb(grid.count.2013.dyad, queen = T)
weighted.neighbors.2013.dyad <- nb2listw(include.self(neighbors.2013.dyad), 
                                         style = "B", 
                                         zero.policy = T)

grid.count.2013.dyad$HOTSPOT <- as.vector(localG_perm(grid.count.2013.dyad$pt_count, 
                                                      weighted.neighbors.2013.dyad))

grid.count.2013.dyad$z_cat <- "Not significant or cold spot"
grid.count.2013.dyad$z_cat[grid.count.2013.dyad$HOTSPOT >= 3.091] <- "Hotspot (P-value: 0.001)"
grid.count.2013.dyad$z_cat[grid.count.2013.dyad$HOTSPOT < 3.091 & grid.count.2013.dyad$HOTSPOT >= 1.96] <- "Hotspot (P-value: 0.05)"
table(grid.count.2013.dyad$z_cat)

grid.count.2013.dyad.sel <- grid.count.2013.dyad[grid.count.2013.dyad$z_cat == "Hotspot (P-value: 0.001)" |
                                                   grid.count.2013.dyad$z_cat == "Hotspot (P-value: 0.05)", ] 

hotspot.2013.dyad <- st_union(grid.count.2013.dyad.sel)

# map
label.2010.dyad <- data.frame(
  x = 42, 
  y = 39,
  text = "2010")

label.2013.dyad <- data.frame(
  x = 51, 
  y = 31, 
  text = "2013")

ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = hotspot.2010.dyad, alpha = 0.6, fill = "blue", color = NA) +
  geom_sf(data = buffer.2010.dyad, alpha = 0.3, fill = "blue", color = NA) +
  geom_sf(data = hotspot.2013.dyad, alpha = 0.6, fill = "red", color = NA) +
  geom_sf(data = buffer.2013.dyad, alpha = 0.3, fill = "red", color = NA) +
  ggtitle("Armed conflict in Syria and Iraq, 2010-2013 - dyad: Spatial change") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data = world_points,
            aes(x = X, y = Y, label = name),
            color = "gray20", 
            fontface = "italic", 
            check_overlap = T,
            size = 3) +
  geom_label(data = label.2010.dyad, 
             aes(x = x, y = y, label = text), 
             color = "blue") +
  geom_label(data = label.2013.dyad, 
             aes(x = x, y = y, label = text), 
             color = "red")
ggsave("figs/syria-iraq-2010-2013-dyad-spatial-change.png", width = 8, height = 4)


# Figure 4.3 --------------------------------------------------------------

points.sf.2013$is_involved <- "IS not involved"
points.sf.2013$is_involved[points.sf.2013$side_a == "IS"] <- "IS involved"
points.sf.2013$is_involved[points.sf.2013$side_b == "IS"] <- "IS involved"

points.sf.2013.is.only <- points.sf.2013 %>% 
  filter(is_involved == "IS involved")

# map
ggplot() +
  geom_sf(data = world, alpha = 0.1) +
  geom_sf(data = points.sf.2013.is.only,
          color = "red", 
          alpha = 0.6, shape = 1) +
  geom_sf(data = buffer.2013, alpha = 0.3,  fill = "grey", color = NA) +

  ggtitle("Armed Conflict in Syria and Iraq, 2013: Only conflict events with IS involved") +
  theme_bw() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(26, 63), ylim = c(22, 43), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)
ggsave("figs/syria-iraq-2013-is-only.png", width = 8, height = 4)
