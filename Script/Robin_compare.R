setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(bipartite)
library(robin)

sf_use_s2(T)
### High res polygons

# Large Marine Ecosystems
LMEs_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/LMEs/", layer = "LMEs66") %>% 
  filter(LME_NAME == "Baltic Sea" | LME_NAME == "North Sea" | LME_NAME ==  "Gulf of Alaska" |
           LME_NAME == "Labrador - Newfoundland" | LME_NAME == "Celtic-Biscay Shelf" | LME_NAME == "Scotian Shelf" |
           LME_NAME == "California Current" | LME_NAME == "Iberian Coastal" | LME_NAME == "Northeast U.S. Continental Shelf" |
           LME_NAME == "Southeast U.S. Continental Shelf" | LME_NAME == "Gulf of Mexico" | LME_NAME == "Norwegian Sea" |
           LME_NAME == "Barents Sea" | LME_NAME == "East Bering Sea" | LME_NAME == "Aleutian Islands") %>% 
  mutate(Coast = case_when(LME_NAME == "Baltic Sea" ~ "Europe",
                           LME_NAME == "North Sea" ~ "Europe",
                           LME_NAME == "Gulf of Alaska" ~ "N_Pacific",
                           LME_NAME == "Labrador - Newfoundland" ~ "N_Atlantic",
                           LME_NAME == "Celtic-Biscay Shelf" ~ "Europe",
                           LME_NAME == "Scotian Shelf" ~ "N_Atlantic",
                           LME_NAME == "California Current" ~ "S_Pacific",
                           LME_NAME == "Iberian Coastal" ~ "Europe",
                           LME_NAME == "Northeast U.S. Continental Shelf" ~ "N_Atlantic",
                           LME_NAME == "Southeast U.S. Continental Shelf" ~ "S_Atlantic",
                           LME_NAME == "Gulf of Mexico" ~ "S_Atlantic",
                           LME_NAME == "Norwegian Sea" ~ "Europe",
                           LME_NAME == "Barents Sea" ~ "Europe",
                           LME_NAME == "East Bering Sea" ~ "N_Pacific",
                           LME_NAME == "Aleutian Islands" ~ "N_Pacific")) %>% 
  mutate(Ocean = case_when(Coast %in% c("S_Pacific", "N_Pacific") ~ "Pacific",
                           Coast %in% c("N_Atlantic", "S_Atlantic") ~ "W_atlantic",
                           Coast %in% c("Europe") ~ "E_atlantic")) %>%
  st_transform(st_crs("EPSG:3995")) %>% 
  st_make_valid()

LME_names <- as.data.frame(LMEs_shelf_seas[,c("LME_NUMBER",
                                              "LME_NAME",
                                              "Coast",
                                              "Ocean")])[,c(1,2,3,4)]


Coastline <- read_sf("../../../Feb2023_Transfer/Oceanographic/Data/GSHHS_shp/l/", layer = "GSHHS_l_L1") %>% 
  st_make_valid() %>% 
  st_transform(st_crs("EPSG:3995")) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))  

# Ecoregions
Ecoregions_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/MEOW/", layer = "meow_ecos")%>% 
  st_transform(st_crs("EPSG:3995")) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))

#####################
# Loading and processing trawl data
#####################

# Fishglob (Maureaud et al. 2024)
trawls <- read_csv("../../../Shelf_Evolution/Data/FishGlob/FishGlob_public_std_clean.csv") %>% 
  rename(taxon = accepted_name) %>% 
  mutate(Latitude = latitude, Longitude = longitude) %>%
  #full_join(LME_names) %>% 
  drop_na(longitude, latitude, Latitude, Longitude) %>% 
  filter(year > 1998) %>%
  filter((survey == "GSL-S") | (depth < 1501 & depth != "NA")) %>% 
  mutate(taxon = gsub(" ", "_", taxon)) %>% 
  filter(grepl("_", taxon))

# Spatial polygons for points (as occurrences)
sf_points <- st_as_sf(trawls, 
                      coords = c("longitude", "latitude"),
                      crs="EPSG:4326") %>% 
  #filter(Ocean == "Pacific") %>% 
  st_transform(st_crs(Ecoregions_shelf_seas))


index <- which(lengths(st_intersects(st_make_valid(Ecoregions_shelf_seas), st_make_valid(sf_points)))>0)
Ecoregions_shelf_seas <- Ecoregions_shelf_seas[index,]

# Cropping Antilles
Coastline_crop <- Coastline %>% 
  st_crop(st_bbox(Ecoregions_shelf_seas))

# Find minimum trawl numbers to subsample down
fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}

grid <- fishnet(sf_points, square = F, cellsize = c(200000, 200000))

#plot(st_geometry(sf_points), add = T, col ="red")

trawl_grid <- grid %>% 
  st_as_sf() %>% 
  mutate(grid_id = row_number()) 

## 100km 
# 92 = 1st quartile, 190 = median including depths <20m 
## 200km 
# 124 = 1st quartile, 337 = median including depths <20m 

World_grid <- st_join(sf_points, trawl_grid) %>% 
  #filter(grid_ID != 14 | grid_ID != 8) %>% 
  group_by(haul_id) %>% 
  mutate(Haul_ID_num = cur_group_id()) %>% 
  ungroup() %>% 
  group_by(grid_id) %>% 
  #count(Haul_ID_num) 
  mutate(unique = n_distinct(Haul_ID_num)) %>% 
  arrange(unique) %>% 
  ungroup() %>% 
  filter(grid_id != "NA") %>% 
  filter(unique > 111) 

############
# Subsample each grid cell down to x
#############

gridcells <- unique(World_grid$grid_id)
grid_ls <- list()  

for(i in gridcells){
  Cell_rarefy <- filter(World_grid, grid_id == i) %>% 
    group_by(unique) %>% 
    sample_n(size = 112) %>%
    ungroup()
  grid_ls[[i]] <- Cell_rarefy
  
}


grid_rarefy <- do.call("rbind", grid_ls)  

# How much is kept after rarefaction?

nrow(grid_rarefy)/nrow(sf_points)

#200 km
# 1st quartile = 1.5%
# median = 2.5%

# Construct graph from edge list
edgelist <- data.frame(as.character(grid_rarefy$grid_id), 
                       as.character(grid_rarefy$taxon),
                       stringsAsFactors = FALSE)
colnames(edgelist) <- c("grid_id", "taxon")

edgelist_u <- unique(edgelist)
nrow(edgelist_u)

# Weights by species occurrences
strength_links <- unlist(apply(edgelist_u, 1,
                               function(x)length(which(edgelist[,1]==x[1] & edgelist[,2]==x[2]))))

## this value must be the same as the total number of rows in the edgelist ---> 29264
sum(strength_links)

graph <- graph_from_edgelist(as.matrix(edgelist_u), directed = FALSE)
E(graph)$weight <- strength_links


#100km / 200km
res <- "200km"

#1stquartile / median
cutoff <- "1stquartile"

############################################################################
# Clustering
############################################################################

set.seed(1234)

# Leiden Algorithm (Traag et al. 2019)
Leiden_map <- cluster_leiden(graph,
                             objective_function = "modularity",
                             resolution_parameter = 2,
                             n_iterations = 20,
                             weights = strength_links)
modularity(graph, membership(Leiden_map), directed = F)

# Infomap algorithm ((Rosvall & Bergstrom, 2008))
Informap_map <- cluster_infomap(graph,
                                nb.trials = 100,
                                e.weights = strength_links)
modularity(graph, membership(Informap_map), directed = F)


####### 
# Statistical comparisons between algorithms with "robin"
# testing stability of clustering (partition) under random perturbations
####### 

graphRandom <- random(graph=graph)
proc_leiden <- robinRobust(graph=graph, graphRandom=graphRandom, method="leiden",  resolution = 2, measure= "vi")

plot(proc_leiden) + theme_bw()

proc_infomap <- robinRobust(graph=graph, graphRandom=graphRandom, 
                            method="infomap", measure ="vi")

plot(proc_infomap)+ theme_bw()

# Interval testing procedure estimating significant differences 
#between curves at each step with multiple corrections
robinFDATest(proc_leiden) 
robinFDATest(proc_infomap)

# Gaussian process test estimating Bayes Factor
robinGPTest(proc_leiden)
robinGPTest(proc_infomap)

# Area under curve against null network clusters
robinAUC(proc_leiden)
robinAUC(proc_infomap)

# Robustness (VI = Measure) of 
comp <- robinCompare(graph=graph, method1="leiden", method2="infomap", 
                     #args1=list(objective_function="modularity"), 
                     args2 = list(nb.trials=100))

# Test directly for significant differences in robustness to perturbations between algorithm
png("../Figures/Robin_FDA.png", width = 17, height = 15, units = "cm", res = 600)
robinFDATest(comp)  
dev.off()

robinGPTest(comp)
robinAUC(comp)
