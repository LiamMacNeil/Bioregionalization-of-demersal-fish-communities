setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(scales)
library(patchwork)
library(brainGraph)
library(RColorBrewer)
library(bipartite)

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

# Remove Antilles
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
# 93 = 1st quartile, 193 = median including depths <20m 
## 200km 
# 112 = 1st quartile, 333 = median including depths <20m 


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
# Nice!
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

# Construct graph fron edgelist
edgelist <- data.frame(as.character(grid_rarefy$grid_id), 
                       as.character(grid_rarefy$taxon),
                       stringsAsFactors = FALSE)
colnames(edgelist) <- c("grid_id", "taxon")

edgelist_u <- unique(edgelist)
nrow(edgelist_u)

# weights by occurrence frequency
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

method <- "Leiden"

############################################################################
# Clustering
############################################################################

set.seed(1234)

# Leiden algorithm (Traag et al. 2019)
map <- cluster_leiden(graph,
                      objective_function = "modularity",
                      resolution_parameter = 2,
                      n_iterations = 100,
                      weights = strength_links)
modularity(graph, membership(map), directed = F)

# How many clusters (biogeographical partitions)?
sizes(map)

colors <- RColorBrewer::brewer.pal(n=length(sizes(map)), name = "Paired")

#png(paste0("../Figures/Networks/Global_network_",cutoff,"_",res,".png"), width=17 , heigh=15, units = "cm", res=600)

plot(graph, 
     #vertex.color=membership(map),
     vertex.color = (rev(colors)[map$membership]), 
     vertex.label=NA, 
     vertex.size=4,
     layout=layout_with_fr(graph))
#dev.off()


##################################################################################
 # Species x site matrix for randomization in null models
web_matrix <- edgelist_u %>% 
  mutate(strength_links = strength_links) %>% 
  pivot_wider(names_from = grid_id, values_from = strength_links) %>%
  replace(is.na(.),0) %>% 
  column_to_rownames("taxon") %>% 
  as.matrix() 

# These can be SLOW   

# randomizes matrix but maintains 
# 1) Number of connections and links are constant (connectence)
# 2) and marginal totals
swapped <- swap.web(web_matrix, N = 999)


# randomizes matrix but maintains dimensionality
# 1) Number of connections and links are constant (connectence) but not link strength (marginal totals)
shuffled <- shuffle.web(web_matrix, N = 999, legacy = F)

# randomizes matrix but maintains dimensionality
# 1) Number of connections and links are constant (connectence)
vaz <- vaznull(web_matrix, N = 999)


##############################################

# Loop clustering calculations and store modularity to compare distribution against observed in Leiden
modularity_shuffled_scores <- list()
modularity_swapped_scores <- list()
modularity_vaz_scores <- list()


for(i in 1:length(shuffled)){
  
  ####################################
  shuffled_data <- shuffled[[i]] %>% 
    as.data.frame() %>% 
    setNames(colnames(web_matrix)) %>% 
    mutate(taxon = rownames(web_matrix)) %>% 
    pivot_longer(cols = c(everything(.), -taxon)) %>% 
    filter(value > 0) %>% 
    rename(grid_ID = name) %>% 
    mutate(grid_ID = as.numeric(grid_ID)) %>% 
    unique()
  
  edgelist_u_shuffled <- shuffled_data %>% 
    dplyr::select(taxon, grid_ID) %>% 
    arrange(grid_ID)
  
  strength_links_shuffled <- shuffled_data %>% 
    dplyr::select(value) %>% 
    unlist() %>% 
    as.vector()
  
  graph_shuffled <- graph_from_edgelist(as.matrix(edgelist_u_shuffled), directed = FALSE)
  E(graph_shuffled)$weight <- strength_links_shuffled
  
  map_shuffled <- cluster_leiden(graph_shuffled,
                        objective_function = "modularity",
                        resolution_parameter = 2,
                        n_iterations = 20,
                        weights = strength_links_shuffled)
  modularity_shuffled_scores[[i]] <- modularity(graph_shuffled, membership =membership(map_shuffled), directed = F)
  
  ####################################
  swapped_data <- swapped[[i]] %>% 
    as.data.frame() %>% 
    setNames(colnames(web_matrix)) %>% 
    mutate(taxon = rownames(web_matrix)) %>% 
    pivot_longer(cols = c(everything(.), -taxon)) %>% 
    filter(value > 0) %>% 
    rename(grid_ID = name) %>% 
    mutate(grid_ID = as.numeric(grid_ID)) %>% 
    unique()
  
  edgelist_u_swapped <- swapped_data %>% 
    dplyr::select(taxon, grid_ID) %>% 
    arrange(grid_ID)
  
  strength_links_swapped <- swapped_data %>% 
    dplyr::select(value) %>% 
    unlist() %>% 
    as.vector()
  
  graph_swapped <- graph_from_edgelist(as.matrix(edgelist_u_swapped), directed = FALSE)
  E(graph_swapped)$weight <- strength_links_swapped
  
  map_swapped <- cluster_leiden(graph_swapped,
                                 objective_function = "modularity",
                                 resolution_parameter = 2,
                                 n_iterations = 20,
                                 weights = strength_links_swapped)
  modularity_swapped_scores[[i]] <- modularity(graph_swapped, membership =membership(map_swapped), directed = F)
  
  ####################################
  vaz_data <- vaz[[i]] %>% 
    as.data.frame() %>% 
    setNames(colnames(web_matrix)) %>% 
    mutate(taxon = rownames(web_matrix)) %>% 
    pivot_longer(cols = c(everything(.), -taxon)) %>% 
    filter(value > 0) %>% 
    rename(grid_ID = name) %>% 
    mutate(grid_ID = as.numeric(grid_ID)) %>% 
    unique()
  
  edgelist_u_vaz <- vaz_data %>% 
    dplyr::select(taxon, grid_ID) %>% 
    arrange(grid_ID)
  
  strength_links_vaz <- vaz_data %>% 
    dplyr::select(value) %>% 
    unlist() %>% 
    as.vector()
  
  graph_vaz <- graph_from_edgelist(as.matrix(edgelist_u_vaz), directed = FALSE)
  E(graph_vaz)$weight <- strength_links_vaz
  
  map_vaz <- cluster_leiden(graph_vaz,
                                objective_function = "modularity",
                                resolution_parameter = 2,
                                n_iterations = 20,
                                weights = strength_links_vaz)
  modularity_vaz_scores[[i]] <- modularity(graph_vaz, membership =membership(map_vaz), directed = F)
  
}

# Bind all null types together from bootstraps
shuffled_mods <- data.frame(do.call("rbind", modularity_shuffled_scores))
swapped_mods <- data.frame(do.call("rbind", modularity_swapped_scores))
vaz_mods <- data.frame(do.call("rbind", modularity_vaz_scores))

Null_modularity <- shuffled_mods %>% 
  rename(Shuffled = do.call..rbind...modularity_shuffled_scores.) %>% 
  mutate(Swapped = swapped_mods$do.call..rbind...modularity_swapped_scores.) %>%
  mutate(Vaznull = vaz_mods$do.call..rbind...modularity_vaz_scores.) %>%
  pivot_longer(cols =c("Shuffled", "Swapped", "Vaznull")) %>% 
  rename(Method = name,
         Modularity = value)

#####################################################
# Comparing Null model modularity
#####################################################

ggplot(Null_modularity, aes(x = Modularity, fill=Method, color = Method))+
  geom_histogram(alpha = 0.5, bins = 1000)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  labs(y = "Frequency")+
  geom_vline(xintercept = modularity(graph, membership(map)), color='#D82632')+
  theme_bw(12)+
  scale_x_continuous(breaks = seq(0,0.7, 0.2))+
  xlim(0, 0.7)+
  coord_cartesian(xlim = c(0, 0.7))
ggsave("../Figures/NullModels_Modularity.png", height = 12, width=15,
       units = "cm", dpi = 600)


