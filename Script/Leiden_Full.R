setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(scales)
library(patchwork)
library(brainGraph)

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
  # After revision, we exclude extreme depths (> 1500 m) or missing depth values, except survey in shallow
  # Gulf of St. Lawrence (south)
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

# Find minimum trawl numbers to subsample down
fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}

# grid size either 200000 or 100000
grid <- fishnet(sf_points, square = F, cellsize = c(100000, 100000))

#plot(st_geometry(sf_points), add = T, col ="red")

trawl_grid <- grid %>% 
  st_as_sf() %>% 
  mutate(grid_id = row_number()) 


# grid cell rarefy determined in the hclust bioregion script 
sums <- st_join(sf_points, trawl_grid) %>% 
  #filter(grid_ID != 14 | grid_ID != 8) %>% 
  group_by(haul_id) %>% 
  mutate(Haul_ID_num = cur_group_id()) %>% 
  ungroup() %>% 
  group_by(grid_id) %>% 
  #count(Haul_ID_num) 
  mutate(unique = n_distinct(Haul_ID_num)) %>% 
  ungroup() %>% 
  group_by(grid_id) %>% 
  #count(Haul_ID_num) 
  mutate(unique = n_distinct(Haul_ID_num)) %>% 
  arrange(unique) 

# What is the min number of trawls per grid cell?
count <- sums[!duplicated(sums$unique),]
summary(count$unique)

## 100km 
# 93 = 1st quartile, 193 = median including depths <20m 
## 200km 
# 111 = 1st quartile, 332 = median including depths <20m 


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
  filter(unique > 92) 

############
# Subsample each grid cell down to x
#############

gridcells <- unique(World_grid$grid_id)
grid_ls <- list()  

for(i in gridcells){
  Cell_rarefy <- filter(World_grid, grid_id == i) %>% 
    group_by(unique) %>% 
    sample_n(size = 93) %>%
    ungroup()
  grid_ls[[i]] <- Cell_rarefy
  
}

grid_rarefy <- do.call("rbind", grid_ls)  

# How much is kept after rarefaction?

nrow(grid_rarefy)/nrow(sf_points)


# Construct graph from edgelist
edgelist <- data.frame(as.character(grid_rarefy$grid_id), 
                       as.character(grid_rarefy$taxon),
                       stringsAsFactors = FALSE)
colnames(edgelist) <- c("grid_id", "taxon")

edgelist_u <- unique(edgelist)
nrow(edgelist_u)

# Weight by occurrence frequency
strength_links <- unlist(apply(edgelist_u, 1,
                               function(x)length(which(edgelist[,1]==x[1] & edgelist[,2]==x[2]))))


graph <- graph_from_edgelist(as.matrix(edgelist_u), directed = FALSE)
E(graph)$weight <- strength_links
V(graph)$type <- V(graph)$name %in% edgelist_u[,1]
V(graph)$shape <- ifelse(V(graph)$type,  "square","circle")


#100km / 200km
res <- "200km"

#1stquartile / median / full
cutoff <- "1stquartile"

# Leiden / Infomap
method <- "Leiden"

############################################################################
#  Leiden clustering
############################################################################

# Searching across the resolution parameter in Leiden
set.seed(Sys.time())
set.seed(1234)

gamma <- seq(0,10,0.5)

nc <- vector("numeric",length(gamma))

for (i in 1:length(gamma)){
  gc <- cluster_leiden(graph, objective_function = "modularity",
                       n_iterations = 20, resolution_parameter = gamma[i],
                       weights = strength_links)
  nc[i] <- length(gc)
}

png(paste0("../Figures/Networks/World_resolution_",cutoff,"_",res,".png"), width = 17, height = 15, units = "cm", res = 600)
plot(gamma,
     nc,
     xlab="Resolution Parameter",
     ylab="Network Clusters",
     col="black",
     bg="maroon",
     main="",
     pch=21,
     cex=1.5,
     cex.main=1.5)

abline(v = 2, lty=2)
dev.off()


############################################
# Final cluster with determined gamma
############################################

set.seed(1234)

# Leiden algorithm (Traag et al. 2019)
map <- cluster_leiden(graph,
                      objective_function = "modularity",
                      resolution_parameter = 2,
                      n_iterations = 100,
                      weights = strength_links
                      )
modularity(graph, membership(map))

# How many clusters (biogeographical partitions)?
sizes(map)


colors <- RColorBrewer::brewer.pal(9, "Spectral")


clusts <- na.omit(as.numeric(map$names))[unique(map$membership)]
bins <- list()

for(i in 1:length(clusts)){
  x <- data.frame(
    grid_id=na.omit(as.numeric(map$names)[which(map$membership==i)]),
    #grid_id=na.omit(as.numeric(map$names))[unique(map$membership==i)],
    #cluster=paste0("Cluster_", i))
    cluster=i)
  bins[[i]] <- x
  
}

mapped_clusters <- do.call("rbind", bins)


map_clusts <- grid_rarefy %>% 
  left_join(mapped_clusters) 

clusts <- grid_rarefy %>% 
  left_join(mapped_clusters) %>%
  # Step 1: Calculate average value per group and arrange
  group_by(cluster) %>%
  summarise(avg_val = median(Latitude), .groups = "drop") %>%
  arrange(desc(avg_val)) %>%
  # Step 2: Assign new group names based on order
  mutate(ID = paste0("Group", row_number())) %>%
  st_drop_geometry()  

map_clusts <- map_clusts %>% 
  left_join(clusts)


# If grid cell is desireable to plot
mapped_grid <- fishnet(map_clusts, 
                       square = F, 
                       cellsize = c(100000,
                                    100000))

colors <- c(brewer.pal(length(unique(map$membership)), "Spectral")
            #"navyblue",
            # "grey", 
            #"black"
)

#colors <- c("#FF7F00", "#FDBF6F","#B2DF8A", "#FB9A99","#33A02C",
 #           "#E31A1C","#A6CEE3","#1F78B4", "grey")

map_grid <- map_clusts %>% 
  mutate(cluster = gsub("_", " ", cluster)) %>% 
  ggplot()+
  geom_sf( aes(color = ID), size=0.01)+
  geom_sf(data=Ecoregions_shelf_seas,color="black", fill = NA, linewidth=0.25)+
  geom_sf(data = Coastline_crop, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(8)+
  theme(
    plot.margin=grid::unit(c(0,0,0,0), "mm"),
    #legend.position = c(0.2, 0.3)
    legend.position = "none"
  )+
  coord_sf(label_axes = "-N",
           expand = F,
           clip = "off")+
  scale_color_manual(values=rev(colors))+
  #scale_color_brewer(palette = "Spectral")+
  guides(color = guide_legend(override.aes = list(size=5)))
map_grid


#ggsave(paste0("../Figures/", method, "_",cutoff,"_",res,".png"), map_grid, width = 17, height = 15, dpi = 600, units = "cm")
#write_sf(map_clusts,paste0("../Data/",method, "_", cutoff, "_Bioregions_",res,".shp"))

# Reorder network vector for color matching (latitude order)
map$cluster <- as.numeric(gsub( "Group", "", 
                                clusts$ID[match(map$membership,clusts$cluster)]))

png(paste0("../Figures/network_",method,"_", cutoff,"_",res,".png"), width=17 , height=15, units = "cm", res=600)
par(bg=NA)
plot((graph), 
     vertex.color = rev(colors)[map$cluster], 
     vertex.label=NA, 
     edge.border = NA,
     vertex.size=2,
     layout=layout_with_lgl(graph))
dev.off()

##################################################################################
# Participation coefficients (PC)
##################################################################################

## Square or not??
# Roger Guimerà & Luís A. Nunes Amaral (2005) says yes
# Bloomfield et al. (2018) says no

# 0.625 are critical apparrently 

selection <- which(is.na((as.numeric(V(graph)$name)))==FALSE)
output <- (part_coeff(graph, membership(map), A = NULL, weighted = T))
output_cells <- output[selection]

# Checking part_coeff source
part_coeff <- function(g, memb, A=NULL, weighted=FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse=FALSE, names=FALSE, attr='weight')
    } else {
      A <- as_adj(g, sparse=FALSE, names=FALSE)
    }
  }
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  1 - ((1 / Ki^2) * rowSums(Kis^2))
}

mapped_clusters$PC <- output_cells

map_clusts_PC <- grid_rarefy %>% 
  left_join(mapped_clusters) %>% 
  arrange(grid_id)

##############################
# Loading bioregionalization with PC values for plotting
map_clusts_PC <- read_sf("../Data/", layer = 'Leiden_1stquartile_Bioregions_200km')

mapped_grid <- fishnet(map_clusts_PC, 
                       square = F, 
                       cellsize = c(200000,
                                    200000))

map_PC <- map_clusts_PC %>% 
  #left_join(mapped_clusters) %>% 
  arrange(grid_id) %>%
  mutate(cluster = gsub("_", " ", cluster)) %>% 
  ggplot()+
  geom_sf(aes(color = PC), size=0.01)+
  #geom_sf(data=cluster_polys,color="black", fill = NA, linewidth=0.1)+
  geom_sf(data=mapped_grid,color="black", fill = NA, linewidth=0.1)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(16)+
  theme(
    plot.margin=grid::unit(c(0,0,0,0), "mm"),
    #legend.position = c(0.2, 0.3)
  )+
  scale_color_steps2(limits = c(0,0.8),
                     breaks = seq(0,0.8, 0.1),
                     midpoint = 0.4,
                     high = muted("red"),
                     mid = "white",
                     low = muted("blue"),
                     name = "Partcipation\nCoefficient")+
  coord_sf(label_axes = "-N",
           expand = F,
           clip = "off")
map_PC
ggsave(paste0("../Figures/PC_", method,"_",cutoff,"_",res,".png"), width = 17, height = 15, dpi = 600, units = "cm")

map_clusts_PC %>% 
  group_by(cluster) %>%
  summarize(PC = mean(PC)) %>% 
  ungroup()

#ggsave(paste0("../Figures/CompleteCases_",cutoff,"_",res,".png"), width = 17, height = 15, dpi = 600, units = "cm")

write_sf(map_clusts_PC, paste0("../Data/",method, "_", cutoff,"_Bioregions_",res,".shp"))

mapped <- map_grid + map_PC + plot_layout(nrow = 1)
ggsave(paste0("../Figures/collage_", method,"_",cutoff,"_",res,".png"), mapped, width = 17, height = 15, dpi = 600, units = "cm")
