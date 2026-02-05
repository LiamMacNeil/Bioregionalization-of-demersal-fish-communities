setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(scales)
library(patchwork)
library(RColorBrewer)

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

# Remove Antilles
Coastline_crop <- Coastline %>% 
  st_crop(st_bbox(Ecoregions_shelf_seas))

# Load in bioregionalizations with participation coefficients
map_clusts_PC <- read_sf("../Data/", layer = 'Leiden_1stquartile_Bioregions_200km')

#####################
# Loading and processing trawl data
#####################

index <- which(lengths(st_intersects(st_make_valid(Ecoregions_shelf_seas), st_make_valid(map_clusts_PC)))>0)
Ecoregions_shelf_seas <- Ecoregions_shelf_seas[index,]

# Find minimum trawl numbers to subsample down
fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}

# grid size either 200000 or 100000
grid <- fishnet(map_clusts_PC, square = F, cellsize = c(200000, 200000))


trawl_grid <- grid %>% 
  st_as_sf() %>% 
  mutate(grid_id = row_number()) 


clusts <- map_clusts_PC %>% 
  # Step 1: Calculate average value per group and arrange
  group_by(cluster) %>%
  summarise(avg_val = median(Latitud), .groups = "drop") %>%
  arrange(desc(avg_val)) %>%
  # Step 2: Assign new group names based on order
  mutate(ID = paste0("Group", row_number())) %>%
  st_drop_geometry()  

map_clusts_PC <- map_clusts_PC %>% 
  left_join(clusts)


# If grid cell is desireable to plot
mapped_grid <- fishnet(map_clusts_PC, 
                       square = F, 
                       cellsize = c(200000,
                                    200000))


colors <- c(brewer.pal(length(unique(map_clusts_PC$ID)), "Spectral"),
            #, "grey", 
            "black"
)

map_grid <- map_clusts_PC %>% 
  mutate(cluster = gsub("_", " ", cluster)) %>% 
  ggplot()+
  geom_sf( aes(color = ID), size=0.01)+
  geom_sf(data=Ecoregions_shelf_seas,color="black", fill = NA, linewidth=0.25)+
  geom_sf(data = Coastline_crop, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(16)+
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


ggsave(paste0("../Figures/", method, "_",cutoff,"_",res,".png"), map_grid, width = 17, height = 15, dpi = 600, units = "cm")
#write_sf(map_clusts,paste0("../Data/",method, "_", cutoff, "_Bioregions_",res,".shp"))

# Reorder network vector for color matching (latitude order)
map$cluster <- as.numeric(gsub( "Group", "", 
                                clusts$ID[match(map$membership,clusts$cluster)]))

png(paste0("../Figures/network_",method,"_", cutoff,"_",res,".png"), width=17 , heigh=15, units = "cm", res=600)
par(bg=NA)
plot((graph), 
     vertex.color = rev(colors)[map$cluster], 
     vertex.label=NA, 
     edge.border = NA,
     vertex.size=2,
     layout=layout_nicely(graph))
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

map_clusts_PC %>% 
  group_by(cluster) %>%
  summarize(PC = mean(PC)) %>% 
  ungroup()
ggsave(paste0("../Figures/PC_", method,"_",cutoff,"_",res,".png"), width = 17, height = 15, dpi = 600, units = "cm")

#map_grid + map_PC + plot_layout(nrow = 1)
#ggsave(paste0("../Figures/CompleteCases_",cutoff,"_",res,".png"), width = 17, height = 15, dpi = 600, units = "cm")

write_sf(map_clusts_PC, paste0("../Data/Rarefied_",method, "_", cutoff,"_Bioregions_",res,".shp"))

