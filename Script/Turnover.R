library(tidyverse)
library(sf)
library(terra)
library(ggpubr)
library(colorspace)
library(patchwork)
library(cluster)
library(ape)
library(betapart)
library(oce)
library(epm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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


# Ecoregions
Ecoregions_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/MEOW/", layer = "meow_ecos")%>% 
  st_transform(st_crs("EPSG:3995")) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))

# Coastline shapefiles
Coastline <- read_sf("../../../Feb2023_Transfer/Oceanographic/Data/GSHHS_shp/l/", layer = "GSHHS_l_L1") %>% 
  st_make_valid() %>% 
  st_transform(st_crs("EPSG:3995")) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))  

# Phylogenetic tree (Rabosky et al. 2018; Chang et al. 2019)
tree <- read.tree("../../../Shelf_Evolution/Data/Phylos/actinopt_12k_raxml.tre")

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
  filter(grepl("_", taxon)) %>% 
  # filter for shallow depths
  filter(taxon %in% tree$tip.label) 

# Spatial polygons for points (as occurrences)
sf_points <- st_as_sf(trawls, 
                      coords = c("longitude", "latitude"), 
                      crs = 4326) %>% 
  st_transform(3995)

###########################################################################

index <- which(lengths(st_intersects(st_make_valid(Ecoregions_shelf_seas), st_make_valid(sf_points)))>0)
Ecoregions_shelf_seas <- Ecoregions_shelf_seas[index,]

##########################################################################
# Phylo dissimilarity in lineages
##########################################################################

# Construct species x grid matrix

fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}

grid <- fishnet(sf_points, square = F, cellsize = c(200000, 200000))

trawl_grid <- grid %>% 
  st_as_sf() %>% 
  mutate(grid_id = row_number())

gridded_trawls <- st_join(sf_points, trawl_grid) %>% 
  #filter(grid_ID != 14 | grid_ID != 8) %>% 
  group_by(haul_id) %>% 
  mutate(Haul_ID_num = cur_group_id()) %>% 
  ungroup() %>% 
  group_by(grid_id) %>% 
  #count(Haul_ID_num) 
  mutate(unique = n_distinct(Haul_ID_num)) %>% 
  arrange(unique) %>% 
  ungroup() %>% 
  filter(grid_id != "NA")

######################
# Phylogenetic tree
######################
# matching underscore for tree
#gridded_trawls$taxon <- gsub(" ", "_", gridded_trawls$taxon)

# shared species between data sources?
length(intersect(gridded_trawls$taxon, (tree$tip.label)))/length(unique(gridded_trawls$taxon))
# 1161/ 1704

species_grid_i <- gridded_trawls %>% 
  st_drop_geometry() %>% 
  #distinct(taxon, .keep_all = T) %>% 
  select(taxon, grid_id) %>% 
  mutate(Presence = 1) %>% 
  # Summarize duplicates
  pivot_wider(names_from = taxon, values_from = Presence,
              values_fn = function(x)sum(x)) %>% 
  replace(is.na(.),0) %>% 
  mutate(across(!grid_id, ~ ifelse(. > 0, 1, 0))) %>% 
  column_to_rownames(var = "grid_id")

pruned.tree <- keep.tip(tree, tree$tip.label[which(tree$tip.label %in% colnames(species_grid_i))])

phylo_diss <- phylo.beta.pair(species_grid_i, pruned.tree, index.family="jaccard")

# Pairwise turnover method
phylo_diss_jtu <- data.frame(as.matrix(phylo_diss[[1]])) %>% 
  rownames_to_column(var = "grid_id") %>% 
  mutate(grid_id = as.numeric(grid_id)) %>% 
  pivot_longer(!grid_id,
               names_to = "grid_pair",
               values_to = "phylo.beta.jtu") %>% 
  filter(phylo.beta.jtu != 0) %>% 
  group_by(grid_id) %>% 
  summarize(SD_jtu = sd(phylo.beta.jtu),
            jtu = mean(phylo.beta.jtu)) %>% 
  ungroup()

phylo_grid <- trawl_grid %>% 
  full_join(phylo_diss_jtu) %>% 
  distinct(grid_id,.keep_all = T) %>% 
  st_as_sf() %>% 
  st_make_valid() %>% 
  st_join(Ecoregions_shelf_seas)

#sf::write_sf(dis_all_df, "../Data/Dissimilairty_outputs/Trait_Phylo_Diss_200km.shp")

###################################
#dis_all_df <- read_sf("../Data/Dissimilairty_outputs/", layer="Trait_Phylo_Diss_200km") 

phylo_grid %>% 
  ggplot()+
  geom_sf(aes(fill = jtu))+
  theme_minimal(12)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme(axis.text = element_text(size=5),
        legend.key.width = unit(0.5, "cm"))+
  coord_sf(label_axes = "",
           expand = F,
           clip = "off")+
  scale_fill_continuous(type = "viridis", 
                        limits = c(0.2,0.9), 
                        breaks = seq(0.2,1,0.3),
                        name = "Turnover")


# Focal radius method (Baselga et al. 2013; Title et al. 2022)

EPM <- createEPMgrid(sf_points, 
                     resolution = 200000, 
                     cellType = 'hex',
                     retainSmallRanges = T,
                     crs(sf_points))

# merge
EPM_phylo <- reduceToCommonTaxa(addPhylo(EPM, tree, verbose = F))


Turnover_phylo <- epm::betadiv_phylogenetic(EPM_phylo, radius = 500000, component = "turnover")

#
ggplot(Turnover_phylo)+
  #geom_sf(aes(color = LME_NAME), fill = NA)+
  geom_sf(data= Turnover_phylo, aes(fill = phylogenetic_turnover))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(16)+
  scale_fill_gradientn(colors = oce.colorsTemperature(20),
                       limits=c(0, 0.6),
                       breaks=seq(0,0.6,0.2),
                       name="Turnover")+
  coord_sf(label_axes = "-N",
           #crs = st_crs(EPM_phylo_trait),
           expand = F,
           clip = "off")+
  guides(color = "none")
ggsave("../Figures/Turnover_200km_500km.png", width = 17, height = 15, units = "cm", dpi = 600)

EPM_100 <- createEPMgrid(sf_points, 
                     resolution = 100000, 
                     cellType = 'hex',
                     retainSmallRanges = T,
                     crs(sf_points))

# merge
EPM_phylo_100 <- reduceToCommonTaxa(addPhylo(EPM_100, tree, verbose = F))


Turnover_phylo_100 <- epm::betadiv_phylogenetic(EPM_phylo_100, radius = 250000, component = "turnover")

#
ggplot(Turnover_phylo_100)+
  #geom_sf(aes(color = LME_NAME), fill = NA)+
  geom_sf(aes(fill = phylogenetic_turnover))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(16)+
  scale_fill_gradientn(colors = oce.colorsTemperature(20),
                       limits=c(0, 0.6),
                       breaks=seq(0,0.6,0.2),
                       name="Turnover")+
  coord_sf(label_axes = "-N",
           #crs = st_crs(EPM_phylo_trait),
           expand = F,
           clip = "off")+
  guides(color = "none")
ggsave("../Figures/Turnover_100km_500km.png", width = 17, height = 15, units = "cm", dpi = 600)

sf::write_sf(st_join(Turnover_phylo, trawl_grid), "../Data/Turnover_200km_complete_500km.shp")
sf::write_sf(st_join(Turnover_phylo_100, trawl_grid), "../Data/Turnover_100km_complete_500km.shp")

##########################################################################
# If bootstrapping the random subsample process
##########################################################################
phylo_diss_metrics <- list()  
boots <- list()
boot <- 99

for(i in 1:boot){
  
  # species x grid matrix
  species_grid_i <- gridded_trawls %>% 
    st_drop_geometry() %>% 
    filter(taxon %in% overall.pruned.tree$tip.label) %>% 
    filter(unique > 50) %>% 
    # Random subsample stratified
    slice(sample(1:n())) %>% 
    # Random subsample stratified?
    group_by(haul_id) %>% 
    mutate(Group = cur_group_id()) %>% 
    ungroup() 
  
  # Index by all rows which match randomized haul ID 
  Match <- species_grid_i$Group[which(species_grid_i$Group %in% unique(species_grid_i$Group)[1:133])]
  
  #Final for matching against global trait space
  species_grid_i <- species_grid_i %>% 
    filter(Group %in% Match) %>%  
    distinct(taxon, .keep_all = T) %>% 
    select(taxon, grid_id) %>% 
    #distinct(taxon, .keep_all = T) %>% 
    mutate(Presence = 1) %>% 
    # Summarize duplicates
    pivot_wider(names_from = taxon, values_from = Presence,
                values_fn = function(x)sum(x)) %>% 
    replace(is.na(.),0) %>% 
    mutate(across(!grid_id, ~ ifelse(. > 0, 1, 0))) %>% 
    column_to_rownames(var = "grid_id")
  
  pruned.tree <- keep.tip(tree, tree$tip.label[which(tree$tip.label %in% colnames(species_grid_i))])
  
  phylo_diss <- phylo.beta.pair(species_grid_i, pruned.tree, index.family="jaccard")
  
  phylo_diss_jtu <- data.frame(as.matrix(phylo_diss[[1]])) %>% 
    rownames_to_column(var = "grid_id") %>% 
    mutate(grid_id = as.numeric(grid_id)) %>% 
    pivot_longer(!grid_id,
                 names_to = "grid_pair",
                 values_to = "phylo.beta.jtu") %>% 
    mutate(bootstrap = i) %>% 
    filter(phylo.beta.jtu != 0)
  
  phylo_diss_jne <- data.frame(as.matrix(phylo_diss[[2]])) %>% 
    rownames_to_column(var = "grid_id") %>% 
    mutate(grid_id = as.numeric(grid_id)) %>% 
    pivot_longer(!grid_id,
                 names_to = "grid_pair",
                 values_to = "phylo.beta.jne")%>% 
    mutate(bootstrap = i) %>% 
    filter(phylo.beta.jne != 0)
  
  phylo_diss_jac <- data.frame(as.matrix(phylo_diss[[3]])) %>% 
    rownames_to_column(var = "grid_id") %>% 
    mutate(grid_id = as.numeric(grid_id)) %>% 
    pivot_longer(!grid_id,
                 names_to = "grid_pair",
                 values_to = "phylo.beta.jac")%>% 
    mutate(bootstrap = i) %>% 
    filter(phylo.beta.jac != 0)
  
  phylo_diss_comb <- full_join(phylo_diss_jtu, full_join(phylo_diss_jne, phylo_diss_jac))
  
  phylo_diss_metrics[[i]] <- phylo_diss_comb
  boots[[i]] <- phylo_diss_metrics
  
}

# need to remove big files above to run this, otherwise memory exhausted
#rm(trawls_traits, gridded_trawls)

phylo_diss_df <- do.call(bind_rows, boots)


#write.csv(phylo_diss_df, "../Data/PCoA_outputs/Turnover_1stQ_200km_6axes_Cont_boot.csv")

phylo_diss_df <- read_csv("../Data/PCoA_outputs/Turnover_1stQ_200km_6axes_Cont_boot.csv")

phylo_diss_df <- phylo_diss_df %>% 
  group_by(grid_id, grid_pair) %>% 
  # median of bootstrap variation
  summarise(phylo.beta.jtu_SD = sd(phylo.beta.jtu),
            phylo.beta.jne_SD = sd(phylo.beta.jne),
            phylo.beta.jac_SD = sd(phylo.beta.jac),
            phylo.beta.jtu = median(phylo.beta.jtu),
            phylo.beta.jne = median(phylo.beta.jne),
            phylo.beta.jac = median(phylo.beta.jac)) %>% 
  ungroup()
