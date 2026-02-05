setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(scales)
library(patchwork)
library(ggalluvial)

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

# Remove Antilles
Coastline_crop <- Coastline %>% 
  st_crop(st_bbox(Ecoregions_shelf_seas))

# Find minimum trawl numbers to subsample down
fishnet <- function(geometry, ...){
  
  grid <- st_make_grid(geometry,...)
  index <- which(lengths(st_intersects(grid, geometry))>0)
  grid[index]
}

# grid size either 200000 or 100000
grid <- fishnet(sf_points, square = F, cellsize = c(200000, 200000))

#plot(st_geometry(sf_points), add = T, col ="red")

trawl_grid <- grid %>% 
  st_as_sf() %>% 
  mutate(grid_id = row_number()) 

# Load existing bioregionalizations with PC values for plotting
map_clusts <- read_sf("../Data/", layer = 'Rarefied_Leiden_1stquartile_Bioregions_200km')

##################################################################################
# Analyzing the community features of each bioregion
##################################################################################

#Unique species per bioregion?

clusts <- map_clusts %>% 
  # Step 1: Calculate average value per group and arrange
  group_by(cluster) %>%
  summarise(avg_val = median(Latitud), .groups = "drop") %>%
  arrange(desc(avg_val)) %>%
  # Step 2: Assign new group names based on order
  mutate(ID = paste0("Group", row_number())) %>%
  st_drop_geometry()  

map_clusts <- map_clusts %>% 
  left_join(clusts)

classes <- map_clusts %>% 
  group_by(ID, class, genus) %>% 
  summarise(freq = n()) %>% 
  mutate(prop = freq/sum(freq)) %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern /nPacific",
                                     ID == "Group2"~"North & n/Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group4"~"Outer /nEuropean Shel",
                                     ID == "Group5"~"(Sub-)Arctic /nAtlantic",
                                     ID == "Group6"~"Temperate /nPacific",
                                     ID == "Group7"~"Northeast /nAtlantic Shelf",
                                     ID == "Group8"~"Southeast n/US Shelf"

  ),
  levels = c("Northern /nPacific",
             "North & n/Celtic Seas",
             "Baltic Sea",
             "Outer /nEuropean Shel",
             "(Sub-)Arctic /nAtlantic",
             "Temperate /nPacific",
             "Northeast /nAtlantic Shelf",
             "Southeast n/US Shelf"
             
  )))

# Vizualize some facets of taxonomy 
ggplot(classes, aes(x = bioregions, y = prop, fill=class))+
  geom_bar(position = "stack", stat = "identity")+
  theme(axis.text.x = element_text(angle=45, hjust=1))


classes %>% 
  group_by(ID) %>% 
  summarize() %>%
  ungroup() %>% 
  mutate(area = st_area(.))


mapped_grid <- fishnet(map_clusts, 
                       square = F, 
                       cellsize = c(200000,
                                    200000))

###################################################################
# Calculating area
###################################################################
grid_bioregions <- st_join(trawl_grid, map_clusts) %>% 
  group_by(ID) %>% 
  summarize()

colors <- c("#FF7F00", "#FDBF6F","#B2DF8A", "#FB9A99","#33A02C",
             "#E31A1C","#A6CEE3","#1F78B4")

# Arrange by latitude
grid_bioregions %>% 
  filter(ID != "NA") %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                     ID == "Group2"~"North & Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group5"~"Outer European Shelf",
                                     ID == "Group4"~"(Sub-)Arctic Atlantic",
                                     ID == "Group6"~"Temperate Pacific",
                                     ID == "Group7"~"Northeast Atlantic Shelf",
                                     ID == "Group8"~"Southeast US Shelf"
                                     
  ),
  levels = c("Northern Pacific",
             "North & Celtic Seas",
             "Baltic Sea",
             "(Sub-)Arctic Atlantic",
             "Outer European Shelf",
             "Temperate Pacific",
             "Northeast Atlantic Shelf",
             "Southeast US Shelf"
             
  ))) %>% 
  ggplot()+
  geom_sf(aes(fill = bioregions))+
  #geom_sf(data = mapped_grid, fill = NA, color = "black", linewidth=0.05)+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(16)+
  theme(
    plot.margin=grid::unit(c(0,0,0,0), "mm")
    #legend.position = c(0.2, 0.3)
    #legend.position = "none"
  )+
  coord_sf(label_axes = "-N",
           expand = F,
           clip = "off")+
  scale_fill_manual(values=rev(colors),
                    name = "Bioregion")+
  #scale_color_brewer(palette = "Spectral")+
  guides(color = guide_legend(override.aes = list(size=5)))

ggsave("../Figures/R1_Bioregions.png", width = 17, height = 15, units = "cm", dpi = 600)

# how many grid cells per?
map_clusts %>% group_by(cluster) %>% summarize(grid_id = sum(n_distinct(grid_id)))


Polygons <- grid_bioregions %>% 
  filter(ID != "NA") %>% 
  group_by(ID) %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_area(.) %>% 
  as.data.frame() %>% 
  mutate(ID = paste0("Group", seq(1, 8, 1))) %>% 
  left_join(grid_bioregions) %>% 
  st_as_sf() %>% 
  rename("Area" = ".",
         geometry = "x") %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                     ID == "Group2"~"North & Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group4"~"Outer European Shelf",
                                     ID == "Group5"~"(Sub-)Arctic Atlantic",
                                     ID == "Group6"~"Temperate Pacific",
                                     ID == "Group7"~"Northeast Atlantic Shelf",
                                     ID == "Group8"~"Southeast US Shelf"
                                     
  ),
  levels = c("Northern Pacific",
             "North & Celtic Seas",
             "Baltic Sea",
             "Outer European Shelf",
             "(Sub-)Arctic Atlantic",
             "Temperate Pacific",
             "Northeast Atlantic Shelf",
             "Southeast US Shelf"
             
  ))) %>% 
  arrange(desc(Area))

Polygons

###################################################################
# Bioregion Richness
###################################################################
map_clusts %>% 
  group_by(cluster) %>% 
  summarise(count = n_distinct(taxon)) %>% 
  ggplot(aes(x = cluster, y = count))+
  geom_col()+
  theme_bw(16)+
  labs(x = "", y = "Count")

map_clusts %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                                    ID == "Group2"~"North & Celtic Seas",
                                                    ID == "Group3"~"Baltic Sea",
                                                    ID == "Group5"~"Outer European Shelf",
                                                    ID == "Group4"~"(Sub-)Arctic Atlantic",
                                                    ID == "Group6"~"Temperate Pacific",
                                                    ID == "Group7"~"Northeast Atlantic Shelf",
                                                    ID == "Group8"~"Southeast US Shelf"),
                                        levels = c("Northern Pacific",
                                                   "North & Celtic Seas",
                                                   "Baltic Sea",
                                                   "Outer European Shelf",
                                                   "(Sub-)Arctic Atlantic",
                                                   "Temperate Pacific",
                                                   "Northeast Atlantic Shelf",
                                                   "Southeast US Shelf"))) %>% 
  group_by(bioregions) %>% 
  summarize(richness = sum(n_distinct(taxon))) %>% 
  ungroup()

###################################################################
# Endemic species count 
###################################################################
map_clusts %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                     ID == "Group2"~"North & Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group4"~"(Sub-)Arctic Atlantic",
                                     ID == "Group5"~"Outer European Shelf",
                                     ID == "Group6"~"Temperate Pacific",
                                     ID == "Group7"~"Northeast Atlantic Shelf",
                                     ID == "Group8"~"Southeast US Shelf"),
                           levels = c("Northern Pacific",
                                      "North & Celtic Seas",
                                      "Baltic Sea",
                                      "Outer European Shelf",
                                      "(Sub-)Arctic Atlantic",
                                      "Temperate Pacific",
                                      "Northeast Atlantic Shelf",
                                      "Southeast US Shelf"))) %>% 
  distinct(taxon, bioregions) %>%                  
  group_by(taxon) %>%
  filter(n_distinct(bioregions) == 1) %>%         # Keep only names found in 1 group
  ungroup() %>%
  count(bioregions, taxon = "unique_name_count")  


###################################################################
# Depth range 
###################################################################
map_clusts %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                     ID == "Group2"~"North & Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group4"~"(Sub-)Arctic Atlantic",
                                     ID == "Group5"~"Outer European Shelf",
                                     ID == "Group6"~"Temperate Pacific",
                                     ID == "Group7"~"Northeast Atlantic Shelf",
                                     ID == "Group8"~"Southeast US Shelf"),
                           levels = c("Northern Pacific",
                                      "North & Celtic Seas",
                                      "Baltic Sea",
                                      "Outer European Shelf",
                                      "(Sub-)Arctic Atlantic",
                                      "Temperate Pacific",
                                      "Northeast Atlantic Shelf",
                                      "Southeast US Shelf"))) %>% 
  group_by(bioregions) %>%
  summarize(min = min(depth, na.rm = T),
          max = max(depth, na.rm = T), 
          mean = mean(depth, na.rm = T)) %>% 
  ungroup()

###################################################################
# Mean PC values  
###################################################################
map_clusts %>% 
  mutate(bioregions=factor(case_when(ID == "Group1"~"Northern Pacific",
                                     ID == "Group2"~"North & Celtic Seas",
                                     ID == "Group3"~"Baltic Sea",
                                     ID == "Group4"~"(Sub-)Arctic Atlantic",
                                     ID == "Group5"~"Outer European Shelf",
                                     ID == "Group6"~"Temperate Pacific",
                                     ID == "Group7"~"Northeast Atlantic Shelf",
                                     ID == "Group8"~"Southeast US Shelf"),
                           levels = c("Northern Pacific",
                                      "North & Celtic Seas",
                                      "Baltic Sea",
                                      "Outer European Shelf",
                                      "(Sub-)Arctic Atlantic",
                                      "Temperate Pacific",
                                      "Northeast Atlantic Shelf",
                                      "Southeast US Shelf"))) %>% 
  group_by(bioregions) %>%
  summarize(PC = mean(PC)) %>% 
  ungroup()

###################################################################
# Alluvial diagram
###################################################################

wide_data <- map_clusts %>%
  st_drop_geometry() %>% 
  select(ID, taxon) %>% 
  group_by(taxon, ID) %>% 
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = ID, values_from = present, values_fill = 0) %>% 
  ungroup()

long_data <- wide_data %>%
  pivot_longer(cols = starts_with("Group"), names_to = "group", values_to = "present") %>% 
  filter(present >0)

library(widyr)
# count the number of times two letters appear together
pairs <- pairwise_count(long_data, group, taxon) %>% 
  mutate(from=factor(case_when(item1 == "Group1" ~"Northern Pacific",
                                     item1 == "Group2"~"North & \nCeltic Seas",
                                     item1 == "Group3"~"Baltic Sea",
                                     item1 == "Group4"~"(Sub-)Arctic \nAtlantic",
                                     item1 == "Group5"~"Outer European \nShelf",
                                     item1 == "Group6"~"Temperate Pacific",
                                     item1 == "Group7"~"Northeast \nAtlantic Shelf",
                                     item1 == "Group8"~"Southeast US Shelf"
                                     
  ),
  levels = c("Northern Pacific",
             "North & \nCeltic Seas",
             "Baltic Sea",
             "(Sub-)Arctic \nAtlantic",
             "Outer European \nShelf",
             "Temperate Pacific",
             "Northeast \nAtlantic Shelf",
             "Southeast US Shelf"
             
  ))) %>% 
  mutate(to=factor(case_when(item2 == "Group1" ~"Northern Pacific",
                                     item2 == "Group2"~"North & \nCeltic Seas",
                                     item2 == "Group3"~"Baltic Sea",
                             item2 == "Group4"~"(Sub-)Arctic \nAtlantic",
                             item2 == "Group5"~"Outer European \nShelf",
                                     item2 == "Group6"~"Temperate Pacific",
                                     item2 == "Group7"~"Northeast \nAtlantic Shelf",
                                     item2 == "Group8"~"Southeast US Shelf"
                                     
  ),
  levels = c("Northern Pacific",
             "North & \nCeltic Seas",
             "Baltic Sea",
             "(Sub-)Arctic \nAtlantic",
             "Outer European \nShelf",
             "Temperate Pacific",
             "Northeast \nAtlantic Shelf",
             "Southeast US Shelf"
             
  )))
  

colors <- c("#FF7F00",
            "#FDBF6F",
            "#B2DF8A",
            "#FB9A99",
            "#33A02C",
            "#E31A1C",
            "#A6CEE3", 
            #"#999933",
            "#1F78B4")

ggplot(pairs,
       aes(y = n, axis1 = from, axis2 = to)) +
  geom_alluvium(width = 1/3 , aes(fill = from)) +
  geom_stratum( width = 1/3 , color = NA, aes(fill = from)) +
  geom_stratum( width = 1/3 , color = NA, aes(fill = to)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("", ""), expand = c(0, 0)) +
  theme_minimal(14)+
  theme(legend.position = "none", axis.title.x = element_text(vjust = 5))+
  labs(x = "Bioregion", y ="Frequency")+
  scale_fill_manual(values=rev(colors),  name = "",na.translate = FALSE)
ggsave("../Figures/Alluvial_revision.png", width = 17, height = 15, units = "cm", dpi = 600)  



