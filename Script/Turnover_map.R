library(tidyverse)
library(sf)
library(terra)
library(ggpubr)
library(colorspace)
library(patchwork)
library(cluster)
library(ape)
library(oce)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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

#Coastline <- read_sf("../../Feb2023_Transfer/Oceanographic/Data/tl_2019_us_coastline/", 
#                     layer = "tl_2019_us_coastline") %>% 
#  st_make_valid() %>% 
#  st_transform(st_crs("EPSG:3005")) %>% 
#  st_crop(st_bbox(LMEs_shelf_seas)) 

Ecoregions_shelf_seas <- read_sf("../../../Shelf_Evolution/Data/MEOW/", layer = "meow_ecos")%>% 
  st_transform(st_crs("EPSG:3995")) %>% 
  st_crop(st_bbox(LMEs_shelf_seas))

# Loading both turnover scales
two <- read_sf("../Data/", layer = "Turnover_200km_complete_500km")
one <- read_sf("../Data/", layer = "Turnover_100km_complete_500km")

x <- ggplot(one)+
  #geom_sf(aes(color = LME_NAME), fill = NA)+
  geom_sf(data= one, aes(fill = phylgn_))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(8)+
  scale_fill_gradientn(colors = oce.colorsTemperature(20),
                       limits=c(0, 0.6),
                       breaks=seq(0,0.6,0.2),
                       name="Turnover")+
  coord_sf(label_axes = "-N",
           #crs = st_crs(EPM_phylo_trait),
           expand = F,
           clip = "off")+
  guides(color = "none")

#
y <- ggplot(two)+
  #geom_sf(aes(color = LME_NAME), fill = NA)+
  geom_sf(aes(fill = phylgn_))+
  geom_sf(data = Coastline, fill = "grey70", color = "black", linewidth=0.05)+
  theme_minimal(8)+
  scale_fill_gradientn(colors = oce.colorsTemperature(20),
                       limits=c(0, 0.6),
                       breaks=seq(0,0.6,0.2),
                       name="Turnover")+
  coord_sf(label_axes = "-N",
           #crs = st_crs(EPM_phylo_trait),
           expand = F,
           clip = "off")+
  guides(color = "none")
 
x+ y + plot_layout(nrow = 1)
ggsave("../Figures/Turnover.png", width = 17, height = 15, units = "cm", dpi = 600) 
