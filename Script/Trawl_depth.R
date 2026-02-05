setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("igraph")
library(tidyverse)
library(sf)
library(oce)
library(scales)
library(colorspace)
library(patchwork)

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
# FishGlob
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

sf_points <- st_as_sf(trawls, 
                      coords = c("longitude", "latitude"),
                      crs="EPSG:4326") %>% 
  st_transform(st_crs(Ecoregions_shelf_seas))

index <- which(lengths(st_intersects(st_make_valid(Ecoregions_shelf_seas), st_make_valid(sf_points)))>0)
Ecoregions_shelf_seas <- Ecoregions_shelf_seas[index,]

map <- sf_points %>% 
  st_join(Ecoregions_shelf_seas) %>% 
  # fix join and filter!
  mutate(ECOREGION = case_when(
    ECOREGION == "Gulf of Maine/Bay of Fundy" ~ "Maine/Fundy",
    ECOREGION == "Gulf of St. Lawrence - Eastern Scotian Shelf" ~ "Gulf St. Lawrence",
    ECOREGION == "Southern Grand Banks - South Newfoundland" ~ "Newfoundland-\nLabrador",
    ECOREGION == "Northern Grand Banks - Southern Labrador" ~ "Newfoundland-\nLabrador",
    ECOREGION == "Southern Norway" ~ "Norwegian Seas",
    ECOREGION == "Northern Norway and Finnmark" ~ "Norwegian Seas",
    ECOREGION == "North American Pacific Fijordland" ~ "Pacific Fijordland",
    ECOREGION == "North and East Barents Sea" ~ "NE. Barents Sea",
    ECOREGION == "Eastern Bering Sea" ~ "E. Bering Sea",
    ECOREGION == "Northern California" ~ "California Current",
    ECOREGION == "Southern California Bight" ~ "California Current",
    ECOREGION == "Northern Gulf of Mexico" ~ "N. Gulf of Mexico",
    ECOREGION == "Oregon, Washington, Vancouver Coast and Shelf" ~ "Vancouver Shelf",
    ECOREGION == "Puget Trough/Georgia Basin" ~ "Vancouver Shelf",
    ECOREGION == "South European Atlantic Shelf" ~ "S. European Shelf",
    .default = ECOREGION)) %>% 
  #filter(ECOREGION != "NA") %>% 
  mutate(ECOREGION = fct_reorder(factor(ECOREGION),Latitude, mean, .desc = T)) %>% 
  ggplot()+
  geom_sf(aes(color = Latitude), fill = NA, size=0.01)+
  geom_sf(data = Ecoregions_shelf_seas, fill = NA, color = "black", linewidth = 0.3)+
  geom_sf(data = Coastline, fill = "grey", color = "black", linewidth = 0.05)+
  theme_minimal(18)+
  #scale_fill_gradientn(colors = oce.colorsViridis(20))+
  coord_sf(label_axes = "",
           expand = F,
           clip = "off")+
  guides(color = "none")+
  scale_color_gradientn(colours = colorspace::diverge_hcl(20, rev=T),
                        name="Latitude")
#ggsave("../Figures/Depth_distributions.png", width = 17, height = 15, dpi = 600, units = "cm")

# Summarize depth variation
Summary <-sf_points %>% 
  group_by(survey, year) %>% 
  summarize(min_depth = min(depth),
            max_depth = max(depth),
    #SD_depth = sd(depth),
            depth = mean(depth),
            Latitude = mean(Latitude)) %>% 
  st_join(Ecoregions_shelf_seas) %>% 
  mutate(ECOREGION = case_when(
    ECOREGION == "Gulf of Maine/Bay of Fundy" ~ "Maine/Fundy",
    ECOREGION == "Gulf of St. Lawrence - Eastern Scotian Shelf" ~ "Gulf \nSt. Lawrence",
    ECOREGION == "Southern Grand Banks - South Newfoundland" ~ "Newfoundland-\nLabrador",
    ECOREGION == "Northern Grand Banks - Southern Labrador" ~ "Newfoundland-\nLabrador",
    ECOREGION == "Gulf of Alaska" ~ "Gulf of\n Alaska",
    ECOREGION == "Southern Norway" ~ "Norwegian \nSeas",
    ECOREGION == "Northern Norway and Finnmark" ~ "Norwegian \nSeas",
    ECOREGION == "North American Pacific Fijordland" ~ "Pacific \nFijordland",
    ECOREGION == "North and East Barents Sea" ~ "NE. \nBarents Sea",
    ECOREGION == "Eastern Bering Sea" ~ "E. \nBering Sea",
    ECOREGION == "Northern California" ~ "California \nCurrent",
    ECOREGION == "Southern California Bight" ~ "California \nCurrent",
    ECOREGION == "Northern Gulf of Mexico" ~ "N. Gulf \nof Mexico",
    ECOREGION == "Oregon, Washington, Vancouver Coast and Shelf" ~ "Vancouver \nShelf",
    ECOREGION == "Puget Trough/Georgia Basin" ~ "Vancouver \nShelf",
    ECOREGION == "South European Atlantic Shelf" ~ "S. European \nShelf",
    .default = ECOREGION)) %>% 
  filter(ECOREGION != "NA") %>% 
  mutate(ECOREGION = fct_reorder(factor(ECOREGION),Latitude, mean, .desc = T)) 
  
# Plot variation
S <- Summary %>% 
  ggplot(aes(x = year, y = depth, color = Latitude))+
  #stat_summary(fun.data = "mean_cl_boot",size=0.5)+
  #geom_pointrange(aes(ymin = depth-SD_depth, ymax = depth+SD_depth), size=0.001)+
  geom_pointrange(aes(ymin = min_depth, ymax = max_depth), size=0.01, fatten = 0)+
  #geom_point(shape = 1, colour = "black")+
  scale_color_gradientn(colours = colorspace::diverge_hcl(20, rev=T),
                       name="Latitude")+
  theme_minimal(6)+
  scale_y_continuous(limits = c(0,1500))+
  scale_y_reverse()+
  #geom_hline(yintercept = 200, color = "yellow", linewidth = 0.5)+
  #geom_hline(yintercept = 1500, color = "black", linewidth = 0.5)+
  labs(x = "Year", y = "Depth (m)")+
  facet_wrap(~ECOREGION)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        strip.text = element_text(size=5))
#ggsave("../Figures/Depth_distributions.png", width = 17, height = 15, dpi = 600, units = "cm")

p <- map + S + plot_layout(nrow = 1)
ggsave("../Figures/Depth_distributions.png", p, width = 21, height = 17, dpi = 600, units = "cm")
