# make a map for Figure 1


# packages
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

# read in dat
dat <- readRDS("Data/FrogID_data.RDS")

# read in spatial dat
urban_areas <- st_read("Data/urban_areas/SUA_2016_AUST.shp")

# Now assign FrogID records to urban areas
urban_areas <- urban_areas %>%
  dplyr::filter(! SUA_NAME16 %in% c("Not in any Significant Urban Area (NSW)",
                                    "Not in any Significant Urban Area (OT)",
                                    "Not in any Significant Urban Area (ACT)",
                                    "Not in any Significant Urban Area (SA)",
                                    "Not in any Significant Urban Area (NT)",
                                    "Not in any Significant Urban Area (Tas.)",
                                    "Not in any Significant Urban Area (WA)",
                                    "Not in any Significant Urban Area (Qld)",
                                    "Not in any Significant Urban Area (Vic.)")) %>%
  mutate(SUA_NAME16=as.character(as.factor(SUA_NAME16))) %>%
  mutate(col.id=1:nrow(.))

aus_ecoregions <- st_read("Data/WWF_ecoregions/WWF_ecoregions.shp")

ggplot()+
  geom_sf(data=aus_ecoregions, aes(fill=as.character(BIOME)))+
  scale_fill_brewer(palette="Set2")+
  geom_sf(data=urban_areas, fill="black", color="black")+
  theme_bw()+
  theme(panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank())+
  guides(fill=FALSE)

ggsave("Figures/map.png", width=4.5, height=4, units="in")
ggsave("Figures/map.eps", width=4.5, height=4, units="in")
ggsave("Figures/map.svg", width=4.5, height=4, units="in")
