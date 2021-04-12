## A script to assign every record the urban area
## and ecoregion for further analyses

# packages
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

# read in dat
dat <- readRDS("Data/FrogID_data.RDS")

# read in spatial dat
urban_areas <- st_read("Data/urban_areas/SUA_2016_AUST.shp")

aus_ecoregions <- st_read("Data/WWF_ecoregions/WWF_ecoregions.shp")

# Assign every FrogID record to an ecoregion
# create metadata df
# used to link the dfs back together later
# the values are associated based on the .html file linked above
# this is done manually
eco_metadata <- data.frame(aus_eco=aus_ecoregions[[8]],
                           col.id=1:375) %>%
  mutate(Ecoregion=case_when(
    aus_eco == 7 ~ "Tropical & Subtropical Grasslands",
    aus_eco == 13 ~ "Deserts & Xeric Shrublands",
    aus_eco == 1 ~ "Tropical & Subtropical Moist Broadleaf Forests",
    aus_eco == 4 ~ "Temperate Broadleaf & Mixed Forests",
    aus_eco == 12 ~ "Mediterranean Forests",
    aus_eco == 14 ~ "Mangroves",
    aus_eco == 8 ~ "Temperate Grasslands", 
    aus_eco == 10 ~ "Montane Grasslands & Shrublands"
  ))

# convert points from FrogID into sf points
# so that they can be used with the polygons
points_sf <- st_as_sf(dat, coords = c("lng", "lat"), 
                      crs = st_crs(aus_ecoregions), agr = "constant")

# use intersect to return a list of points and the associated
# polygon they belong to - uses integer ids for rows and columns
assigned_points <- as.data.frame(st_intersects(points_sf, aus_ecoregions))

# join the original dataframe with the intersected df
# then remove the row and col ids
dat2 <- dat %>%
  mutate(row.id = 1:nrow(.)) %>%
  left_join(assigned_points, by="row.id") %>%
  left_join(., eco_metadata, by="col.id") %>%
  dplyr::select(-row.id, -col.id, -aus_eco)

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


# convert points from FrogID into sf points
# so that they can be used with the polygons
points_sf <- st_as_sf(dat2, coords = c("lng", "lat"), 
                      crs = st_crs(urban_areas), agr = "constant")

# use intersect to return a list of points and the associated
# polygon they belong to - uses integer ids for rows and columns
assigned_points <- as.data.frame(st_intersects(points_sf, urban_areas))

# join the original dataframe with the intersected df
# then remove the row and col ids
dat3 <- dat2 %>%
  mutate(row.id = 1:nrow(.)) %>%
  left_join(assigned_points, by="row.id") %>%
  left_join(., urban_areas, by="col.id") %>%
  dplyr::select(-row.id, -col.id, -geometry)

# now export this data as analysis dat
saveRDS(dat3, "Data/dat_for_urban_level_analysis.RDS")



