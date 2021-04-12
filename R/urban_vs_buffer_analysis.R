# This is an analysis where each city/SUA receives a buffer
# and then frog biodiversity is calculated within the SUA and the buffer
# to create a paired analysis of SUA and associated buffers

# packages
library(dplyr)
library(vegan)
library(ggplot2)
library(sf)
library(tidyr)
library(tibble)
library(scales)
library(ape)
library(picante)
library(phylobase)
library(patchwork)
library(readr)
library(lubridate)

# read in data for analysis
# add a column whether a record is within a SUA or not
# filter out records which were not appropriately assigned an Ecoregion
# remove the Montane Grasslands & Shrublands as it has a few records and no urban areas within it apparently
dat <- readRDS("Data/dat_for_urban_level_analysis.RDS") %>%
  mutate(urban=ifelse(is.na(SUA_NAME16)==TRUE, "Non-urban", "Urban")) %>%
  dplyr::filter(complete.cases(Ecoregion)) %>%
  dplyr::filter(Ecoregion != "Montane Grasslands & Shrublands")

# turn into spatial points for later
spatial_points <- dat %>%
  st_as_sf(coords=c("lng", "lat"), crs=4283)

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

# summarize the number of obs per SUA
SUA_summary <- dat %>%
  group_by(SUA_NAME16) %>%
  summarize(N=n())

ggplot(SUA_summary, aes(x=N))+
  geom_histogram(color="black", fill="gray30", bins=40)+
  scale_x_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Number of observations")+
  ylab("Number of SUAs")+
  geom_vline(xintercept=150, color="red", linetype="dashed")

ggsave("Figures/hist_of_SUA_vs_obs.png", width=4.8, height=3.8, units="in")

# quite a few SUAs have very few observations!
# so we will use a cutoff where we only investigate
# those SUAs >= 150 observations
# this will limit the number of SUAs used in the analysis
SUAs <- SUA_summary %>%
  dplyr::filter(N>=150) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  .$SUA_NAME16

# read in phylo tree
# and prepare the tree for analyses
tree <- ape::read.tree(file="Data/phylo_data/Australia/aus_phylo_tree.tre")

## Make a list of spp in our dataset
our_species <- dat %>%
  mutate(species=gsub(" ", "_", .$species)) %>%
  dplyr::select(species) %>%
  distinct() %>%
  .$species

our_sp_list <- data.frame(species = our_species)

## Make a list of spp in the phylo tree
spp_in_tree <- data.frame(species = tree[["tip.label"]])

## Filter the tree list to spp not in our dataset 
no_data_sp <- tree$tip.label[!tree$tip.label %in% our_species]

## Drop the spp in the tree that are not in our dataset so that 
## you're left with spp in the tree that are also in our list
our_tree <- drop.tip(tree, no_data_sp)
our_tree_sp <- data.frame(species = our_tree[["tip.label"]]) %>% 
  arrange(.$species)

## View spp that are in our list that are not in our tree 
anti_join(our_sp_list, our_tree_sp) 

## NOTE that there are mismatches b/w the tree and our taxonomy:


## rename species in our dataset
## so that they match with the phylogenetic tree
dat <- dat %>%
  mutate(species2=gsub("Litoria verreauxii", "Cyclorana verreauxii", species)) %>%
  mutate(species2=gsub("Litoria nudidigitus", "Litoria nudidigita", species2)) %>%
  mutate(species2=gsub("Cyclorana verrucosa", "Litoria verrucosa", species2)) %>%
  mutate(species2=gsub("Uperoleia saxatilis", "Uperoleia russelli", species2)) %>%
  mutate(species2=gsub("Cyclorana occidentalis", "Litoria verrucosa", species2)) %>%
  mutate(species2=gsub("Litoria burrowsae", "Litoria burrowsi", species2)) %>%
  mutate(species2=gsub("Litoria bella", "Litoria chloris", species2)) %>%
  mutate(species2=gsub("Uperoleia mahonyi", "Uperoleia laevigata", species2))

# now get our tree again
# using the renamed species
## Make a list of spp in our dataset
our_species2 <- dat %>%
  mutate(species2=gsub(" ", "_", .$species2)) %>%
  dplyr::select(species2) %>%
  distinct() %>%
  .$species2

our_sp_list2 <- data.frame(species = our_species2)

## Make a list of spp in the phylo tree
spp_in_tree2 <- data.frame(species = tree[["tip.label"]])

## Filter the tree list to spp not in our dataset 
no_data_sp2 <- tree$tip.label[!tree$tip.label %in% our_species2]

## Drop the spp in the tree that are not in our dataset so that 
## you're left with spp in the tree that are also in our list
our_tree2 <- drop.tip(tree, no_data_sp2)
our_tree_sp2 <- data.frame(species = our_tree2[["tip.label"]]) %>% 
  arrange(.$species)

## View spp that are in our list that are not in our tree 
anti_join(our_sp_list2, our_tree_sp2) 

# Now write a function that takes a given SUA
# and calculates the necessary data for both that SUA 
# and for the buffer surrounding that SUA
# for different buffer sizes as well

SUA_specific_summary_function <- function(SUA_name, buffer_size){
  
  # get the specific city spatial components from
  # the urban areas shapefile
  city_spatial <- urban_areas %>%
    dplyr::filter(SUA_NAME16==SUA_name) 
  
  # first get summary of data from within the SUA
  within_summary <- dat %>%
    dplyr::filter(SUA_NAME16==SUA_name) %>%
    group_by(SUA_NAME16) %>%
    summarize(Number_obs=n(),
              SR=length(unique(species)))
  
  SD <- dat %>%
    dplyr::filter(SUA_NAME16==SUA_name) %>%
    group_by(SUA_NAME16, species) %>%
    summarize(species_count=n()) %>%
    group_by(SUA_NAME16) %>%
    summarize(SD=diversity(species_count))
  
  PD <- dat %>%
    dplyr::filter(SUA_NAME16==SUA_name) %>%
    mutate(species2=gsub(" ", "_", .$species2)) %>%
    group_by(SUA_NAME16, species2) %>%
    summarize(species_count=n()) %>%
    pivot_wider(names_from=species2, values_from=species_count, values_fill=list(species_count=0)) %>%
    as.matrix() %>%
    pd(., tree=our_tree2) %>%
    as.data.frame() %>%
    dplyr::select(PD) %>%
    .$PD
  
  within_summary <- within_summary %>%
    left_join(., SD) %>%
    mutate(PD=PD) %>%
    mutate(area=as.numeric(st_area(city_spatial)))
  
  
  # now create buffer and get all records within buffer
  # first need to convert the spatial data to UTM
  lonlat2UTM <- function(lonlat) {
    utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
    if(lonlat[2] > 0) {
      utm + 32600
    } else{
      utm + 32700
    }
  }
  
  # get coordinates 
  coords <- st_centroid(city_spatial) %>%
    st_coordinates() %>%
    as.data.frame()
  
  # get UTM
  EPSG_2_UTM <- lonlat2UTM(c(coords$X, coords$Y))
  
  # transform city to UTM
  city_projected <- st_transform(st_as_sf(city_spatial), EPSG_2_UTM)
  
  # create buffer and then transform back to EPSG of city
  buffer <- city_projected %>%
    st_buffer(dist=buffer_size) %>%
    st_transform(crs=st_crs(city_spatial))
  
  # get the difference between the city and buffer
  difference <- st_difference(buffer, city_spatial)
  
  # now get all the observations that are within this 'difference' polygon
  points_in_buffer <- difference %>%
    st_intersects(spatial_points) %>%
    as.data.frame() %>%
    left_join(., dat %>%
                mutate(col.id=1:nrow(.)), by="col.id")
  
  buffer_PD <- points_in_buffer %>%
    mutate(species2=gsub(" ", "_", .$species2)) %>%
    group_by(species2) %>%
    summarize(species_count=n()) %>%
    pivot_wider(names_from=species2, values_from=species_count, values_fill=list(species_count=0)) %>%
    as.matrix() %>%
    pd(., tree=our_tree2) %>%
    as.data.frame() %>%
    dplyr::select(PD) %>%
    .$PD
  
  # summarize stats for buffer
  buffer_summary <- data.frame(SUA_NAME16=SUA_name) %>%
    mutate(buffer_Number_obs=nrow(points_in_buffer)) %>%
    mutate(buffer_SR=length(unique(points_in_buffer$species))) %>%
    mutate(buffer_SD=points_in_buffer %>%
             group_by(species) %>%
             summarize(species_count=n()) %>%
             .$species_count %>%
             diversity()) %>%
    mutate(buffer_PD=buffer_PD) %>%
    mutate(buffer_area=as.numeric(st_area(difference)))
  
  # create a final summary
  final_summary <- within_summary %>%
    left_join(., buffer_summary)
  
}

# now apply the function to each SUA
SUA_buffer_comparison_50_km <- bind_rows(lapply(SUAs, function(x) {SUA_specific_summary_function(x, 50000)}))

SUA_buffer_comparison_100_km <- bind_rows(lapply(SUAs, function(x) {SUA_specific_summary_function(x, 100000)}))

SUA_buffer_comparison_150_km <- bind_rows(lapply(SUAs, function(x) {SUA_specific_summary_function(x, 150000)}))

# make a plot of these data at 100 km
# buffer sizes
# first get rid of any cities which had small number of obs in buffer
SR_plot <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, SR, buffer_SR) %>%
  dplyr::rename(Urban=SR) %>%
  dplyr::rename(`Non-urban`=buffer_SR) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to="value") %>%
  ggplot(., aes(x=comparison, y=value, group=SUA_NAME16))+
  geom_point()+
  geom_line()+
  xlab("")+
  ylab("Species Richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

SR_plot

ggsave("Figures/SR_buffer_comp.png", width=4.5, height=4, units="in")
ggsave("Figures/SR_buffer_comp.eps", width=4.5, height=4, units="in")
ggsave("Figures/SR_buffer_comp.svg", width=4.5, height=4, units="in")

SD_plot <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, SD, buffer_SD) %>%
  dplyr::rename(Urban=SD) %>%
  dplyr::rename(`Non-urban`=buffer_SD) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to="value") %>%
  ggplot(., aes(x=comparison, y=value, group=SUA_NAME16))+
  geom_point()+
  geom_line()+
  xlab("")+
  ylab("Shannon Diversity")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

SD_plot

ggsave("Figures/SD_buffer_comp.png", width=4.5, height=4, units="in")
ggsave("Figures/SD_buffer_comp.eps", width=4.5, height=4, units="in")
ggsave("Figures/SD_buffer_comp.svg", width=4.5, height=4, units="in")

PD_plot <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, PD, buffer_PD) %>%
  dplyr::rename(Urban=PD) %>%
  dplyr::rename(`Non-urban`=buffer_PD) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to="value") %>%
  ggplot(., aes(x=comparison, y=value, group=SUA_NAME16))+
  geom_point()+
  geom_line()+
  xlab("")+
  ylab("Phylogenetic Diversity")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

PD_plot

ggsave("Figures/PD_buffer_comp.png", width=4.5, height=4, units="in")
ggsave("Figures/PD_buffer_comp.eps", width=4.5, height=4, units="in")
ggsave("Figures/PD_buffer_comp.svg", width=4.5, height=4, units="in")

SR_plot + SD_plot + PD_plot + plot_layout(ncol=1)

ggsave("Figures/Figure_3.png", height=8, width=4, units="in")

SR_mod <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, SR, buffer_SR) %>%
  dplyr::rename(Urban=SR) %>%
  dplyr::rename(`Non-urban`=buffer_SR) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to=c("SR")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, Number_obs, buffer_Number_obs) %>%
              dplyr::rename(Urban=Number_obs) %>%
              dplyr::rename(`Non-urban`=buffer_Number_obs) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="N")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, area, buffer_area) %>%
              dplyr::rename(Urban=area) %>%
              dplyr::rename(`Non-urban`=buffer_area) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="Area")) %>%
  lm(SR ~ N + comparison, data=.)

summary(SR_mod)

SD_mod <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, SD, buffer_SD) %>%
  dplyr::rename(Urban=SD) %>%
  dplyr::rename(`Non-urban`=buffer_SD) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to=c("SD")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, Number_obs, buffer_Number_obs) %>%
              dplyr::rename(Urban=Number_obs) %>%
              dplyr::rename(`Non-urban`=buffer_Number_obs) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="N")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, area, buffer_area) %>%
              dplyr::rename(Urban=area) %>%
              dplyr::rename(`Non-urban`=buffer_area) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="Area")) %>%
  lm(SD ~ N + comparison, data=.)

summary(SD_mod)


PD_mod <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  dplyr::select(SUA_NAME16, PD, buffer_PD) %>%
  dplyr::rename(Urban=PD) %>%
  dplyr::rename(`Non-urban`=buffer_PD) %>%
  pivot_longer(-SUA_NAME16, names_to="comparison", values_to=c("PD")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, Number_obs, buffer_Number_obs) %>%
              dplyr::rename(Urban=Number_obs) %>%
              dplyr::rename(`Non-urban`=buffer_Number_obs) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="N")) %>%
  left_join(., SUA_buffer_comparison_100_km %>%
              dplyr::filter(buffer_Number_obs >=150) %>%
              dplyr::select(SUA_NAME16, area, buffer_area) %>%
              dplyr::rename(Urban=area) %>%
              dplyr::rename(`Non-urban`=buffer_area) %>%
              pivot_longer(-SUA_NAME16, names_to="comparison", values_to="Area")) %>%
  lm(PD ~ N + comparison, data=.)

summary(PD_mod)


# make a figure showing the SUAs and the species richness
# as a function of area of a SUA
dat %>%
  dplyr::filter(SUA_NAME16 %in% SUAs) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  group_by(SUA_NAME16) %>%
  summarize(SR=length(unique(species)),
            area=mean(AREASQKM16)) %>%
  ggplot(., aes(x=area, y=SR))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Area (square km)")+
  ylab("Species Richness")

ggsave("Figures/area_richness_curve.png", width=4.8, height=3.8, units="in")

dat %>%
  dplyr::filter(SUA_NAME16 %in% SUAs) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  group_by(SUA_NAME16) %>%
  summarize(SR=length(unique(species)),
            area=mean(AREASQKM16)) %>%
  lm(log10(SR) ~ log10(area), data=.) %>%
  summary()


# summary stuff
comparison_summary <- SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  mutate(SR_comp=((SR/buffer_SR)*100)) %>%
  mutate(SD_comp=((SD/buffer_SD)*100)) %>%
  mutate(PD_comp=((PD/buffer_PD)*100))
  
mean(comparison_summary$SR_comp)
mean(comparison_summary$SD_comp)
mean(comparison_summary$PD_comp)

# write out table for supplementary table
SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  write_csv(., "Results/buffer_analysis_results.csv")


# How does the area relate to each other
SUA_buffer_comparison_100_km %>%
  dplyr::filter(buffer_Number_obs >=150) %>%
  ggplot(., aes(x=area/1000000, y=buffer_area/1000000))+
  geom_point()+
  scale_x_log10(labels=comma)+
  scale_y_log10(labels=comma)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Urban area (km2)")+
  ylab("Buffer area (km2)")+
  geom_smooth(method="lm")

ggsave("Figures/buffer_vs_urban_area_size.png", width=4.8, height=3.8, units="in")



###################################################################################
###################################################################################
###################################################################################
############################# BUT TEST THE INFLUENCE OF TEMPORAL BIASES ON THESE RESULTS
########### It is possible that the 150 cutoff could be from different months for different seasons
########### So let's check that first
# summarize the number of obs per SUA
# First make a figure of the overall monthly bias
dat %>%
  mutate(MONTH=month(date, abbr=TRUE, label=TRUE)) %>%
  group_by(MONTH) %>%
  summarize(N=n()) %>%
  ggplot(., aes(x=MONTH, y=N))+
  geom_col()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Number of observations")

# get a summary of the number of obs
# per month per SUA
SUA_summary2 <- dat %>%
  group_by(SUA_NAME16) %>%
  mutate(N=n()) %>%
  mutate(MONTH=month(date, abbr=TRUE, label=TRUE)) %>%
  group_by(SUA_NAME16, MONTH) %>%
  summarize(number_obs=n(),
            total_obs=mean(N))

# now get rid of any SUA that has <150 total obs
# and make a plot of ten of them
SUA_summary2 %>%
  dplyr::filter(total_obs>=150) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  ggplot(., aes(x=MONTH, y=number_obs))+
  geom_col()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Number of observations")+
  facet_wrap(~SUA_NAME16, scales="free_y")
  

# try another way to look at this
SUA_summary2 %>%
  dplyr::filter(total_obs>=150) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  group_by(SUA_NAME16) %>%
  mutate(number_obs_scaled=scales::rescale(number_obs)) %>%
  ggplot(., aes(x=MONTH, y=number_obs_scaled, group=SUA_NAME16))+
  geom_line(alpha=0.9, color="gray80")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Number of observations (scaled)")


# now let's plot it by ecoregion
# because if the temporal biases are equalish within ecoregions
# then this implies the bias is equal among SUAs
# therefore, we need to get the 'mode' for each SUA
# mode is necessary because a couple SUAs have span multiple ecoregions
test <- dat %>% 
  group_by(SUA_NAME16) %>% 
  summarize(N=length(unique(Ecoregion)))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

SUA_summary2 %>%
  dplyr::filter(total_obs>=150) %>%
  dplyr::filter(complete.cases(SUA_NAME16)) %>%
  group_by(SUA_NAME16) %>%
  mutate(number_obs_scaled=scales::rescale(number_obs)) %>%
  left_join(., dat %>% 
              group_by(SUA_NAME16) %>%
              summarize(Ecoregion=Mode(Ecoregion))) %>%
  ggplot(., aes(x=MONTH, y=number_obs_scaled, group=SUA_NAME16, color=Ecoregion))+
  geom_line()+
  scale_color_brewer(palette = "Set2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Number of observations (scaled)")+
  facet_wrap(~Ecoregion)+
  guides(color=FALSE)+
  theme(axis.text.x=element_text(size=6))+
  theme(strip.text.x = element_text(size = 7.5))

ggsave("Figures/ecoregion_temporal_comparison.png", width=7.4, height=5.8, units="in")

# Look at Perth and Darwin as examples
SUA_summary2 %>%
  dplyr::filter(SUA_NAME16 %in% c("Darwin", "Perth")) %>%
  ggplot(., aes(x=MONTH, y=number_obs, group=SUA_NAME16, fill=SUA_NAME16))+
  geom_col()+
  scale_fill_brewer(palette="Set2")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Number of observations")+
  facet_wrap(~SUA_NAME16, scales="free_y")+
  guides(fill=FALSE)

ggsave("Figures/perth_darwin_temporal_comparison.png", width=6.8, height=4.8, units="in")





