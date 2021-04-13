# This is a script to do an analysis of the urban areas
# versus non urban areas
# broad-level analyses showing the differences in frog species richness and diversity of these two regions

# packages
library(dplyr)
library(vegan)
library(ggplot2)
library(sf)
library(tidyr)
library(tibble)
library(readr)

# read in data for analysis
# add a column whether a record is within a SUA or not
# filter out records which were not appropriately assigned an Ecoregion
# remove the Montane Grasslands & Shrublands as it has a few records and no urban areas within it apparently
dat <- readRDS("Data/dat_for_urban_level_analysis.RDS") %>%
  mutate(urban=ifelse(is.na(SUA_NAME16)==TRUE, "Non-urban", "Urban")) %>%
  dplyr::filter(complete.cases(Ecoregion)) %>%
  dplyr::filter(Ecoregion != "Montane Grasslands & Shrublands")

family <- read_csv("Data/Australian_Frog_List_Plus_ConsStat_Family_April_2020.csv") %>%
  dplyr::rename(IUCN=`IUCN Status`)

# summary stuff
length(unique(dat$species))
nrow(dat)

# first look at overall the bias in records from SUAs or not
dat %>%
  .$urban %>%
  table() %>%
  as.data.frame() %>%
  mutate(total_records=sum(Freq)) %>%
  mutate(percent=(Freq/total_records)*100)

# species richness in urban versus non-urban areas
dat %>%
  group_by(urban) %>%
  summarize(species_richness=length(unique(species)),
            number_records=length(unique(id)))

# find the number of species unique to urban areas and unique to non-urban areas
# and the species found in both
temp <- dat %>%
  group_by(urban, species) %>%
  summarize(species_richness=length(unique(species))) %>%
  group_by(species) %>%
  mutate(N=n())

# number of species found in both urban and non urban
both_urban_and_nonurban <- temp %>%
  dplyr::select(species, N) %>%
  distinct() %>%
  dplyr::filter(N==2) 

nrow(both_urban_and_nonurban)

# number of species found only in non-urban
only_non_urban <- temp %>%
  dplyr::filter(N==1) %>%
  dplyr::filter(urban=="Non-urban") 

nrow(only_non_urban)

# number of species found only in urban
only_urban <- temp %>%
  dplyr::filter(N==1) %>%
  dplyr::filter(urban=="Urban")

nrow(only_urban)

# Now look at total species richness in urban versus non-urban areas by ecoregion
ecoregion_summary <- dat %>%
  group_by(Ecoregion, urban) %>%
  summarize(species_richness=length(unique(species)),
            number_records=length(unique(id)),
            total_urban_area=sum(AREASQKM16))

shared_species <- dat %>%
  group_by(Ecoregion, urban, species) %>%
  summarize(N=n()) %>%
  group_by(species, Ecoregion) %>%
  mutate(N=n()) %>%
  dplyr::filter(N==2) %>%
  group_by(Ecoregion, urban) %>%
  summarize(shared_species=n())

# threatened species found in urban areas
dat %>%
  left_join(., family) %>%
  group_by(urban, IUCN) %>%
  summarize(species_richness=length(unique(species)),
            number_records=length(unique(id)))



ecoregion_summary %>%
  dplyr::select(-total_urban_area) %>%
  write_csv(., "Results/ecoregion_urban_richness_summary.csv")


# NOW REPEAT ALL OF THE ABOVE BUT FOR FAMILY
# find the number of species unique to urban areas and unique to non-urban areas
# and the species found in both
temp2 <- dat %>%
  left_join(., family) %>%
  group_by(urban, Family) %>%
  summarize(family_richness=length(unique(Family)),
            family_records=n()) %>%
  group_by(Family) %>%
  mutate(N=n())

# number of species found in both urban and non urban
temp2 %>%
  dplyr::select(Family, N) %>%
  distinct() %>%
  dplyr::filter(N==2) %>%
  nrow(.)

# number of species found only in non-urban
temp2 %>%
  dplyr::filter(N==1) %>%
  dplyr::filter(urban=="Non-urban") %>%
  nrow(.)

# number of species found only in urban
temp2 %>%
  dplyr::filter(N==1) %>%
  dplyr::filter(urban=="Urban") %>%
  nrow(.)

# Now look at total species richness in urban versus non-urban areas by ecoregion
ecoregion_summary2 <- dat %>%
  left_join(., family) %>%
  group_by(Ecoregion, urban) %>%
  summarize(family_richness=length(unique(Family)),
            number_records=length(unique(id)),
            total_urban_area=sum(AREASQKM16))

shared_species <- dat %>%
  group_by(Ecoregion, urban, species) %>%
  summarize(N=n()) %>%
  group_by(species, Ecoregion) %>%
  mutate(N=n()) %>%
  dplyr::filter(N==2) %>%
  group_by(Ecoregion, urban) %>%
  summarize(shared_species=n())

# threatened species found in urban areas
dat %>%
  left_join(., family) %>%
  group_by(urban, IUCN) %>%
  summarize(species_richness=length(unique(species)),
            number_records=length(unique(id)))

# get a few example species that are found only in urban and
# non-urban areas
ex_sp <- dat %>%
  left_join(., family) %>%
  group_by(urban, IUCN, species) %>%
  summarize(number_records=length(unique(id)))

both <- ex_sp %>%
  dplyr::filter(species %in% both_urban_and_nonurban$species)

non_urban <- ex_sp %>%
  dplyr::filter(species %in% only_non_urban$species)

urban <- ex_sp %>%
  dplyr::filter(species %in% only_urban$species)
