# An analysis to assess urbanization gradient influences on
# biodiversity of frogs

# packages
library(dplyr)
library(vegan)
library(ggplot2)
library(sf)
library(tidyr)
library(tibble)
library(ape)
library(picante)
library(phylobase)
library(purrr)

# read in data for analysis
# add a column whether a record is within a SUA or not
# filter out records which were not appropriately assigned an Ecoregion
# remove the Montane Grasslands & Shrublands as it has a few records and no urban areas within it apparently
dat <- readRDS("Data/dat_for_urban_level_analysis.RDS") %>%
  mutate(urban=ifelse(is.na(SUA_NAME16)==TRUE, "Non-urban", "Urban")) %>%
  dplyr::filter(complete.cases(Ecoregion)) %>%
  dplyr::filter(Ecoregion != "Montane Grasslands & Shrublands")

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

# read in data from GEE
# from Gracie
urbanization_dat <- readRDS("Data/FrogIDdata_wModScores.RDS") %>%
  dplyr::filter(id %in% dat$id) %>%
  dplyr::filter(complete.cases(zviirs_dnb_avg_rad.1km)) %>%
  left_join(., dat %>%
              dplyr::select(id, Ecoregion) %>%
              distinct()) %>%
  left_join(., dat %>%
              dplyr::select(species, species2) %>%
              distinct())

# perform a resampling analysis that selects
# a random quantile value of the distribution of urban values
# and then selects all observations (within an ecoregion)
# that are within a quantile of that distribution
# and calculates the diversity metrics of those records
# so some records will get used many times
# and there should be less samples at the top end of the
# distribution

# first get a list of unique urban values in the analysis
urban_dist <- urbanization_dat %>%
  dplyr::select(zviirs_dnb_avg_rad.1km)

# make a plot of this distribution
ggplot(urbanization_dat, aes(x=zviirs_dnb_avg_rad.1km))+
  geom_density()+
  geom_rug(alpha = 1/2)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(strip.text=element_text(size=7))+
  scale_x_log10()+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  ylab("Density of observations")+
  facet_wrap(~Ecoregion, scales="free")

ggsave("Figures/urbanization_observation_density_distributions.png", width=7.5, height=5, units="in")

gradient_analysis_function <- function(ecoregion_name){
  
  temp <- urbanization_dat %>%
    dplyr::filter(Ecoregion==ecoregion_name)
  
  # get a random quantile to sample
  random_value <- sample(seq(0.01, 0.99, 0.001), 1)
  
  upper_quantile <- random_value + 0.05
  upper_quantile <- ifelse(upper_quantile>=1, 1, upper_quantile)
  
  lower_quantile <- random_value - 0.05
  lower_quantile <- ifelse(lower_quantile<=0, 0, lower_quantile)
  
  # get the corresponding value of the distribution of urban gradient
  urban_value <- quantile(temp$zviirs_dnb_avg_rad.1km, random_value)
  urban_value_lower <- quantile(temp$zviirs_dnb_avg_rad.1km, lower_quantile)
  urban_value_upper <- quantile(temp$zviirs_dnb_avg_rad.1km, upper_quantile)
  
  filtered_dat <- temp %>%
    dplyr::filter(zviirs_dnb_avg_rad.1km >= urban_value_lower) %>%
    dplyr::filter(zviirs_dnb_avg_rad.1km <= urban_value_upper) 
  
  filtered_dat <- if (nrow(filtered_dat)>=100) {
    filtered_dat %>%
      sample_n(100)
  } else {
    filtered_dat
  }
  
  SR=length(unique(filtered_dat$species))
  
  SD=filtered_dat %>%
    group_by(species) %>%
    summarize(species_count=n()) %>%
    .$species_count %>%
    diversity()
  
  PD=filtered_dat %>%
    mutate(species2=gsub(" ", "_", .$species2)) %>%
    group_by(species2) %>%
    summarize(species_count=n()) %>%
    pivot_wider(names_from=species2, values_from=species_count, values_fill=list(species_count=0)) %>%
    as.matrix() %>%
    pd(., tree=our_tree2) %>%
    as.data.frame() %>%
    dplyr::select(PD) %>%
    .$PD
  
  summary <- data.frame(SR=SR,
                        SD=SD,
                        PD=PD,
                        VIIRS_lights_1km=urban_value,
                        Number_obs=nrow(filtered_dat)) %>%
    mutate(Ecoregion=ecoregion_name)
  
  return(summary)
}

# Note that deserts throws an error because there is not enough data present for this ecoregion
# due to eliminating some records for this code
temp_grasslands <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Temperate Grasslands")))
deserts <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Deserts & Xeric Shrublands")))
med_forests <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Mediterranean Forests")))
temp_forests <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Temperate Broadleaf & Mixed Forests")))
tropical_grasslands <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Tropical & Subtropical Grasslands")))
tropical_forests <- do.call("rbind", purrr::rerun(10000, gradient_analysis_function("Tropical & Subtropical Moist Broadleaf Forests")))

urban_gradients <- temp_grasslands %>%
  #bind_rows(deserts) %>%
  bind_rows(med_forests) %>%
  bind_rows(temp_forests) %>%
  bind_rows(tropical_grasslands) %>%
  bind_rows(tropical_forests) %>%
  mutate(Ecoregion=as.factor(as.character(Ecoregion)))

SR_plot <- ggplot(urban_gradients, aes(x=VIIRS_lights_1km, y=SR, color=Ecoregion))+
  #geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))+
  #geom_smooth(method="lm")+
  scale_color_brewer(palette="Set2")+
  scale_x_log10()+
  #scale_y_log10()+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  ylab("Species Richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)
  #facet_wrap(~Ecoregion, scales="free")

SR_plot

ggsave("Figures/SR_gradient.png", width=4.5, height=4, units="in")
ggsave("Figures/SR_gradient.eps", width=4.5, height=4, units="in")
ggsave("Figures/SR_gradient.svg", width=4.5, height=4, units="in")

SD_plot <- ggplot(urban_gradients, aes(x=VIIRS_lights_1km, y=SD, color=Ecoregion))+
  #geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))+
  scale_color_brewer(palette="Set2")+
  scale_x_log10()+
  #scale_y_log10()+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  ylab("Shannon Diversity")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)
#facet_wrap(~Ecoregion, scales="free")

SD_plot

ggsave("Figures/SD_gradient.png", width=4.5, height=4, units="in")
ggsave("Figures/SD_gradient.eps", width=4.5, height=4, units="in")
ggsave("Figures/SD_gradient.svg", width=4.5, height=4, units="in")

PD_plot <- ggplot(urban_gradients, aes(x=VIIRS_lights_1km, y=PD, color=Ecoregion))+
  #geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))+
  scale_color_brewer(palette="Set2")+
  scale_x_log10()+
  #scale_y_log10()+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  ylab("Phylogenetic Diversity")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)
#facet_wrap(~Ecoregion, scales="free")

PD_plot

ggsave("Figures/PD_gradient.png", width=4.5, height=4, units="in")
ggsave("Figures/PD_gradient.eps", width=4.5, height=4, units="in")
ggsave("Figures/PD_gradient.svg", width=4.5, height=4, units="in")

gradient_legend <- ggplot(urban_gradients, aes(x=VIIRS_lights_1km, y=SR, color=Ecoregion))+
  #geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))+
  #geom_smooth(method="lm")+
  scale_color_brewer(palette="Set2")+
  scale_x_log10()+
  #scale_y_log10()+
  xlab(bquote('Average radiance ('* 'nW' ~cm^-2~sr^-1*')'))+
  ylab("Species Richness")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

gradient_legend

ggsave("Figures/gradient_legend.png", width=4.5, height=4, units="in")
ggsave("Figures/gradient_legend.eps", width=4.5, height=4, units="in")
ggsave("Figures/gradient_legend.svg", width=4.5, height=4, units="in")

library(mgcv)

SR_mod <- gam(SR ~ s(log10(VIIRS_lights_1km), k=6) + s(Ecoregion, bs="re"), data=urban_gradients)
plot(SR_mod)
summary(SR_mod)

SR_mod.2 <- gam(SR ~ s(log10(VIIRS_lights_1km), k=6, by=Ecoregion), data=urban_gradients)
plot(SR_mod.2)




