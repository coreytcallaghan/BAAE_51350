# This is a script to get resampled summary statistics of
# inside and outside urban areas

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
library(patchwork)
library(lme4)
library(lmerTest)
library(broom)

# source function for plotting later
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

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


#############################
#############################
# Species Richness resampling
SR_function <- function(draw){
  
  apply_over_sample_sizes <- function(N){
    
    temp <- dat %>%
      group_by(Ecoregion, urban) %>%
      sample_n(N)
    
    SR_summary <- temp %>%
      group_by(Ecoregion, urban) %>%
      summarize(SR=length(unique(species))) %>%
      mutate(sample_size=N)
    
    return(SR_summary)
  }
  
  SR_results_intermediate <- bind_rows(lapply(seq(50, 100, by=10), function(x) {apply_over_sample_sizes(x)})) %>%
    mutate(resample_number=draw)
  
  return(SR_results_intermediate)
}

SR_results <- bind_rows(lapply(c(1:1000), function(x) {SR_function(x)}))

SR_plot <- SR_results %>%
  dplyr::filter(sample_size==100) %>%
  ggplot(., aes(x=stringr::str_wrap(Ecoregion, 15), y=SR, fill=urban))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
  geom_point(aes(x=stringr::str_wrap(Ecoregion, 15), y=SR, color=urban), 
             position = position_jitter(width = .15), size = .5, alpha = 0.8)+
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_fill_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  scale_color_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  xlab("")+
  ylab("Species Richness")+
  guides(fill=FALSE)+
  guides(color=FALSE)+
  theme(axis.text.y=element_text(size=6))+
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.x=element_blank())+
  #ggtitle("A")+
  theme(plot.title = element_text(vjust=-6.5, hjust=0.01, size=7))

SR_plot

##############################
##############################
# Shannon diversity resampling
SD_function <- function(draw){
  
  apply_over_sample_sizes <- function(N){
    
    temp <- dat %>%
      group_by(Ecoregion, urban) %>%
      sample_n(N) %>%
      group_by(Ecoregion, urban, species) %>%
      summarize(species_count=n())
    
    SD_summary <- temp %>%
      group_by(Ecoregion, urban) %>%
      summarize(SD=diversity(species_count)) %>%
      mutate(sample_size=N)
    
    return(SD_summary)
  }
  
  SD_results_intermediate <- bind_rows(lapply(seq(50, 100, by=10), function(x) {apply_over_sample_sizes(x)})) %>%
    mutate(resample_number=draw)
  
  return(SD_results_intermediate)
}

SD_results <- bind_rows(lapply(c(1:1000), function(x) {SD_function(x)}))

SD_plot <- SD_results %>%
  dplyr::filter(sample_size==100) %>%
  ggplot(., aes(x=stringr::str_wrap(Ecoregion, 15), y=SD, fill=urban))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
  geom_point(aes(x=stringr::str_wrap(Ecoregion, 15), y=SD, color=urban), 
             position = position_jitter(width = .15), size = .5, alpha = 0.8)+
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_fill_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  scale_color_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  xlab("")+
  ylab("Shannon Diversity")+
  guides(fill=FALSE)+
  guides(color=FALSE)+
  theme(axis.text.y=element_text(size=6))+
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.x=element_blank())+
  #ggtitle("B")+
  theme(plot.title = element_text(vjust=-6.5, hjust=0.01, size=7))

SD_plot


##############################
##############################
# Phylogenetic diversity resampling
PD_function <- function(draw){
  
  apply_over_sample_sizes <- function(N){
    
    temp <- dat %>%
      group_by(Ecoregion, urban) %>%
      sample_n(N) %>%
      group_by(Ecoregion, urban, species2) %>%
      summarize(species_count=n())
    
    apply_over_ecoregions <- function(ecoregion){
      
      temp2 <- temp %>%
        dplyr::filter(Ecoregion==ecoregion)
      
      urban_nonurban_matrix <- temp2 %>%
        ungroup() %>%
        mutate(species2=gsub(" ", "_", .$species2)) %>%
        pivot_wider(names_from=species2, values_from=species_count, values_fill=list(species_count=0)) %>%
        dplyr::select(-Ecoregion) %>%
        column_to_rownames(var="urban") %>%
        as.matrix()
      
      pd_results <- pd(urban_nonurban_matrix, tree=our_tree2) %>%
        as.data.frame() %>%
        rownames_to_column(var="urban") %>%
        mutate(Ecoregion=ecoregion)
      
    }
    
    PD_summary <- bind_rows(lapply(unique(temp$Ecoregion), function(x) {apply_over_ecoregions(x)})) %>%
      mutate(sample_size=N)
    
    return(PD_summary)
  }
  
  PD_results_intermediate <- bind_rows(lapply(seq(90, 100, by=10), function(x) {apply_over_sample_sizes(x)})) %>%
    mutate(resample_number=draw)
  
  return(PD_results_intermediate)
}

PD_results <- bind_rows(lapply(c(1:1000), function(x) {PD_function(x)})) %>%
  dplyr::select(Ecoregion, urban, PD, sample_size, resample_number)

PD_plot <- PD_results %>%
  dplyr::filter(sample_size==100) %>%
  ggplot(., aes(x=stringr::str_wrap(Ecoregion, 15), y=PD, fill=urban))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8)+
  geom_point(aes(x=stringr::str_wrap(Ecoregion, 15), y=PD, color=urban), 
             position = position_jitter(width = .15), size = .5, alpha = 0.8)+
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5)+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_fill_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  scale_color_manual(breaks=c("Urban", "Non-urban"), values=c("gray20", "#006400"))+
  xlab("")+
  ylab("Phylogenetic Diversity")+
  #guides(fill=FALSE)+
  guides(color=FALSE)+
  theme(axis.text.y=element_text(size=6))+
  theme(panel.grid.major.x=element_blank())+
  theme(panel.grid.minor.x=element_blank())+
  #ggtitle("C")+
  theme(plot.title = element_text(vjust=-6.5, hjust=0.01, size=7))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())

PD_plot

SR_plot

ggsave("Figures/SR_resampled.png", width=4.5, height=4, units="in")
ggsave("Figures/SR_resampled.eps", width=4.5, height=4, units="in")
ggsave("Figures/SR_resampled.svg", width=4.5, height=4, units="in")

SD_plot

ggsave("Figures/SD_resampled.png", width=4.5, height=4, units="in")
ggsave("Figures/SD_resampled.eps", width=4.5, height=4, units="in")
ggsave("Figures/SD_resampled.svg", width=4.5, height=4, units="in")

PD_plot

ggsave("Figures/PD_resampled.png", width=4.5, height=4, units="in")
ggsave("Figures/PD_resampled.eps", width=4.5, height=4, units="in")
ggsave("Figures/PD_resampled.svg", width=4.5, height=4, units="in")

#SR_plot + SD_plot + PD_plot + plot_layout(ncol=1)

#ggsave("Figures/Figure_2.png", height=8, width=4, units="in")

# perform some simple models to see whether there is statistical difference between
# each urban and non urban form
SR_dat <- SR_results %>%
  dplyr::filter(sample_size==100)
SR_mod <- lmer(SR ~ urban + (1|Ecoregion), data=SR_dat)
summary(SR_mod)
anova(SR_mod)

SR_eco_models <- SR_results %>%
  dplyr::filter(sample_size==100) %>%
  group_by(Ecoregion) %>%
  nest() %>%
  mutate(mod=map(data, ~lm(SR ~ urban, data=.x))) %>%
  mutate(tidied=map(mod, tidy),
         glanced=map(mod, glance),
         augmented=map(mod, augment))

SR_mod_summary <- SR_eco_models %>%
  unnest(tidied) %>%
  dplyr::select(1, 4:8)

SR_basic_summary <- SR_dat %>%
  group_by(Ecoregion, urban) %>%
  summarize(mean=mean(SR),
            sd=sd(SR),
            N=n()) %>%
  mutate(se=sd/sqrt(N))



SD_dat <- SD_results %>%
  dplyr::filter(sample_size==100)
SD_mod <- lmer(SD ~ urban + (1|Ecoregion), data=SD_dat)
summary(SD_mod)
anova(SD_mod)


SD_eco_models <- SD_results %>%
  dplyr::filter(sample_size==100) %>%
  group_by(Ecoregion) %>%
  nest() %>%
  mutate(mod=map(data, ~lm(SD ~ urban, data=.x))) %>%
  mutate(tidied=map(mod, tidy),
         glanced=map(mod, glance),
         augmented=map(mod, augment))

SD_mod_summary <- SD_eco_models %>%
  unnest(tidied) %>%
  dplyr::select(1, 4:8)

SD_basic_summary <- SD_dat %>%
  group_by(Ecoregion, urban) %>%
  summarize(mean=mean(SD),
            sd=sd(SD),
            N=n()) %>%
  mutate(se=sd/sqrt(N))



PD_dat <- PD_results %>%
  dplyr::filter(sample_size==100)
PD_mod <- lmer(PD ~ urban + (1|Ecoregion), data=PD_dat)
summary(PD_mod)
anova(PD_mod)


PD_eco_models <- PD_results %>%
  dplyr::filter(sample_size==100) %>%
  group_by(Ecoregion) %>%
  nest() %>%
  mutate(mod=map(data, ~lm(PD ~ urban, data=.x))) %>%
  mutate(tidied=map(mod, tidy),
         glanced=map(mod, glance),
         augmented=map(mod, augment))

PD_mod_summary <- PD_eco_models %>%
  unnest(tidied) %>%
  dplyr::select(1, 4:8)

PD_basic_summary <- PD_dat %>%
  group_by(Ecoregion, urban) %>%
  summarize(mean=mean(PD),
            sd=sd(PD),
            N=n()) %>%
  mutate(se=sd/sqrt(N))
