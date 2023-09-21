###########Compute and Interpret Quality of Functional Spaces ##########
#https://cmlmagneville.github.io/mFD/articles/Compute_and_interpret_quality_of_functional_spaces.html#tutorials-data
library(mFD)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(ade4)

SpeciesTraits <- read.csv("SpeciesTraitsDataframe.csv", header=T, row.names =1)
View(SpeciesTraits)


# plot the table:
knitr::kable(head(SpeciesTraits, 
             caption = "Species x traits dataframe based on *fish* dataset"))


TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)

# Remove fuzzy traits for this example and thus remove lat column:
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.

#Error as factor 

SpeciesTraits$Size <- as.factor(SpeciesTraits$Size)
SpeciesTraits$Mobility <- as.factor(SpeciesTraits$Mobility)
SpeciesTraits$Activity <- as.factor(SpeciesTraits$Activity)
SpeciesTraits$Schooling <- as.factor(SpeciesTraits$Schooling)
SpeciesTraits$Position <- as.factor(SpeciesTraits$Position)
SpeciesTraits$Diet <- as.factor(SpeciesTraits$Diet)

# Compute Functional Distance
sp_dist_fish <- mFD::funct.dist(sp_tr         = SpeciesTraits,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

sp_dist_fish

# sum up the distance matrix:
summary(as.matrix(sp_dist_fish))
#The Gower distances range from < 0.01 to 0.790. For instance, Gower distances between blackberry and 3 other fruits are:
  

#3. Compute functional space, quality metrics and plot them
#3.1. Compute functional spaces and associated quality metrics
#We now compute varying number of functional space from 1 to 9 dimensions based on a 
#PCoA as well as an UPGMA dendrogram using mFD::quality.fspaces() function. 
#We also compute 4 quality metrics (= all combinations of deviation weighting and distance scaling) 
#(details: mFD General Workflow tutorial, step 4.1).


# Compute Functional Spaces Quality (to retrieve species coordinates)
fspaces_quality_fish <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")
fspaces_quality_fish



#a data frame gathering for each space (in rows), values of quality metric(s) (in columns)
round(fspaces_quality_fish$"quality_fspaces", 3)            # Quality metrics of spaces




#Lists with details required for other tasks in step 4 to plot functional space quality and in step 5 to plot functional space.
#NOTE The space with the best quality has the lowest 
#quality metric. Here, thanks to mad values, we can see that the 
#4D space is the best one. That is why the following of this tutorial
# will use this multidimensional space.

#4.2. Illustrating the quality of the selected functional spaces#

mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fish,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average","pcoa_1d", "pcoa_2d", "pcoa_3d",
                                 "pcoa_4d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#####5. Test correlation between functional axes and traits
sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

fish_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = SpeciesTraits, 
  sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# Print traits with significant effect:
fish_tr_faxes$"tr_faxes_stat"[which(fish_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# Return plots:
fish_tr_faxes$"tr_faxes_plot"

sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"
sp_faxes_coord_fish

#We can thus see that PC1 is mostly driven by Climate (temperate on the left and tropical
#on the right) and Plant Type (forb & shrub on the left vs tree & vine on the right) and 
#Size (large fruits on the right) with weaker influence of Seed (eta2 < 0.25). Then, PC2 
#is mostly driven by Seed (no seed on the left and pit seed on the right) with weaker infl-
#uence of Plant Type. PC3 is driven by only one trait, Size. And finally PC4 is mostly driven by Sugar 
#(high sugar content on the right and low sugar content on the left) with a weaker influence of Plant Type.

#6 plot functional space

sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fish,
  faxes           = NULL,
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot$patchwork

##7. Compute functional diversity indices & plot them ####
####7.1. Functional alpha diversity indices in a multidimensional space ####

sp_dist_fish <- mFD::funct.dist(sp_tr         = SpeciesTraits,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fish <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish, 
  maxdim_pcoa         = 10,
  deviation_weighting = 'absolute',
  fdist_scaling       = FALSE,
  fdendro             = 'average')
############################ Plot species position in a functional space Individuals graphs ####
sp_dist_fish1 <- mFD::funct.dist(sp_tr         = SpeciesTraits,
                                  tr_cat         = TraitDescription.,
                                  metric         = "gower",
                                  scale_euclid   = "scale_center",
                                  ordinal_var    = "classic",
                                  weight_type    = "equal",
                                  stop_if_NA     = TRUE)

fspaces_quality_fish1 <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")


sp_faxes_coord_fish1 <- fspaces_quality_fish1$details_fspaces$sp_pc_coord

mFD::funct.space.plot(
  sp_faxes_coord    = sp_faxes_coord_fish1[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")],
  faxes             = NULL,
  name_file         = NULL,
  faxes_nm          = NULL,
  range_faxes       = c(NA, NA),
  color_bg          = "grey95",
  color_pool          = "darkturquoise",
  fill_pool           = "white",
  shape_pool          = 21,
  size_pool           = 1,
  plot_ch           = TRUE,
  color_ch          = "darkblue",
  fill_ch           = "white",
  alpha_ch          = 1,
  plot_vertices     = TRUE,
  color_vert        = "darkturquoise",
  fill_vert         = "darkturquoise",
  shape_vert        = 22,
  size_vert         = 1,
  plot_sp_nm         = NULL,
  nm_size            = 3,
  nm_color           = "black",
  nm_fontface        = "plain",
  check_input        = TRUE)

#OJO CORRER #
# Retrieve species coordinates matrix:
sp_faxes_coord_fish <- fspaces_quality_fish$details_fspaces$sp_pc_coord
sp_faxes_coord_fish

#OJO código de especies 
# Compute alpha diversity indices

SpeciesAsbSitio <- read.csv("Suma x Sitio (Sum1) .csv", header=T, row.names = 1)
SpeciesAsbSitio <- as.matrix(SpeciesAsbSitio)
SpeciesAsbSitio
View(SpeciesAsbSitio)


alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = SpeciesAsbSitioTodo, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

fd_ind_values_fish <- alpha_fd_indices_fish$functional_diversity_indices
fd_ind_values_fish

