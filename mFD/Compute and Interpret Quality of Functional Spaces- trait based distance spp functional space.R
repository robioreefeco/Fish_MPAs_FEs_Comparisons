###########Compute and Interpret Quality of Functional Spaces ##########
#https://cmlmagneville.github.io/mFD/articles/Compute_and_interpret_quality_of_functional_spaces.html#tutorials-data
library(mFD)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(ade4)
library(ggthemes)
library(ggpubr)
library(reshape)
library(zoo)
library(gghalves)
library(ggridges)

SpeciesTraits <- read.csv("SpeciesTraitsDataframe.csv", header=T, row.names =1)
View(SpeciesTraits)
SpeciesTraits

# plot the table:
knitr::kable(head(SpeciesTraits, 
             caption = "Species x traits dataframe based on *fish* dataset"))

TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.

CASpeciesAsb <- read.csv("CASpeciesAsb.csv", header=T, row.names = 1)
CASpeciesAsb <- as.matrix(CASpeciesAsb)
CASpeciesAsb

PASpeciesAsb <- read.csv("PASpeciesAsb.csv", header=T, row.names = 1)
PASpeciesAsb <- as.matrix(PASpeciesAsb)
PASpeciesAsb

FRSpeciesAsb <- read.csv("FRSpeciesAsb.csv", header=T, row.names = 1)
FRSpeciesAsb <- as.matrix(FRSpeciesAsb)
FRSpeciesAsb
# Remove fuzzy traits for this example and thus remove lat column:


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

fish_tr_faxes

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
sp_faxes_coord_fish1

mFD::funct.space.plot(
  sp_faxes_coord    = sp_faxes_coord_fish1[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")],
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

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
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

big_plot


big_plot2 <- mFD::funct.space.plot(
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

#OJO CORRER #
# Retrieve species coordinates matrix:
sp_faxes_coord_fish <- fspaces_quality_fish$details_fspaces$sp_pc_coord
sp_faxes_coord_fish

# Compute alpha diversity indices

#All transects for Boxplot Functional Indices







#All Sites normal general analysis
alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = CASpeciesAsb, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

#### The function has two main outputs:
###a data frame gathering indices values in each assemblage (for FIde values,
###there are as many columns as there are axes to the studied functional space).

fd_ind_values_fish <- alpha_fd_indices_fish$functional_diversity_indices
fd_ind_values_fish



#####a details list of data frames and sublists gathering information 
#####such as coordinates of centroids, distances and identity of the nearest neighbour,
#####distances to the centroid, etc. The user does not have to directly use it but 
###it will be useful if FD indices are then plotted. It can be retrieved through:

details_list_fish <- alpha_fd_indices_fish$"details"
details_list_fish

###### FRic ######

#####Then, you can plot functional indices using the mFD::alpha.multidim.plot() 
###function as follows:

#CA
plots_alpha_fishPA <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("PA"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "pink", asb2 = "pink"),
  color_vert               = c(pool = "grey50", asb1 = "pink", asb2 = "pink"),
  fill_sp                  = c(pool = NA, asb1 = "pink", asb2 = "pink"),
  fill_vert                = c(pool = NA, asb1 = "pink", asb2 = "pink"),
  color_ch                 = c(pool = NA, asb1 = "pink", asb2 = "dpink"),
  fill_ch                  = c(pool = "white", asb1 = "pink", asb2 = "pink"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness  PA
plots_alpha_fishPA$"fric"$"patchwork"

plots_alpha_fishPA <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("PA"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "pink", asb2 = "pink"),
  color_vert               = c(pool = "grey50", asb1 = "pink", asb2 = "pink"),
  fill_sp                  = c(pool = NA, asb1 = "pink", asb2 = "pink"),
  fill_vert                = c(pool = NA, asb1 = "pink", asb2 = "pink"),
  color_ch                 = c(pool = NA, asb1 = "pink", asb2 = "pink"),
  fill_ch                  = c(pool = "white", asb1 = "pink", asb2 = "pink"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness between CA and FR
plots_alpha_fishPA$"fric"$"patchwork"

plots_alpha_fishPAFR <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("FR"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "darkpink", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "darkpink", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness between PA and FR
plots_alpha_fishPAFR$"fric"$"patchwork"



######## Functional Beta Diversity: 7.2. Functional beta diversity indices based 
#####on multidimensional space #####

asb_sp_fish_occ <- asb_sp_fish_summ$"asb_sp_occ"
asb_sp_fish_occ

beta_fd_indices_fish <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = asb_sp_fish_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)


beta_fd_indices_fish


#######a dist object with beta indices values for each pair of assemblages:
head(beta_fd_indices_fish$"pairasb_fbd_indices", 10)

####a list containing details such as inputs, vertices of the global pool and 
#of each assemblage and FRic values for each assemblage

beta_fd_indices_fish$"details"

###a vector containing the FRic value for each assemblage retrieved through the details_beta list:
beta_fd_indices_fish$"details"$"asb_FRic"

######a list of vectors containing names of species being vertices of the convex hull for each assemblage retrieved through the details_beta list:
beta_fd_indices_fish$"details"$"asb_vertices"

########## Checar #https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html#functional-beta-diversity-indices-based-on-multidimensional-space



#########Compute Functional alpha-Diversity indices based on Hill Numbers #######

# Compute alpha fd hill indices:
alpha_hill_fish <- alpha.fd.hill(
  asb_sp_w         = SpeciesAsbSitio, 
  sp_dist          = sp_dist_fish, 
  q                = c(0, 1, 2),
  tau              = 'mean', 
  check_input      = TRUE, 
  details_returned = TRUE)

alpha_hill_fish 


# Compute beta functional hill indices:
beta_hill_fish <- beta.fd.hill(
  asb_sp_w         = SpeciesAsbSitio, 
  sp_dist          = sp_dist_fish, 
  q                = c(0,1,2), 
  tau              = 'mean',
  beta_type        = 'Jaccard', 
  check_input      = TRUE, 
  details_returned = TRUE)

beta_hill_fish

# Then use the mFD::dist.to.df function to ease visualizing result:
## for q = 0:
mFD::dist.to.df(list_dist = list(FDq2 = beta_hill_fish$beta_fd_q$q0))

mFD::dist.to.df(list_dist = list(FDq2 = beta_hill_fish$beta_fd_q$q1))

mFD::dist.to.df(list_dist = list(FDq2 = beta_hill_fish$beta_fd_q$q2))

#### PA FD #####

#All Sites normal general analysis sustituir PASpeciesAsb a CASpeciesAsb y FRSpeciesAsb
alpha_fd_indices_fish <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fish[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = FRSpeciesAsb, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

#### The function has two main outputs:
###a data frame gathering indices values in each assemblage (for FIde values,
###there are as many columns as there are axes to the studied functional space).

fd_ind_values_fish <- alpha_fd_indices_fish$functional_diversity_indices
fd_ind_values_fish



#####a details list of data frames and sublists gathering information 
#####such as coordinates of centroids, distances and identity of the nearest neighbour,
#####distances to the centroid, etc. The user does not have to directly use it but 
###it will be useful if FD indices are then plotted. It can be retrieved through:

details_list_fish <- alpha_fd_indices_fish$"details"
details_list_fish

###### FRic ######

#####Then, you can plot functional indices using the mFD::alpha.multidim.plot() 
###function as follows:

#CA
plots_alpha_fishFR <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("FR"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  color_vert               = c(pool = "grey50", asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  fill_sp                  = c(pool = NA, asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  fill_vert                = c(pool = NA, asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  color_ch                 = c(pool = NA, asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  fill_ch                  = c(pool = "white", asb1 = "#5B9BD5", asb2 = "#2E75B6"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness  PA
plots_alpha_fishFR$"fric"$"patchwork"

plots_alpha_fishPA <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("PA"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#FF99FF", asb2 = "#FF66FF"),
  color_vert               = c(pool = "grey50", asb1 = "#FF99FF", asb2 = "#FF66FF"),
  fill_sp                  = c(pool = NA, asb1 = "#FF99FF", asb2 = "#FF66FF"),
  fill_vert                = c(pool = NA, asb1 = "#FF99FF", asb2 = "#FF66FF"),
  color_ch                 = c(pool = NA, asb1 = "#FF99FF", asb2 = "#FF66FF"),
  fill_ch                  = c(pool = "white", asb1 = "#FF99FF", asb2 = "#FF66FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness between CA and FR
plots_alpha_fishPA$"fric"$"patchwork"

plots_alpha_fishPAFR <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_fish,
  plot_asb_nm              = c("FR"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "darkpink", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "darkpink", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "darkpink", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Ploting functional richness between PA and FR
plots_alpha_fishPA$"fric"$"patchwork"
