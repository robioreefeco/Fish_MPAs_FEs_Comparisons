library(mFD)
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)


#############CA###########

speciestraitsCA <- read.csv("CASpeciesTraits.csv", header = TRUE, row.names = 1)
View(speciestraitsCA)

SpeciesAsbCA <- read.csv("CASpeciesAsb.csv", header = TRUE, row.names = 1)
View(SpeciesAsbCA)

TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)

# Remove fuzzy traits for this example and thus remove lat column:
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.


#Error es factor

speciestraitsCA$Size <- as.factor(speciestraitsCA$Size)
speciestraitsCA$Mobility <- as.factor(speciestraitsCA$Mobility)
speciestraitsCA$Activity <- as.factor(speciestraitsCA$Activity)
speciestraitsCA$Schooling <- as.factor(speciestraitsCA$Schooling)
speciestraitsCA$Position <- as.factor(speciestraitsCA$Position)
speciestraitsCA$Diet <- as.factor(speciestraitsCA$Diet)

# Summarize Species x Traits data
mFD::sp.tr.summary(tr_cat = TraitDescription, sp_tr = speciestraitsCA)

# Summarize Assemblages Data

matrixespeciesCA <- as.matrix(SpeciesAsbCA)
matrixespeciesCA
mFD::asb.sp.summary(asb_sp_w = matrixespeciesCA)



# Summary of the assemblages * species dataframe:
asb_sp_fishCA_summ <- mFD::asb.sp.summary(asb_sp_w = matrixespeciesCA)
asb_sp_fishCA_summ

#transformacion de matriz de biomasa a matriz de ocurrencia peces 2000
asb_sp_fishCA_occ <- asb_sp_fishCA_summ$"asb_sp_occ"
asb_sp_fishCA_occ

# Computing distances between species based on functional traits
sp_dist_fishCA <- mFD::funct.dist(sp_tr         = speciestraitsCA,
                                    tr_cat        = TraitDescription.,
                                    metric        = "gower",
                                    scale_euclid  = "scale_center",
                                    ordinal_var   = "classic",
                                    weight_type   = "equal",
                                    stop_if_NA    = TRUE)

sp_dist_fishCA

#entidades funcionale (FL)

sp_to_fe_fishCA <- mFD::sp.to.fe(
  sp_tr       = speciestraitsCA, 
  tr_cat      = TraitDescription., 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE) 

sp_to_fe_fishCA


#a vector containing FEs names:
sp_to_fe_fishCA$"fe_nm"

#a vector containing for each species, the FE it belongs to:
sp_feCA <- sp_to_fe_fishCA$"sp_fe"
sp_feCA

#a data frame containing for FEs, the values of traits for each FE:
fe_tr <- sp_to_fe_fishCA$"fe_tr"
fe_tr

#a vector containing the number of species per FE:

fe_nb_spCA <- sp_to_fe_fishCA$"fe_nb_sp"
fe_nb_spCA

#a detailed list containing vectors or list with supplementary information about FEs:

sp_to_fe_fishCA$"details_fe"     

#indices de moulliot
#Compute the set of indices based on number of species in Functional Entities
mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishCA_occ, 
  sp_to_fe         = sp_to_fe_fishCA,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE,
  details_returned = TRUE) 

#Illustrate Functional Diversity indices based on Functional Entities

# Compute alpha fd indices
alpha_fd_fe_fishCA <- mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishCA_occ, 
  sp_to_fe         = sp_to_fe_fishCA,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE, 
  details_returned = TRUE)

alpha_fd_fe_fishCA

# compute trait-based distances:
dist_fishCA <- mFD::funct.dist(
  sp_tr         = speciestraitsCA,
  tr_cat        = TraitDescription.,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# sum up the distance matrix:
summary(as.matrix(dist_fishCA))

# use quality.fpscaes function to compute quality metrics:
fspaces_quality_fishCA <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishCA,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")


round(fspaces_quality_fishCA$"quality_fspaces", 3) 

#4.2. Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishCA,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")




mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishCA,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d",
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#5. Test correlation between functional axes and traits

sp_faxes_coord_fishCA <- fspaces_quality_fishCA$"details_fspaces"$"sp_pc_coord"


fishCA_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = speciestraitsCA, 
  sp_faxes_coord = sp_faxes_coord_fishCA[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

fishCA_tr_faxes$"tr_faxes_stat"[which(fishCA_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]


fishCA_tr_faxes$"tr_faxes_plot"


#6. Plot functional space (POLIGONOS)

sp_faxes_coord_fishCA <- fspaces_quality_fishCA$"details_fspaces"$"sp_pc_coord"



big_plotCA <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fishCA[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkorange",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "orange",
  fill_vert       = "orange",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


big_plotCA$patchwork



#Compute a set of alpha functional indices for a set of assemblages


# Compute functional distance 
sp_dist_fishCA <- mFD::funct.dist(sp_tr         = speciestraitsCA,
                                    tr_cat        = TraitDescription.,
                                    metric        = "gower",
                                    scale_euclid  = "scale_center",
                                    ordinal_var   = "classic",
                                    weight_type   = "equal",
                                    stop_if_NA    = TRUE)

# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fishCA <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishCA, 
  maxdim_pcoa         = 10,
  deviation_weighting = 'absolute',
  fdist_scaling       = FALSE,
  fdendro             = 'average')


# Retrieve species coordinates matrix:
SpeciesAsbCA  <- as.matrix(SpeciesAsbCA)
SpeciesAsbCA 

sp_faxes_coord_fishCA <- fspaces_quality_fishCA$details_fspaces$sp_pc_coord

# Compute alpha diversity indices
alpha_fd_indices_fishCA <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fishCA[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = SpeciesAsbCA, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

# Retrieve alpha diversity indices table (INDICES)
fd_ind_values_fishCA <- alpha_fd_indices_fishCA$functional_diversity_indices
fd_ind_values_fishCA


#############PA###########

speciestraitsPA <- read.csv("PASpeciesTraits.csv", header = TRUE, row.names = 1)
View(speciestraitsPA)

SpeciesAsbPA <- read.csv("PASpeciesAsb.csv", header = TRUE, row.names = 1)
View(SpeciesAsbCA)

TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)

# Remove fuzzy traits for this example and thus remove lat column:
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.


#Error es factor

speciestraitsPA$Size <- as.factor(speciestraitsPA$Size)
speciestraitsPA$Mobility <- as.factor(speciestraitsPA$Mobility)
speciestraitsPA$Activity <- as.factor(speciestraitsPA$Activity)
speciestraitsPA$Schooling <- as.factor(speciestraitsPA$Schooling)
speciestraitsPA$Position <- as.factor(speciestraitsPA$Position)
speciestraitsPA$Diet <- as.factor(speciestraitsPA$Diet)

# Summarize Species x Traits data
mFD::sp.tr.summary(tr_cat = TraitDescription, sp_tr = speciestraitsPA)

# Summarize Assemblages Data

matrixespeciesPA <- as.matrix(SpeciesAsbPA)
matrixespeciesPA
mFD::asb.sp.summary(asb_sp_w = matrixespeciesPA)



# Summary of the assemblages * species dataframe:
asb_sp_fishPA_summ <- mFD::asb.sp.summary(asb_sp_w = matrixespeciesPA)
asb_sp_fishPA_summ

#transformacion de matriz de biomasa a matriz de ocurrencia peces 2000
asb_sp_fishPA_occ <- asb_sp_fishPA_summ$"asb_sp_occ"
asb_sp_fishPA_occ

# Computing distances between species based on functional traits
sp_dist_fishPA <- mFD::funct.dist(sp_tr         = speciestraitsPA,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

sp_dist_fishPA

#entidades funcionale (FL)

sp_to_fe_fishPA <- mFD::sp.to.fe(
  sp_tr       = speciestraitsPA, 
  tr_cat      = TraitDescription., 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE) 

sp_to_fe_fishPA


#a vector containing FEs names:
sp_to_fe_fishPA$"fe_nm"

#a vector containing for each species, the FE it belongs to:
sp_fePA <- sp_to_fe_fishPA$"sp_fe"
sp_fePA

#a data frame containing for FEs, the values of traits for each FE:
PAfe_tr <- sp_to_fe_fishPA$"fe_tr"
PAfe_tr

#a vector containing the number of species per FE:

fe_nb_spPA <- sp_to_fe_fishPA$"fe_nb_sp"
fe_nb_spPA

#a detailed list containing vectors or list with supplementary information about FEs:

sp_to_fe_fishPA$"details_fe"     

#indices de moulliot
#Compute the set of indices based on number of species in Functional Entities
mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishPA_occ, 
  sp_to_fe         = sp_to_fe_fishPA,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE,
  details_returned = TRUE) 

#Illustrate Functional Diversity indices based on Functional Entities

# Compute alpha fd indices
alpha_fd_fe_fishPA <- mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishPA_occ, 
  sp_to_fe         = sp_to_fe_fishPA,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE, 
  details_returned = TRUE)

alpha_fd_fe_fishPA

# compute trait-based distances:
dist_fishPA <- mFD::funct.dist(
  sp_tr         = speciestraitsPA,
  tr_cat        = TraitDescription.,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# sum up the distance matrix:
summary(as.matrix(dist_fishPA))

# use quality.fpscaes function to compute quality metrics:
fspaces_quality_fishPA <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishPA,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")


round(fspaces_quality_fishPA$"quality_fspaces", 3) 

#4.2. Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishPA,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")




mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishPA,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d",
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#5. Test correlation between functional axes and traits

sp_faxes_coord_fishPA <- fspaces_quality_fishPA$"details_fspaces"$"sp_pc_coord"


fishPA_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = speciestraitsCA, 
  sp_faxes_coord = sp_faxes_coord_fishCA[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

fishPA_tr_faxes$"tr_faxes_stat"[which(fishPA_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]


fishPA_tr_faxes$"tr_faxes_plot"


#6. Plot functional space (POLIGONOS)

sp_faxes_coord_fishPA <- fspaces_quality_fishPA$"details_fspaces"$"sp_pc_coord"



big_plotPA <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fishPA[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "pink",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "pink",
  fill_vert       = "pink",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


big_plotPA$patchwork



#Compute a set of alpha functional indices for a set of assemblages


# Compute functional distance 
sp_dist_fishPA <- mFD::funct.dist(sp_tr         = speciestraitsPA,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fishPA <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishCA, 
  maxdim_pcoa         = 10,
  deviation_weighting = 'absolute',
  fdist_scaling       = FALSE,
  fdendro             = 'average')


# Retrieve species coordinates matrix:
SpeciesAsbPA  <- as.matrix(SpeciesAsbPA)
SpeciesAsbPA 

sp_faxes_coord_fishPA <- fspaces_quality_fishPA$details_fspaces$sp_pc_coord

# Compute alpha diversity indices
alpha_fd_indices_fishPA <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fishPA[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = SpeciesAsbPA, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

# Retrieve alpha diversity indices table (INDICES)
fd_ind_values_fishPA <- alpha_fd_indices_fishPA$functional_diversity_indices
fd_ind_values_fishPA



#############FR###########

speciestraitsFR <- read.csv("FRSpeciesTraits.csv", header = TRUE, row.names = 1)
View(speciestraitsPA)

SpeciesAsbFR <- read.csv("FRSpeciesAsb.csv", header = TRUE, row.names = 1)
View(SpeciesAsbCA)

TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)

# Remove fuzzy traits for this example and thus remove lat column:
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.


#Error es factor

speciestraitsFR$Size <- as.factor(speciestraitsFR$Size)
speciestraitsFR$Mobility <- as.factor(speciestraitsFR$Mobility)
speciestraitsFR$Activity <- as.factor(speciestraitsFR$Activity)
speciestraitsFR$Schooling <- as.factor(speciestraitsFR$Schooling)
speciestraitsFR$Position <- as.factor(speciestraitsFR$Position)
speciestraitsFR$Diet <- as.factor(speciestraitsFR$Diet)

# Summarize Species x Traits data
mFD::sp.tr.summary(tr_cat = TraitDescription, sp_tr = speciestraitsFR)

# Summarize Assemblages Data

matrixespeciesFR <- as.matrix(SpeciesAsbFR)
matrixespeciesFR
mFD::asb.sp.summary(asb_sp_w = matrixespeciesFR)



# Summary of the assemblages * species dataframe:
asb_sp_fishFR_summ <- mFD::asb.sp.summary(asb_sp_w = matrixespeciesFR)
asb_sp_fishFR_summ

#transformacion de matriz de biomasa a matriz de ocurrencia peces 2000
asb_sp_fishFR_occ <- asb_sp_fishFR_summ$"asb_sp_occ"
asb_sp_fishFR_occ

# Computing distances between species based on functional traits
sp_dist_fishFR <- mFD::funct.dist(sp_tr         = speciestraitsFR,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

sp_dist_fishFR

#entidades funcionale (FL)

sp_to_fe_fishFR <- mFD::sp.to.fe(
  sp_tr       = speciestraitsFR, 
  tr_cat      = TraitDescription., 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE) 

sp_to_fe_fishFR


#a vector containing FEs names:
sp_to_fe_fishFR$"fe_nm"

#a vector containing for each species, the FE it belongs to:
sp_feFR <- sp_to_fe_fishFR$"sp_fe"
sp_feFR

#a data frame containing for FEs, the values of traits for each FE:
FRfe_tr <- sp_to_fe_fishFR$"fe_tr"
FRfe_tr

#a vector containing the number of species per FE:

fe_nb_spFR <- sp_to_fe_fishFR$"fe_nb_sp"
fe_nb_spFR

#a detailed list containing vectors or list with supplementary information about FEs:

sp_to_fe_fishFR$"details_fe"     

#indices de moulliot
#Compute the set of indices based on number of species in Functional Entities
mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishFR_occ, 
  sp_to_fe         = sp_to_fe_fishFR,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE,
  details_returned = TRUE) 

#Illustrate Functional Diversity indices based on Functional Entities

# Compute alpha fd indices
alpha_fd_fe_fishFR <- mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fishFR_occ, 
  sp_to_fe         = sp_to_fe_fishFR,
  ind_nm           = c("fred", "fored", "fvuln"),
  check_input      = TRUE, 
  details_returned = TRUE)

alpha_fd_fe_fishFR

# compute trait-based distances:
dist_fishFR <- mFD::funct.dist(
  sp_tr         = speciestraitsFR,
  tr_cat        = TraitDescription.,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# sum up the distance matrix:
summary(as.matrix(dist_fishFR))

# use quality.fpscaes function to compute quality metrics:
fspaces_quality_fishFR <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishFR,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")


round(fspaces_quality_fishFR$"quality_fspaces", 3) 

#4.2. Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishFR,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")




mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fishFR,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d",
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#5. Test correlation between functional axes and traits

sp_faxes_coord_fishFR <- fspaces_quality_fishFR$"details_fspaces"$"sp_pc_coord"


fishFR_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = speciestraitsFR, 
  sp_faxes_coord = sp_faxes_coord_fishFR[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

fishFR_tr_faxes$"tr_faxes_stat"[which(fishFR_tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]


fishFR_tr_faxes$"tr_faxes_plot"


#6. Plot functional space (POLIGONOS)

sp_faxes_coord_fishFR <- fspaces_quality_fishFR$"details_fspaces"$"sp_pc_coord"



big_plotFR <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fishFR[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkblue",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blue",
  fill_vert       = "blue",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


big_plotFR$patchwork



#Compute a set of alpha functional indices for a set of assemblages


# Compute functional distance 
sp_dist_fishFR <- mFD::funct.dist(sp_tr         = speciestraitsFR,
                                  tr_cat        = TraitDescription.,
                                  metric        = "gower",
                                  scale_euclid  = "scale_center",
                                  ordinal_var   = "classic",
                                  weight_type   = "equal",
                                  stop_if_NA    = TRUE)

# Compute functional spaces quality to retrieve species coordinates matrix:
fspaces_quality_fishFR <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fishFR, 
  maxdim_pcoa         = 10,
  deviation_weighting = 'absolute',
  fdist_scaling       = FALSE,
  fdendro             = 'average')


# Retrieve species coordinates matrix:
SpeciesAsbFR  <- as.matrix(SpeciesAsbFR)
SpeciesAsbFR 

sp_faxes_coord_fishFR <- fspaces_quality_fishFR$details_fspaces$sp_pc_coord

# Compute alpha diversity indices
alpha_fd_indices_fishFR <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_fishFR[, c('PC1', 'PC2', 'PC3', 'PC4')],
  asb_sp_w         = SpeciesAsbFR, 
  ind_vect         = c('fdis', 'fmpd', 'fnnd', 'feve', 'fric', 'fdiv', 
                       'fori', 'fspe'),
  scaling          = TRUE, 
  check_input      = TRUE, 
  details_returned = TRUE)

# Retrieve alpha diversity indices table (INDICES)
fd_ind_values_fishFR <- alpha_fd_indices_fishFR$functional_diversity_indices
fd_ind_values_fishFR


