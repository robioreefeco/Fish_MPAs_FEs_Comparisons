library(mFD)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(ade4)

SpeciesAsb <- read.csv("SpeciesAsb.csv", header=T, row.names = 1)
SpeciesAsb <- as.matrix(SpeciesAsb)
SpeciesAsb
View(SpeciesAsb)


SpeciesTraits <- read.csv("SpeciesTraitsDataframe.csv", header=T, row.names =1)
View(SpeciesTraits)

TraitDescription <- read.csv("TraitDescriptiondataframe.csv")
View(TraitDescription)

# Remove fuzzyweight traits for this example and thus remove lat column:
TraitDescription. <- TraitDescription[ , -3]
TraitDescription.


##############Compute and Interpret Quality of Functional Spaces ##########
####https://cmlmagneville.github.io/mFD/articles/Compute_and_interpret_quality_of_functional_spaces.html ####

# compute trait-based distances:
dist_fish <- mFD::funct.dist(
  sp_tr         = SpeciesTraits,
  tr_cat        = TraitDescription,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

###############Summarize Species x Traits data frame FISHES###############

SpeciesTraits$Size <- as.factor(SpeciesTraits$Size)
SpeciesTraits$Mobility <- as.factor(SpeciesTraits$Mobility)
SpeciesTraits$Activity <- as.factor(SpeciesTraits$Activity)
SpeciesTraits$Schooling <- as.factor(SpeciesTraits$Schooling)
SpeciesTraits$Position <- as.factor(SpeciesTraits$Position)
SpeciesTraits$Diet <- as.factor(SpeciesTraits$Diet)

# Summarize Species x Traits data
mFD::sp.tr.summary(tr_cat = TraitDescription, sp_tr = SpeciesTraits)
View(TraitDescription)
summary(TraitDescription)
summary(SpeciesTraits)
str(TraitDescription)

# Summarize Assemblages Data
mFD::asb.sp.summary(asb_sp_w = SpeciesAsb)




########### How to Deal With Functional Entities ########### 
####### 1. Why Functional Entities (FEs)? #######
### https://cmlmagneville.github.io/mFD/articles/How_to_deal_with_Functional_Entities.html ###

# Size grouped into only 6 categories:
SpeciesTraits[ , "Size"] <- as.character(SpeciesTraits[ , "Size"])

SpeciesTraits[which(SpeciesTraits[ , "Size"] %in% "0-7cm", "Size")] <- "S1"
SpeciesTraits[which(SpeciesTraits[ , "Size"] == "7.1-15cm"), "Size"]  <- "S2"
SpeciesTraits[which(SpeciesTraits[ , "Size"] == "15.1-30cm"), "Size"] <- "S3"
SpeciesTraits[which(SpeciesTraits[ , "Size"] == "30.1-50cm"), "Size"] <- "S4"
SpeciesTraits[which(SpeciesTraits[ , "Size"] == "50.1-80cm"), "Size"] <- "S5"
SpeciesTraits[which(SpeciesTraits[ , "Size"] == ">80cm"), "Size"] <- "S6"

SpeciesTraits[ , "Size"] <- factor(SpeciesTraits[, "Size"], levels = c("S1", "S2", "S3","S4", "S5", "S6"), ordered = TRUE)

#Mobility grouped into only 3 categories:

SpeciesTraits[ , "Mobility"] <- as.character(SpeciesTraits[, "Mobility"])

SpeciesTraits[which(SpeciesTraits[ , "Mobility"] != "Sed"), "Mobility"] <- "Sed"
SpeciesTraits[which(SpeciesTraits[ , "Mobility"] != "Mob"), "Mobility"] <- "Mob"
SpeciesTraits[which(SpeciesTraits[ , "Mobility"] != "VMob"), "Mobility"] <- "VMob"

SpeciesTraits[ , "Mobility"] <- factor(SpeciesTraits[ , "Mobility"], levels = c("Sed", "Mob", "VMob"), ordered = TRUE)

#Activity grouped into only 3 categories:
SpeciesTraits[ , "Activity"] <- as.character(SpeciesTraits[, "Activity"])

SpeciesTraits[which(SpeciesTraits[ , "Activity"] != "Day"), "Activity"] <- "Day"
SpeciesTraits[which(SpeciesTraits[ , "Activity"] != "Both"), "Activity"] <- "Both"
SpeciesTraits[which(SpeciesTraits[ , "Activity"] != "Night"), "Activity"] <- "Night"

SpeciesTraits[ , "Activity"] <- factor(SpeciesTraits[ , "Activity"], levels = c("Day", "Both", "Night"), ordered = TRUE)


#Schooling grouped into only 5 categories:
SpeciesTraits[ , "Schooling"] <- as.character(SpeciesTraits[, "Schooling"])

SpeciesTraits[which(SpeciesTraits[ , "Schooling"] != "Sol"), "Schooling"] <- "Sol"
SpeciesTraits[which(SpeciesTraits[ , "Schooling"] != "Pair"), "Schooling"] <- "Pair"
SpeciesTraits[which(SpeciesTraits[ , "Schooling"] != "SmallG"), "Schooling"] <- "SmallG"
SpeciesTraits[which(SpeciesTraits[ , "Schooling"] != "MedG"), "Schooling"] <- "MedG"
SpeciesTraits[which(SpeciesTraits[ , "Schooling"] != "LargeG"), "Schooling"] <- "LargeG"

SpeciesTraits[ , "Schooling"] <- factor(SpeciesTraits[ , "Schooling"], levels = c("Sol", "Pair", "SmallG","MedG", "LargeG"), ordered = TRUE)


#Position grouped into only 3 categories:
SpeciesTraits[ , "Position"] <- as.character(SpeciesTraits[, "Position"])

SpeciesTraits[which(SpeciesTraits[ , "Position"] != "Bottom"), "Position"] <- "Bottom"
SpeciesTraits[which(SpeciesTraits[ , "Position"] != "Low"), "Position"] <- "Low"
SpeciesTraits[which(SpeciesTraits[ , "Position"] != "High"), "Position"] <- "High"

SpeciesTraits[ , "Position"] <- factor(SpeciesTraits[ , "Position"], levels = c("Bottom", "Low", "High"), ordered = TRUE)


#Diet grouped into only 7 categories:
SpeciesTraits[ , "Diet"] <- as.character(SpeciesTraits[, "Diet"])

SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "HD"), "Diet"] <- "HD"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "HM"), "Diet"] <- "HM"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "IS"), "Diet"] <- "IS"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "IM"), "Diet"] <- "IM"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "PK"), "Diet"] <- "PK"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "FC"), "Diet"] <- "FC"
SpeciesTraits[which(SpeciesTraits[ , "Diet"] != "OM"), "Diet"] <- "OM"

SpeciesTraits[ , "Diet"] <- factor(SpeciesTraits[ , "Diet"], levels = c("HD", "HM", "IS", "IM", "PK", "FC", "OM"), ordered = TRUE)

knitr::kable(head(SpeciesTraits), caption = "Species x traits dataframe based on *Fish* dataset")

####### Species x assemblages dataframe based on fish dataset############
SpeciesAsb
knitr::kable(as.data.frame(SpeciesAsb[1:6, 1:6]), 
             caption = "Species x assemblages dataframe based on *fish* dataset")



######Traits types based on fish & samples dataset ####

# only keep traits 1 - 4:
TraitDescription <- TraitDescription[1:4, ]

knitr::kable(head(TraitDescription), 
             caption = "Traits types based on *Fish & Sites* dataset")



########## Summarize species assemblages: ############
asb_sp_fish_summ <- mFD::asb.sp.summary(SpeciesAsb)
asb_sp_fish_summ
View(asb_sp_fish_summ)

# retrieve species occurrences for the first 3 assemblages (fish videos):
head(asb_sp_fish_summ$asb_sp_occ, 3)

asb_sp_fish_occ <- asb_sp_fish_summ$"asb_sp_occ"
asb_sp_fish_occ
head(asb_sp_fish_summ$asb_sp_occ)



View(asb_sp_fish_summ)


############# Gather Spp into FE's ############

#########FE##########

# Compute gathering species into FEs:

sp_to_fe_fish <- mFD::sp.to.fe(
  sp_tr       = SpeciesTraits, 
  tr_cat      = TraitDescription, 
  fe_nm_type  = "fe_rank", 
  check_input = TRUE)
sp_to_fe_fish #CORRER ESTO SOLAMENTE

View(sp_to_fe_fish)

#a vector containing FEs names:

sp_to_fe_fish$"fe_nm"

#a vector containing for each species, the FE it belongs to:
sp_fe <- sp_to_fe_fish$"sp_fe"
sp_fe

#a data frame containing for FEs, the values of traits for each FE:
fe_tr <- sp_to_fe_fish$"fe_tr"
fe_tr

#a vector containing the number of species per FE:
fe_nb_sp <- sp_to_fe_fish$"fe_nb_sp"
fe_nb_sp

#a detailed list containing vectors or list with supplementary information about FEs:
sp_to_fe_fish$"details_fe"


View(sp_to_fe_fish)
str(sp_to_fe_fish)


View(asb_sp_fish_occ)
str(asb_sp_fish_occ)

sp_to_fe_fish


# Get the occurrence dataframe:
asb_sp_fish_summ <- mFD::asb.sp.summary(asb_sp_w = SpeciesAsb) 
asb_sp_fish_occ <- asb_sp_fish_summ$'asb_sp_occ'
asb_sp_fish_occ

############function FE ############
alpha.fd.fe <- function(asb_sp_occ, sp_to_fe_, ind_nm = c("fred", "fored", 
                                                           "fvuln"), 
                        check_input = TRUE, details_returned = TRUE) {}

# Compute alpha fd indices:
alpha_fd_fe_fish <- mFD::alpha.fd.fe(
  asb_sp_occ       = asb_sp_fish_occ, 
  sp_to_fe         = sp_to_fe_fish,
  ind_nm           = c('fred', 'fored', 'fvuln'),
  check_input      = TRUE, 
  details_returned = TRUE)

# dataframe with indices values for each assemblage:
alpha_fd_fe_fish$"asb_fdfe"
options(max.print = 99999999)

summary(alpha_fd_fe_fish$"asb_fdfe")
str(alpha_fd_fe_fish$"asb_fdfe")
###5. Plot functional indices based on FEs

alpha_fd_fe_fish_plot <- mFD::alpha.fd.fe.plot(
  alpha_fd_fe       = alpha_fd_fe_fish,
  plot_asb_nm       = c("CAAgosto2021T1C1V1"),
  plot_ind_nm       = c("fred", "fored", "fvuln"),
  name_file         = NULL,
  color_fill_fored  = "darkolivegreen2",
  color_line_fred   = "darkolivegreen4",
  color_fill_bar    = "grey80",
  color_fill_fvuln  = "lightcoral",
  color_arrow_fvuln = "indianred4",
  size_line_fred    = 3,
  size_arrow_fvuln  = 1,
  check_input       = TRUE) 
  

alpha_fd_fe_fish_plot 
