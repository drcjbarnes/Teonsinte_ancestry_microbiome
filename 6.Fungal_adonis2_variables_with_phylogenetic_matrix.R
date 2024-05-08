library("geodist")
library(dplyr)
library(stringr)
library(funrar) #Converting into relative abundances 
library(vegan)
library(RColorBrewer)
library(ggplot2); packageVersion("ggplot2")######added this


setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Fungi/Final_fungal_analyses/")

#PERMANOVA on phylogenetic matrix
phylo.matrix <- readRDS(file = "phylo.matrix.RDS")
phylo_metadata <- readRDS(file = "phylo_metadata.RDS")
phylo_data <- readRDS(file = "phylo_data.RDS")
phylo_standardized_meta <- readRDS(file = "phylo_standardized_meta.RDS")


##################################################### RDA soil ##############################################
data_metadata_soil<-phylo_standardized_meta[,c("clay_1m","silt_1m","sand_1m", "pH_H20_top","pH_KCl_top")]#getting the soil metadata to perform a PCA to subsequently filter for "redundant " parameters

PCA_soil<-rda(data_metadata_soil ~ phylo_standardized_meta$clay_1m + phylo_standardized_meta$silt_1m + phylo_standardized_meta$sand_1m +
                phylo_standardized_meta$pH_H20_top + phylo_standardized_meta$pH_KCl_top, scale=FALSE)
summary(PCA_soil)
screeplot(PCA_soil)

PCA_soil_1<-vegan::scores(PCA_soil,choices=1)#finally no error, the vegan package has to be specified, there might be another function "scores" loaded from another package.
PCA_soil_1<-PCA_soil_1$sites#extract only sites from the scores
PCA_soil_2<-vegan::scores(PCA_soil,choices=2)#and from RDA2
PCA_soil_2<-PCA_soil_2$sites
PCA_soil_3<-vegan::scores(PCA_soil,choices=3)
PCA_soil_3<-PCA_soil_3$sites


######################################################## Temperature ##################################################
data_metadata_clim<-phylo_standardized_meta[,c("Ann.Mean.Tmp","Mean.Tmp.Wet.Q","Mean.Tmp.Dry.Q","Mean.Tmp.Wrm.Q","Mean.Tmp.Cld.Q")]#getting the soil metadata to perform a PCA to subsequently filter for "redundant " parameters
head(data_metadata_clim)
PCA<-rda(data_metadata_clim,scale=FALSE)
plot (PCA)
summary(PCA)
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
biplot(PCA, choices = c(1,2), type = c("text", "points"), xlim = c(-5,10)) 
PCA$CA$eig[1:2]

PCA_clim<-rda(data_metadata_clim ~ phylo_standardized_meta$Ann.Mean.Tmp + phylo_standardized_meta$Mean.Tmp.Wet.Q + phylo_standardized_meta$Mean.Tmp.Dry.Q + phylo_standardized_meta$Mean.Tmp.Wrm.Q + phylo_standardized_meta$Mean.Tmp.Cld.Q,scale=FALSE)
summary(PCA_clim)

#scores(PCA_clim,display = 'species')#did not work, error "Error in UseMethod("scores") : no applicable method for 'scores' applied to an object of class "c('rda', 'cca')""
PCA_clim_1<-vegan::scores(PCA_clim,choices=1)#finally no error, the vegan package has to be specified, there might be another function "scores" loaded from another package.
PCA_clim_1<-PCA_clim_1$sites#extract only sites from the scores
PCA_clim_2<-vegan::scores(PCA_clim,choices=2)#and from RDA2
PCA_clim_2<-PCA_clim_2$sites
PCA_clim_2


############################################################# PCNM ###################################################
coordinates<-phylo_standardized_meta[,c("LON","LAT")]#extract only coordinate columns from metadatatable
list(coordinates)

#library("geodist")
#distancematrix <- geodist(c(data_metadata_early$LON[1],data_metadata_early$LAT[1]),c(data_metadata_early$LON[2],data_metadata_early$LAT[2]),measure="haversine")
distancematrix <- geodist(coordinates,coordinates,measure="haversine") #all coordinates against all coordinates distance matrix
PCA_coord<-rda(coordinates ~ distancematrix,scale=FALSE)
plot(PCA_coord)
summary(PCA_coord)

PCNM_coord_1<-vegan::scores(PCA_coord,choices=1)
PCNM_coord_1<-PCNM_coord_1$sites#extract only sites from the scores
PCNM_coord_2<-vegan::scores(PCA_coord,choices=2)#and from RDA2
PCNM_coord_2<-PCNM_coord_2$sites
PCNM_coord_2


########################################### with controlled for mexicana and parviglumis ##################################
adonis2(phylo.matrix ~ phylo_standardized_meta$ALT,group=data_metadata_early$Taxonomia)
phylo_standardized_meta$ALT
adonis2(phylo.matrix ~ PCA_clim_1,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCA_clim_2,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCA_soil_1,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCA_soil_2,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCA_soil_chem_1,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCA_soil_chem_2,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ phylo_standardized_meta$Ann.Prc,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCNM_coord_1,group=data_metadata_early$Taxonomia)
adonis2(phylo.matrix ~ PCNM_coord_2,group=data_metadata_early$Taxonomia)

adonis2(phylo.matrix ~ phylo_standardized_meta$Ann.Mean.Tmp,group=data_metadata_early$Taxonomia)










###############Individual parviglumis and mexicana analyses

early_data <- readRDS(file="early_data.RDS")
metadata_early <- readRDS(file="early_metadata.RDS")
metadata_standardized <- readRDS(file="Fungal_metadata_standardized.RDS")
Taxonomy <- readRDS(file="Final_fungal_taxonomy.RDS")
early_richness <- readRDS(file = "early_richness.RDS")
early_shannons <- readRDS(file = "fungal_early_shannons.RDS")

########Very important to select the sub-species you want to analyse
#a <- which(metadata_early$Taxonomia == "Zea_mays_subsp._parviglumis")
#a <- which(metadata_early$Taxonomia == "Zea_mays_subsp._mexicana")
#a <- which(metadata_early$Taxonomia != "Undetermined_spp.")

early_data <- early_data[a , ]
dim(early_data)
metadata_early <- metadata_early[a , ]
dim(metadata_early)
metadata_standardized <- metadata_standardized[a , ]
dim(metadata_standardized)
early_richness <- early_richness[a]
early_shannons <- early_shannons[a]


####################################################################################################
###########################  with standardised numerical metadata  #################################
####################################################################################################

metadata_numerical <- metadata_early[,c("ALT","pH_H20_top","pH_KCl_top","clay_1m","silt_1m","sand_1m","Ann.Mean.Tmp","Mean.Tmp.Wet.Q","Mean.Tmp.Dry.Q","Mean.Tmp.Wrm.Q","Mean.Tmp.Cld.Q","Ann.Prc","LAT","LON")]
metadata_standardized <- metadata_numerical %>% mutate_all(~(scale(.) %>% as.vector))#scale scales all variables by subtracting the mean and dividing by the standard deviation ("z-score")

data_metadata_soil <- metadata_standardized[,c("clay_1m","silt_1m","sand_1m", "pH_H20_top","pH_KCl_top")]#getting the soil metadata to perform a PCA to subsequently filter for "redundant " parameters
head(data_metadata_soil)
PCA <- rda(data_metadata_soil,scale=FALSE)
plot (PCA)
summary(PCA)
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
biplot(PCA, choices = c(1,2), type = c("text", "points"), xlim = c(-5,10)) 
PCA$CA$eig[1:2]

PCA_soil<-rda(data_metadata_soil ~ metadata_standardized$clay_1m + metadata_standardized$silt_1m + metadata_standardized$sand_1m +
                metadata_standardized$pH_H20_top + metadata_standardized$pH_KCl_top,scale=FALSE)
summary(PCA_soil)
screeplot(PCA_soil)

PCA_soil_1 <- vegan::scores(PCA_soil,choices=1)#finally no error, the vegan package has to be specified, there might be another function "scores" loaded from another package.
PCA_soil_1 <- PCA_soil_1$sites#extract only sites from the scores
PCA_soil_2 <- vegan::scores(PCA_soil,choices=2)#and from RDA2
PCA_soil_2 <- PCA_soil_2$sites
PCA_soil_3 <-vegan::scores(PCA_soil,choices=3)#and from RDA2
PCA_soil_3 <-PCA_soil_3$sites


####################bray curtis multidimensional scaling just in the early data###########
EarlyDataBray<-vegdist(early_data, Type = "bray", binary = FALSE)#calculate distance ,matrix (early lines only)
EarlyDataMDS <- monoMDS(EarlyDataBray) # Performs multidimensional scaling for the data
plot(EarlyDataMDS)

Early_Nmds1 <- EarlyDataMDS$points[,1] #Extracts dimension 1 as a vector
Early_Nmds2 <- EarlyDataMDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

Early_Nmds=cbind(Early_Nmds1, Early_Nmds2) #Makes a new sheet with dimension1 and 2
Early_Nmds<-as.data.frame(Early_Nmds) #Creates a new dataframe with dimensions 1 and 2

######################adonis on the early data################################
adonis2(EarlyDataBray ~ PCA_soil_1)
adonis2(EarlyDataBray ~ PCA_soil_2)
adonis2(EarlyDataBray ~ PCA_soil_3)
##############################################################################
##############################################################################

#PCA of the environment climate variables(metadata) (Mean Temperatures)

data_metadata_clim<-metadata_standardized[,c("Ann.Mean.Tmp","Mean.Tmp.Wet.Q","Mean.Tmp.Dry.Q","Mean.Tmp.Wrm.Q","Mean.Tmp.Cld.Q", "Ann.Prc")]#getting the soil metadata to perform a PCA to subsequently filter for "redundant " parameters
head(data_metadata_clim)
PCA<-rda(data_metadata_clim,scale=FALSE)
plot (PCA)
summary(PCA)
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
biplot(PCA, choices = c(1,2), type = c("text", "points"), xlim = c(-5,10)) 
PCA$CA$eig[1:2]


#prcomp.PCA.soil<-prcomp(data_metadata_clim, center = TRUE,scale. = TRUE)#different method to compute (proper) PCA
#summary(prcomp.PCA.soil)
#scores(prcomp.PCA.soil)

#biplot(prcomp.PCA.soil)
PCA_clim<-rda(data_metadata_clim ~ metadata_standardized$Ann.Mean.Tmp + metadata_standardized$Mean.Tmp.Wet.Q + metadata_standardized$Mean.Tmp.Dry.Q + metadata_standardized$Mean.Tmp.Wrm.Q + metadata_standardized$Mean.Tmp.Cld.Q,scale=FALSE)
summary(PCA_clim)

#scores(PCA_clim,display = 'species')#did not work, error "Error in UseMethod("scores") : no applicable method for 'scores' applied to an object of class "c('rda', 'cca')""
PCA_clim_1<-vegan::scores(PCA_clim,choices=1)#finally no error, the vegan package has to be specified, there might be another function "scores" loaded from another package.
PCA_clim_1<-PCA_clim_1$sites#extract only sites from the scores
PCA_clim_2<-vegan::scores(PCA_clim,choices=2)#and from RDA2
PCA_clim_2<-PCA_clim_2$sites
PCA_clim_2


######################adonis T for the early data################################
adonis2(EarlyDataBray ~ PCA_clim_1)
adonis2(EarlyDataBray ~ PCA_clim_2)

##############################################################################
##############################################################################
coordinates <- metadata_early[,c("LON","LAT")]#extract only coordinate columns from metadatatable
list(coordinates)

#library("geodist")
#distancematrix <- geodist(c(data_metadata_early$LON[1],data_metadata_early$LAT[1]),c(data_metadata_early$LON[2],data_metadata_early$LAT[2]),measure="haversine")
distancematrix <- geodist(coordinates, coordinates,measure="haversine") #all coordinates against all coordinates distance matrix
PCA_coord<-rda(coordinates ~ distancematrix, scale=FALSE)
plot(PCA_coord)
summary(PCA_coord)

PCNM_coord_1 <- vegan::scores(PCA_coord,choices=1)
PCNM_coord_1 <- PCNM_coord_1$sites#extract only sites from the scores
PCNM_coord_2 <- vegan::scores(PCA_coord,choices=2)#and from RDA2
PCNM_coord_2 <- PCNM_coord_2$sites
PCNM_coord_2

adonis2(EarlyDataBray ~ PCNM_coord_1)
adonis2(EarlyDataBray ~ PCNM_coord_2)

#Combined model
adonis2(EarlyDataBray ~ PCA_clim_1 + PCNM_coord_1)

adonis2(EarlyDataBray ~  metadata_standardized$ALT)
adonis2(EarlyDataBray ~  PCA_clim_1)


adonis2(EarlyDataBray ~  metadata_early$ALT)

###Final parviglumis model
#None

###Final mexicana model
adonis2(EarlyDataBray ~  metadata_early$ALT)
kruskal.test(early_richness, metadata_early$ALT)


#####Subspecies effect 
adonis2(EarlyDataBray ~ metadata_early$Taxonomia)
kruskal.test(early_richness, metadata_early$Taxonomia)

#Altitude effect while controlling for taxonomia
adonis2(EarlyDataBray ~  metadata_standardized$ALT, strata = metadata_early$Taxonomia)



cor.test(early_shannons, metadata_standardized$ALT, method = "spearman")
cor.test(early_shannons, metadata_standardized$Ann.Prc, method = "spearman")
cor.test(early_shannons, metadata_standardized$Ann.Mean.Tmp, method = "spearman")
cor.test(early_shannons, PCA_clim_1, method = "spearman")
cor.test(early_shannons, PCA_clim_2, method = "spearman")
cor.test(early_shannons, PCA_soil_1, method = "spearman")
cor.test(early_shannons, PCA_soil_2, method = "spearman")
cor.test(early_shannons, PCA_soil_3, method = "spearman")
cor.test(early_shannons, PCNM_coord_1, method = "spearman")
cor.test(early_shannons, PCNM_coord_2, method = "spearman")



