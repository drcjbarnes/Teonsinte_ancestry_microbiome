library(dplyr)
library(stringr)
library(funrar) #Converting into relative abundances 
library(vegan)
library(RColorBrewer)
library(ggplot2); packageVersion("ggplot2")######added this

setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Fungi/Final_fungal_analyses/")

metadata <- read.table("final_metadata_subspecies_relabelled.txt", header = TRUE,fill = TRUE)
#Richness <- readRDS(file="Fungal_Richness.RDS")
Richness <- readRDS(file="PhyloDiv.RDS")
rel.data <- readRDS(file="Fungal_relative_data_maria.RDS")

metadata <- metadata[order(metadata$PlantID),]
Richness <- Richness[order(rownames(rel.data))]
rel.data <- rel.data[order(rownames(rel.data)),]

fungal_shannons_div <- diversity(rel.data, index = 'shannon')
#saveRDS(fungal_shannons_div, file = "fungal_shannons_div.RDS")

a <- which(metadata$Type == "Landrace")
early_data <- rel.data[a , ]
early_richness <- Richness[a]
early_metadata <- metadata[a , ]
early_shannons <- fungal_shannons_div[a]
#saveRDS(early_shannons, file = "fungal_early_shannons.RDS")


a <- which(metadata$GENOTYPE == "Soil" | metadata$GENOTYPE == "B73" | metadata$GENOTYPE == "CML312")
metadata$group <- metadata$GENOTYPE
metadata$group[-a] <- "Teosinte"

aggregate(fungal_shannons_div, by = list(Category = metadata$group), FUN = mean)
aggregate(fungal_shannons_div, by = list(Category = metadata$group), FUN = sd)
kruskal.test(fungal_shannons_div ~ metadata$group)
dunnTest(fungal_shannons_div ~ metadata$group, data = metadata)


####################################################################################################
###########################  with standardised numerical metadata  #################################
####################################################################################################

metadata_numerical <- early_metadata[,c("ALT","pH_H20_top","pH_KCl_top","clay_1m","silt_1m","sand_1m","Ann.Mean.Tmp","Mean.Tmp.Wet.Q","Mean.Tmp.Dry.Q","Mean.Tmp.Wrm.Q","Mean.Tmp.Cld.Q","Ann.Prc","LAT","LON")]
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
#saveRDS(Early_Nmds, file = "Fungal_Early_Nmds.RDS")

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


coordinates<-early_metadata[,c("LON","LAT")]#extract only coordinate columns from metadatatable
list(coordinates)

library("geodist")
#distancematrix <- geodist(c(data_metadata_early$LON[1],data_metadata_early$LAT[1]),c(data_metadata_early$LON[2],data_metadata_early$LAT[2]),measure="haversine")
distancematrix <- geodist(coordinates, coordinates,measure="haversine") #all coordinates against all coordinates distance matrix
PCA_coord<-rda(coordinates ~ distancematrix,scale=FALSE)
plot(PCA_coord)
summary(PCA_coord)

PCNM_coord_1<-vegan::scores(PCA_coord,choices=1)
PCNM_coord_1<-PCNM_coord_1$sites#extract only sites from the scores
PCNM_coord_2<-vegan::scores(PCA_coord,choices=2)#and from RDA2
PCNM_coord_2<-PCNM_coord_2$sites
PCNM_coord_2

##############################################################################
##############################################################################
adonis2(EarlyDataBray ~ PCNM_coord_1)
adonis2(EarlyDataBray ~ PCNM_coord_2)




#Final model
adonis2(EarlyDataBray ~  metadata_standardized$ALT)

#Taxonomia exploration
adonis2(EarlyDataBray ~ data_metadata_early$Taxonomia)
adonis2(EarlyDataBray ~ data_metadata_early$Taxonomia + PCA_clim_1 + PCNM_coord_1)















#################################################################################################################
#########################"correlate" richness against av.T and precipitation (scatterplots)#####################################
#################################################################################################################
#metadata_standardized$Richness

#saveRDS(metadata_standardized, file="Bacterial_metadata_standardized.RDS")
#saveRDS(early_data , file="early_data.RDS")
#saveRDS(early_metadata , file="early_metadata.RDS")
#saveRDS(early_richness, file = "early_richness.RDS")

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





############################################### for subspecies separately #########################################

only_mexicana <- subset(early_metadata, Taxonomia== "Zea_mays_subsp._mexicana")
only_mexicana_data <- only_mexicana[,-(1:56)]
only_mexicana_data <- only_mexicana_data[,(1:(length(only_mexicana_data)-1))]

only_parviglumis <- subset(early_metadata, Taxonomia=="Zea_mays_subsp._parviglumis")
only_parviglumis_data <- only_parviglumis[,-(1:56)]
only_parviglumis_data <- only_parviglumis_data[,(1:(length(only_parviglumis_data)-1))]

EarlyDataBray_mexicana <- vegdist(only_mexicana_data, Type = "bray", binary = FALSE)
EarlyDataBray_parviglumis <- vegdist(only_parviglumis_data, Type = "bray", binary = FALSE)

adonis2(EarlyDataBray_mexicana ~ only_mexicana$ALT,permutations=100000)
adonis2(EarlyDataBray_parviglumis ~ only_parviglumis$ALT,permutations=100000)

#saveRDS(only_mexicana,file="Bacterial_data_metadata_early_mexicana.RDS")
#saveRDS(only_parviglumis,file="Bacterial_data_metadata_early_parviglumis.RDS")



