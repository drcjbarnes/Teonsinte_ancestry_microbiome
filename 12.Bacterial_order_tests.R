library(dplyr)
library(stringr)
library(funrar) #Converting into relative abundances 
library(vegan)
library(RColorBrewer)
library(ggplot2); packageVersion("ggplot2")######added this
library(FSA)
#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Bacteria/Final_bacterial_analyses/")

metadata <- read.table("final_metadata_subspecies_relabelled.txt", header = TRUE,fill = TRUE)
#Richness <- readRDS(file="Bacterial_Richness.RDS")
Richness <- readRDS(file="PhyloDiv.RDS")
rel.data <- readRDS(file="Bacterial_relative_data.RDS")
table(metadata$GENOTYPE)
metadata <- metadata[order(metadata$PlantID),]
Richness <- Richness[order(rownames(rel.data))]
rel.data <- rel.data[order(rownames(rel.data)),]
bacterial_shannons_div <- diversity(rel.data, index = 'shannon')

a <- which(metadata$GENOTYPE == "Soil" | metadata$GENOTYPE == "B73" | metadata$GENOTYPE == "CML312")
metadata$group <- metadata$GENOTYPE
metadata$group[-a] <- "Teosinte"


a <- which(metadata$GENOTYPE == "B73")
B73_data <- rel.data[a , ]
B73_richness <- Richness[a]
B73_metadata <- metadata[a , ]
B73_shannons <- bacterial_shannons_div[a]

B73_metadata$B73_order <- c(1, 10, 2, 3, 4, 5, 6, 7, 8, 9)

B73_metadata$B73_order <- factor(B73_order, 
  levels =  c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  labels =  c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))


B73DataBray<-vegdist(B73_data, Type = "bray", binary = FALSE) #calculate distance ,matrix (B73 lines only)
B73DataMDS <- monoMDS(B73DataBray) #Performs multidimensional scaling for the data
plot(B73DataMDS)

B73_Nmds1 <- B73DataMDS$points[,1] #Extracts dimension 1 as a vector
B73_Nmds2 <- B73DataMDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

B73_Nmds=cbind(B73_Nmds1, B73_Nmds2) #Makes a new sheet with dimension1 and 2
B73_Nmds<-as.data.frame(B73_Nmds) #Creates a new dataframe with dimensions 1 and 2
#saveRDS(B73_Nmds, file = "Bacterial_B73_Nmds.RDS")

ggplot(data = B73_Nmds, aes( y = B73_Nmds2, x = B73_Nmds1, Type="p", colour = B73_order)) + 
  geom_point(size = 4) + labs(x="Dimension 1", y = "Dimension 2") + 
  theme(plot.title=element_text(hjust=0)) +
  scale_color_brewer(palette = "RdYlBu",name=NULL)

plot(B73_shannons, B73_metadata$B73_order)








a <- which(metadata$GENOTYPE == "CML312")
CML312_data <- rel.data[a , ]
CML312_richness <- Richness[a]
CML312_metadata <- metadata[a , ]
CML312_shannons <- bacterial_shannons_div[a]

CML312_metadata$CML312_order <- c(1, 10, 2, 4, 5, 6, 7, 8, 9)

CML312_metadata$CML312_order <- factor(CML312_metadata$CML312_order, 
  levels =  c(1, 10, 2, 4, 5, 6, 7, 8, 9),
  labels =  c(1, 10, 2, 4, 5, 6, 7, 8, 9))


CML312DataBray<-vegdist(CML312_data, Type = "bray", binary = FALSE) #calculate distance ,matrix (CML312 lines only)
CML312DataMDS <- monoMDS(CML312DataBray) #Performs multidimensional scaling for the data
plot(CML312DataMDS)

CML312_Nmds1 <- CML312DataMDS$points[,1] #Extracts dimension 1 as a vector
CML312_Nmds2 <- CML312DataMDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

CML312_Nmds=cbind(CML312_Nmds1, CML312_Nmds2) #Makes a new sheet with dimension1 and 2
CML312_Nmds<-as.data.frame(CML312_Nmds) #Creates a new dataframe with dimensions 1 and 2
#saveRDS(CML312_Nmds, file = "Bacterial_CML312_Nmds.RDS")

ggplot(data = CML312_Nmds, aes( y = CML312_Nmds2, x = CML312_Nmds1, Type="p", colour = CML312_order)) + 
  geom_point(size = 4) + labs(x="Dimension 1", y = "Dimension 2") + 
  theme(plot.title=element_text(hjust=0)) +
  scale_color_brewer(palette = "RdYlBu",name=NULL)

plot(CML312_shannons, CML312_metadata$CML312_order)




