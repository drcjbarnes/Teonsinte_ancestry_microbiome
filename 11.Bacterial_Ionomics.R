

library(dplyr)
library(stringr)
library(funrar) #Converting into relative abundances 
library(vegan)
library(RColorBrewer)
library(ggplot2); packageVersion("ggplot2")######added this
#library(devtools)
#install_github("thomasp85/patchwork")
library(patchwork)

setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Ruairidh_maize_followup_and_Maria_Bunner/Sept23_finishing_files/Bacteria/Final_bacterial_analyses/")


ion.data <- read.table(file = "Root_DNA_Teosinte_2017.txt")
colnames(ion.data) <- ion.data[1,]
ion.data <- ion.data[-1,]

early_data <- readRDS("early_data.RDS")
early_metadata <- readRDS("early_metadata.RDS")



early_metadata$GENOTYPE
a <- which(early_metadata$GENOTYPE %in%  ion.data$CIMMYTMA)
early_data <- early_data[a, ]
early_metadata <- early_metadata[a, ]

early_metadata$GENOTYPE
ion.data$CIMMYTMA

ions <- ion.data[, 5:25]
ions <- mutate_all(ions, function(x) as.numeric(as.character(x)))
ion_standardised <- ions %>% mutate_all(~(scale(.) %>% as.vector))#scale scales all variables by subtracting the mean and dividing by the standard deviation ("z-score")


PCA <- princomp(ion_standardised)
screeplot(PCA)

PC1 <- PCA$scores[,1]
PC2 <- PCA$scores[,2]
PC3 <- PCA$scores[,3]

EarlyDataBray<-vegdist(early_data, Type = "bray", binary = FALSE) 
EarlyDataMDS <- monoMDS(EarlyDataBray) 
plot(EarlyDataMDS)

Early_Nmds1 <- EarlyDataMDS$points[,1] 
Early_Nmds2 <- EarlyDataMDS$points[,2] 

Early_Nmds=cbind(Early_Nmds1, Early_Nmds2)
Early_Nmds<-as.data.frame(Early_Nmds) 

adonis2(EarlyDataBray ~ PC1)
adonis2(EarlyDataBray ~ PC2)
adonis2(EarlyDataBray ~ PC3)






