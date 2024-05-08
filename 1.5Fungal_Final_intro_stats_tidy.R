
library(tidyr)
library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)
#library(taRifx)#Error in library(taRifx) : es gibt kein Paket namens ‘taRifx’

library(vegan)
library(mvabund)
library(class) 
library(lme4)
library(MASS)
library(data.table)
library(funrar) #Converting into relative abundances 

library(ggplot2); theme_set(theme_bw(base_size = 14))
library(ggrepel)
library(RColorBrewer)
library("ggsci")
library(corrplot)
library(grid)
library(gridExtra)
library(ggpubr)

library(ape)
library(phytools)
library(caper)
library(picante)

library(mclust)
library(cluster)

library(mapdata)
library(maps)
library(googleway)
library(ggmap)

#library(dada2)
library(metacoder)
library(taxa)
#install_github("thomasp85/patchwork")
library(patchwork)


setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Fungi/Final_fungal_analyses/")

rel.data <- readRDS(file = "Fungal_relative_data_maria.RDS") #changed to "_maria", #uses the rel.data from the final_fungal:process_triplicates_tidy.R script
dim(rel.data)

#metadata <- read.table("final_metadata.txt", header = TRUE,fill = TRUE)
metadata <- read.table("final_metadata_subspecies_relabelled.txt", header = TRUE,fill = TRUE)#in this file the some of the teosintes have relabeled subspecies metadata
dim(metadata)

fungal_taxonomy <- readRDS(file = "Final_fungal_taxonomy.RDS")
#Richness <- readRDS(file="Fungal_Richness.RDS")#uses the Richness from the final_fungal_process_triplicates_tidy.R script
Richness <- readRDS(file="PhyloDiv.RDS")#uses the Richness from the final_fungal_process_triplicates_tidy.R script

BrayData <- vegdist(rel.data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
mds <- monoMDS(BrayData) # Performs multidimensional scaling for the data
plot(mds)

mds.Nmds1<-mds$points[,1] #Extracts dimension 1 as a vector
mds.Nmds2<-mds$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd
Nmds <- cbind(mds.Nmds1, mds.Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

#colourCount = length(unique(metadata$GENOTYPE))##########  object 'metadata_aggregated' not found
colourCount = length(unique(metadata$GENOTYPE))
#getPalette = colorRampPalette(brewer.pal(3, "Set1"))

metadata$GENOTYPE[11:57] <-"Teosintes"
metadata$GENOTYPE <- factor(metadata$GENOTYPE, 
  levels = c("B73", "CML312", "Teosintes", "Soil"),
  labels = c("B73", "CML312", "Teosintes", "Soil"))


Fig1d <- ggplot(data = Nmds, aes( y = mds.Nmds2, x = mds.Nmds1, Type="p", colour = metadata$GENOTYPE)) + 
  geom_point(size = 4) + labs(x="Dimension 1", y = "Dimension 2") + #geom_point(size = 4, alpha = 0.5) to make dots more  translucent
  theme(plot.title=element_text(hjust=0))+
  scale_color_brewer(palette = "RdYlBu",name=NULL)
Fig1d

#metadata$GENOTYPE_simplified <- metadata$GENOTYPE
#x <- which(metadata$GENOTYPE_simplified == "CML312" | metadata$GENOTYPE_simplified == "B73" |
#             metadata$GENOTYPE_simplified == "Soil")
#metadata$GENOTYPE_simplified[-x] <- "Teosintes"
#metadata$GENOTYPE_simplified

Fig1b<- ggplot(metadata, aes(x = GENOTYPE, y = Richness, 
  fill = GENOTYPE)) + geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = "RdYlBu", name=NULL)+
  labs(x=NULL)
Fig1b


saveRDS(Richness, file = "fungal_richness.RDS")
saveRDS(metadata, file = "fungal_metadata_aggregated.RDS")
saveRDS(Nmds, file = "fungal_NMDS.RDS")



ord <- BrayData
groups <- metadata$GENOTYPE
mod <- betadisper(ord, groups)
permutest(mod)
#if the dispersion is different between groups, then examine
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

library(FSA)
kruskal.test(Richness ~ metadata$GENOTYPE)
dunnTest(Richness ~ metadata$GENOTYPE, data = metadata)





(Fig1a+ Fig1b)/(Fig1c+Fig1d)+plot_annotation(tag_levels = list(c('a - Bacteria', 'b - Fungi','c - Bacteria', 'd - Fungi')))+ plot_layout(guides = "collect")#for the plot first run the bacterial "data"Data_analysis" script

dim(rel.data)
dim(fungal_taxonomy)

AMF_taxa <- which(fungal_taxonomy$Order == "Glomerales")
AMF_data <- rel.data[ , AMF_taxa]
dim(AMF_data)
summary(colSums(AMF_data))


aggregate(Richness, by = list(Category = metadata$GENOTYPE), FUN = mean)
aggregate(Richness, by = list(Category = metadata$GENOTYPE), FUN = sd)

