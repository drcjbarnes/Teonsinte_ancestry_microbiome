library(vegan)
library(dplyr)
library(stringr)
library(funrar)
library(RColorBrewer)
library(ggplot2); theme_set(theme_bw(base_size = 14))
library(ggrepel)
library("ggsci")
library(corrplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(ape)

setwd("Final_fungal_analyses/")

#integration of phylogenetic matrix of plants into bacterial dataset
#run the data_analysis and adonis scripts before!

phylo.matrix <- readRDS(file = "phylo.matrix.RDS")
Taxonomy <- readRDS(file="Final_fungal_taxonomy.RDS")
phylo.matrix <- readRDS(file = "phylo.matrix.RDS")
phylo_metadata <- readRDS(file = "phylo_metadata.RDS")
phylo_data <- readRDS(file = "phylo_data.RDS")
phylo_standardized_meta <- readRDS(file = "phylo_standardized_meta.RDS")

#a <- which(phylo_metadata$Taxonomia == "Zea_mays_subsp._parviglumis")
#a <- which(phylo_metadata$Taxonomia == "Zea_mays_subsp._mexicana")

phylo_data <- phylo_data[a , ]
dim(phylo_data)
phylo_metadata <- phylo_metadata[a , ]
dim(phylo_metadata)
phylo_standardized_meta <- phylo_standardized_meta[a , ]
dim(phylo_standardized_meta)
phylo.matrix <- phylo.matrix[a , a ]

EarlyDataBray_sorted <- vegdist(phylo_data, Type = "bray", binary = FALSE)#calculate distance ,matrix (early lines only)
EarlyDataMDS <- monoMDS(EarlyDataBray_sorted) # Performs multidimensional scaling for the data
plot(EarlyDataMDS)

phyl_dist_mantel <- mantel(EarlyDataBray_sorted, phylo.matrix) #correlates two matrices. So bray curtis and phylogenetic distance matrix.
phyl_dist_mantel
phyl_dist_protest<-protest(EarlyDataBray_sorted, phylo.matrix) #Does the same, but with Procrustes ordination test.
phyl_dist_protest

#-----------------------------------------------------------PCOA function --------------------------------------------------------
PCOA_phylo<- pcoa(phylo.matrix)
summary(PCOA_phylo)
adonis2(EarlyDataBray_sorted ~ PCOA_phylo$vectors[,1])
adonis2(EarlyDataBray_sorted ~ PCOA_phylo$vectors[,2])
adonis2(EarlyDataBray_sorted ~ PCOA_phylo$vectors[,3])
adonis2(EarlyDataBray_sorted ~ PCOA_phylo$vectors[,4])



