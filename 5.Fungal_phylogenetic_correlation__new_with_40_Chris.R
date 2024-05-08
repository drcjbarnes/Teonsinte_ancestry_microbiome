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


setwd("Final_fungal_analyses/")

#integration of phylogenetic matrix of plants into bacterial dataset
#run the data_analysis and adonis scripts before!

phylo.matrix <- read.table(file = ("TeosinteAccessions_maxMiss30_DArTSeq_fst_matrix_final_without_spaces_chris.txt"), header=TRUE, row.names = 1) #This is the phylogenetic distance matrix, the rows are in the same order as the columns.
dim(phylo.matrix)
row.names(phylo.matrix)
colnames(phylo.matrix)
#remove outlier
phylo.matrix <- phylo.matrix[-c(which(colnames(phylo.matrix) == "CIMMYTMA9478")), -c(which(rownames(phylo.matrix) == "CIMMYTMA9478"))]
dim(phylo.matrix)

a <- (which(early_metadata$GENOTYPE %in% rownames(phylo.matrix)))
phylo_data <- early_data[a, ]
phylo_metadata <- early_metadata[a, ]
phylo_standardized_meta <- metadata_standardized[a, ]

phylo_metadata$GENOTYPE
phylo.matrix <- phylo.matrix[order(row.names(phylo.matrix)), order(colnames(phylo.matrix))]
row.names(phylo.matrix)
row.names(phylo.matrix)

phylo_data <- phylo_data[order(phylo_metadata$GENOTYPE),]
phylo_standardized_meta <- phylo_standardized_meta[order(phylo_metadata$GENOTYPE),]
phylo_metadata <- phylo_metadata[order(phylo_metadata$GENOTYPE),]


EarlyDataBray_sorted<-vegdist(phylo_data, Type = "bray", binary = FALSE)#calculate distance ,matrix (early lines only)
EarlyDataMDS <- monoMDS(EarlyDataBray_sorted) # Performs multidimensional scaling for the data
plot(EarlyDataMDS)

phylo.matrix[is.na(phylo.matrix)] <- 0 #overwrite the "NA" with 0,because the mantel and procrust test can't handle NA


phyl_dist_mantel <- mantel(EarlyDataBray_sorted, phylo.matrix) #correlates two matrices. So bray curtis and phylogenetic distance matrix.
phyl_dist_mantel
#phyl_dist_protest<-protest(EarlyDataBray_sorted,phylo.matrix,permutations = 100000)
phyl_dist_protest <- protest(EarlyDataBray_sorted,phylo.matrix) #Does the same, but with Procrustes ordination test.
phyl_dist_protest


#-----------------------------------------------------------PCOA function --------------------------------------------------------
PCOA_phylo<- pcoa(phylo.matrix)
summary(PCOA_phylo)
adonis2(EarlyDataBray_sorted~PCOA_phylo$vectors[,1])
adonis2(EarlyDataBray_sorted~PCOA_phylo$vectors[,2])
adonis2(EarlyDataBray_sorted~PCOA_phylo$vectors[,3])
adonis2(EarlyDataBray_sorted~PCOA_phylo$vectors[,4])


#saveRDS(phylo.matrix, file = "phylo.matrix.RDS")
#saveRDS(phylo_metadata, file = "phylo_metadata.RDS")
#saveRDS(phylo_data, file = "phylo_data.RDS")
#saveRDS(phylo_standardized_meta, file = "phylo_standardized_meta.RDS")


