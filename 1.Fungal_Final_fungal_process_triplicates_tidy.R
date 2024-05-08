
library(tidyr)
library(tidyverse)
library(reshape2)
library(dplyr)
library(tibble)
#library(taRifx)

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
#register_google(key = "AIzaSyDJgKddMxuJwi5t6hlU-n60fAdhWd8tv_M")

#library(dada2)
library(metacoder)
library(taxa)



setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Fungi/Final_fungal_analyses/")

fungal_taxonomy <- readRDS(file = "fungal_taxonomy.RDS")
class(fungal_taxonomy)
dim(fungal_taxonomy)

fungal_data <- readRDS(file = "fungal_data.RDS")
fungal_data <- as.data.frame(fungal_data)
dim(fungal_data)

fungal_taxonomy <- readRDS(file = "fungal_taxonomy.RDS")
dim(fungal_taxonomy)

metadata <- read.table("Fungal_metadata.txt", header = TRUE)
dim(metadata)
metadata <- metadata[ which(metadata$SampleID %in% rownames(fungal_data)), ]
dim(metadata)

fungal_data <- fungal_data[ c(order(rownames(fungal_data))), ]
rownames(fungal_data)
metadata$SampleID

fungal_data <- fungal_data[ -c(which(metadata$SampleNegative == "Negative")) , ]
dim(fungal_data)
metadata <- metadata[ -c(which(metadata$SampleNegative == "Negative")) , ]
dim(metadata)


frequency <- as.data.frame(metadata %>% group_by(PlantID) %>% tally())
z <- which(frequency$n > "1")
z <- frequency[c(z),1]

data_binary <- (fungal_data > 0) # logical, or
data_binary <- (data_binary > 0)*1L # integer 01

data_binary_aggregated <- aggregate(data_binary, 
  by = list(Category = metadata$PlantID, metadata$GENOTYPE), FUN = sum)
metadata_aggregated <- data_binary_aggregated[,1:2]
colnames(metadata_aggregated) <- c("PlantID", "GENOTYPE")
data_binary_aggregated <- data_binary_aggregated[,-c(1:2)]
data_binary_aggregated <- (data_binary_aggregated > 1) # logical, or
data_binary_aggregated <- (data_binary_aggregated > 0)*1L # integer 01

data_aggregated <- aggregate(fungal_data, 
  by = list(Category = metadata$PlantID, metadata$GENOTYPE), FUN = sum)
data_aggregated <- data_aggregated[,-c(1:2)]
summed_data <- data_aggregated*data_binary_aggregated

#colnames(summed_data) <- paste("seq", (1:ncol(summed_data)), sep = "_") #Makes the colnames much simpler and shorter!

fungal_taxonomy["kingdom"][fungal_taxonomy["kingdom"] == ""]<-"unknown"#replace the empty fields in the fungal_taxonomy table with "unknown, to filter later out

plastids <- which(fungal_taxonomy == "Chloroplast" | fungal_taxonomy == "Mitochondria"| fungal_taxonomy == "unknown")#this was commented out before, added "| fungal_taxonomy == "unknown, because in table there were a lot of species with empty kingdom 
summed_data <- summed_data[, -c(plastids)]#
dim(summed_data)#
fungal_taxonomy <- fungal_taxonomy[-c(plastids), ]#Error: object 'Taxonomy' not found -> changed to "fungal_taxonomy"
dim(fungal_taxonomy)#

empty_ASVs <- which(colSums(summed_data) == 0)
summed_data <- summed_data[, -c(empty_ASVs)]
dim(summed_data)
fungal_taxonomy <- fungal_taxonomy[-c(empty_ASVs), ]
dim(fungal_taxonomy)

summary(rowSums(summed_data))
sum(summed_data)
#persistences <- colSums(summed_data > 0)
#x <- which(persistences > 1)
#length(x)
#summed_data <- summed_data[, c(x)] #Removes taxa that are all 0s
#dim(summed_data)
#colnames(summed_data[235])


#fungal_taxonomy <- fungal_taxonomy[c(x),]
#dim(fungal_taxonomy)

colnames(summed_data) <- paste("OTU", (1:ncol(summed_data)), sep = "_") #Makes the colnames much simpler and shorter!
rownames(fungal_taxonomy) <- colnames(summed_data)

Richness <- rarefy(summed_data, (min(rowSums(summed_data))), se = FALSE, MARGIN = 1) #Richness will vary based on number of reads. Therefore I calculated richness based on the lowest samples reads.

library(pracma)
summed_data <- as.matrix(summed_data)
rel_data <- nthroot(as.matrix(summed_data), 4)
rel.data <- as.data.frame((make_relative(rel_data)*100))

colnames(fungal_taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

a <- which(metadata_aggregated$GENOTYPE == "CIMMYTMA9478")
rel.data <- rel.data[-(a) , ]
Richness <- Richness[-(a)]
metadata_aggregated <- metadata_aggregated[-(a) , ]

Richness <- as.data.frame(Richness)
rownames(Richness) <- metadata_aggregated$PlantID
rownames(rel.data) <- metadata_aggregated$PlantID

metadata_aggregated$GENOTYPE
rownames(Richness)
rownames(rel.data)
metadata_aggregated$PlantID

Richness <- Richness[order(rownames(rel.data)),]
rel.data <- rel.data[order(rownames(rel.data)),]

#saveRDS(rel.data, file = "Fungal_relative_data_maria.RDS")
######saveRDS(metadata_aggregated, file = "Fungal_metadata_initial.RDS")
#saveRDS(fungal_taxonomy, file = "Final_fungal_taxonomy.RDS")
#saveRDS(Richness, file = "Fungal_Richness.RDS")


#seqs <- colnames(fungal_data)
#otab <- otu_table(fungal_data, taxa_are_rows = FALSE)
#colnames(otab) <- paste0("seq", seq(ncol(otab)))
#otab = t(otab)
#write.table(seqs, "dada_seqs.txt",quote=FALSE)
#write.table(otab, "dada_table.txt",quote=FALSE,sep="\t")

#Run this in Terminal to create a fasta file with all the sequences and OTU in.
#grep -v '^x' dada_seqs.txt | awk '{print ">seq"$1"\n"$2}' > dada_seqs.fa; rm dada_seqs.txt
#big_fasta <- read.fasta(file = "dada_seqs.fa")
#colnames(fungal_data) <- rownames(otab)

#saveRDS(big_fasta, file = "big_fasta.RDS")
readRDS(file = "big_fasta.RDS")

#write.fasta(sequences = big_fasta, names = names(big_fasta), file.out = "big_fasta.fasta")
big_fasta <- readDNAStringSet(file = "big_fasta.fasta") #Reading as fasta file corrects the format for now. Need to convert in script, probsbly with MSA. 
big_fasta_aligned <- msa(big_fasta, method = "ClustalW", type = "dna", order = "input")


alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}


alignment2Fasta(big_fasta_aligned, 'big_fasta_aligned.fasta')
sink() 



tree <- read.newick(file = "big_fasta_aligned_RAxML_Tree.newick")
plot(tree)
colnames(rel.data) <- str_replace(colnames(rel.data), "_", "")


remove_tips <- which(!(tree$tip.label %in% colnames(rel.data)))
#remove_tips <- !(tree$tip.label %in% colnames(rel.data))
tree <- drop.tip(tree, remove_tips)
plot(tree)  
length(tree$tip.label)

length(sample.order <- (tree$tip.label))
colnames(rel.data)

sample.order <- tree$tip.label
phylo.rel.data <- rel.data[,c(sample.order)]
colnames(phylo.rel.data)

PD <- pd(phylo.rel.data, tree, include.root = FALSE)
PD <- aggregate(PD$PD, by = list(Category = metadata_aggregated$PlantID), FUN = mean)
PD <- PD$x 
#saveRDS(PD, file = "PhyloDiv.RDS")
PD <- readRDS("PhyloDiv.RDS")



aggregate(PD, by = list(Category = metadata$GENOTYPE), FUN = mean)
aggregate(PD, by = list(Category = metadata$GENOTYPE), FUN = sd)










