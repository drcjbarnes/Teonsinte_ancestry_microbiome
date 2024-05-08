
#install.packages("pracma")#had to install this for nthroot function


library(dplyr)
library(stringr)
library(funrar) #Converting into relative abundances 
library(vegan)
library(RColorBrewer)
library(ggplot2); packageVersion("ggplot2")######added this
#library(devtools)
#install_github("thomasp85/patchwork")
library(patchwork)


library(msa)
library(ape)
library(ips)
library(phangorn)
library(DECIPHER)
library("seqinr")
library(phytools)
library(picante)




setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Bacteria/Final_bacterial_analyses/")

big_fasta <- read.fasta(file = "/Users/christopherbarnes/Library/CloudStorage/Dropbox/Sequencing_data/Maize_Ruairidh_followup_data/dada_seqs.fa") #FASTA file with all seqs in.
metadata <- read.table("metadata_outlier_removed.txt", header = TRUE,fill = TRUE)

data <- read.table("data.txt", header = TRUE, row.names = 1)
dim(data)
data<-data[,-11]
data<-data[,-83]
data<-data[,-155]#remove outlier columns
dim(data)
data <- as.data.frame(t(data))

Taxonomy <- read.table("Taxonomy.txt", header = TRUE)
Taxonomy <- str_split_fixed(Taxonomy$Taxonomy, ";", 7)
rownames(Taxonomy) <- colnames(data)
colnames(Taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Taxonomy <- as.data.frame(Taxonomy)
dim(Taxonomy)

#NegData <- data[c(which(metadata$SampleNegative == "Negative")),]
#dim(NegData)
#NegMetadata <- metadata[c(which(metadata$Individual == "Negative")),]
#dim(NegMetadata)
#NegASVPersistences <- colSums(NegData > 0)
#contaminants_data <- rel.data[, c(which(NegASVPersistences > 10))]
#write.table(contaminants_data, file = "contaminants_data.txt")
#rel.data <- rel.data[, -c(which(NegASVPersistences > 10))]
#dim(rel.data)
#mean(rowSums(rel.data))
#contaminants_taxonomy <- Taxonomy[c(which(NegASVPersistences > 10)), ]
#write.table(contaminants_taxonomy, file = "contaminants_taxonomy")
#Taxonomy <- Taxonomy[-c(which(NegASVPersistences > 10)), ]
#dim(Taxonomy)
#Taxonomy <- as.data.frame(Taxonomy)

#####There were no negatives that passed the quality checks.

sum(data)

metadata$PLANTID
frequency <- as.data.frame(metadata %>% group_by(PLANTID) %>% tally())
z <- which(frequency$n > "1")
z <- frequency[c(z),1]

data_binary <- (data > 0) # logical, or
data_binary <- (data_binary > 0)*1L # integer 01

data_binary_aggregated <- aggregate(data_binary, 
  by = list(Category = metadata$PLANTID, metadata$GENOTYPE), FUN = sum)
metadata_aggregated <- data_binary_aggregated[,1:2]
colnames(metadata_aggregated) <- c("PlantID", "GENOTYPE")
data_binary_aggregated <- data_binary_aggregated[,-c(1:2)]
data_binary_aggregated <- (data_binary_aggregated > 1) # logical, or
data_binary_aggregated <- (data_binary_aggregated > 0)*1L # integer 01

data_aggregated <- aggregate(data, 
  by = list(Category = metadata$PLANTID, metadata$GENOTYPE), FUN = sum)

data_aggregated <- data_aggregated[,-c(1:2)]
summed_data <- data_aggregated*data_binary_aggregated

#low.read.samples <- which(rowSums(summed_data) < 50)
#summed_data <- summed_data[-c(low.read.samples),]
#dim(summed_data)
#metadata <- metadata[-c(low.read.samples),]
#dim(metadata)

Archaea <-  which(Taxonomy$Kingdom == "Archaea")
mean(rowSums(summed_data[,c(Archaea)]))
Eukaryota <- which(Taxonomy$Kingdom == "Eukaryota")
mean(rowSums(summed_data[,c(Eukaryota)]))
NA.taxa <- which(Taxonomy$Kingdom == "NA")
mean(rowSums(summed_data[, NA.taxa]))
Bacteria.taxa <- which(Taxonomy$Kingdom == "Bacteria")
mean(rowSums(summed_data[, Bacteria.taxa]))

Not.bacteria.taxa <- which(Taxonomy$Kingdom != "Bacteria")
mean(rowSums(summed_data[, Not.bacteria.taxa]))

summed_data <- summed_data[, -c(Not.bacteria.taxa)]
Taxonomy <- Taxonomy[-c(Not.bacteria.taxa), ]
dim(Taxonomy)
big_fasta <- big_fasta[-c(Not.bacteria.taxa)]
length(big_fasta)

plastids <- which(Taxonomy$Order == "Chloroplast" | Taxonomy$Family == "Mitochondria")
summed_data <- summed_data[, -c(plastids)]
dim(summed_data)
Taxonomy <- Taxonomy[-c(plastids), ]
dim(Taxonomy)
big_fasta <- big_fasta[-c(plastids)]
length(big_fasta)

empty_ASVs <- which(colSums(summed_data) == 0)
summed_data <- summed_data[, -c(empty_ASVs)]
dim(summed_data)
Taxonomy <- Taxonomy[-c(empty_ASVs), ]
dim(Taxonomy)
big_fasta <- big_fasta[-c(empty_ASVs)]
length(big_fasta)


#persistences <- colSums(summed_data > 0)
#x <- which(persistences > 5)
#length(x)
#summed_data <- summed_data[, c(x)] #Removes taxa that are all 0s
#dim(summed_data)

min(rowSums(summed_data))
Richness <- rarefy(summed_data, 303, se = FALSE, MARGIN = 1)
Richness <- as.data.frame(Richness)
summary(rowSums(summed_data))
sum(summed_data)
summed_data <- as.matrix(summed_data)
rel_data <- pracma::nthroot(as.matrix(summed_data), 4) #######throws ERROR "Error in nthroot(as.matrix(summed_data), 4) : could not find function "nthroot"" -> added "pracma:"
rel.data <- as.data.frame((make_relative(rel_data)*100))

BrayData <- vegdist(rel.data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
mds <- monoMDS(BrayData) # Performs multidimensional scaling for the data
plot(mds)

mds.Nmds1<-mds$points[,1] #Extracts dimension 1 as a vector
mds.Nmds2<-mds$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd
Nmds <- cbind(mds.Nmds1, mds.Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

colourCount = length(unique(metadata_aggregated$GENOTYPE))
#getPalette = colorRampPalette(brewer.pal(3, "Set1"))

metadata_aggregated$GENOTYPE[11:57] <- "Teosintes"
metadata_aggregated$GENOTYPE <- factor(metadata_aggregated$GENOTYPE, 
  levels = c("B73", "CML312", "Teosintes", "Soil"),
  labels = c("B73", "CML312", "Teosintes", "Soil"))

Fig1c <- ggplot(data = Nmds, aes( y = mds.Nmds2, x = mds.Nmds1, Type="p", colour = metadata_aggregated$GENOTYPE)) + 
  geom_point(size = 3) + labs(x = "Dimension 1", y = "Dimension 2", colour = NULL) + 
  theme(plot.title=element_text(hjust=0)) +
  scale_colour_manual(values = c("red", "orange", "#a9e37c", "brown")) +
  geom_point(shape = 1, size = 3,colour = "black") + ggtitle("C")
Fig1c

Richness_bact <- Richness#rename it for plotting, because the bacterial and fungal plots are going to be plotted in one figure
metadata_aggregated$GENOTYPE[11:57] <- "Teosintes"
Fig1a <- ggplot(data = metadata_aggregated, aes(x = GENOTYPE, y = Richness_bact$Richness, 
  fill = GENOTYPE)) + geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("red", "orange", "#a9e37c", "brown")) +
  labs(x=NULL, y="Richness") + ggtitle("A")
Fig1a #for the figure first run this script before the fungal script

saveRDS(Richness_bact$Richness, file = "bacterial_richness.RDS")
saveRDS(metadata_aggregated, file = "bacterial_metadata_aggregated.RDS")
saveRDS(Nmds, file = "bacterial_NMDS.RDS")

rownames(rel.data) <- metadata_aggregated$PlantID
rownames(Richness) <- metadata_aggregated$PlantID

#saveRDS(Richness, file = "Bacterial_Richness.RDS")
#saveRDS(Taxonomy, file="Final_bacterial_taxonomy.RDS")
#saveRDS(rel.data, file = "Bacterial_relative_data.RDS")
#saveRDS(metadata_aggregated, file = "metadata_aggregated.RDS)

aggregate(Richness, by = list(Category = metadata_aggregated$GENOTYPE), FUN = mean)
aggregate(Richness, by = list(Category = metadata_aggregated$GENOTYPE), FUN = sd)


#saveRDS(big_fasta, "big_fasta.RDS")
#big_fasta <- readRDS("big_fasta.RDS")

library(vegan)
class(rel.data)
metadata$shannons <- diversity(rel_data, index = "shannon")


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

length(sample.order <- (tree$tip.label))
dim(rel.data)
colnames(rel.data) <- str_replace(colnames(rel.data), "SeqVar", "seq")
#length(which(x %in% tree$tip.label))

sample.order <- tree$tip.label
phylo.rel.data <- rel.data[,c(sample.order)]
colnames(phylo.rel.data)


PD <- pd(phylo.rel.data, tree, include.root = FALSE)
PD <- aggregate(PD$PD, by = list(Category = metadata_aggregated$PlantID), FUN = mean)
#saveRDS(PD, file = "PhyloDiv.RDS")
PD <- readRDS(file = "PhyloDiv.RDS")

aggregate(PD, by = list(Category = metadata_aggregated$GENOTYPE), FUN = mean)
aggregate(PD, by = list(Category = metadata_aggregated$GENOTYPE), FUN = sd)










ord <- BrayData
groups <- metadata_aggregated$GENOTYPE
mod <- betadisper(ord, groups)
permutest(mod)
#if the dispersion is different between groups, then examine
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

library(FSA)
kruskal.test(Richness ~ metadata_aggregated$GENOTYPE)
dunnTest(Richness ~ metadata_aggregated$GENOTYPE, data = metadata)












