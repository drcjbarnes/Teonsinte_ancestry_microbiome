library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(patchwork)
library(ggpubr)
#bacterial

setwd("Final_fungal_analyses/")

early_data <- readRDS(file = "early_data.RDS")
dim(early_data)
Taxonomy <- readRDS(file = "Final_fungal_taxonomy.RDS")
dim(Taxonomy)
metadata_early <- readRDS(file = "early_metadata.RDS")
dim(metadata_early)
metadata_standardized <- readRDS(file="Fungal_metadata_standardized.RDS")
dim(metadata_standardized)


#######To look at the changes in relative abundances at family level, can be changed to different taxonomic levels. 
t.rel.data<- as.data.frame(t(early_data)) 
dim(t.rel.data)

t.rel.data.summed <- aggregate(t.rel.data, by = list(Category = Taxonomy$Family), FUN = sum) #This uses taxonomy to sum abundances by family. 
dim(t.rel.data.summed)
rel.data.summed <- as.data.frame(t(t.rel.data.summed))
names(rel.data.summed) <- lapply(rel.data.summed[1, ], as.character)
rel.data.summed <- rel.data.summed[-1,]
rel.data.summed <- mutate_all(rel.data.summed, function(x) as.numeric(as.character(x)))
dim(t.rel.data.summed)

wide.data <- data.frame("Genotype" = early_metadata$GENOTYPE, "Altitude" = metadata_standardized$ALT, rel.data.summed) #This merges families with selected metadata parameters. Allows for barcharts to be made later.

##########################################################
a <- which(colMeans(rel.data.summed) > 1)
rel.data.summed1 <- rel.data.summed[, a] #This is a threshold percentage to remove families that are less than X% on average.
colSums(rel.data.summed1 > 0)


##########################################################
Minor <- (100 - rowSums(rel.data.summed1, na.rm = FALSE, dims = 1)) #This is useful to create a minor group, so bars add up to 100%.
Newbie <- data.frame(rel.data.summed1)
Wide.data.filtered <- data.frame(wide.data[,1:2], Newbie) #Note this is 1:2 because I had 2 metadata parameters. Change to be the number you have.


############Filter data###############
taxa_df <- rel.data.summed1#write for family grouped data in taxa_df for filtering and for correlations
occurrences <- colSums(taxa_df> 0)
taxa_df <- taxa_df[, which(occurrences > 4)] #Then filters based on the threshold value
#taxonomy[which(occurrences > 5), ] #Then filters based on the threshold value


##############Running correlations for whole datasets##################

cor.test.r2 <- list()
cor.test.p <- list()

cor.test(wide.data$Altitude, taxa_df[,1], method = 'spearman')#test correlation for the for loop, to test the unscaled metadata just use "metadata$ALT"

for (i in 1:ncol(taxa_df)){
  cor_test_output <- cor.test(metadata_standardized$ALT, taxa_df[,i], method = 'spearman',)#omits NAs in altiturde metadata, so no worries about them
  cor.test.r2[[i]] <- cbind(cor_test_output$estimate)
  cor.test.p[[i]] <- cbind(cor_test_output$p.value)
}

cor.test.r2 <- unlist(cor.test.r2, use.names = FALSE)
cor.test.p <- unlist(cor.test.p, use.names = FALSE)

cor_test_outputs <- data.frame("R2" = cor.test.r2, "P-value" = cor.test.p)
rownames(cor_test_outputs) <- colnames(taxa_df)

#########
cor_test_outputs$p.adj <- p.adjust(cor_test_outputs$P.value, method = "fdr")#"fdr", "none", "BH", "bonferroni"....
#########

#########Chris addons from here###########
cor_test_outputs_sorted <- cor_test_outputs %>% arrange(P.value)

a <- rownames(cor_test_outputs_sorted)[1:5]
b <- which(colnames(taxa_df) %in% a)
colnames(taxa_df)
top10_data_alt <- taxa_df [ , b]
dim(top10_data_alt)
colnames(top10_data_alt)

wide.data <- data.frame("ALT" = early_metadata$ALT, top10_data_alt)
long_taxa_data <- melt(wide.data, id.vars = c("ALT"), value.name = "Abundance")

my_palette <- brewer.pal(name="Greens", n=9)[5:9]

ggplot(long_taxa_data, aes(ALT, Abundance, colour = variable)) +
  geom_point(size = 3) + 
  labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  facet_wrap(. ~ variable, scales = 'free', nrow = 1) + 
  stat_smooth(geom = "smooth",method="lm", formula = y ~ x, se = FALSE, 
  aes(color=variable), linewidth=2) + 
  scale_colour_manual(values = c(my_palette)) +
  geom_point(shape = 1, size = 3,colour = "black")

#saveRDS(long_taxa_data, file = "long_family_data_fungi.RDS")




#######To look at richness of a family, but can be changed to be various taxonomic levels.

#Convert data to binary.
t.rel.data <- (t.rel.data > 0) # logical, or
t.rel.data <- (t.rel.data > 0)*1L # integer 01

#t.rel.data<- as.data.frame(t(data)) #Transforms your data (samples and columns)
t.rel.data<- as.data.frame(t(early_data))
dim(t.rel.data)

t.rel.data.summed <- aggregate(t.rel.data, by = list(Category = Taxonomy$Family), FUN = sum) #This uses taxonomy to sum abundances by family. 
dim(t.rel.data.summed)
rel.data.summed <- as.data.frame(t(t.rel.data.summed))
names(rel.data.summed) <- lapply(rel.data.summed[1, ], as.character)
rel.data.summed <- rel.data.summed[-1,]
rel.data.summed <- mutate_all(rel.data.summed, function(x) as.numeric(as.character(x)))

wide.data <- data.frame("Genotype" = early_metadata$GENOTYPE, "Altitude" = early_metadata$ALT, rel.data.summed) #This merges families with selected metadata parameters. Allows for barcharts to be made later.

#rel.data.summed1 <- data[,colMeans(rel.data.summed) > 3] #This is a threshold percentage to remove families that are less than X% on average.
#Minor <- (100 - rowSums(rel.data.summed1, na.rm = FALSE, dims = 1)) #This is useful to create a minor group, so bars add up to 100%.
#Newbie <- data.frame(rel.data.summed1)
Newbie <- data.frame(rel.data.summed)

Wide.data.filtered <- data.frame(Wide.data[,1:2], Newbie) #Note this is 1:2 because I had 2 metadata parameters. Change to be the number you have.
long.data <- melt(Wide.data.filtered, id.vars = c("Genotype", "Altitude"),  variable.name = "Family", value.name = "Abundance") #Converts from wide to long format.

colourCount = length(unique(long.data$Family))
getPalette = colorRampPalette(brewer.pal(4, "Accent"))

long.data<-subset(x=long.data,Altitude != "NA")#filter to get only early lines



############Filter data###############
taxa_df<-rel.data.summed1#write for family grouped data in taxa_df for filtering and for correlations
occurrences <- colSums(taxa_df> 0)
taxa_df <- taxa_df[, which(occurrences > 4)] #Then filters based on the threshold value
#taxonomy[which(occurrences > 5), ] #Then filters based on the threshold value

##############Running correlations for whole datasets##################

cor.test.r2 <- list()
cor.test.p <- list()

cor.test(Wide.data$Altitude, taxa_df[,1], method = 'spearman')#test correlation for the for loop, to test the unscaled metadata just use "metadata$ALT"

for (i in 1:ncol(taxa_df)){
  cor_test_output <- cor.test(Wide.data$Altitude, taxa_df[,i], method = 'spearman',)#omits NAs in altiturde metadata, so no worries about them
  cor.test.r2[[i]] <- cbind(cor_test_output$estimate)
  cor.test.p[[i]] <- cbind(cor_test_output$p.value)
}

cor.test.r2 <- unlist(cor.test.r2, use.names = FALSE)
cor.test.p <- unlist(cor.test.p, use.names = FALSE)

cor_test_outputs <- data.frame("R2" = cor.test.r2, "P-value" = cor.test.p)
rownames(cor_test_outputs) <- colnames(taxa_df)

#########
cor_test_outputs$p.adj <- p.adjust(cor_test_outputs$P.value, method = "fdr")#"fdr", "none", "BH", "bonferroni"....
#########

#significant_cor_taxa <- cor_test_outputs[which(cor_test_outputs$P.value < 0.05),] #Then filters based on the threshold value
significant_cor_taxa <- cor_test_outputs[which(cor_test_outputs$p.adj < 0.05),]#take the adjusted p-value
#sig_var_taxa <- rownames(significant_cor_taxa) %in% colnames(taxa_df) #####Give you a data table of only your significant taxa. 
sig_var_taxa_alt <- taxa_df[,which(colnames(taxa_df)%in% rownames(significant_cor_taxa))]#Give you a data table of only your significant taxa.
sig_var_taxa_alt

cor_test_outputs_sorted<-cor_test_outputs %>% arrange(p.adj)
top10_cor_taxa_alt <- cor_test_outputs_sorted[c(1:10),]#takes the p-vlue (not adjusted)
top10_cor_taxa_alt


