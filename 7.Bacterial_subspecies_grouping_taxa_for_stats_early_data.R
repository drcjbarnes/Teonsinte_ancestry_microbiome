library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(patchwork)
library(ggpubr)
#bacterial with subspecies 


setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Teosintes_microbiome_ancestry/Sept23_finishing_files/Bacteria/Final_bacterial_analyses/")

early_data <- readRDS(file="early_data.RDS")
metadata_early <- readRDS(file="early_metadata.RDS")
metadata_standardized <- readRDS(file="Bacterial_metadata_standardized.RDS")
Taxonomy <- readRDS(file="Final_bacterial_taxonomy.RDS")
early_shannons <- readRDS(file = "bacterial_early_shannons.RDS")

#a <- which(metadata_early$Taxonomia == "Zea_mays_subsp._parviglumis")
#a <- which(metadata_early$Taxonomia == "Zea_mays_subsp._mexicana")

early_data <- early_data[a , ]
dim(early_data)
metadata_early <- metadata_early[a , ]
dim(metadata_early)
metadata_standardized <- metadata_standardized[a , ]
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

wide.data <- data.frame("Genotype" = metadata_early$GENOTYPE, "Altitude" = metadata_standardized$ALT, rel.data.summed) #This merges families with selected metadata parameters. Allows for barcharts to be made later.


##########################################################
a <- which(colMeans(rel.data.summed) > 1)
rel.data.summed1 <- rel.data.summed[, a]
colSums(rel.data.summed1 > 0)


############Filter data###############
taxa_df <- rel.data.summed1 #write for family grouped data in taxa_df for filtering and for correlations
occurrences <- colSums(taxa_df> 0)
taxa_df <- taxa_df[, which(occurrences > 4)] #Then filters based on the threshold value


##############Running correlations for whole datasets##################

cor.test.r2 <- list()
cor.test.p <- list()

cor.test(wide.data$Altitude, taxa_df[,1], method = 'spearman')#test correlation for the for loop, to test the unscaled metadata just use "metadata$ALT"

for (i in 1:ncol(taxa_df)){
  cor_test_output <- cor.test(wide.data$Altitude, taxa_df[,i], method = 'spearman',)#omits NAs in altiturde metadata, so no worries about them
  cor.test.r2[[i]] <- cbind(cor_test_output$estimate)
  cor.test.p[[i]] <- cbind(cor_test_output$p.value)
}

cor.test.r2 <- unlist(cor.test.r2, use.names = FALSE)
cor.test.p <- unlist(cor.test.p, use.names = FALSE)

cor_test_outputs <- data.frame("R2" = cor.test.r2, "P-value" = cor.test.p)
rownames(cor_test_outputs) <- colnames(taxa_df)

#########
cor_test_outputs$p.adj <- p.adjust(cor_test_outputs$P.value, method = "fdr")#"fdr", "none", "BH", "bonferroni"....
cor_test_outputs_sorted<-cor_test_outputs %>% arrange(P.value)
top10_cor_taxa_alt <- cor_test_outputs_sorted[c(1:10),]#takes the p-vlue (not adjusted)
top10_cor_taxa_alt

top10_cor_taxa_alt[1]<-rownames(top10_cor_taxa_alt)
colnames(top10_cor_taxa_alt)[1]<-"Family"


#########Chris addons from here###########
cor_test_outputs_sorted <- cor_test_outputs %>% arrange(P.value)
cor_test_outputs_sorted <- cor_test_outputs_sorted[-(which(rownames(cor_test_outputs_sorted) == "NA")) , ]

a <- rownames(cor_test_outputs_sorted)[1:5]
b <- which(colnames(taxa_df) %in% a)
colnames(taxa_df)
top10_data_alt <- taxa_df [ , b]
dim(top10_data_alt)
colnames(top10_data_alt)


wide.data <- data.frame("ALT" = metadata_early$ALT, top10_data_alt)
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



#cor_tests_purvi <- cor_test_outputs_sorted
#cor_tests_mexa <- cor_test_outputs_sorted




######Exclude none determined 

early_data <- readRDS(file="early_data.RDS")
metadata_early <- readRDS(file="early_metadata.RDS")
metadata_standardized <- readRDS(file="Bacterial_metadata_standardized.RDS")
Taxonomy <- readRDS(file="Final_bacterial_taxonomy.RDS")

a <- which(metadata_early$Taxonomia == "Undetermined_spp.")
early_data <- early_data[-a , ]
dim(early_data)
metadata_early <- metadata_early[-a , ]
dim(metadata_early)
metadata_standardized <- metadata_standardized[-a , ]
dim(metadata_standardized)

early_shannons <- readRDS(file="bacterial_early_shannons.RDS")
length(early_shannons)
early_shannons <- early_shannons[-a]
length(early_shannons)
kruskal.test(early_shannons, metadata_early$Taxonomia)

t.rel.data<- as.data.frame(t(early_data)) 
dim(t.rel.data)

t.rel.data.summed <- aggregate(t.rel.data, by = list(Category = Taxonomy$Family), FUN = sum) #This uses taxonomy to sum abundances by family. 
dim(t.rel.data.summed)
rel.data.summed <- as.data.frame(t(t.rel.data.summed))
names(rel.data.summed) <- lapply(rel.data.summed[1, ], as.character)
rel.data.summed <- rel.data.summed[-1,]
rel.data.summed <- mutate_all(rel.data.summed, function(x) as.numeric(as.character(x)))
rel.data.summed <- rel.data.summed[ , -c(which(colnames(rel.data.summed) == "NA"))]
rel.data.summed <- rel.data.summed[ , -c(which(colnames(rel.data.summed) == "Unknown_Family"))]
rel.data.summed <- rel.data.summed [ , names(sort(colMeans(rel.data.summed), decreasing = TRUE))[1:5]]

wide.data <- data.frame("Genotype" = metadata_early$GENOTYPE, "ALT" = metadata_standardized$ALT, 
  "Taxonomia" = metadata_early$Taxonomia, rel.data.summed) #This merges families with selected metadata parameters. Allows for barcharts to be made later.

long_taxa_data <- melt(wide.data, id.vars = c("Genotype", "ALT", "Taxonomia"), value.name = "Abundance")

long_taxa_data$Taxonomia

long_taxa_data$Taxonomia <- str_replace(long_taxa_data$Taxonomia, "Zea_mays_subsp._parviglumis", "Z. parviglumis")
long_taxa_data$Taxonomia <- str_replace(long_taxa_data$Taxonomia, "Zea_mays_subsp._mexicana", "Z. mexicana")

my_palette <- brewer.pal(name="Set1", n=4)[3:4]

ggplot(long_taxa_data, aes(ALT, Abundance, colour = Taxonomia)) +
  geom_point(size = 3) + 
  labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  facet_wrap(. ~ variable, scales = 'free', nrow = 1) + 
  stat_smooth(geom = "smooth",method="lm", formula = y ~ x, se = FALSE, 
              aes(color=Taxonomia), linewidth=2) + 
  scale_colour_manual(values = c(my_palette)) +
  geom_point(shape = 1, size = 3,colour = "black")

saveRDS(long_taxa_data, file = "sub_species_family_lineplot_data.RDS")












