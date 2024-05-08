#Fungi
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)

setwd("Final_fungal_analyses/")

early_data <- readRDS(file="early_data.RDS")
Taxonomy <- readRDS(file="Final_fungal_taxonomy.RDS")
data_metadata_early <- readRDS(file="early_metadata.RDS")
metadata_standardized <- readRDS(file="Bacterial_metadata_standardized.RDS")

############Filter data###############
occurrences <- colSums(early_data > 0) #Counts in how many samples the taxa is present within.
taxa_df <- early_data[, which(occurrences > 4)] #Then filters based on the threshold value
#taxonomy[which(occurrences > 5), ] #Then filters based on the threshold value

##############Running correlations for whole datasets##################

cor.test.r2 <- list()
cor.test.p <- list()

cor.test(metadata_standardized$ALT, taxa_df[,1], method = 'spearman')#test correlation for the for loop, to test the unscaled metadata just use "metadata$ALT"

for (i in 1:ncol(taxa_df)){
  cor_test_output <- cor.test(metadata_standardized$ALT, taxa_df[,i], method = 'spearman',)#omits NAs in altiturde metadata, so no worries about them
  cor.test.r2[[i]] <- cbind(cor_test_output$estimate)
  cor.test.p[[i]] <- cbind(cor_test_output$p.value)
}

cor.test.r2 <- unlist(cor.test.r2, use.names = FALSE)
cor.test.p <- unlist(cor.test.p, use.names = FALSE)

cor_test_outputs <- data.frame("R2" = cor.test.r2, "P-value" = cor.test.p)
rownames(cor_test_outputs) <- colnames(taxa_df)
cor_test_outputs$p.adj <- p.adjust(cor_test_outputs$P.value, method = "fdr")#"fdr", "none", "BH", "bonferroni"....
cor_test_outputs
cor_test_outputs_taxonomy<-Taxonomy[which(rownames(Taxonomy)%in% rownames(cor_test_outputs)),]
cor_test_outputs_taxonomy
#significant_cor_taxa <- cor_test_outputs[which(cor_test_outputs$P.value < 0.05),] #Then filters based on the threshold value
significant_cor_taxa <- cor_test_outputs[which(cor_test_outputs$p.adj < 0.05),]#take the adjusted p-value
#sig_var_taxa <- rownames(significant_cor_taxa) %in% colnames(taxa_df) #####Give you a data table of only your significant taxa. 
sig_var_taxa_alt <- taxa_df[,which(colnames(taxa_df)%in% rownames(significant_cor_taxa))]#Give you a data table of only your significant taxa.
sig_var_taxa_alt

cor_test_outputs_sorted <- cor_test_outputs %>% arrange(p.adj)
top10_cor_taxa_alt <- cor_test_outputs_sorted[c(1:16),]
length(which(cor_test_outputs_sorted$P.value < 0.0222555))

########Chris addons#########

alt_output <- cbind(cor_test_outputs, cor_test_outputs_taxonomy)
alt_output_sorted <- alt_output %>% arrange(P.value)
a <- rownames(alt_output_sorted)[1:16]
b <- which(colnames(taxa_df) %in% a)
colnames(taxa_df)
top10_data_alt <- taxa_df [ , b]
dim(top10_data_alt)
colnames(top10_data_alt)

rownames(Taxonomy) 
b <- which(rownames(Taxonomy) %in% a)
top10_taxonomy_alt <- Taxonomy [b , ]
labels <- paste(top10_taxonomy_alt$Genus, rownames(top10_taxonomy_alt), sep = "_")
colnames(top10_data_alt) <- labels

wide.data <- data.frame("ALT" = early_metadata$ALT, top10_data_alt)
colnames(wide.data)

#Taxonomy[which(rownames(Taxonomy) == "OTU_505"),]

colnames(wide.data)[3] <- "Sordariales_OTU_44"
colnames(wide.data)[4] <- "Chaetomiaceae_OTU_53"
colnames(wide.data)[9] <- "Sporormiaceae_OTU_265"
colnames(wide.data)[10] <- "Phaeosphaeriaceae_OTU_314"
colnames(wide.data)[15] <- "Ceratobasidiaceae_OTU_505"

long_taxa_data <- melt(wide.data, id.vars = c("ALT"), value.name = "Abundance")

colourCount = length(unique(long_taxa_data$variable))
getPalette = colorRampPalette(brewer.pal(16, "Set2"))

ggplot(long_taxa_data, aes(ALT, Abundance, colour = variable)) +
  geom_point(size = 3) + 
  labs(x="Altitude (m)", y = "Rel. Ab.", colour = NULL, shape = NULL) + 
  facet_wrap(. ~ variable, scales = 'free', nrow = 2) + 
  stat_smooth(geom = "smooth",method="lm", formula = y ~ x, se = FALSE, 
  aes(color=variable), linewidth=2) + scale_colour_manual(values = getPalette(colourCount)) +
  geom_point(shape = 1,size = 3,colour = "black")

#saveRDS(long_taxa_data, file = "long_taxa_data_fungi.RDS")


