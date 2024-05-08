#bact
library(dplyr)
library(stringr)

setwd("Final_bacterial_analyses/")

early_data <- readRDS(file="early_data.RDS")
metadata_early <- readRDS(file="early_metadata.RDS")
metadata_standardized <- readRDS(file="Bacterial_metadata_standardized.RDS")
Taxonomy <- readRDS(file="Final_bacterial_taxonomy.RDS")

#a <- which(early_metadata$Taxonomia == "Zea_mays_subsp._parviglumis")
#a <- which(early_metadata$Taxonomia == "Zea_mays_subsp._mexicana")

early_data <- early_data[a , ]
dim(early_data)
metadata_early <- metadata_early[a , ]
dim(metadata_early)
metadata_standardized <- metadata_standardized[a , ]
dim(metadata_standardized)

taxa_df <- early_data

############Filter data###############
occurrences <- colSums(taxa_df> 0) #Counts in how many samples the taxa is present within.
taxa_df <- taxa_df[, which(occurrences > 4)] #Then filters based on the threshold value

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
cor_test_outputs <- cbind(cor_test_outputs, cor_test_outputs_taxonomy)

cor_test_outputs_sorted <- cor_test_outputs %>% arrange(P.value)
top10_cor_taxa_alt <- cor_test_outputs_sorted[c(1:10),]

#cor_tests_purvi <- cor_test_outputs_sorted
#cor_tests_mexa <- cor_test_outputs_sorted


