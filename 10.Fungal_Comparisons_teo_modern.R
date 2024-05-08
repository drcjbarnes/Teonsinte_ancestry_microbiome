

setwd("Final_fungal_analyses/")

metadata <- read.table("final_metadata_subspecies_relabelled.txt", header = TRUE,fill = TRUE)
Richness <- readRDS(file="Fungal_Richness.RDS")
rel.data <- readRDS(file="Fungal_relative_data_maria.RDS")
Taxonomy <- readRDS(file="Final_fungal_taxonomy.RDS")

metadata <- metadata[order(metadata$PlantID),]
Richness <- Richness[order(rownames(rel.data))] 
rel.data <- rel.data[order(rownames(rel.data)),]

soil_samples <- which(metadata$GENOTYPE == "Soil")

metadata <- metadata[-(soil_samples) , ]
dim(metadata)
Richness <- Richness[-(soil_samples)]
length(Richness)
rel.data <- rel.data[-(soil_samples) , ]
dim(rel.data)

metadata$compar <- cut(metadata$ALT, breaks=c(0,1000,2000,3000),labels=c("0-1000","1000-2000","2000-3000"))
metadata$compar <- c(metadata$GENOTYPE[1:19], metadata$compar[20:66])
metadata$compar <- factor(metadata$compar, levels = c("1", "2", "3", "B73", "CML312"), labels = c("0-1000", "1001-2000", "2001-3000", "B73", "CML312"))

#######To look at the changes in relative abundances at family level, can be changed to different taxonomic levels. 
t.rel.data<- as.data.frame(t(rel.data)) 
dim(t.rel.data)

t.rel.data.summed <- aggregate(t.rel.data, by = list(Category = Taxonomy$Family), FUN = sum) #This uses taxonomy to sum abundances by family. 
dim(t.rel.data.summed)
rel.data.summed <- as.data.frame(t(t.rel.data.summed))
names(rel.data.summed) <- lapply(rel.data.summed[1, ], as.character)
rel.data.summed <- rel.data.summed[-1,]
rel.data.summed <- mutate_all(rel.data.summed, function(x) as.numeric(as.character(x)))
dim(t.rel.data.summed)


##########################################################
a <- which(colMeans(rel.data.summed) > 1)
rel.data.summed1 <- rel.data.summed[, a] #This is a threshold percentage to remove families that are less than X% on average.
colSums(rel.data.summed1 > 0)


##########################################################
Minor <- (100 - rowSums(rel.data.summed1, na.rm = FALSE, dims = 1)) #This is useful to create a minor group, so bars add up to 100%.
Newbie <- data.frame(rel.data.summed1)
#Wide.data.filtered <- data.frame(wide.data[,1:2], Newbie) #Note this is 1:2 because I had 2 metadata parameters. Change to be the number you have.
Wide.data.filtered <- aggregate(Newbie, by = list(Category = metadata$compar), FUN = mean) #This uses taxonomy to sum abundances by family. 
long_taxa_data <- melt(Wide.data.filtered, id.vars = c("Category"), value.name = "Abundance")


#my_palette <- brewer.pal(name="Greens", n=9)[5:9]

FigureS3 <- ggplot(long_taxa_data, aes(x = Category, y = Abundance, fill = variable)) + #FigureS3 <- ggplot(long.data, aes(x = Type, y = Abundance, fill = Family)) + #FigureS3 <- ggplot(long.data_binned_alt, aes(x = new_bin, y = value, fill = variable)) 
  geom_bar(stat = "identity", colour="black") + #facet_grid( ~factor(Genotype,levels=genotype_order), space = "free", scales = "free") +#+ facet_grid( ~ Genotype, space = "free", scales = "free") +
  theme(strip.text.x = element_text(angle = 90)) + theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  xlab("Altitude in m")+ylab("Abundance")
FigureS3


key_families <- c("Diatrypaceae", "Glomeraceae", "Lentitheciaceae", "Mortierellaceae", "Pleosporaceae")
Newbie <- Newbie[ , key_families]
dim(Newbie)

Wide.data.filtered <- data.frame("compar" = metadata$compar, "Altitude" = metadata$ALT, Newbie) #Note this is 1:2 because I had 2 metadata parameters. Change to be the number you have.
long_taxa_data <- melt(Wide.data.filtered, id.vars = c("compar", "Altitude"), value.name = "Abundance")

colourCount = length(unique(long_taxa_data$compar))
getPalette = colorRampPalette(brewer.pal(5, "Greens"))

colourCount = length(unique(long_taxa_data$compar))
getPalette2 = colorRampPalette(brewer.pal(5, "Blues"))

x <- c(getPalette(colourCount)[2:4], getPalette2(colourCount)[4:5])

ggplot(long_taxa_data, aes(x = variable, y = Abundance, 
  fill = compar)) + geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette = "Set1") + ggtitle("B") +
  labs(y = "Rel. Ab. (%)", x = NULL, fill = NULL) + scale_fill_manual(values = c(x))

#saveRDS(long_taxa_data, file = "long_family_data_fungal_maize_teo_comp.RDS")



#aggregate(Wide.data.filtered$Caulobacteraceae, by = list(Category = metadata$compar), FUN = mean)
#aggregate(Wide.data.filtered$Caulobacteraceae, by = list(Category = metadata$compar), FUN = sd)

kruskal.test(Wide.data.filtered$Diatrypaceae, metadata$compar)
kruskal.test(Wide.data.filtered$Glomeraceae, metadata$compar)
kruskal.test(Wide.data.filtered$Lentitheciaceae, metadata$compar)
kruskal.test(Wide.data.filtered$Mortierellaceae, metadata$compar)
kruskal.test(Wide.data.filtered$Pleosporaceae, metadata$compar)





####################bray curtis multidimensional scaling just in the early data###########

library(ecodist)
DataBray<-vegdist(rel.data, Type = "bray", binary = FALSE)#calculate distance ,matrix (early lines only)
pcoaVS <- pco(DataBray, negvals = "zero", dround = 0)

Nmds <-cbind(pcoaVS$vectors[,1], pcoaVS$vectors[,2]) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

scrs <- data.frame(cbind(Nmds, "compar" = metadata$compar))
cent <- aggregate(cbind(Nmds1, Nmds2) ~ metadata$compar, data = Nmds, FUN = mean)
segs <- merge(scrs, setNames(cent, c('compar','oNMDS1','oNMDS2')), by = 'compar', sort = FALSE)

ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
  colour = metadata$compar)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("B") +
  theme(legend.position = "right") +
  #theme(legend.text=element_text(size=8)) +
  #scale_shape_manual(name = "Skin status", values = c(20, 18, 15)) + #theme(legend.position = "none") +
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2)) + 
  geom_point(shape = 1,size = 4,colour = "black") + 
  labs(x = "PCoA1", y = "PCoA2") + scale_color_manual(values = c(x))
#+ facet_grid(metadata$Plant ~ metadata$Site) 





setwd("/Users/christopherbarnes/Library/CloudStorage/Dropbox/Maize_projects/Ruairidh_maize_followup_and_Maria_Bunner/Sept23_finishing_files/Bacteria/Chris_complete_rewrite/")


key_taxa <- rel.data[ , ("SeqVar81")]

aggregate(key_taxa, by = list(Category = metadata$compar), FUN = mean) #This uses taxonomy to sum abundances by family. 

aggplot(rel.data, aes(x = metadata$compar, y = key_taxa, 
                      fill = metadata$compar)) + geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("A") +
  labs(y = "Rel. Ab. (%)", x = NULL, fill = NULL) + scale_fill_manual(values = c(x))



