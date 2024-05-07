#bring in the ASV table
abund_table <- read.csv2("ASVs.csv", row.names = 1, check.names = FALSE)

#Transpose the OTU table
abund_table <- t(abund_table)

#bring in the taxonomy table table
taxa <- read.csv2("taxonomy.csv", row.names = 1, check.names = FALSE)

#Bring the meta table in (This table alredy have the order level proportions of different dietary items found in each sample)
meta_table <- read.csv2("STAT.csv", row.names = 1, check.names = FALSE)
meta_table <- meta_table[rownames(abund_table),]

#Load following packages
library(vegan)
library(RColorBrewer)
library(microbiome)
library(ape)
library(viridis)
library(phyloseq)
library(ggplot2)
library(grid)
library(gridExtra)
library(pairwiseAdonis)
library(picante)
library(adephylo)
library(tidyverse)
library(dplyr)
library(ggordiplots)


OTU1 <- otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX1 <- phyloseq::tax_table(as.matrix(taxa))
SAM <- sample_data(meta_table)

physeq <- merge_phyloseq(phyloseq(OTU1, TAX1), SAM)


###Generate PCOA plots to visualize the data
#PCOA
pc <- ordinate(physeq, "PCoA", "jaccard")
p2 <- plot_ordination(physeq, pc, type="samples", color="Location",shape = "Host",  axes = c(1,2))
p2 + geom_point(size=6)+ scale_color_manual(values=c("#0926b5","#9e280b","#d19a19","#7676e8", "#c2c2f0","#990000","#3838f5","#d19a19","#050705","#006666","#1f2cd1",  "#d19a19")) + theme_classic() +  scale_shape_manual(values=c(16,17)) + stat_ellipse(aes(group =Host)) 

#Download the vectors for cordinates
write.csv2(pc$vectors, file = "Diet_cord.csv")

#Now using the downloaded coordinates (we have also added meta data into the downloaded file: Diet_cord2.csv) we generate a 3d plot
library(scatterplot3d)
Cord <- read.csv2("Diet_cord2.csv", row.names = 1, check.names = FALSE)

#with scatterplot3d
scatterplot3d(Cord[,4:6])
shapes = c(16,17) 
shapes <- shapes[as.factor(Cord$Host)]
scatterplot3d(Cord[,4:6], pch = Cord$Location)

colors <- c("#0926b5","#9e280b","#d19a19")
colors <- colors[as.factor(Cord$Location)]
scatterplot3d(Cord[,4:6], pch = shapes,cex.symbols=2.5, color=colors, angle = 35, grid=TRUE, box=FALSE)


###Conduct PERMANOVAS using Adonis
bDist <- distance(physeq, method = "bray")

AD1 <- adonis2(bDist ~ physeq@sam_data$Host + physeq@sam_data$Location + physeq@sam_data$Host*physeq@sam_data$Location, by = "margin", permutations = 10000)

AD1


### Calculate alpha diversity indexes of microbiomes
Div <- alpha(physeq, index = c("observed",	"chao1", "diversity_inverse_simpson",	"diversity_gini_simpson",	"diversity_shannon","dominance_core_abundance" ))
write.csv2(Div, file = "Diversity.csv")

#### Generate correlations between Diet proportions and microbiome genus level relative abundance usinng trans_env function in microeco package).

#subsample to seperate the two host species
Genet <-subset_samples(physeq, Host == "Genetta_spp")
Civet <-subset_samples(physeq, Host == "C_civetta")

#env correlations
library(microeco)
library(file2meco)

data <- phyloseq2meco(Civet)

#calculate abundances
data$cal_abund()

#corelation matrix
#for verterbrate dietary items
t1 <- trans_env$new(dataset = data, add_data = data$sample_table[, 4:16])
t1$cal_cor (use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$plot_cor()
Genet_order <- c("Hyracoidea",	"Macroscelidea",	"Primates",	"Rodentia",	"Eulipotyphla",	"Anura",	"Coliiformes",	"Galliformes",	"Squamata")
Civet_order <- c("Primates",	"Rodentia","H_Artiodactyla",	"Anura",	"Galliformes",	"Passeriformes",	"Squamata")

t1$plot_cor(filter_feature = c("","*", "**"), text_x_order = Civet_order)
write.csv2(t1$res_cor, file = "Vert_cor_Civet.csv")

#for inverterbrate dietary items
t1 <- trans_env$new(dataset = data, add_data = data$sample_table[, 17:27])
t1$cal_cor (use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$plot_cor()

Genet_order <- c( "Aranae","Blattodea", "Coleoptera",	"Diptera", "Orthoptera",	"Lepidoptera",	"Hemiptera", "Hymenoptera")
Civet_order <- c( "Spirostreptida", "Coleoptera",	"Diptera", "Orthoptera",	"Lepidoptera",	"Hemiptera", "Hymenoptera")

t1$plot_cor(filter_feature = c("","*", "**"), text_x_order = Civet_order)
t1$plot_cor(filter_feature = c("","*", "**"))

write.csv2(t1$res_cor, file = "Invert_cor_Cevet.csv")

#for plants
t1 <- trans_env$new(dataset = data, add_data = data$sample_table[, 28:54])
t1$cal_cor (use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$plot_cor()
t1$plot_cor(filter_feature = c("","*","**"))
write.csv2(t1$res_cor, file = "Plant_cor_Cevet.csv")



