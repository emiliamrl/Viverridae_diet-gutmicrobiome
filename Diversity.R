###Diversity analysis###

##Beta diversity or PCoA##
#First load the packages
library(phyloseq)
library(Biostrings)
library(ggplot2)
library("vegan")

#Load the ASV dataset amd format it
relative_gutMB <- read.csv2("ASV_relativeabundance.csv", row.names = 1)
asv <- relative_gutMB[, 3:8363] #Just the ASV data
metadata <- relative_gutMB[, 1:2] #Just the metadata
#Load in the taxonomy information 
taxonomy <- read.csv2("ASVs_taxonomy.csv", row.names = 1) 
taxa <- as.matrix(taxonomy) #Reformat to matrix

#Permanova _Gut mb
library(dplyr)
adonis2(formula = asv~Host+Location, data = relative_gutMB, method = "bray", permutations = 999, by="terms")

#Create a phyloseq object annd transform it
ps <- phyloseq(otu_table(asv, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(taxa))  
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

#Make the PCoA ordination
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

#Fix legend text 
genetta_nam1 <- expression(paste(italic("Genetta"),plain(" spp.")))
civetta_nam1 <- expression(paste(italic("C. civetta      ")))


#Plot the ordination
pcoa_plot <- plot_ordination(ps.prop, ord.nmds.bray, color="Location", shape="Host", axes = 2:3)
pcoa_plot
pcoa_plot2 <- pcoa_plot + 
  scale_color_manual(values = c("#00665d", "#993404", "#dec27d")) + 
  geom_point(size = 3)+
  labs(y= "PCo 3 (7.4%)", x= "PCo 2 (9.5%)")+   #Check the percentage of the original graph to specify
  theme_classic()+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(legend.text = element_text(size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+ 
  theme(legend.position="right")+
  ggtitle(label = "(b)", subtitle="Gut microbiome")+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  stat_ellipse(type="t", show.legend = FALSE, aes(group = Host))+
  scale_shape_discrete(labels=c(civetta_nam1, genetta_nam1))
pcoa_plot2

#Diet
presence_absence <- read.csv2("OTU_presence_absence.csv", row.names = 1)
otu <- presence_absence[, 3:975]
metadata_d <- presence_absence[, 1:2]
taxonomy_d <- read.csv2("OTU_taxonomy.csv", row.names = 1) 
taxa_d <- as.matrix(taxonomy_d) #Reformat to matrix

#Permanova - Diet
adonis2(formula = otu~Host+Location, data = presence_absence, method = "jaccard", permutations = 999, by="terms")

#Create a phyloseq object annd transform it
ps_d <- phyloseq(otu_table(otu, taxa_are_rows=FALSE),
               sample_data(metadata_d),
               tax_table(taxa_d))  
ps.prop_d <- transform_sample_counts(ps_d, function(otu) otu/sum(otu))

#Make the PCoA ordination
ord.nmds.jaccard <- ordinate(ps.prop_d, method="PCoA", distance="jaccard")

#Plot the ordination
pcoa_plot_2 <- plot_ordination(ps.prop_d, ord.nmds.jaccard, color="Location", shape ="Host", axes = 1:2)
pcoa_plot_2
pcoa_plot3 <- pcoa_plot_2 + 
  scale_color_manual(values = c("#00665d", "#993404", "#dec27d")) + 
  geom_point(size = 3)+
  labs(y= "PCo 2 (3.4%)", x= "PCo 1 (4%)")+   #Check the percentage of the original graph to specify
  theme_classic()+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(legend.text = element_text(size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+ 
  theme(legend.position="right")+
  ggtitle(label = "(a)", subtitle="Diet")+
  theme(plot.subtitle = element_text(hjust = 0.5))+
  scale_shape_discrete(labels=c(civetta_nam1, genetta_nam1))+
  stat_ellipse(type="t", show.legend = FALSE, aes(group = Host))
pcoa_plot3

##Alpha diversity or species richness##
#First load the packages
library(magrittr)
library("ggpubr")
library("dplyr")

#Load in the data
#Gut microbiome 
data_gutMB <- read.csv2("ASVs_counts.csv", row.names = 1)
asv_gutMB <- data_gutMB[, 3:8363] #Sort ASVs 
meta_gutMB <- data_gutMB[, 1:2] #Sort metadata
#The diet
data_diet <- read.csv2("OTU_presence_absence.csv", row.names = 1)
otu_diet <- data_diet[, 3:975] #Sort OTUs
meta_diet <- data_diet[, 1:2] #Sort metadata

#Compute the richness
alpha_gutMB <- specnumber(asv_gutMB)  #gut microbiome
alpha_diet <- specnumber(otu_diet)  #diet

#Test if the data is normally distriibuted with a shapiro test
#gut microbiome
barplot(alpha_gutMB, ylab="no. of species", main="Richness across all samples")
shapiro <- shapiro.test(alpha_gutMB)
shapiro$p.value
#diet
barplot(alpha_diet, ylab="no. of species", main="Richness across all samples")
shapiro <- shapiro.test(alpha_diet)
shapiro$p.value
#If p<0.05, data is normally distributed and we need to check for homogeneity of variance. 
#If data is not normally distributed we can use a Kruskal Wallis. 

#Bartlett's test, to check for homogeneity. 
#gut microbiome
bartlett.test(alpha_gutMB, meta_gutMB$Location)
bartlett.test(alpha_gutMB, meta_gutMB$Host)
#diet
bartlett.test(alpha_diet, meta_diet$Location)
bartlett.test(alpha_diet, meta_diet$Host)
#if the test is significant we can continue with the anova, if not we should do Welch one-way test.

#For the Welch one-way test we need to first make a dataframe including richness pr sample and metadata.
#Gut MB
alpha_gutMB_test <- cbind(alpha_gutMB, meta_gutMB)
#Diet
alpha_diet_test <- cbind(alpha_diet, meta_diet)

#Welch one-way test
#GutMB
oneway.test(alpha_gutMB ~ Host, data = alpha_gutMB_test)
#Split the gut MB into the two host species
sorted_gutmb <- alpha_gutMB_test[order(alpha_gutMB_test$Host), ]
civet_gutmb <- filter(sorted_gutmb, Host == "C. civetta")
genet_gutmb <- filter(sorted_gutmb, Host == "Genetta spp.")

oneway.test(alpha_gutMB ~ Location, data = civet_gutmb)
oneway.test(alpha_gutMB ~ Location, data = genet_gutmb)

#Diet - Host - Anova since the homogenity test was significant
anova <- aov(alpha_diet ~ Host, data = alpha_diet_test)
summary(anova)
#Diet - Location
sorted_diet <- alpha_diet_test[order(alpha_diet_test$Host), ]
civet_diet <- filter(sorted_diet, Host == "C. civetta")
genet_diet <- filter(sorted_diet, Host == "Genetta spp.")

aov_civ <- aov(alpha_diet ~ Location, data = civet_diet)
summary(aov_civ)
aov_gen <- aov(alpha_diet ~ Location, data = genet_diet)
summary(aov_gen)

#Now we want to visualize the data through nice boxplots. 
genetta_nam <- expression(paste(italic("Genetta"),plain(" spp.")))
civetta_nam <- expression(paste(italic("C. civetta")))

#GutMB
alpha_gutmb_plot <- ggplot(data = alpha_gutMB_test, aes(x=Host, y=alpha_gutMB)) + 
  geom_boxplot(aes(fill=Location))+
  scale_fill_manual(values=c("#00665d", "#993404", "#dec27d"))+
  scale_x_discrete(labels= c("C. civetta" = civetta_nam, "Genetta spp." = genetta_nam))+
  labs(y= "Gut microbiome richness", x = " ")+
  theme_classic()+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  geom_signif(comparisons = list(c("C. civetta", "Genetta spp.")), 
              map_signif_level=TRUE)+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 14))+
  ylim(0, 650)+ 
  theme(legend.position="right")+
  ggtitle(label = "(d)", subtitle="Gut microbiome")+
  theme(plot.subtitle = element_text(hjust = 0.5))
alpha_gutmb_plot
#Diet
alpha_diet_plot <- ggplot(data = alpha_diet_test, aes(x=Host, y=alpha_diet)) + 
  geom_boxplot(aes(fill=Location))+
  scale_fill_manual(values=c("#00665d", "#993404", "#dec27d"))+
  scale_x_discrete(labels= c("C. civetta" = civetta_nam, "Genetta spp." = genetta_nam))+
  labs(y= "Diet richness", x = " ")+
  theme_classic()+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  geom_signif(comparisons = list(c("C. civetta", "Genetta spp.")), 
              map_signif_level=TRUE)+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 14))+
  ylim(0, 70)+ 
  theme(legend.position="right")+
  ggtitle(label = "(c)", subtitle="Diet")+
  theme(plot.subtitle = element_text(hjust = 0.5))
alpha_diet_plot

#Combine the plots
row1 <- ggarrange(pcoa_plot3, pcoa_plot2, ncol = 2, nrow = 1)
row1

row2 <- ggarrange(alpha_diet_plot, alpha_gutmb_plot, ncol = 2, nrow = 1)
row2 

#Combine the doubble boxplot with the PCoA
plot <- ggarrange(row1, row2, ncol = 1, nrow = 2)
plot

