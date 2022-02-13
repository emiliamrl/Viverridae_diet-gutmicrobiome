#Species Accumulation curve 
#Load in the data 
ASVs_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/ASVs_counts_civetta.csv", row.names=1, sep=";")
ASVs_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/ASVs_counts_genetta.csv", row.names=1, sep=";")
coleop_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Coleop_Civetta.csv", row.names=1, sep=";")
coleop_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Coleop_Genetta.csv", row.names=1, sep=";")
mam_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Mam16S_Civetta.csv", row.names=1, sep=";")
mam_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Mam16S_Genetta.csv", row.names=1, sep=";")
trac_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Trac01_Civetta.csv", row.names=1, sep=";")
trac_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Trac01_Genetta.csv", row.names=1, sep=";")
vert_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Vert01_Civetta.csv", row.names=1, sep=";")
vert_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Vert01_Genetta.csv", row.names=1, sep=";")
zeale_civ <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Zeale_Civetta.csv", row.names=1, sep=";")
zeale_gen <- read.csv("~/Desktop/MSc_article/Species_accumulation_curve/OTU_97_Zeale_Genetta.csv", row.names=1, sep=";")

#Create the species accumulation curve
library("vegan")
acc_civ_16S <- specaccum(ASVs_civ, method = "exact", permutations = 100)
acc_gen_16S <- specaccum(ASVs_gen, method = "exact", permutations = 100)
acc_civ_coleop <- specaccum(coleop_civ, method = "exact", permutations = 100)
acc_gen_coleop <- specaccum(coleop_gen, method = "exact", permutations = 100)
acc_civ_mam <- specaccum(mam_civ, method = "exact", permutations = 100)
acc_gen_mam <- specaccum(mam_gen, method = "exact", permutations = 100)
acc_civ_trac <- specaccum(trac_civ, method = "exact", permutations = 100)
acc_gen_trac <- specaccum(trac_gen, method = "exact", permutations = 100)
acc_civ_vert <- specaccum(vert_civ, method = "exact", permutations = 100)
acc_gen_vert <- specaccum(vert_gen, method = "exact", permutations = 100)
acc_civ_zeale <- specaccum(zeale_civ, method = "exact", permutations = 100)
acc_gen_zeale <- specaccum(zeale_gen, method = "exact", permutations = 100)

data_civ_16S <- data.frame(Samples=acc_civ_16S$sites, Richness=acc_civ_16S$richness, SD=acc_civ_16S$sd)
data_gen_16S <- data.frame(Samples=acc_gen_16S$sites, Richness=acc_gen_16S$richness, SD=acc_gen_16S$sd)
data_civ_coleop <- data.frame(Samples=acc_civ_coleop$sites, Richness=acc_civ_coleop$richness, SD=acc_civ_coleop$sd)
data_gen_coleop <- data.frame(Samples=acc_gen_coleop$sites, Richness=acc_gen_coleop$richness, SD=acc_gen_coleop$sd)
data_civ_mam <- data.frame(Samples=acc_civ_mam$sites, Richness=acc_civ_mam$richness, SD=acc_civ_mam$sd)
data_gen_mam <- data.frame(Samples=acc_gen_mam$sites, Richness=acc_gen_mam$richness, SD=acc_gen_mam$sd)
data_civ_trac <- data.frame(Samples=acc_civ_trac$sites, Richness=acc_civ_trac$richness, SD=acc_civ_trac$sd)
data_gen_trac <- data.frame(Samples=acc_gen_trac$sites, Richness=acc_gen_trac$richness, SD=acc_gen_trac$sd)
data_civ_vert <- data.frame(Samples=acc_civ_vert$sites, Richness=acc_civ_vert$richness, SD=acc_civ_vert$sd)
data_gen_vert <- data.frame(Samples=acc_gen_vert$sites, Richness=acc_gen_vert$richness, SD=acc_gen_vert$sd)
data_civ_zeale <- data.frame(Samples=acc_civ_zeale$sites, Richness=acc_civ_zeale$richness, SD=acc_civ_zeale$sd)
data_gen_zeale <- data.frame(Samples=acc_gen_zeale$sites, Richness=acc_gen_zeale$richness, SD=acc_gen_zeale$sd)


#Plot with ggplot
library(ggplot2)

civetta_16S <- ggplot() +
  geom_point(data=data_civ_16S, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_16S, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_16S ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "ASVs")+
  ggtitle("A")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_16S

genetta_16S <- ggplot() +
  geom_point(data=data_gen_16S, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_16S, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_16S ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "ASVs")+
  ggtitle("B")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_16S

civetta_coleop <- ggplot() +
  geom_point(data=data_civ_coleop, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_coleop, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_coleop ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("C")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_coleop

genetta_coleop <- ggplot() +
  geom_point(data=data_gen_coleop, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_coleop, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_coleop ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("D")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_coleop

civetta_mam <- ggplot() +
  geom_point(data=data_civ_mam, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_mam, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_mam ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("E")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_mam

genetta_mam <- ggplot() +
  geom_point(data=data_gen_mam, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_mam, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_mam ,aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("F")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_mam

civetta_trac <- ggplot() +
  geom_point(data=data_civ_trac, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_trac, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_trac, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("G")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_trac

genetta_trac <- ggplot() +
  geom_point(data=data_gen_trac, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_trac, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_trac, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("H")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_trac

civetta_vert <- ggplot() +
  geom_point(data=data_civ_vert, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_vert, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_vert, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("I")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_vert

genetta_vert <- ggplot() +
  geom_point(data=data_gen_vert, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_vert, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_vert, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("J")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_vert

civetta_zeale <- ggplot() +
  geom_point(data=data_civ_zeale, aes(x=Samples, y=Richness)) +
  geom_line(data=data_civ_zeale, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_civ_zeale, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("K")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
civetta_zeale

genetta_zeale <- ggplot() +
  geom_point(data=data_gen_zeale, aes(x=Samples, y=Richness)) +
  geom_line(data=data_gen_zeale, aes(x=Samples, y=Richness)) +
  geom_ribbon(data=data_gen_zeale, aes(x=Samples, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+
  labs(y= "OTUs")+
  ggtitle("L")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))
genetta_zeale

#Combine plots
library("ggpubr")
ABC <- ggarrange(civetta_16S, genetta_16S, civetta_coleop, ncol = 3, nrow = 1)
DEF <- ggarrange(genetta_coleop, civetta_mam, genetta_mam, ncol = 3, nrow = 1)
GHI <- ggarrange(civetta_trac, genetta_trac, civetta_vert, ncol = 3, nrow = 1)
JKL <- ggarrange(genetta_vert, civetta_zeale, genetta_zeale, ncol = 3, nrow = 1)


plot <- ggarrange(ABC, DEF, GHI, JKL, ncol = 1, nrow = 4)
plot
