#With vegan
# Bring in the ASV table
civetta_mb <- read.csv("ASV_civetta.csv", row.names=1, sep=";")

#Bring in diet table
civetta_diet <- read.csv("civetta_otu.csv", row.names=1, sep=";")

library(vegan)

# Make distances (I think Jaccard would be better since diets stuff only have precens absence data)
#first distances for microbiomes
mbc_dist <- vegdist(civetta_mb, method="jaccard")

#distance for diets
dc_dist <-vegdist(civetta_diet, method="jaccard")

shapiro <- shapiro.test(mbc_dist)
shapiro$p.value
shapiro <- shapiro.test(dc_dist)
shapiro$p.value

#normal mantel test
Man<- mantel(mbc_dist, dc_dist, method="pearson", permutations=10000, strata = NULL,na.rm = FALSE)
Man

#We need a data.frame/matrix format
library(reshape2)
df_mbc <- melt(as.matrix(mbc_dist),varnames = c("row", "col"))
df_mbc <- as.data.frame(df_mbc)
df_dc <- melt(as.matrix(dc_dist),varnames = c("row", "col"))
df_dc <- as.data.frame(df_dc)

#Change column name for the distance value
library(dplyr)
df_mbc %>% rename(value_mb = value) -> df_mbc
df_dc %>% rename(value_d = value) -> df_dc
#Join the two dataframes
df_mbc %>% inner_join(df_dc) -> dfc

library(ggplot2)
head(dfc)
civetta_plot <- ggplot(dfc, aes(x=value_mb, y=value_d)) +
  geom_bin2d() +
  scale_fill_steps(n.breaks = 15, low = "#c7eae5",high = "#003b2f", name = "No. of comparisons")+
  xlab("Jaccard similarity of the gut microbiome") +
  ylab("Jaccard similarity of the diet")+
  geom_smooth(method = "lm", colour = "black")+
  theme_classic()+
  ggtitle("(a)")+
  theme(axis.title = element_text(face="bold"))+
  theme(title = element_text(face="bold", size=14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.key.size = unit(1.5, 'cm'))+
  annotate(geom="text", x=0.2, y=3, label = "p>0.05", size=7)
civetta_plot

#Repeat for genetta spp. 
genetta_mb <- read.csv("ASV_genetta.csv", row.names=1, sep=";")

genetta_diet <- read.csv("genetta_otu.csv", row.names=1, sep=";")

mbg_dist <- vegdist(genetta_mb, method="jaccard")
dg_dist <-vegdist(genetta_diet, method="jaccard")

shapiro <- shapiro.test(mbg_dist)
shapiro$p.value
shapiro <- shapiro.test(dg_dist)
shapiro$p.value

Man<- mantel(mbg_dist, dg_dist, method="pearson", permutations=10000, strata = NULL,na.rm = FALSE)
Man

df_mbg <- melt(as.matrix(mbg_dist),varnames = c("row", "col"))
df_mbg <- as.data.frame(df_mbg)
df_dg <- melt(as.matrix(dg_dist),varnames = c("row", "col"))
df_dg <- as.data.frame(df_dg)

df_mbg %>% rename(value_mb = value) -> df_mbg
df_dg %>% rename(value_d = value) -> df_dg

df_mbg %>% inner_join(df_dg) -> dfg

head(dfg)
genetta_plot <- ggplot(dfg, aes(x=value_mb, y=value_d)) +
  geom_bin2d() +
  scale_fill_steps(n.breaks = 15, low = "#c7eae5",high = "#003b2f", name = "No. of comparisons")+
  xlab("Jaccard similarity of the gut microbiome") +
  ylab("Jaccard similarity of the diet")+
  geom_smooth(method = "lm", colour = "black")+
  theme_classic()+
  ggtitle("(b)")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.key.size = unit(1.5, 'cm'))+
  annotate(geom="text", x=0.2, y=1.23, label = "p<0.05", size=7)
genetta_plot


##Correlation
alpha_gutMB_civ <- specnumber(civetta_mb)
alpha_diet_civ <- specnumber(civetta_diet) 

alpha_gutMB_gen <- specnumber(genetta_mb)
alpha_diet_gen <- specnumber(genetta_diet) 

civ_test <- cor.test(alpha_gutMB_civ, alpha_diet_civ, method="spearman", exact=FALSE)
civ_test

gen_test <- cor.test(alpha_gutMB_gen, alpha_diet_gen, method="spearman", exact=FALSE)
gen_test

civetta_df <- data.frame(alpha_gutMB_civ, alpha_diet_civ)
genetta_df <- data.frame(alpha_gutMB_gen, alpha_diet_gen)

civ_plot <- ggplot(civetta_df, aes(x=alpha_gutMB_civ, y=alpha_diet_civ))+
  geom_point()+
  geom_smooth(method = "auto", colour = "black")+
  xlab("Gut microbiome richness")+
  ylab("Diet richness")+
  theme_classic()+
  ggtitle("(c)")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 12))+
  annotate(geom="text", x=250, y=60, label = "p>0.05", size=7)
civ_plot 

gen_plot <- ggplot(genetta_df, aes(x=alpha_gutMB_gen, y=alpha_diet_gen))+
  geom_point()+
  geom_smooth(method = "auto", colour = "black")+
  xlab("Gut microbiome richness")+
  ylab("Diet richness")+
  theme_classic()+
  ggtitle("(d)")+
  theme(axis.title = element_text(face="bold", size = 14))+
  theme(title = element_text(face="bold", size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(axis.text = element_text(size = 14))+
  theme(plot.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 12))+
  annotate(geom="text", x=100, y=45, label = "p>0.05", size=7)
gen_plot 

##combine plots
library(ggpubr)


man_plot <- ggarrange(civetta_plot, genetta_plot, ncol = 2, nrow = 1)
man_plot

cor_plot <- ggarrange(civ_plot, gen_plot, ncol = 2, nrow = 1)
cor_plot

plot <- ggarrange(man_plot, cor_plot, ncol = 1, nrow = 2)
plot

