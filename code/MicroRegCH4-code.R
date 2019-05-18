rm(list=ls())


# Load packages and set theme for plotting ----
library(vegan); library(lmtest); library(tidyverse)
library(cowplot); library(ggplot2); library(ggthemes); library(RColorBrewer)

# Uses theme_base() from the ggthemes library
# Adds theme arguments to make the plot transparent and change axis spacing and sizes 
# Sets to default using theme_set()

theme_set(theme_bw()+theme(panel.grid = element_blank(), 
                           plot.margin=unit(c(0.5,0.1,0.5,0.5),"cm"),
                           plot.title = element_text(size=10, hjust = 0.5),
                           axis.text.y = element_text(size=8, angle=90, hjust = c(0.5,0.5)),
                           axis.text.x = element_text(size=8),
                           axis.title.y = element_text(size=10),
                           axis.title.x = element_text(size=10),
                           legend.title=element_text(size=10), 
                           legend.text=element_text(size=8),
                           legend.key.height = unit(0.5, "cm")))

# Read in data and set up for MCC and non-MCC analyses ----
#Read in OTU table
otu=read.csv("MicrobeRegCH4/data/otuTable.csv", stringsAsFactors = F) 
#Environmental metadata
env=read.csv("MicrobeRegCH4/data/envTable.csv", stringsAsFactors = F)

#Subset only the methanogen community 
otuM=otu[grep("Methano", otu$Consensus.Lineage),]
#Subset only OTUs that are not in the methanogen community 
otuB=otu[!(otu$OTU.ID%in%otuM$OTU.ID),]

#Calculate methanogen relative abundance 
env$methanogenRelAbund=colSums(otuM[,2:15])/colSums(otu[,2:15])

# PCOA analyses ----
#MCC
#Transform to relative abundance data
otuM.ra=decostand(otuM[,2:15], method="total", MARGIN = 2) 
otuM.ra$OTU.ID=otuM$OTU.ID; otuM.ra$tax=otuM$Consensus.Lineage #add OTU.ID and taxonomy info

#Generate distance matrix
distM=vegdist(t(otuM.ra[,1:14]), method="bray")
#Generate principle coordinate scaling with eigenvalues
pcoaM=cmdscale(distM, eig=T)
#Create new dataframe 
pcoaM.DF=merge(data.frame(LakeID=rownames(pcoaM$points),
                    axis1=pcoaM$points[,1],
                    axis2=pcoaM$points[,2]),env)

#Generate plot of MCC PCoA
mccPCOA=ggplot(pcoaM.DF, aes(x=axis1, y=axis2))+
  ggtitle("Lake sediment MCC")+
  geom_hline(yintercept = 0, linetype="dotted", color="grey80")+
  geom_vline(xintercept = 0, linetype="dotted", color="grey80")+
  geom_point(aes(size=OMperc, color=pH))+
  geom_point(aes(size=OMperc), shape=1)+
  geom_text(aes(label=LakeID), hjust=c(1.2), vjust=c(-0.5), size=2)+
  xlab("PCoA Axis 1 [29.6%]")+ylab("PCoA Axis 2 [20.4%]")+
  scale_color_gradient(low="grey90", high="black")+
  scale_size_continuous(range=c(1.5,4))+
  scale_x_continuous(limits=c(-0.3, 0.6))+
  scale_y_continuous(limits=c(-0.4, 0.5))+
  theme(legend.title = element_text(size=8),
        legend.text=element_text(size=6),
        legend.key.height = unit(0.3, "cm"), 
        legend.justification = c(1, 0), 
        legend.position = c(1, 0),
        legend.box = "horizontal",
        legend.background = element_blank())

#non-MCC
#Transform to relative abundance data
otuB.ra=decostand(otuB[,2:15], method="total", MARGIN = 2) 
otuB.ra$OTU.ID=otuB$OTU.ID; otuB.ra$tax=otuB$Consensus.Lineage #add OTU.ID and taxonomy info

#Generate distance matrix
distB=vegdist(t(otuB.ra[,1:14]), method="bray")
#Generate principle coordinate scaling with eigenvalues
pcoaB=cmdscale(distB, eig=T)
#Create new dataframe 
pcoaB.DF=merge(data.frame(LakeID=rownames(pcoaB$points), 
                          axis1=pcoaB$points[,1],
                          axis2=pcoaB$points[,2]), env)
#Generate plot of microbial community pcoa
nonPCOA=ggplot(pcoaB.DF, aes(x=axis1, y=axis2))+
  ggtitle("Lake sediment non-MCC")+
  geom_hline(yintercept = 0, linetype="dotted", color="grey80")+
  geom_vline(xintercept = 0, linetype="dotted", color="grey80")+
  geom_point(aes(size=OMperc, color=pH))+
  geom_point(aes(size=OMperc), shape=1)+
  geom_text(aes(label=LakeID), hjust=c(1.2), vjust=c(-0.5), size=2)+
  xlab("PCoA Axis 1 [21.1%]")+ylab("PCoA Axis 2 [15.8%]")+
  scale_color_gradient(low="grey90", high="black")+
  scale_size_continuous(range=c(1.5,4))+
  scale_x_continuous(limits=c(-0.5, 0.5))+
  scale_y_continuous(limits=c(-0.4, 0.4))+
  theme(legend.title = element_text(size=8),
        legend.text=element_text(size=6),
        legend.key.height = unit(0.3, "cm"), 
        legend.justification = c(1, 0), 
        legend.position = c(1, 0),
        legend.box = "horizontal",
        legend.background = element_blank())

#NMDS 
nmdsB=metaMDS(t(otuB.ra[,1:14]), distance="bray")
nmdsB

# Statistics ----
# MCC analysis
#Is community composition explained by pH or chlA?
adonis(distM~env$pH) 
adonis(distM~env$chlA)
#What about other environmental variables?
adonis(distM~env$OMperc)
adonis(distM~env$TN)
adonis(distM~env$TP)

#Linear regressions: PCoA1 ~ environmental variables
summary(lm(pcoaM.DF$axis1~env$pH))
summary(lm(pcoaM.DF$axis1~env$chlA)) 
summary(lm(pcoaM.DF$axis1~env$OMperc)) 
summary(lm(pcoaM.DF$axis1~env$TN))
summary(lm(pcoaM.DF$axis1~env$TP))

#Linear regression: PCoA2  ~ environmental variables
summary(lm(pcoaM.DF$axis2~env$pH))
summary(lm(pcoaM.DF$axis2~env$chlA)) 
summary(lm(pcoaM.DF$axis2~env$OMperc)) 
summary(lm(pcoaM.DF$axis2~env$TN))
summary(lm(pcoaM.DF$axis2~env$TP))

# Non-MCC analyses
#Is community composition explained by pH or chlA?
adonis(distB~env$pH) 
adonis(distB~env$chlA)
#What about other environmental variables?
adonis(distB~env$OMperc)
adonis(distB~env$TN)
adonis(distB~env$TP)

#Linear regressions: PCoA1 ~ environmental variables
summary(lm(pcoaB.DF$axis1~env$pH))
summary(lm(pcoaB.DF$axis1~env$chlA)) 
summary(lm(pcoaB.DF$axis1~env$OMperc)) 
summary(lm(pcoaB.DF$axis1~env$TN))
summary(lm(pcoaB.DF$axis1~env$TP))

#Linear regression: PCoA2  ~ environmental variables
summary(lm(pcoaB.DF$axis2~env$pH))
summary(lm(pcoaB.DF$axis2~env$chlA)) 
summary(lm(pcoaB.DF$axis2~env$OMperc)) 
summary(lm(pcoaB.DF$axis2~env$TN))
summary(lm(pcoaB.DF$axis2~env$TP))

#Create figure 2
legend=get_legend(nmccPCOA)
plot_grid(mccPCOA+guides(color=F, size=F), 
          nmccPCOA+guides(color=F, size=F), 
          legend,
          nrow=1, rel_widths = c(1,1,.3), labels=c("a", "b"))
ggsave("MicrobeRegCH4/figures/Fig2.pdf", width=7, height=3.5)


# Are the two distance matrices correlated? 
mantel(distM, distB)

# What is the best predictor of methanogenesis? ----
#Environmental 
summary(lm(env$CH4prod~env$chlA))
summary(lm(env$CH4prod~env$pH))
summary(lm(env$CH4prod~env$OMperc))
summary(lm(env$CH4prod~env$TN))
summary(lm(env$CH4prod~env$TP))
fit11=lm(env$CH4prod~env$TP+env$chlA)

#Microbial 
summary(lm(env$CH4prod~env$methanogenRelAbund))
summary(lm(env$CH4prod~pcoaM.DF$axis1))
summary(lm(env$CH4prod~pcoaB.DF$axis1))

#Environmental + Microbial 
fit1=lm(env$CH4prod~env$chlA)
fit2=lm(env$CH4prod~env$chlA+env$methanogenRelAbund)
fit3=lm(env$CH4prod~env$chlA+pcoaM.DF$axis1)
fit4=lm(env$CH4prod~env$chlA+pcoaB.DF$axis1)

fit12=lm(env$CH4prod~env$TP+env$chlA+pcoaB.DF$axis1)

lrtest(fit1, fit2)
lrtest(fit1, fit3)
lrtest(fit1, fit4)

fit5=lm(env$CH4prod~env$TP)
fit6=lm(env$CH4prod~env$TP+env$methanogenRelAbund)
fit7=lm(env$CH4prod~env$TP+pcoaM.DF$axis1)
fit8=lm(env$CH4prod~env$TP+pcoaB.DF$axis1)

lrtest(fit5, fit6)
lrtest(fit5, fit7)
lrtest(fit5, fit8)


# Plot observed vs predicted methanogenesis
pcoaB.DF$predCH4.chlA=coef(fit1)[1]+coef(fit1)[2]*pcoaB.DF$chlA
pcoaB.DF$predCH4.chlA.pcoaB1=coef(fit3)[1]+coef(fit3)[2]*pcoaB.DF$chlA+coef(fit3)[3]*pcoaB.DF$axis1


ggplot(pcoaB.DF, aes(x=predCH4.chlA, y=CH4prod))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.text.y = element_text(angle=90, hjust=c(0.5,0.5)))+
  xlim(2750, 19000)+ylim(2750, 19000)+
  geom_abline(slope = 1, linetype="dashed", color="grey50")+geom_point(size=2)+
  ylab("Observed methanogenesis")+
  xlab("Predicted methanogenesis\nPredictors: chl a")
  
ggplot(pcoaB.DF, aes(x=predCH4.chlA.pcoaB1, y=CH4prod))+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.text.y = element_text(angle=90, hjust=c(0.5,0.5)))+
  xlim(2750, 19000)+ylim(2750, 19000)+
  geom_abline(slope = 1, linetype="dashed", color="grey50")+geom_point(size=3)+
  ylab("Observed methanogenesis")+
  xlab("Predicted methanogenesis\nPredictors: chl a + non-MCC PCoA 1")+
  annotate("text", x=5000, y=17000, label=paste("R^2 ==", 0.54), parse=T)
ggsave("figures/observed-predicted-full.pdf", height=4, width=4)

# Create order-level barplot of non-MCC ----
#Use tidyverse to create order column
otuB.ra$phylum=otuB.ra$tax %>% str_split(";") %>% lapply(. %>% {.[2]}) %>% unlist() %>% str_replace(" p__", "")
otuB.ra$phylum[otuB.ra$phylum==""]=NA
unique(otuB.ra$phylum)
#Summarise based on order
nmccBar=otuB.ra[,c(1:14,17)] %>% group_by(phylum) %>% summarise_all(funs(sum))
#Create row sums 
nmccBar$sum=rowSums(nmccBar[,2:15])/sum(nmccBar[,2:15])

#Clean up naming conventions 
nmccBar$phylum[nmccBar$phylum=="[Caldithrix]"]="Caldithrix"
nmccBar$phylum[nmccBar$phylum=="[Parvarchaeota]"]="Parvarchaeota"

#Get Other group
otherNMCC=nmccBar[nmccBar$sum<0.02 | (is.na(nmccBar$phylum)),]

#Throw out all the others 
nmccBar[nmccBar$sum<0.02 | (is.na(nmccBar$phylum)),1]=NA
#Summarise 
nmccBar=nmccBar[,c(1:15)] %>% group_by(phylum) %>% summarise_all(funs(sum)) %>% gather("lakeID", "relAbund", 2:15)

#Create bar plot

nmccBar$lakeID=factor(nmccBar$lakeID, levels = c("BR", "MO", "BO", "WL", "PA", "PE",
                                                 "FO", "TU", "CR", "BE","HB", "BA","CB", "NG")) #order lakeID


nMCC.bar=ggplot(nmccBar, aes(x=lakeID, y=relAbund))+geom_col(aes(fill=phylum))+
  scale_fill_brewer(palette = "Spectral", na.value="grey30")+
  xlab("Lake ID")+ylab("Relative Abundance")+
  theme(legend.title = element_text(size=8), legend.text=element_text(size=6),
        legend.key.height = unit(0.3, "cm"), legend.position = "top")+
  scale_y_continuous(limits = c(0,1), expand = c(0.01,0.01))+
  guides(fill=guide_legend(nrow=4,byrow=F))
nMCC.bar


non.MCC.plot=plot_grid(nonPCOA, nMCC.bar, labels=c("a", "b"), align = "hv")
non.MCC.plot
#ggsave("figures/non-MCC.pdf", width=7, height=4)

#Summarise other ----
otherNMCC=otherNMCC[,c(1:15)] %>% group_by(phylum) %>% summarise_all(funs(sum)) 
otherNMCC=otherNMCC[,c(1:15)] %>% gather("lakeID", "relAbund", 2:15)

#
colourCount = length(unique(otherNMCC$phylum))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

nMCC.other.bar=ggplot(otherNMCC, aes(x=lakeID, y=relAbund))+geom_col(aes(fill=phylum))+
  scale_fill_manual(values = getPalette(colourCount), na.value="grey10")+
  xlab("Lake ID")+ylab("Relative Abundance")+
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=8),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom")+
  guides(fill=guide_legend(nrow=10,byrow=F))
nMCC.other.bar
ggsave("figures/non-MCC-barOther.pdf", width=6.5, height=5)


# Create order-level barplot of MCC ----
#Use tidyverse to create order column
otuM.ra$order=otuM.ra$tax %>% str_split(";") %>% lapply(. %>% {.[c(4,6)]}) %>% 
  lapply(. %>% {paste(.[1], .[2], sep=";")}) %>% unlist() %>% 
  str_replace(" o__", "") %>% str_replace(" g__", "")
unique(otuM.ra$order)
otuM.ra$order[otuM.ra$order=="NA;NA"]="Methanocellales;"
otuM.ra$order[otuM.ra$order=="E2;"]="Methanomassiliicoccales;"

unique(otuM.ra$order)

otuM.ra=otuM.ra[,c(1:14,17)] %>% group_by(order) %>% summarise_all(funs(sum))
#Create row sums 
otuM.ra$sum=rowSums(otuM.ra[,2:15])/sum(otuM.ra[,2:15])

#Get Other group
otherM=otuM.ra[otuM.ra$sum<0.01 | (is.na(otuM.ra$order)),]

#Throw out all the others 
otuM.ra[otuM.ra$sum<0.02 | (is.na(otuM.ra$order)),1]=NA
#Summarise 
otuM.ra=otuM.ra[,c(1:15)] %>% group_by(order) %>% summarise_all(funs(sum)) %>% gather("lakeID", "relAbund", 2:15)

#Create bar plot
otuM.ra$lakeID=factor(otuM.ra$lakeID, levels = c("BR", "MO", "BO", "WL", "PA", "PE",
                                                 "FO", "TU", "CR", "BE","HB", "BA","CB", "NG")) #order lakeID

MCC.bar=ggplot(otuM.ra, aes(x=lakeID, y=relAbund))+geom_col(aes(fill=order))+
  scale_fill_brewer(palette = "Spectral", na.value="grey30")+
  xlab("Lake ID")+ylab("Relative Abundance")+
  theme(legend.title = element_text(size=8), legend.text=element_text(size=6),
        legend.key.height = unit(0.3, "cm"), legend.position = "top")+
  scale_y_continuous(limits = c(0,1), expand = c(0.01,0.01))+
  guides(fill=guide_legend(nrow=4,byrow=F))
MCC.bar


MCC.plot=plot_grid(mccPCOA, MCC.bar, labels=c("a", "b"), align = "hv")
MCC.plot
#ggsave("figures/MCC.pdf", width=7, height=4)

#Summarise other ----
otherM=otherM[,c(1:15)] %>% group_by(order) %>% summarise_all(funs(sum)) 
otherM=otherM[,c(1:15)] %>% gather("lakeID", "relAbund", 2:15)

#
colourCount = length(unique(otherM$order))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

MCC.other.bar=ggplot(otherM, aes(x=lakeID, y=relAbund))+geom_col(aes(fill=order))+
  scale_fill_manual(values = getPalette(colourCount), na.value="grey10")+
  xlab("Lake ID")+ylab("Relative Abundance")+
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=8),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom")+
  guides(fill=guide_legend(nrow=5,byrow=F))
MCC.other.bar
ggsave("figures/MCC-barOther.pdf", width=6.5, height=5)







