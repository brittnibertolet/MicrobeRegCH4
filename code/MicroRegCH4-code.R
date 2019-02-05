rm(list=ls())

########################################
#### Load packages and read in data ####
########################################
library(vegan)
library(ggplot2)
library(lmtest)

#Read in OTU table
otu=read.csv("git-repo/otuTable.csv", stringsAsFactors = F) 
#Environmental metadata
env=read.csv("git-repo/envTable.csv", stringsAsFactors = F)

#Subset only the methanogen community 
otuM=otu[grep("Methano", otu$Consensus.Lineage),]
#Subset only OTUs that are not in the methanogen community 
otuB=otu[!(otu$OTU.ID%in%otuM$OTU.ID),]


#############################################################
#### Methanogen microbial community composition analyses ####
#############################################################
#Transform to relative abundance data
otuM.ra=decostand(otuM[,2:15], method="total", MARGIN = 2) 
otuM.ra$OTU.ID=otuM$OTU.ID; otuM.ra$tax=otuM$Consensus.Lineage #add OTU.ID and taxonomy info

#Generate distance matrix
distM=vegdist(t(otuM.ra[,1:14]), method="bray")
#Generate principle coordinate scaling with eigenvalues
pcoaM=cmdscale(distM, eig=T)
#Create new dataframe 
pcoaM.DF=data.frame(LakeID=rownames(pcoaM$points),axis1=pcoaM$points[,1],axis2=pcoaM$points[,2])

#Generate plot of MCC PCoA
ggplot(pcoaM.DF, aes(x=axis1, y=axis2))+
  theme_bw()+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                   axis.text.y = element_text(angle=90, hjust = c(0.5,0.5)))+
  ggtitle("Lake sediment MCC")+
  geom_hline(yintercept = 0, linetype="dotted", color="grey80")+
  geom_vline(xintercept = 0, linetype="dotted", color="grey80")+
  geom_point(size=2)+
  geom_text(aes(label=LakeID), hjust=c(1.2), vjust=c(-0.5), size=3)+
  xlab("PCoA Axis 1 [28.9%]")+ylab("PCoA Axis 2 [20.1%]")+
  scale_x_continuous(limits=c(-0.3, 0.6))+
  scale_y_continuous(limits=c(-0.3, 0.45))
#ggsave("figures/MCC-PCOA.pdf", height=3.5, width=3.5)

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

#################################################################
#### Non-methanogen microbial community composition analyses ####
#################################################################
#Transform to relative abundance data
otuB.ra=decostand(otuB[,2:15], method="total", MARGIN = 2) 
otuB.ra$OTU.ID=otuB$OTU.ID; otuB.ra$tax=otuB$Consensus.Lineage #add OTU.ID and taxonomy info

#Generate distance matrix
distB=vegdist(t(otuB.ra[,1:14]), method="bray")
#Generate principle coordinate scaling with eigenvalues
pcoaB=cmdscale(distB, eig=T)
#Create new dataframe 
pcoaB.DF=data.frame(LakeID=rownames(pcoaB$points), 
                        axis1=pcoaB$points[,1],
                        axis2=pcoaB$points[,2])
#Generate plot of microbial community pcoa
ggplot(pcoaB.DF, aes(x=axis1, y=axis2))+
  theme_bw()+theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
                   axis.text.y = element_text(angle=90, hjust = c(0.5,0.5)))+
  ggtitle("Lake sediment non-MCC")+
  geom_hline(yintercept = 0, linetype="dotted", color="grey80")+
  geom_vline(xintercept = 0, linetype="dotted", color="grey80")+
  geom_point(size=2)+geom_text(aes(label=LakeID), hjust=c(1.2), vjust=c(-0.5), size=3)+
  xlab("PCoA Axis 1 [20.4%]")+ylab("PCoA Axis 2 [15.4%]")+
  scale_x_continuous(limits=c(-0.5, 0.35))+
  scale_y_continuous(limits=c(-0.35, 0.35))
ggsave("figures/N-MCC-PCOA.pdf", height=3.5, width=3.5)

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


###################################################
#### Are the two distance matrices correlated? #### 
###################################################
mantel(distM, distB)

############################################################################
#### Do metrics of community composition correlate with methanogenesis? ####
############################################################################

#Is non-MCC PCoA axes related to methanogenesis
fit1=lm(env$CH4prod~env$chlA) 
fit2=lm(env$CH4prod~env$chlA+pcoaM.DF$axis1)
fit3=lm(env$CH4prod~env$chlA+pcoaB.DF$axis1)
fit4=lm(env$CH4prod~env$chlA+env$pH)
fit5=lm(env$CH4prod~env$chlA+env$OMperc)

lrtest(fit1, fit2)
lrtest(fit1, fit3)
lrtest(fit1, fit4)
lrtest(fit1, fit5)




