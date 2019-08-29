#R script to summarize sense results in a table for Herring MSE, June 2016, based on Beer Belly script

#setwd("C:/0_Data/BeerBelly/BeerBelly_MS")
setwd("~/Data/Projects/MSE/HerringMSE/GOM")
#setwd("~/Dropbox/MSE/HerringMSE/GOM")

GOMsensedat<-read.csv("results/GOM_Diffs_novul_allquants.csv")
GOMsensedatvul<-read.csv("results/GOM_Diffs_vul_allquants.csv")

#need to change column names coming out of results folder to have Actor and Acted and possibly other sensible names
names(GOMsensedat)[1:2]<-c("Acted", "Actor")
names(GOMsensedatvul)[1:2]<-c("Acted", "Actor")


levels(GOMsensedat$Actor)
#[1] "Bacteria"                             
#[2] "Baleen Whales"                        
#[3] "Demersals- benthivores"               
#[4] "Demersals- omnivores"                 
#[5] "Demersals- piscivores"                
#[6] "Detritus-POC"                         
#[7] "Discard"                              
#[8] "Fishery"                              
#[9] "Gelatinous Zooplankton"               
#[10] "HMS"                                  
#[11] "Large Copepods"                       
#[12] "Larval-juv fish- all"                 
#[13] "Macrobenthos- crustaceans"            
#[14] "Macrobenthos- molluscs"               
#[15] "Macrobenthos- other"                  
#[16] "Macrobenthos- polychaetes"            
#[17] "Medium Pelagics- (piscivores & other)"
#[18] "Megabenthos- filterers"               
#[19] "Megabenthos- other"                   
#[20] "Micronekton"                          
#[21] "Microzooplankton"                     
#[22] "Odontocetes"                          
#[23] "Phytoplankton- Primary Producers"     
#[24] "Pinnipeds"                            
#[25] "Sea Birds"                            
#[26] "Sharks- pelagics"                     
#[27] "Shrimp et al."                        
#[28] "Small copepods"                       
#[29] "Small Pelagics- anadromous"           
#[30] "Small Pelagics- commercial"           
#[31] "Small Pelagics- other"                
#[32] "Small Pelagics- squid"                

levels(GOMsensedat$Acted)
#[1] "avP_Bacteria"                             
#[2] "avP_Baleen.Whales"                        
#[3] "avP_Demersals..benthivores"               
#[4] "avP_Demersals..omnivores"                 
#[5] "avP_Demersals..piscivores"                
#[6] "avP_Detritus.POC"                         
#[7] "avP_Discard"                              
#[8] "avP_Gelatinous.Zooplankton"               
#[9] "avP_HMS"                                  
#[10] "avP_Large.Copepods"                       
#[11] "avP_Larval.juv.fish..all"                 
#[12] "avP_Macrobenthos..crustaceans"            
#[13] "avP_Macrobenthos..molluscs"               
#[14] "avP_Macrobenthos..other"                  
#[15] "avP_Macrobenthos..polychaetes"            
#[16] "avP_Medium.Pelagics...piscivores...other."
#[17] "avP_Megabenthos..filterers"               
#[18] "avP_Megabenthos..other"                   
#[19] "avP_Micronekton"                          
#[20] "avP_Microzooplankton"                     
#[21] "avP_Odontocetes"                          
#[22] "avP_Phytoplankton..Primary.Producers"     
#[23] "avP_Pinnipeds"                            
#[24] "avP_Sea.Birds"                            
#[25] "avP_Sharks..pelagics"                     
#[26] "avP_Shrimp.et.al."                        
#[27] "avP_Small.copepods"                       
#[28] "avP_Small.Pelagics..anadromous"           
#[29] "avP_Small.Pelagics..commercial"           
#[30] "avP_Small.Pelagics..other"                
#[31] "avP_Small.Pelagics..squid"                
#[32] "X"                 

GOMpel.actor.sensedat<-subset(GOMsensedat, Actor=="Small Pelagics- commercial")
GOMpel.acted.sensedat<-subset(GOMsensedat, Acted=="avP_Small.Pelagics..commercial")

GOMphy.actor.sensedat<-subset(GOMsensedat, Actor=="Phytoplankton- Primary Producers")

GOMbird.acted.sensedat<-subset(GOMsensedat, Acted=="avP_Sea.Birds")
GOMHMS.acted.sensedat<-subset(GOMsensedat, Acted=="avP_HMS")
GOMbal.acted.sensedat<-subset(GOMsensedat, Acted=="avP_Baleen.Whales")

#counts of total significant changes by system for table
with(GOMpel.actor.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant poll 90% quantile
with(GOMpel.acted.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant spp impacts 90% quantile
with(GOMphy.actor.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant lgphy 90% quantile

GOMpelv.actor.sensedat<-subset(GOMsensedatvul, Actor=="Small Pelagics- commercial")
GOMpelv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Small.Pelagics..commercial")

GOMphyv.actor.sensedat<-subset(GOMsensedatvul, Actor=="Phytoplankton- Primary Producers")

GOMbirdv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Sea.Birds")
GOMHMSv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_HMS")
GOMbalv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Baleen.Whales")
GOMpinnv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Pinnipeds")
GOMgfishv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Demersals..piscivores")
GOModontv.acted.sensedat<-subset(GOMsensedatvul, Acted=="avP_Odontocetes")


#counts of total significant changes by system for table
with(GOMpelv.actor.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant poll 90% quantile
with(GOMpelv.acted.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant spp impacts 90% quantile
with(GOMphyv.actor.sensedat, sum((X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), na.rm=T)) #total significant lgphy 90% quantile

sigtest<-function(dat)
{
  ((dat$X5.<0. & dat$X95.<0.)|(dat$X5.>0. & dat$X95.>0.))
}

#choose one of the below for boxplots

boxdat<-GOMpel.actor.sensedat #impacts of 10% increase in herring survival, vuls at defaults
#boxdat<-GOMpel.acted.sensedat
#boxdat<-GOMphy.actor.sensedat

#boxdat<-GOMbird.acted.sensedat
#boxdat<-GOMHMS.acted.sensedat
#boxdat<-GOMbal.acted.sensedat


boxdat<-GOMpelv.actor.sensedat #impacts of 10% increase in herring survival, vuls varying
#boxdat<-GOMpelv.acted.sensedat
#boxdat<-GOMphyv.actor.sensedat

#boxdat<-GOMbirdv.acted.sensedat
#boxdat<-GOMHMSv.acted.sensedat
#boxdat<-GOMbalv.acted.sensedat
#boxdat<-GOMpinnv.acted.sensedat
#boxdat<-GOMgfishv.acted.sensedat
#boxdat<-GOModontv.acted.sensedat


boxdat<-boxdat[complete.cases(boxdat),] #gets rid of NA rows

#pdf(paste(filename, "_box.pdf", sep=""), width=10.4, height=9.5)
par(mfrow=c(1,1))
par(mar=c(14,4,4,2)+0.1)
par(oma=c(2,0,0,0))
plot((1:(dim(boxdat)[1])), seq(min(boxdat$X5.), max(boxdat$X95.), length.out=dim(boxdat)[1]), 
     type="n", xlab="", ylab="Proportional Difference from Base", axes=F)
abline(h=0.0, lty=2)
#mtext(paste(outname))
#use for who does Actor matter to
#axis(1, at=1:(dim(boxdat)[1]), lab=boxdat$Acted[1:dim(boxdat)[1]], las=2)
#use for who matters to Actor
axis(1, at=1:(dim(boxdat)[1]), lab=boxdat$Actor[1:dim(boxdat)[1]], las=2)

mtext("Species", side=1, outer=T)
axis(2)
for (j in 1:dim(boxdat)[1]){
  rect((j-0.2),boxdat$X25.[j],(j+0.2),boxdat$X75.[j], col="grey")
  arrows(j, boxdat$X25.[j], j, boxdat$X5.[j], angle=90, length=0.1)
  arrows(j, boxdat$X75.[j], j, boxdat$X95.[j], angle=90, length=0.1)
}
#dev.off()


with(GOArecdown, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant recdown 90% quantile
with(GOAfishup, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant fishup 90% quantile
with(GOAboth, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant both 90% quantile



#GOA increasing/decreasing counts for table
#increasing
with(GOApol.actor.sensedat, sum((X5.>0. & X95.>0.))) #increasing significant poll 90% quantile
with(GOApol.acted.sensedat, sum((X5.>0. & X95.>0.))) #increasing significant spp impacts 90% quantile
with(GOAlgphy.actor.sensedat, sum((X5.>0. & X95.>0.))) #increasing significant lgphy 90% quantile
with(GOAsmphy.actor.sensedat, sum((X5.>0. & X95.>0.))) #increasing significant smphy 90% quantile

with(GOArecdown, sum((CI_5>0. & CI_95>0.))) #increasing significant recdown 90% quantile
with(GOAfishup, sum((CI_5>0. & CI_95>0.))) #increasing significant fishup 90% quantile
with(GOAboth, sum((CI_5>0. & CI_95>0.))) #increasing significant both 90% quantile

#decreasing
with(GOApol.actor.sensedat, sum((X5.<0. & X95.<0.))) #decreasing significant poll 90% quantile
with(GOApol.acted.sensedat, sum((X5.<0. & X95.<0.))) #decreasing significant spp impacts 90% quantile
with(GOAlgphy.actor.sensedat, sum((X5.<0. & X95.<0.))) #decreasing significant lgphy 90% quantile
with(GOAsmphy.actor.sensedat, sum((X5.<0. & X95.<0.))) #decreasing significant smphy 90% quantile

with(GOArecdown, sum((CI_5<0. & CI_95<0.))) #decreasing significant recdown 90% quantile
with(GOAfishup, sum((CI_5<0. & CI_95<0.))) #decreasing significant fishup 90% quantile
with(GOAboth, sum((CI_5<0. & CI_95<0.))) #decreasing significant both 90% quantile

#names of significantly changing groups
subset(GOApol.actor.sensedat, (X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), select=Acted)

#names and 50% range, decreasing spp
subset(GOApol.actor.sensedat, (X5.<0. & X95.<0.), select=c(Acted, X25., X75.))

#names and 50% range, increasing spp
subset(GOApol.actor.sensedat, (X5.>0. & X95.>0.), select=c(Acted, X25., X75.))

#who matters to pollock? who brings them down?
subset(GOApol.acted.sensedat, (X5.<0. & X95.<0.), select=c(Actor, X25., X75.))

#who matters to pollock? who props them up?
subset(GOApol.acted.sensedat, (X5.>0. & X95.>0.), select=c(Actor, X25., X75.))


# try a table with bars to the side?
library(ggplot2)
library(gridExtra)
library(plotrix)
library(grid)

########  Fig 6  ######################

sigtest<-function(dat)
{
  ((dat$X5.<0. & dat$X95.<0.)|(dat$X5.>0. & dat$X95.>0.))
}

GOAgroups<-c(
"W. Pollock", 
"W. Pollock_Juv", 
"Arrowtooth",
"P. Halibut",
"P. Cod",
"Longnose skate",
"Other skates",
"Steller Sea Lion",
"Steller Sea Lion_Juv",
"Resident seals",
"N. Fur Seal",
"N. Fur Seal_Juv",
"Resident Killers",
"Copepods",
"Euphausiids",
"Pandalidae",
"NP shrimp"
)

GOAgroupdat<-data.frame(subset(GOApol.actor.sensedat, Acted %in% GOAgroups, select=c(Acted, X5., X25., X50., X75., X95.)))
GOAgroupdat$sig <- ifelse(sigtest(GOAgroupdat), "yes", "no")
GOAgroupdat$Acted <- factor(GOAgroupdat$Acted, levels = rev(c("W. Pollock", 
                                                          "W. Pollock_Juv", 
                                                          "Arrowtooth",
                                                          "P. Halibut",
                                                          "P. Cod",
                                                          "Longnose skate",
                                                          "Other skates",
                                                          "Steller Sea Lion",
                                                          "Steller Sea Lion_Juv",
                                                          "Resident seals",
                                                          "N. Fur Seal",
                                                          "N. Fur Seal_Juv",
                                                          "Resident Killers",
                                                          "Copepods",
                                                          "Euphausiids",
                                                          "Pandalidae",
                                                          "NP shrimp")))

group1GOA <- ggplot(GOAgroupdat, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","green")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-0.1, 1.99), breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  #axis.break(2, .5, style="zigzag")
  coord_flip()

EBSgroups<-c(
  "W..Pollock", 
  "W..Pollock_Juv", 
  "Arrowtooth",
  "P..Halibut",
  "P..Cod",
  "Alaska.skate",
  "Other.skates",
  "Steller.Sea.Lion",
  "Steller.Sea.Lion_Juv",
  "Resident.seals",
  "N..Fur.Seal",
  "N..Fur.Seal_Juv",
  "Resident.Killers",
  "Copepods",
  "Euphausiids",
  "Pandalidae",
  "NP.shrimp"
)


EBSgroupdat<-data.frame(subset(EBSpol.actor.sensedat, Acted %in% EBSgroups, select=c(Acted, X5., X25., X50., X75., X95.)))
EBSgroupdat$sig <- ifelse(sigtest(EBSgroupdat), "yes", "no")
EBSgroupdat$Acted <- factor(EBSgroupdat$Acted, levels = rev(c("W..Pollock", 
                                                              "W..Pollock_Juv", 
                                                              "Arrowtooth",
                                                              "P..Halibut",
                                                              "P..Cod",
                                                              "Alaska.skate",
                                                              "Other.skates",
                                                              "Steller.Sea.Lion",
                                                              "Steller.Sea.Lion_Juv",
                                                              "Resident.seals",
                                                              "N..Fur.Seal",
                                                              "N..Fur.Seal_Juv",
                                                              "Resident.Killers",
                                                              "Copepods",
                                                              "Euphausiids",
                                                              "Pandalidae",
                                                              "NP.shrimp")))

group1EBS <- ggplot(EBSgroupdat, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","blue")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-0.2, .5), breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()

group1title<-paste("Experiment 1: Adult pollock increase 10%")

#grid.arrange(group1GOA, group1EBS, ncol=2)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("Fig6.pdf", width=13.6, height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(20, 46))) # 20 rows, 45 columns
grid.text(group1title, vp = vplayout(1:1, 1:46))
print(group1GOA, vp = vplayout(2:20, 1:30))  
print(group1EBS, vp = vplayout(2:20, 31:46))
dev.off()

########  Fig 7  ###################### same as 6 but with GOApol.acted.sensedat and EBS

GOAgroupdat2<-data.frame(subset(GOApol.acted.sensedat, Actor %in% GOAgroups, select=c(Actor, X5., X25., X50., X75., X95.)))
GOAgroupdat2$sig <- ifelse(sigtest(GOAgroupdat2), "yes", "no")
GOAgroupdat2$Actor <- factor(GOAgroupdat2$Actor, levels = rev(c("W. Pollock", 
                                                              "W. Pollock_Juv", 
                                                              "Arrowtooth",
                                                              "P. Halibut",
                                                              "P. Cod",
                                                              "Longnose skate",
                                                              "Other skates",
                                                              "Steller Sea Lion",
                                                              "Steller Sea Lion_Juv",
                                                              "Resident seals",
                                                              "N. Fur Seal",
                                                              "N. Fur Seal_Juv",
                                                              "Resident Killers",
                                                              "Copepods",
                                                              "Euphausiids",
                                                              "Pandalidae",
                                                              "NP shrimp")))

group2GOA <- ggplot(GOAgroupdat2, aes(Actor, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","green")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-0.75, 4.3), breaks = seq(-.8,5,.2), labels = seq(-.8,5,.2)) +
  labs(y="GOA Pollock proportional change from base", x="") +
  #axis.break(2, .5, style="zigzag")
  coord_flip()

EBSgroups2<-c(
  "W. Pollock", 
  "W. Pollock_Juv", 
  "Arrowtooth",
  "P. Halibut",
  "P. Cod",
  "Alaska skate",
  "Other skates",
  "Steller Sea Lion",
  "Steller Sea Lion_Juv",
  "Resident seals",
  "N. Fur Seal",
  "N. Fur Seal_Juv",
  "Resident Killers",
  "Copepods",
  "Euphausiids",
  "Pandalidae",
  "NP shrimp"
)


EBSgroupdat2<-data.frame(subset(EBSpol.acted.sensedat, Actor %in% EBSgroups2, select=c(Actor, X5., X25., X50., X75., X95.)))
EBSgroupdat2$sig <- ifelse(sigtest(EBSgroupdat2), "yes", "no")
EBSgroupdat2$Actor <- factor(EBSgroupdat2$Actor, levels = rev(c("W. Pollock", 
                                                              "W. Pollock_Juv", 
                                                              "Arrowtooth",
                                                              "P. Halibut",
                                                              "P. Cod",
                                                              "Alaska skate",
                                                              "Other skates",
                                                              "Steller Sea Lion",
                                                              "Steller Sea Lion_Juv",
                                                              "Resident seals",
                                                              "N. Fur Seal",
                                                              "N. Fur Seal_Juv",
                                                              "Resident Killers",
                                                              "Copepods",
                                                              "Euphausiids",
                                                              "Pandalidae",
                                                              "NP shrimp")))

group2EBS <- ggplot(EBSgroupdat2, aes(Actor, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","blue")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-0.1, .6), breaks = seq(-1,1,.2), labels = seq(-1,1,.2)) +
  labs(y="EBS Pollock proportional change from base", x="") +
  coord_flip()

group2title<-paste("Experiment 2: Other groups individually increase 10%")

#grid.arrange(group1GOA, group1EBS, ncol=2)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("Fig7.pdf", width=13.6, height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(20, 75))) # 20 rows, 45 columns
grid.text(group2title, vp = vplayout(1:1, 1:75))
print(group2GOA, vp = vplayout(2:20, 1:54))  # the big plot covers rows 1 and cols 1:4
print(group2EBS, vp = vplayout(2:20, 55:75))
dev.off()

########  Fig 8  ###################### same as 6 but with GOAlgphy.actor.sensedat and EBS

GOAgroupdat3<-data.frame(subset(GOAlgphy.actor.sensedat, Acted %in% GOAgroups, select=c(Acted, X5., X25., X50., X75., X95.)))
GOAgroupdat3$sig <- ifelse(sigtest(GOAgroupdat3), "yes", "no")
GOAgroupdat3$Acted <- factor(GOAgroupdat3$Acted, levels = rev(c("W. Pollock", 
                                                              "W. Pollock_Juv", 
                                                              "Arrowtooth",
                                                              "P. Halibut",
                                                              "P. Cod",
                                                              "Longnose skate",
                                                              "Other skates",
                                                              "Steller Sea Lion",
                                                              "Steller Sea Lion_Juv",
                                                              "Resident seals",
                                                              "N. Fur Seal",
                                                              "N. Fur Seal_Juv",
                                                              "Resident Killers",
                                                              "Copepods",
                                                              "Euphausiids",
                                                              "Pandalidae",
                                                              "NP shrimp")))

group3GOA <- ggplot(GOAgroupdat3, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","green")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-.2, .6), breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  #axis.break(2, .5, style="zigzag")
  coord_flip()


EBSgroupdat3<-data.frame(subset(EBSlgphy.actor.sensedat, Acted %in% EBSgroups, select=c(Acted, X5., X25., X50., X75., X95.)))
EBSgroupdat3$sig <- ifelse(sigtest(EBSgroupdat3), "yes", "no")
EBSgroupdat3$Acted <- factor(EBSgroupdat3$Acted, levels = rev(c("W..Pollock", 
                                                              "W..Pollock_Juv", 
                                                              "Arrowtooth",
                                                              "P..Halibut",
                                                              "P..Cod",
                                                              "Alaska.skate",
                                                              "Other.skates",
                                                              "Steller.Sea.Lion",
                                                              "Steller.Sea.Lion_Juv",
                                                              "Resident.seals",
                                                              "N..Fur.Seal",
                                                              "N..Fur.Seal_Juv",
                                                              "Resident.Killers",
                                                              "Copepods",
                                                              "Euphausiids",
                                                              "Pandalidae",
                                                              "NP.shrimp")))

group3EBS <- ggplot(EBSgroupdat3, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25., fill=sig), fatten = 0.1, width = 0.5) +
  scale_fill_manual(values=c("white","blue")) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(limits = c(-0.2, .6), breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()

pdf("Fig8.pdf", width=13.6, height=8)
grid.arrange(group3GOA, group3EBS, ncol=2, top="Experiment 3: Large phytoplankton increase 10%")
dev.off()

#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

#pdf("Fig8.pdf", width=13.6, height=8)
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(20, 46))) # 20 rows, 45 columns
#print(group3GOA, vp = vplayout(1:20, 1:30))  # the big plot covers rows 1 and cols 1:4
#print(group3EBS, vp = vplayout(1:20, 31:46))
#dev.off()

########  Fig A1   #################
GOAsigdat<-data.frame(subset(GOApol.actor.sensedat, (X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), select=c(Acted, X5., X25., X50., X75., X95.)))
#GOAsigdat<-GOAsigdat[with(GOAsigdat, order(-X50.)), ]
GOAsigdat$Acted <- factor(GOAsigdat$Acted, levels = GOAsigdat$Acted[order(abs(GOAsigdat$X50.))])

EBSsigdat<-data.frame(subset(EBSpol.actor.sensedat, (X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), select=c(Acted, X5., X25., X50., X75., X95.)))
#EBSsigdat<-EBSsigdat[with(EBSsigdat, order(-X50.)), ]
EBSsigdat$Acted <- factor(EBSsigdat$Acted, levels = EBSsigdat$Acted[order(abs(EBSsigdat$X50.))])

pGOA <- ggplot(GOAsigdat, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25.), fill="green", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  coord_flip()

pEBS <- ggplot(EBSsigdat, aes(Acted, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25.), fill="blue", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()


#grid.arrange(pGOA, pEBS, ncol=2)    

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("FigS1.pdf", width=13.6, height=8)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(43, 32))) # 43 rows, 31 columns
  grid.text(group1title, vp = vplayout(1:1, 1:32))
  print(pGOA, vp = vplayout(2:17, 1:23))  # the big plot covers rows 1 and cols 1:4
  print(pEBS, vp = vplayout(2:43, 24:32))
dev.off()
######################################################

##############  Fig A2  #################################
GOApollimpact<-data.frame(subset(GOApol.acted.sensedat, (X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), select=c(Actor, X5., X25., X50., X75., X95.)))
#GOAsigdat<-GOAsigdat[with(GOAsigdat, order(-X50.)), ]
GOApollimpact$Actor <- factor(GOApollimpact$Actor, levels = GOApollimpact$Actor[order(abs(GOApollimpact$X50.))])

EBSpollimpact<-data.frame(subset(EBSpol.acted.sensedat, (X5.<0. & X95.<0.)|(X5.>0. & X95.>0.), select=c(Actor, X5., X25., X50., X75., X95.)))
#EBSsigdat<-EBSsigdat[with(EBSsigdat, order(-X50.)), ]
EBSpollimpact$Actor <- factor(EBSpollimpact$Actor, levels = EBSpollimpact$Actor[order(abs(EBSpollimpact$X50.))])

pollGOA <- ggplot(GOApollimpact, aes(Actor, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25.), fill="green", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-.8,4.3,.2), labels = seq(-.8,4.3,.2)) +
  labs(y="GOA Pollock proportional change from base", x="") +
  coord_flip()

pollEBS <- ggplot(EBSpollimpact, aes(Actor, X50.)) +
  geom_errorbar(aes(ymax = X95., ymin = X5.), height = 0.15) +
  geom_crossbar(aes(ymax = X75., ymin = X25.), fill="blue", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,1,.2), labels = seq(-1,1,.2)) +
  labs(y="EBS Pollock proportional change from base", x="") +
  coord_flip()


#grid.arrange(pollGOA, pollEBS, ncol=2)    

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("FigS2.pdf", width=13.6, height=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(28, 80))) # 43 rows, 31 columns
grid.text(group2title, vp = vplayout(1:1, 1:80))
print(pollGOA, vp = vplayout(2:15, 1:63))  # the big plot covers rows 1 and cols 1:4
print(pollEBS, vp = vplayout(2:28, 64:80))
dev.off()
######################################################

##############  Fig A3  #################################
#with(GOArecdown, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant recdown 90% quantile


GOAsigdat3<-data.frame(subset(GOArecdown, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#GOAsigdat<-GOAsigdat[with(GOAsigdat, order(-X50.)), ]
GOAsigdat3$Acted <- factor(GOAsigdat3$Acted, levels = GOAsigdat3$Acted[order(abs(GOAsigdat3$CI_50))])

EBSsigdat3<-data.frame(subset(EBSrecdown, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#EBSsigdat<-EBSsigdat[with(EBSsigdat, order(-X50.)), ]
EBSsigdat3$Acted <- factor(EBSsigdat3$Acted, levels = EBSsigdat3$Acted[order(abs(EBSsigdat3$CI_50))])

pGOA3 <- ggplot(GOAsigdat3, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="green", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  coord_flip()

pEBS3 <- ggplot(EBSsigdat3, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="blue", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()

groupS3title<-paste("Experiment 5: Climate-driven decreased pollock recruitment")

#grid.arrange(pGOA3, pEBS3, ncol=2)    

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("FigS3.pdf", width=10, height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(54, 19))) # 54 rows, 19 columns
print(pGOA3, vp = vplayout(2:19, 1:6))  # the big plot covers rows 1 and cols 1:4
print(pEBS3, vp = vplayout(2:54, 7:19))
grid.text(groupS3title, vp = vplayout(1:1, 1:19))
dev.off()


##############  Fig A4  #################################
#with(GOAfishup, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant recdown 90% quantile


GOAsigdat4<-data.frame(subset(GOAfishup, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#GOAsigdat<-GOAsigdat[with(GOAsigdat, order(-X50.)), ]
GOAsigdat4$Acted <- factor(GOAsigdat4$Acted, levels = GOAsigdat4$Acted[order(abs(GOAsigdat4$CI_50))])

EBSsigdat4<-data.frame(subset(EBSfishup, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#EBSsigdat<-EBSsigdat[with(EBSsigdat, order(-X50.)), ]
EBSsigdat4$Acted <- factor(EBSsigdat4$Acted, levels = EBSsigdat4$Acted[order(abs(EBSsigdat4$CI_50))])

pGOA4 <- ggplot(GOAsigdat4, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="green", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  coord_flip()

pEBS4 <- ggplot(EBSsigdat4, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="blue", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()

groupS4title<-paste("Experiment 6: Increased pollock fishing")

#grid.arrange(pGOA4, pEBS4, ncol=2)    

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("FigS4.pdf", width=10, height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(45, 19))) # 54 rows, 19 columns
print(pGOA4, vp = vplayout(2:7, 1:6))  # the big plot covers rows 1 and cols 1:4
print(pEBS4, vp = vplayout(2:45, 7:19))
grid.text(groupS4title, vp = vplayout(1:1, 1:19))
dev.off()

##############  Fig A5  #################################
#with(GOAboth, sum((CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.))) #total significant recdown 90% quantile


GOAsigdat5<-data.frame(subset(GOAboth, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#GOAsigdat<-GOAsigdat[with(GOAsigdat, order(-X50.)), ]
GOAsigdat5$Acted <- factor(GOAsigdat5$Acted, levels = GOAsigdat5$Acted[order(abs(GOAsigdat5$CI_50))])

EBSsigdat5<-data.frame(subset(EBSboth, (CI_5<0. & CI_95<0.)|(CI_5>0. & CI_95>0.), select=c(Acted, CI_5, CI_25, CI_50, CI_75, CI_95)))
#EBSsigdat<-EBSsigdat[with(EBSsigdat, order(-X50.)), ]
EBSsigdat5$Acted <- factor(EBSsigdat5$Acted, levels = EBSsigdat5$Acted[order(abs(EBSsigdat5$CI_50))])

pGOA5 <- ggplot(GOAsigdat5, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="green", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,2,.1), labels = seq(-1,2,.1)) +
  labs(y="GOA Proportional change from base", x="") +
  coord_flip()

pEBS5 <- ggplot(EBSsigdat5, aes(Acted, CI_50)) +
  geom_errorbar(aes(ymax = CI_95, ymin = CI_5), height = 0.15) +
  geom_crossbar(aes(ymax = CI_75, ymin = CI_25), fill="blue", fatten = 0.1, width = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash") +  
  scale_y_continuous(breaks = seq(-1,1,.1), labels = seq(-1,1,.1)) +
  labs(y="EBS Proportional change from base", x="") +
  coord_flip()

groupS5title<-paste("Experiment 7: Both decreased pollock recruitment and increased pollock fishing")
#grid.arrange(pGOA5, pEBS5, ncol=2)    

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("FigS5.pdf", width=10, height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(33, 19))) # 54 rows, 19 columns
print(pGOA5, vp = vplayout(2:6, 1:6))  # the big plot covers rows 1 and cols 1:4
print(pEBS5, vp = vplayout(2:33, 7:19))
grid.text(groupS5title, vp = vplayout(1:1, 1:19))
dev.off()

###############################################################


#########################   not using below here
#testing
with(pol.actor.sensedat, sum((X5.<0. )))
with(pol.actor.sensedat, sum((X95.<0.)))
with(pol.actor.sensedat, sum((X5.>0. )))
with(pol.actor.sensedat, sum((X95.>0.)))

with(pol.actor.sensedat, sum((X5.<0. & X95.<0.))) #total significant decrease with pollock increase
with(pol.actor.sensedat, sum((X5.>0. & X95.>0.))) #total significant increase with pollock increase


#not working below
def.sig = function(dat, lowsig, highsig)
{
  (dat$lowsig<0. & dat$highsig.<0.)|(dat$lowsig>0. & dat$highsig>0.)
}

find.sig = function(dat)
{
  nsig<-sum((dat$X5.<0. & dat$X95.<0.)|(dat$X5.>0. & dat$X95.>0.))
  return(nsig)
}

find.sig(pol.actor.sensedat)
find.sig(lgphy.actor.sensedat)

aggregate(sensedat, by=list(levels(sensedat$Actor)), find.sig)

tapply(sensedat, sensedat$Actor, find.sig)

