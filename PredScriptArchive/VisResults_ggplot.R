rm(list=ls())

# visualize predator results
# November 2016

library(tools)
library(dplyr)
library(ggplot2)

# read in saved rds files from CompileAllPredResults.R

allres <- readRDS("allres.rds")

# example self organizing map SOM visualization from
# https://en.wikibooks.org/wiki/Data_Mining_Algorithms_In_R/Clustering/Self-Organizing_Maps_(SOM)

#library("kohonen")
#data("wines")
#wines.sc <- scale(wines)
#set.seed(7)
#wine.som <- som(data = wines.sc, grid = somgrid(5, 4, "hexagonal"))
#plot(wine.som, main = "Wine data")

# possible workflow:
# boxplots across all metrics to see which are sensitive
# put sensitive ones into SOM?

# tradeoff plots with herring and other preds from allres
names(allres)

ternherrSSB <- ggplot(allres, aes(x=MedianSSB, y=MedPredProdtern, colour=CR)) +
  geom_point()

ternherrSSB + facet_wrap(c("Bias", "M", "Wt"))

ternherrYield <-ggplot(allres, aes(x=Yield, y=MedPredProdtern, colour=CR)) +
  geom_point()

ternherrYield + facet_wrap(c("Bias", "M", "Wt"))

#what is up with the increasing productivity in the loM BB3rPerc runs? find crnums to look at tern time series
highprod <- allres %>%
  filter(M=="LoM", CR=="BB3yrPerc", MedPredProdtern>1.1)

highprod$CRnum
# looks like the terns have been reduced very low and are recovering during the last 50 yrs of many time series
# may be exacerbated by herring B swings but can't tell; anyway legit because higher prod possible at low pop size

#plot tern pop size, prod, and herring pop size, yield?

ternpopherrSSB <- ggplot(allres, aes(x=MedianSSB, y=MedPredNtern, colour=CR)) +
  geom_point()

ternpopherrSSB + facet_wrap(c("Bias", "M", "Wt"))

ternpopherrSSBstat <- ggplot(allres, aes(x=MedSSBrelSSBmsy, y=MedPredNtern, colour=CR)) +
  geom_point()

ternpopherrSSBstat + facet_wrap(c("Bias", "M", "Wt"))
ggsave("ternN_herrSSBrelSSBmsy.png", scale=1.3)

ternpopherrB0stat <- ggplot(allres, aes(x=MedSSBrelSSBzero, y=MedPredNtern, colour=CR)) +
  geom_point()

ternpopherrB0stat + facet_wrap(c("Bias", "M", "Wt"))
ggsave("ternN_herrSSBrelSSB0.png", scale=1.3)


ternpopYield <- ggplot(allres, aes(x=Yield, y=MedPredNtern, colour=CR))+
  geom_point()

ternpopYield + facet_wrap(c("Bias", "M", "Wt"))

#after adding groundfish

gfpopherrSSB <- ggplot(allres, aes(x=MedianSSB, y=MedPredNgf, colour=CR)) +
  geom_point()

gfpopherrSSB + facet_wrap(c("Bias", "M", "Wt"))

gfstatherrSSB <- ggplot(allres, aes(x=MedianSSB, y=MedPredB_statusgf, colour=CR)) +
  geom_point()

gfstatherrSSB + facet_wrap(c("Bias", "M", "Wt"))

gfstatherrYield <- ggplot(allres, aes(x=Yield, y=MedPredB_statusgf, colour=CR)) +
  geom_point()

gfstatherrYield + facet_wrap(c("Bias", "M", "Wt"))

gf25statherrYield <- ggplot(allres, aes(x=Yield, y=Q25PredB_statusgf, colour=CR)) +
  geom_point()

gf25statherrYield + facet_wrap(c("Bias", "M", "Wt"))

gfstatternprod <- ggplot(allres, aes(x=MedPredProdtern, y=MedPredB_statusgf, colour=CR)) +
  geom_point()

gfstatternprod + facet_wrap(c("Bias", "M", "Wt"))

################### General herring results

herrrange <- ggplot(allres, aes(x=MedianSSB, y=SurpProd, colour=CR))+
  geom_point()

herrrange + facet_wrap(c("Bias", "M", "Wt"))

herryield <- ggplot(allres, aes(x=MedianSSB, y=Yield, colour=CR))+
  geom_point()

herryield + facet_wrap(c("Bias", "M", "Wt"))

herryieldrel <- ggplot(allres, aes(x=MedSSBrelSSBzero, y=YieldrelMSY, colour=CR))+
  geom_point()

herryieldrel + facet_wrap(c("Bias", "M", "Wt"))

herryieldterns <-ggplot(allres, aes(x=YieldrelMSY, y=MedPredProdtern, colour=CR))+
  geom_point()

herryieldterns + facet_wrap(c("Bias", "M", "Wt"))

MedPropYrs_goodProd_Targplustern

herryieldternprod <-ggplot(allres, aes(x=YieldrelMSY, y=MedPropYrs_goodProd_Targplustern, colour=CR))+
  geom_point()

herryieldternprod + facet_wrap(c("Bias", "M", "Wt"))

#things to plot for paper:
#stationary
stouffer
fisherchi
p50_NR
PropClosure
MedSSBrelSSBmsy
PropSSBrelhalfSSBmsy
YieldrelMSY
Yvar
MedPropYrs_goodProd_Targplustern #or Av?
MedPropYrs_okBstatusgf #or Av?
MedPropYrs_goodAvWtstatus #or Av?

metrics <- c('stouffer',
             'fisherchi',
             'p50_NR',
             'PropClosure',
             'MedSSBrelSSBmsy',
             'PropSSBrelhalfSSBmsy',
             'YieldrelMSY',
             'Yvar',
             'MedPropYrs_goodProd_Targplustern', 
             'MedPropYrs_okBstatus', 
             'MedPropYrs_goodAvWtstatus'
             )

metlabels <- c('stationary',
               'stableeq',
               'revenue',
               'closure',
               'relSSB',
               'overfished',
               'relyield',
               'yieldvar',
               'ternprod',
               'dogstatus',
               'tunawt'
               )

OMs <- unique(allres$OM)

OMlabels <- c('LowFastBiased', 
              'LowSlowBiased', 
              'LowFastCorrect', 
              'LowSlowCorrect',  
              'HighFastBiased',  
              'HighSlowBiased', 
              'HighFastCorrect',
              'HighSlowCorrect'
              )

#OM (8) x metric (11) boxplot--

plotting <- function(d, xvar, yvars, xlabs, ylabs){
  P <- list()
  for (i in 1:length(yvars)){
    if(i < length(yvars)){
      p <- ggplot(d, aes_string(x=xvar, y=yvars[i])) + 
      geom_boxplot() + labs(y=ylabs[i]) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
    if(i == length(yvars)){
      p <- ggplot(d, aes_string(x=xvar, y=yvars[i])) + 
        geom_boxplot() + labs(y=ylabs[i]) +
        scale_x_discrete(labels = xlabs)
    }
    P <- c(P, list(p))
  }
  return(list(plots=P, num=length(yvars)))
}

plotting <- function(d, xvar, yvars, xlabs, ylabs){
  P <- list()
  for (i in 1:length(yvars)){
    if(i < length(yvars)){
      p <- ggplot(d, aes_string(x=xvar, y=yvars[i])) + 
        geom_boxplot() + labs(y=ylabs[i]) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
    if(i == length(yvars)){
      p <- ggplot(d, aes_string(x=xvar, y=yvars[i])) + 
        geom_boxplot() + labs(y=ylabs[i]) +
        scale_x_discrete(labels = xlabs)
    }
    gp <- ggplotGrob(p)
    P <- c(P, list(gp))
  }
  return(list(plots=P, num=length(yvars)))
}



#PLOTS <- plotting(d)
#do.call(grid.arrange, c(PLOTS$plots, nrow = PLOTS$num))
library(gridExtra)
library(grid)

metOM <- plotting(allres, "OM", metrics, OMlabels, metlabels)

png("OMex2.png", width = 480, height = 980)
do.call(grid.arrange, c(metOM$plots, nrow = metOM$num))
dev.off()

# look at egg package and this page https://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends
gtable_rbind(metOM$plots, size="min")

statOM <- ggplot(allres, aes(x=OM, y=stouffer)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
fischrOM <- ggplot(allres, aes(x=OM, y=fisherchi)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p50NROM <- ggplot(allres, aes(x=OM, y=p50_NR)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
propclOM <- ggplot(allres, aes(x=OM, y=PropClosure)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
relSSBOM <- ggplot(allres, aes(x=OM, y=MedSSBrelSSBmsy)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
probOFOM <- ggplot(allres, aes(x=OM, y=PropSSBrelhalfSSBmsy)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
relyieldOM <- ggplot(allres, aes(x=OM, y=YieldrelMSY)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
yvarOM <- ggplot(allres, aes(x=OM, y=Yvar)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
probternOM <- ggplot(allres, aes(x=OM, y=MedPropYrs_goodProd_Targplustern)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
probdogOM <- ggplot(allres, aes(x=OM, y=MedPropYrs_okBstatus)) + geom_boxplot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
probtunaOM <- ggplot(allres, aes(x=OM, y=MedPropYrs_goodAvWtstatus)) + geom_boxplot() 


png("OMex.png", width = 480, height = 980)
grid.arrange(
statOM,
fischrOM,
p50NROM, 
propclOM, 
relSSBOM, 
probOFOM, 
relyieldOM, 
yvarOM,
probternOM, 
probdogOM, 
probtunaOM, nrow=11
)
dev.off()

#metric x CR boxplot--> select CR

#this one is still across all OMs, interesting that CR makes a difference for some
statCR <- ggplot(allres, aes(x=CR, y=stouffer)) + geom_boxplot() 
fischrCR <- ggplot(allres, aes(x=CR, y=fisherchi)) + geom_boxplot() 
p50NRCR <- ggplot(allres, aes(x=CR, y=p50_NR)) + geom_boxplot() 
propclCR <- ggplot(allres, aes(x=CR, y=PropClosure)) + geom_boxplot() 
relSSBCR <- ggplot(allres, aes(x=CR, y=MedSSBrelSSBmsy)) + geom_boxplot() 
probOFCR <- ggplot(allres, aes(x=CR, y=PropSSBrelhalfSSBmsy)) + geom_boxplot() 
relyieldCR <- ggplot(allres, aes(x=CR, y=YieldrelMSY)) + geom_boxplot() 
yvarCR <- ggplot(allres, aes(x=CR, y=Yvar)) + geom_boxplot() 
probternCR <- ggplot(allres, aes(x=CR, y=MedPropYrs_goodProd_Targplustern)) + geom_boxplot() 
probdogCR <- ggplot(allres, aes(x=CR, y=MedPropYrs_okBstatus)) + geom_boxplot() 
probtunaCR <- ggplot(allres, aes(x=CR, y=MedPropYrs_goodAvWtstatus)) + geom_boxplot() 

png("CRex.png", width = 480, height = 980)
grid.arrange(
  statCR,
  fischrCR,
  p50NRCR, 
  propclCR, 
  relSSBCR, 
  probOFCR, 
  relyieldCR, 
  yvarCR,
  probternCR, 
  probdogCR, 
  probtunaCR, nrow=11
)
dev.off()
#x-y tradeoff plots for 1-3 BB CRs: 


#moneystabternprob <- ggplot(allres, aes(x=AvPropYrs_goodProd_Targplustern, y=stationary, colour=CR))+
#  geom_point()

#moneystabternprob + facet_wrap(c("Bias", "M", "Wt"))

#moneyternprob <- ggplot(allres, aes(x=AvPropYrs_goodProd_Targplustern, y=nr, colour=CR))+
#  geom_point() +
#  #geom_hline(yintercept=) +
#  geom_vline(xintercept=0.9)

#moneyternprob + facet_wrap(c("Bias", "M", "Wt"))

#moneytunawtprob <- ggplot(allres, aes(x=AvPropYrs_goodAvWtstatus, y=nr, colour=CR))+
#  geom_point()

#moneytunawtprob + facet_wrap(c("Bias", "M", "Wt"))

ternstunaprob <- ggplot(allres, aes(x=AvPropYrs_goodAvWtstatus, y=AvPropYrs_goodProd_Targplustern, colour=CR))+
  geom_point()

ternstunaprob + facet_wrap(c("Bias", "M", "Wt"))

terns30unfished <- ggplot(allres, aes(x=MedSSBrelSSBzero, y=AvPropYrs_goodProd_Targplustern, colour=CR))+
  geom_point() +
  geom_hline(yintercept=0.9) +
  geom_vline(xintercept=0.3)

terns30unfished + facet_wrap(c("Bias", "M", "Wt"))

# how often does ABC fall below a threshold: 90,000 suggested
# how often does B fall below 40% SSB0 (or a numeric treshold like 400K Total B)
# need to look at time series to do this for real
# but try the 25%ile which I dont have for herring--dang

Fvsyield <- ggplot(allres, aes(x=YieldrelMSY, y=PropFrelFmsy, colour=CR))+
  geom_point() +
  geom_hline(yintercept=0.5) #+
  #geom_vline(xintercept=0.3)

Fvsyield + facet_wrap(c("Bias", "M", "Wt"))

BB_1_3 <- allres %>%
  filter(CR %in% c("BB", "BB3yr"))

Nclosures <- BB_1_3 %>%
  group_by(OM,CR) %>%
  filter(PropClosure<0.01) %>%
  summarize(under1p = n())

YieldvsClosure <- ggplot(BB_1_3, aes(x=YieldrelMSY, y=PropClosure, colour=CR))+
  geom_point() +
  geom_hline(yintercept=0.2) #+

YieldvsClosure + facet_wrap(c("Bias", "M", "Wt"))

