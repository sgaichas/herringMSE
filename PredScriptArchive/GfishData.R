rm(list=ls())
# groundfish data analysis for herring MSE
# November 2016

library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)

setwd("~/Documents/0_Data/MSE/HerringMSE")
#setwd("~/Data/Projects/MSE/HerringMSE")

codGB<-read.csv("GBcodpop.csv")
codGOM<-read.csv("GOMcodpop.csv")
dogfish<-read.csv("dogfishpop.csv")
gfishdiet<-read.csv("Gfishdiets.csv")
herring <- read.csv("herringassessout2015.csv")
herring <- herring[,1:6]
dogfishrec <- read.csv("dogfishrec.csv")

# add species column, tidy and combine datasets
codGB <- cbind(codGB, species=rep("GBcod", dim(codGB)[1]))
codGOM <- cbind(codGOM, species=rep("GOMcod", dim(codGOM)[1]))
dogfish <- cbind(dogfish, species=rep("dogfish", dim(dogfish)[1]))
dogfish <- full_join(dogfish, dogfishrec)

gfishdiet <- gfishdiet %>%
  gather(preypred, dietprop, 2:9) %>%
  separate(preypred, c("prey", "species"), sep="_")

codpops <- full_join(codGB, codGOM)
names(codpops)[1]<-"Yr"
names(dogfish)[2:12] <- c("LargeB", "MedB", "RecAge1t",
                          "TotJan1B","survSSB", "survSSB_MA","stocExpB",
                          "stocTotB", "stocSSB", "SSBkgtow", "SSB")
#put all in t not 1000t
dogfish <- dogfish %>%
  mutate_at(vars(LargeB, MedB, RecAge1t,TotJan1B, survSSB, survSSB_MA,
                 stocExpB, stocTotB, stocSSB, SSB),
            funs(.*1000))

dogfish <- dogfish %>%
  mutate(RecAge1 = RecAge1t/PupAvgWtkg) %>% #rec in 1000 fish matches cod units
  mutate(TotJan1B = replace(TotJan1B, TotJan1B==0, NA)) #2014 missed survey

gfishpops <- full_join(codpops, dogfish)

#plots to see if this worked
gfishSSB <- ggplot(gfishpops, aes(x=Yr, y=SSB, colour=species)) + geom_point()
gfishTotB <- ggplot(gfishpops, aes(x=Yr, y=TotJan1B, colour=species)) + geom_point()
gfishRec <- ggplot(gfishpops, aes(x=Yr, y=RecAge1, colour=species)) + geom_point()

gfishSSB #+ geom_point(aes(x=Yr, y=survSSB, colour="survSSB dogfish")) +
  #geom_point(aes(x=Yr, y=stocSSB, colour="stocSSB dogfish"))
gfishTotB
gfishRec

gfishfood <- ggplot(gfishdiet, aes(x=Year, y=dietprop, colour=species)) + geom_point()
gfishfood + facet_wrap(~prey)
ggsave("gfishfood.png")

#merge fishpops and fish diets
gfishpopfood <- left_join(gfishpops, gfishdiet, by=c("Yr"="Year", "species"="species"))

#is gfish rec or B related to herring prop in diet?

gfishrecdiet <- ggplot(gfishpopfood, aes(x=dietprop, y=RecAge1, colour=prey)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) 
gfishrecdiet + facet_wrap(~species)

gfishSSBdiet <- ggplot(gfishpopfood, aes(x=dietprop, y=SSB, colour=prey)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) 
gfishSSBdiet + facet_wrap(~species)

gfishtotBdiet <- ggplot(gfishpopfood, aes(x=dietprop, y=TotJan1B, colour=prey)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) 
gfishtotBdiet + facet_wrap(~species)

# merge herring assessment out with gfish
gfishpopfoodherr <- left_join(gfishpopfood, herring, by=c("Yr"="yr"))

#is herring in the ecosystem related to herring in diet?
gfishdiet_herrtotB <- ggplot(gfishpopfoodherr, aes(x=HerringTotalB, y=dietprop, colour=prey))+ 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Herring Total B and diet proportion")
gfishdiet_herrtotB + facet_wrap(~species)

gfishdiet_herrSSB <- ggplot(gfishpopfoodherr, aes(x=HerringSSB, y=dietprop, colour=prey))+ 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Herring SSB and diet proportion")
gfishdiet_herrSSB + facet_wrap(~species)

gfishdiet_herrrec <- ggplot(gfishpopfoodherr, aes(x=HerringAge1Rec, y=dietprop, colour=prey))+ 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Herring Recruitment and diet proportion")
gfishdiet_herrrec + facet_wrap(~species)

# which relationships are signficant?
# herring in diet and SSB, TotB, Rec
sppdat <- gfishpopfoodherr %>%
  group_by(species) %>%
  filter(prey == "CLUHAR") %>%
  #filter(prey == "CLUALL") %>%
  #do(tidy(cor.test(.$dietprop, .$RecAge1,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$SSB,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$TotJan1B,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$RecAge1,alternative="greater", method="pearson", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$SSB,alternative="greater", method="pearson", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$TotJan1B,alternative="greater", method="pearson", na.action=na.rm)))
  #sppdat
  #do(fitHerDiet = lm(RecAge1 ~ dietprop, data = .))
  do(fitHerDiet = lm(SSB ~ dietprop, data = .))
  #do(fitHerDiet = lm(TotJan1B ~ dietprop, data = .))
sppcoef <- tidy(sppdat, fitHerDiet)
filter(sppcoef, p.value<0.1)

#correlation results:
# only dogfish spearmans rank correlations with TotJan1B (survB)
# are significant at p<0.05 (both CLUALL and CLUHAR respectively below), 
#species   estimate statistic    p.value                          method
#<chr>      <dbl>     <dbl>      <dbl>                          <fctr>
#1 dogfish  0.3638663  5813.626 0.01236172 Spearman's rank correlation rho
#2   GBcod -0.5043053 12690.320 0.99927029 Spearman's rank correlation rho
#3  GOMcod -0.4616158  8746.309 0.99657598 Spearman's rank correlation rho
#species   estimate statistic    p.value                          method
#<chr>      <dbl>     <dbl>      <dbl>                          <fctr>
#1 dogfish  0.3086940  6317.846 0.02966801 Spearman's rank correlation rho
#2   GBcod -0.4292952 12057.534 0.99599262 Spearman's rank correlation rho
#3  GOMcod -0.1888369  7114.000 0.85477167 Spearman's rank correlation rho

# pearsons correlation significant at p<0.05 for CLUHAR, P<0.1 for CLUALL
#species   estimate statistic    p.value parameter   conf.low conf.high
#<chr>      <dbl>     <dbl>      <dbl>     <int>      <dbl>     <dbl>
#1 dogfish  0.3182772  2.014417 0.02574168        36  0.0516519         1
#2   GBcod -0.4895355 -3.321315 0.99894775        35 -0.6737284         1
#3  GOMcod -0.3377073 -1.997634 0.97270176        31 -0.5728874         1
#species   estimate  statistic    p.value parameter    conf.low conf.high
#<chr>      <dbl>      <dbl>      <dbl>     <int>       <dbl>     <dbl>
#1 dogfish  0.2285908  1.4088473 0.08372977        36 -0.04529794         1
#2   GBcod -0.5096583 -3.5044867 0.99936354        35 -0.68811055         1
#3  GOMcod -0.1069803 -0.5990794 0.72326342        31 -0.38651727         1


#lm results:
# all intercepts significant at p<0.05, 
# only slope for GB cod significant at p<0.05 but negative
# for CLUHAR with TotJan1B
# all slopes significant at p<0.1, dogfish positive, GOM cod negative
#species        term    estimate   std.error statistic      p.value
#<chr>       <chr>       <dbl>       <dbl>     <dbl>        <dbl>
#1 dogfish (Intercept) 450760.7294 50505.78790  8.924932 1.183714e-10
#2 dogfish    dietprop  10887.6686  5404.87223  2.014417 5.148337e-02
#3   GBcod (Intercept)  71227.9447  8527.44377  8.352790 7.524678e-10
#4   GBcod    dietprop  -2847.7154   857.40604 -3.321315 2.104508e-03
#5  GOMcod (Intercept)  18985.0974  1848.41522 10.271013 1.685662e-11
#6  GOMcod    dietprop   -185.9602    93.09024 -1.997634 5.459648e-02

#does herring in diet reflect herring assessed pop?
sppdat <- gfishpopfoodherr %>%
  group_by(species) %>%
  filter(prey == "CLUHAR") %>%
  #filter(prey == "CLUALL") %>%
  #do(tidy(cor.test(.$dietprop, .$HerringAge1Rec,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$HerringSSB,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$HerringTotalB,alternative="greater", method="spearman", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$HerringAge1Rec,alternative="greater", method="pearson", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$HerringSSB,alternative="greater", method="pearson", na.action=na.rm)))
  #do(tidy(cor.test(.$dietprop, .$HerringTotalB,alternative="greater", method="pearson", na.action=na.rm)))
  #sppdat
  #do(fitHerDiet = lm(dietprop ~ HerringAge1Rec, data = .))
  #do(fitHerDiet = lm(dietprop ~ HerringSSB, data = .))
  do(fitHerDiet = lm(dietprop ~ HerringTotalB, data = .))
sppcoef <- tidy(sppdat, fitHerDiet)
filter(sppcoef, p.value<0.05)

#correlation results:
# all correlations significantly positive except for diet/herringrec pearsons

#lm results:
# CLUHAR
# intercepts only significant for herringrec predicting diet prop
# slopes only significant with herringSSB predicting diet comp
# slopes for dogfish and GBcod significant for herringTotB predicting diet comp
# CLUALL
# intercepts all and slope for GBcod herring rec predicting diet comp
# dogfish and GBcod slopes significant for herringSSB predicting diet comp
# GBcod slope only significant for herringTotalB predictig diet comp

# Conclusion: Herring eaten in proportion to B. 
# all preds eat herring in proportion to herring pop by some measure
# Conculsion 2: NOT USING COD. 
# cod relationships with herring in diet are negative, probably spurious
# Conclusion 3: USE DOGFISH
# dogfish total B has positive relationship to herring in diet, SSB weak positive
# dogfish recruitment not related to herring diet prop
# will herring affect growth or survival?

#Dogfish population model
dogfishfoodherr <- gfishpopfoodherr %>%
  filter(species=="dogfish") 

dogfishSR <- ggplot(dogfishfoodherr, aes(x=SSB, y=RecAge1, colour="kalmanSSB")) + geom_point()
dogfishSR + geom_point(aes(x=survSSB, y=RecAge1, colour="survSSB")) +
  geom_point(aes(x=stocSSB, y=RecAge1, colour="stocSSB"))

dogfishpoptrend <- ggplot(dogfishfoodherr, aes(x=Yr, y=TotJan1B)) + geom_point()
dogfishdiet <- ggplot(dogfishfoodherr, aes(x=Yr, y=dietprop, colour=prey)) + geom_point()
pupwtdiet <- ggplot(dogfishfoodherr, aes(x=dietprop, y=PupAvgWtkg, colour=prey)) + geom_point()
Fwtdiet <- ggplot(dogfishfoodherr, aes(x=dietprop, y=FAvgWtkg, colour=prey)) + geom_point()

dogfishdiet
pupwtdiet
Fwtdiet

dogfishpoptrend + geom_point(aes(x=Yr, y=SSB, colour="SSB")) +
  geom_point(aes(x=Yr, y=stocSSB, colour="stocSSB")) + 
  geom_point(aes(x=Yr, y=survSSB, colour="survSSB"))  
  
dogfishlow <- dogfishfoodherr %>%
  filter(Yr>1994 & Yr<2008)

lowdogfishdiet <- ggplot(dogfishlow, aes(x=Yr, y=dietprop, colour=prey)) + geom_point()
lowpupwtdiet <- ggplot(dogfishlow, aes(x=dietprop, y=PupAvgWtkg, colour=prey)) + geom_point()
lowFwtdiet <- ggplot(dogfishlow, aes(x=dietprop, y=FAvgWtkg, colour=prey)) + geom_point()

lowdogfishdiet + stat_smooth(method = "lm", formula = y ~ x, size = 1)
lowpupwtdiet + stat_smooth(method = "lm", formula = y ~ x, size = 1)
pupwtdiet + stat_smooth(method = "lm", formula = y ~ x, size = 1)

pupwttrend <- ggplot(dogfishfoodherr, aes(x=Yr, y=PupAvgWtkg)) + geom_point() +
  ggtitle("Dogfish average pup weight from NEFSC surveys")
pupwttrend

##### Dogfish SR fit to full dataset 
#try a different parameterization with only 2 pars
BHfun <- function(x,a,b){
  (x/a)/(1+(1/a*b)*x)
}

Rfun <- function(x,a,b){
 exp(a)*x*exp(-b*x)
}

#Rsteep <- function(Sdat, steep, Smax){
#    Smax*exp(1.25*log(5*steep)*(1-(Sdat/Smax)))/1000
#}

BHsteep <- function(Sdat, steep, Rmax, Smax){
  4*steep*Rmax*Sdat/
    (Smax*(1-steep)+Sdat*(5*steep-1))
}


alldatBH <- dogfishfoodherr %>%
  filter(!is.na(SSB)) %>%
  do(mod2 = nls((1/(RecAge1)) ~ betastar + alphastar*(1/SSB), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

coefBH <- tidy(alldatBH, mod2)

coefBH
#term     estimate    std.error statistic    p.value
#1 alphastar 1.372062e+01 5.543059e+00 2.4752800 0.01815507
#2  betastar 3.536811e-05 9.537265e-05 0.3708412 0.71292915

alldatBH2 <- dogfishfoodherr %>%
  filter(!is.na(stocSSB)) %>%
  do(mod2 = nls((1/(RecAge1)) ~ betastar + alphastar*(1/stocSSB), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

coefBH2 <- tidy(alldatBH2, mod2)

coefBH2
#term      estimate    std.error  statistic      p.value
#1 alphastar  2.568637e+01 7.136382e+00  3.5993548 0.0009526774
#2  betastar -3.740767e-05 8.850028e-05 -0.4226842 0.6750398224

alldatBH3 <- dogfishfoodherr %>%
  filter(!is.na(survSSB)) %>%
  do(mod2 = nls((1/(RecAge1)) ~ betastar + alphastar*(1/survSSB), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

coefBH3 <- tidy(alldatBH3, mod2)

coefBH3
#term     estimate    std.error statistic     p.value
#1 alphastar 1.412076e+01 4.0939551618 3.4491739 0.001098166
#2  betastar 4.447005e-05 0.0000509195 0.8733404 0.386344803

dogfishSR + geom_point(aes(x=survSSB, y=RecAge1, colour="survSSB")) +
  geom_point(aes(x=stocSSB, y=RecAge1, colour="stocSSB")) +
  stat_function(fun=BHfun, args=list(a= 1.372062e+01,b=3.536811e-05), col="red")+
  stat_function(fun=BHfun, args=list(a= 2.568637e+01,b=-3.740767e-05), col="green")+
  stat_function(fun=BHfun, args=list(a= 1.412076e+01,b= 4.447005e-05), col="blue")

# Ricker fits which Paul used for ref pts
alldatR <- dogfishfoodherr %>%
  filter(!is.na(SSB)) %>%
  do(mod2=nls(log(RecAge1/SSB) ~ a - beta * SSB,
                    start=list(a = 1.0 , beta = 1.0), data=.))
coefR <- tidy(alldatR, mod2)

coefR
#term      estimate    std.error statistic      p.value
#1    a -2.859753e+00 4.510571e-01 -6.340113 2.446160e-07
#2 beta -8.908488e-06 4.932608e-06 -1.806040 7.927872e-02

alldatR2 <- dogfishfoodherr %>%
  filter(!is.na(SSB)) %>%
  do(mod2=nls(log(RecAge1/stocSSB) ~ a - beta * stocSSB,
              start=list(a = 1.0 , beta = 1.0), data=.))
coefR2 <- tidy(alldatR2, mod2)

coefR2
#term      estimate    std.error statistic      p.value
#1    a -2.982499e+00 3.714417e-01 -8.029521 1.539837e-09
#2 beta -4.113589e-06 2.619958e-06 -1.570098 1.251413e-01

alldatR3 <- dogfishfoodherr %>%
  filter(!is.na(SSB)) %>%
  do(mod2=nls(log(RecAge1/survSSB) ~ a - beta * survSSB,
              start=list(a = 1.0 , beta = 1.0), data=.))
coefR3 <- tidy(alldatR3, mod2)

coefR3
#term      estimate    std.error  statistic      p.value
#1    a -2.390590e+00 3.834186e-01 -6.2349328 3.376910e-07
#2 beta  6.005462e-07 2.494492e-06  0.2407489 8.111155e-01

dogfishSR + geom_point(aes(x=survSSB, y=RecAge1, colour="survSSB")) +
  geom_point(aes(x=stocSSB, y=RecAge1, colour="stocSSB")) +
  stat_function(fun=Rfun, args=list(a= -2.859753e+00, b=-8.908488e-06), col="red")+
  stat_function(fun=Rfun, args=list(a= -2.982499e+00, b=-4.113589e-06), col="green")+
  stat_function(fun=Rfun, args=list(a= -2.390590e+00, b=6.005462e-07), col="blue")

#invent based on Pauls ests a=0.005 or 5 in pups/SSB in tons
# see big preds spreadsheet--this translates to steepness of ~0.97 
# this is modeling only adult females as whole population
# probably ok because also most of the catch!
# visually not wonderful comparison but not increasing matches
# ricker used by Paul, eq recruitment not terrible,
# and this model stinks for dogfish doesn't it

dogfishSR + #geom_point(aes(x=survSSB, y=RecAge1, colour="survSSB")) +
  #geom_point(aes(x=stocSSB, y=RecAge1, colour="stocSSB")) +
  stat_function(fun=BHsteep, args=list(steep=0.97, Rmax=12705, Smax=300000), col="black") +
  stat_function(fun=BHfun, args=list(a= 1.372062e+01,b=3.536811e-05), col="red")+
  stat_function(fun=Rfun, args=list(a= -2.859753e+00, b=-8.908488e-06), col="red", lty=3)
  
# what relationship to herring?
# dogfish total B reacted + to herring in diet (a coincidence?)
# recruitment not at all
# tiny weak evidence of pup avg wt increasing with herring in diet
# but only during recent period of crappy recruitment!
# no relationship with adult female avg wt
# pup avg wt lowest since early 1990s--energy density change?

# weak impact on growth combined with weak impact on survival?

dogherrB <- ggplot(dogfishfoodherr, aes(x=HerringTotalB, y=TotJan1B)) +
  geom_point() +
  stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
  stat_smooth(method = "lm", formula = y ~ x, size = 1)

dogherrSSB <- ggplot(dogfishfoodherr, aes(x=HerringSSB, y=SSB)) +
  geom_point() +
  stat_smooth(method = "loess", formula = y ~ x, size = 1) + 
  stat_smooth(method = "lm", formula = y ~ x, size = 1)

#compare sim with and without M effect
# warning only works with output of MSE_GroundfishLink in memory
dogsim <- data.frame(simB1 = PredBNoHerrNoFish, simB2 = PredBHerrNoFish)
SSB <- c(dogfish$SSB[44:49],rep(0,144))
  
dogsim <-cbind(dogsim, SSB, simyr = seq(1:150))

#dogsim <- dogsim %>%
# mutate(SSB = replace(SSB, SSB==0, NA)) #2014 missed survey

comparesims <- ggplot(dogsim, aes(x=simyr-5, y=SSB)) + geom_point()

comparesims + geom_line(aes(x=simyr, y=simB1, colour="noHerr")) + 
  geom_line(aes(x=simyr, y=simB2, colour="HerrM"))

