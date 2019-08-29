rm(list=ls())

# Tern data
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
#library(minpack.lm)

#setwd("~/Documents/0_Data/MSE/HerringMSE")
setwd("~/Data/Projects/MSE/HerringMSE")
ternprod<-read.csv("TernProdHerring.csv")

ternprod_long <- ternprod %>%
  gather(spp_site, chickprod, Common.1:Arctic.13) %>%
  separate(spp_site, into=c("species", "site"), sep="\\.")

#basic plotting
ptest1 <- ggplot(ternprod_long, aes(x=HerringSSB, y=chickprod, colour=species, group=site)) +
  geom_point() +
  ggtitle("Herring SSB and tern chick productivity")

ptest2 <- ggplot(subset(ternprod_long, species=="Common"), aes(x=HerringSSB, y=chickprod)) +
  geom_point() +
  ggtitle("Herring SSB and Common tern chick productivity")

ptest3 <- ggplot(subset(ternprod_long, species=="Arctic"), aes(x=HerringSSB, y=chickprod)) +
  geom_point() +
  ggtitle("Herring SSB and Arctic tern chick productivity")

#######################################################################
#is there a relationship between herring SSB and tern chick production?
#######################################################################
p1 <- ggplot(ternprod_long, aes(x=HerringSSB, y=chickprod)) +
  geom_point() 

p1 + facet_wrap(~species)

p1 + facet_wrap(~site)

sitedat <- ternprod_long %>%
  group_by(site) %>%
  do(fit1 = lm(chickprod ~ HerringSSB, data = .))

  
sitecoef <- tidy(sitedat, fit1)

filter(sitecoef, p.value<0.05)
#Source: local data frame [11 x 6]
#Groups: site [10]
#
#site        term     estimate    std.error statistic      p.value
#(chr)       (chr)        (dbl)        (dbl)     (dbl)        (dbl)
#1     10 (Intercept) 1.872729e+00 5.520519e-01  3.392306 6.859192e-03
#2     11 (Intercept) 9.664599e-01 4.287808e-01  2.253972 3.402887e-02
#3     12 (Intercept) 9.450428e-01 3.101058e-01  3.047485 5.714768e-03
#4     13 (Intercept) 7.581351e-01 2.842887e-01  2.666779 1.841548e-02
#5      2 (Intercept) 9.796870e-01 2.233022e-01  4.387270 1.304906e-04
#6      4 (Intercept) 7.622121e-01 1.150552e-01  6.624753 2.467679e-07
#7      5 (Intercept) 3.231874e-01 1.498014e-01  2.157439 3.970093e-02
#8      5  HerringSSB 7.874285e-07 2.650429e-07  2.970947 6.035239e-03
#9      7 (Intercept) 6.836047e-01 1.624374e-01  4.208420 2.265601e-04
#10     8 (Intercept) 9.529669e-01 3.619677e-01  2.632741 1.882927e-02
#11     9 (Intercept) 1.120223e+00 3.686809e-01  3.038463 9.508725e-03

sitecoef1<- sitecoef %>%
  group_by(site) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

#names(sitecoef1)[3] <- "Intercept" #platform specific?s
names(sitecoef1) <- gsub("(", "", names(sitecoef1), fixed = TRUE)
names(sitecoef1) <- gsub(")", "", names(sitecoef1), fixed = TRUE)


p1 + facet_wrap(~site) + geom_abline(aes(slope=HerringSSB, intercept=Intercept), sitecoef1)

p1.1 <- ggplot(ternprod_long, aes(x=Yr, y=HerringSSB/mean(HerringSSB, na.rm=T))) + geom_line() + geom_point(aes(y=chickprod/mean(chickprod, na.rm=T), colour=species))

p1.1 + facet_wrap(~site)

p1.1

#######################################################################
# again based on age1 recruitment
#######################################################################
p2 <- ggplot(ternprod_long, aes(x=HerringAge1Rec, y=chickprod)) +
  geom_point() 

p2 + facet_wrap(~species)

p2 + facet_wrap(~site)

sitedat2 <- ternprod_long %>%
  group_by(site) %>%
  do(fit2 = lm(chickprod ~ HerringAge1Rec, data = .))


sitecoef2 <- tidy(sitedat2, fit2)

filter(sitecoef2, p.value<0.05)
#Source: local data frame [11 x 6]
#Groups: site [11]
#
#site        term  estimate  std.error statistic      p.value
#(chr)       (chr)     (dbl)      (dbl)     (dbl)        (dbl)
#1     10 (Intercept) 1.3046335 0.21812329  5.981175 1.354706e-04
#2     11 (Intercept) 0.7807566 0.22406143  3.484565 2.001946e-03
#3     12 (Intercept) 1.2989558 0.17164683  7.567607 1.097802e-07
#4     13 (Intercept) 1.3199950 0.14886180  8.867251 4.052137e-07
#5      2 (Intercept) 0.7160657 0.11648973  6.147028 9.250733e-07
#6      4 (Intercept) 0.7600061 0.05911340 12.856748 9.771569e-14
#7      5 (Intercept) 0.7325215 0.09359484  7.826516 1.589399e-08
#8      6 (Intercept) 0.7731875 0.21359615  3.619857 1.606211e-03
#9      7 (Intercept) 0.8093000 0.08501915  9.519032 1.995033e-10
#10     8 (Intercept) 1.3416607 0.22971966  5.840426 3.251024e-05
#11     9 (Intercept) 1.1769952 0.20486455  5.745236 6.763109e-05


sitecoef2.1 <- sitecoef2 %>%
  group_by(site) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

#names(sitecoef2.1)[3] <- "Intercept"
names(sitecoef2.1) <- gsub("(", "", names(sitecoef2.1), fixed = TRUE)
names(sitecoef2.1) <- gsub(")", "", names(sitecoef2.1), fixed = TRUE)


p2 + facet_wrap(~site) + geom_abline(aes(slope=HerringAge1Rec, intercept=Intercept), sitecoef2.1)

p2.1 <- ggplot(ternprod_long, aes(x=Yr, y=HerringAge1Rec/mean(HerringAge1Rec, na.rm=T))) + geom_line() + geom_point(aes(y=chickprod/mean(chickprod, na.rm=T), colour=species))

p2.1 + facet_wrap(~site)

p2.1

#######################################################################
# again based on total biomass 
#######################################################################
p3 <- ggplot(ternprod_long, aes(x=HerringTotB, y=chickprod)) +
  geom_point() 

p3 + facet_wrap(~species)

p3 + facet_wrap(~site)

sitedat3 <- ternprod_long %>%
  group_by(site) %>%
  do(fit3 = lm(chickprod ~ HerringTotB, data = .))


sitecoef3 <- tidy(sitedat3, fit3)

filter(sitecoef3, p.value<0.05)
#Source: local data frame [6 x 6]
#Groups: site [6]

#site        term  estimate std.error statistic      p.value
#(chr)       (chr)     (dbl)     (dbl)     (dbl)        (dbl)
#1    12 (Intercept) 0.8505712 0.4037098  2.106888 4.623928e-02
#2    13 (Intercept) 0.8466264 0.3772976  2.243922 4.152692e-02
#3     2 (Intercept) 1.1014752 0.2829083  3.893399 5.113907e-04
#4     4 (Intercept) 0.9142651 0.1464078  6.244648 7.050952e-07
#5     5 (Intercept) 0.4771900 0.2125313  2.245270 3.283104e-02
#6     7 (Intercept) 0.7837919 0.2084806  3.759543 7.653867e-04


sitecoef3.1 <- sitecoef3 %>%
  group_by(site) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

#names(sitecoef3.1)[3] <- "Intercept"
names(sitecoef3.1) <- gsub("(", "", names(sitecoef3.1), fixed = TRUE)
names(sitecoef3.1) <- gsub(")", "", names(sitecoef3.1), fixed = TRUE)


p3 + facet_wrap(~site) + geom_abline(aes(slope=HerringTotB, intercept=Intercept), sitecoef3.1)

p3.1 <- ggplot(ternprod_long, aes(x=Yr, y=HerringTotB/mean(HerringTotB, na.rm=T))) + geom_line() + geom_point(aes(y=chickprod/mean(chickprod, na.rm=T), colour=species))

p3.1 + facet_wrap(~site)

p3.1

######################################################################
# can we do any nonlinear fits to Bev Holt function?
######################################################################


sitedatBH <- ternprod_long %>%
  group_by(site) %>%
  filter(!is.na(chickprod)) %>%
  filter(site > 1, site != 11) %>%
  #select(chickprod, HerringAge1Rec) %>%
  do(mod2 = nls((1/chickprod) ~ betastar + alphastar*(1/HerringAge1Rec), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

sitecoefBH <- tidy(sitedatBH, mod2)

filter(sitecoefBH, p.value<0.05)
#Source: local data frame [10 x 6]
#Groups: site [9]
#
#site      term     estimate    std.error statistic      p.value
#<chr>     <chr>        <dbl>        <dbl>     <dbl>        <dbl>
#1     10  betastar 7.093499e-01 1.843022e-01  3.848840 3.218306e-03
#2     12  betastar 1.173062e+00 2.758501e-01  4.252532 3.000652e-04
#3     13  betastar 9.245447e-01 1.694027e-01  5.457673 8.439150e-05
#4      2  betastar 1.274005e+00 3.431008e-01  3.713209 8.343201e-04
#5      4  betastar 1.164004e+00 1.822350e-01  6.387378 4.746859e-07
#6      5 alphastar 5.766466e+06 2.807976e+06  2.053602 4.946019e-02
#7      5  betastar 1.161232e+00 2.389693e-01  4.859335 4.074194e-05
#8      7  betastar 1.640927e+00 2.825325e-01  5.807922 2.701184e-06
#9      8  betastar 1.298682e+00 2.543989e-01  5.104904 1.293251e-04
#10     9  betastar 1.039144e+00 4.095681e-01  2.537170 2.478674e-02

sitecoefBH.1 <- sitecoefBH %>%
  group_by(site) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

names(sitecoefBH.1) <- gsub("(", "", names(sitecoefBH.1), fixed = TRUE)
names(sitecoefBH.1) <- gsub(")", "", names(sitecoefBH.1), fixed = TRUE)

spdatBH <- ternprod_long %>%
  group_by(species) %>%
  filter(!is.na(chickprod)) %>%
  filter(site > 1, site != 11) %>%
  do(mod3 = nls((1/chickprod) ~ betastar + alphastar*(1/HerringAge1Rec), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

spcoefBH <- tidy(spdatBH, mod3)

spcoefBH.1 <- spcoefBH %>%
  group_by(species) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

names(spcoefBH.1) <- gsub("(", "", names(spcoefBH.1), fixed = TRUE)
names(spcoefBH.1) <- gsub(")", "", names(spcoefBH.1), fixed = TRUE)


BHfun <- function(x,a,b){
  (x/a)/(1+(1/a*b)*x)
}

pBHdat <- left_join(ternprod_long, sitecoefBH.1) %>%
  mutate(BHfit = BHfun(HerringAge1Rec, alphastar, betastar))

pBHspdat <- left_join(ternprod_long, spcoefBH.1) %>%
  mutate(BHfit = BHfun(HerringAge1Rec, alphastar, betastar))
  

pBH <- ggplot(pBHdat, aes(x=HerringAge1Rec, y=chickprod, colour=species)) +
  geom_point() +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(limits = c(0,NA))
             
pBH + facet_wrap(~site) + 
  #geom_line(aes(x=HerringAge1Rec, y=BHfit, col="BHfit")) +
  #geom_smooth(aes(x=HerringAge1Rec, y=BHfit, col="BHfit"), fullrange=T)
  stat_smooth(aes(x=HerringAge1Rec, y=BHfit, col="BHfit"), fullrange=TRUE)


pBHsp <- ggplot(pBHspdat, aes(x=HerringAge1Rec, y=chickprod, colour=species)) +
  geom_point() +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(limits = c(0,NA))

pBHsp + facet_wrap(~species) + 
  #geom_line(aes(x=HerringAge1Rec, y=BHfit, col="BHfit"))
  geom_smooth(aes(x=HerringAge1Rec, y=BHfit, col="BHfit"), fullrange=T)
  #stat_smooth(aes(col="BHfit"), formula = y ~ BHfun(HerringAge1Rec, alphastar, betastar), fullrange=TRUE)


pBH.1 <- ggplot(ternprod_long, aes(x=Yr, y=HerringAge1Rec/mean(HerringAge1Rec, na.rm=T))) + geom_line() + geom_point(aes(y=chickprod/mean(chickprod, na.rm=T), colour=species))

pBH.1 + facet_wrap(~site)



#brute force loop testing below

#testBH <- ternprod_long %>%
#  filter(site==5) %>%
#  filter(!is.na(chickprod)) #%>%
#  #select(chickprod, HerringAge1Rec)
#
#mod2<-with(testBH, nls((1/chickprod) ~ betastar + alphastar*(1/HerringAge1Rec), 
#                    start=list(alphastar = 1.0 , betastar = 1.0)))
#
#with(testBH, plot(HerringAge1Rec, chickprod, xlim=c(0,max(HerringAge1Rec)), ylim=c(0,max(chickprod, na.rm=T)), main="BH test site 5"))
#
#alphaBH<-1/coef(mod2)[1]
#betaBH<-1/coef(mod2)[1]*coef(mod2)[2]
#with(testBH, lines(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100), alphaBH*(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))/(1+(betaBH*seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))), col="blue"))
#

#for(i in unique(ternprod_long$site)) {
#  dat<-filter(ternprod_long, site==i, !is.na(chickprod))
#  mod2<-with(dat, nls((1/chickprod) ~ betastar + alphastar*(1/HerringAge1Rec), 
#                         start=list(alphastar = 1.0 , betastar = 1.0)))
#  
#  with(dat, plot(HerringAge1Rec, chickprod, xlim=c(0,max(HerringAge1Rec)), ylim=c(0,max(chickprod, na.rm=T)), main=paste("BH test site ", i)))
#  
#  alphaBH<-1/coef(mod2)[1]
#  betaBH<-1/coef(mod2)[1]*coef(mod2)[2]
#  with(dat, lines(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100), alphaBH*(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))/(1+(betaBH*seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))), col="blue"))
#
#}

#HerringSSB
#function breaks entirely for sites 1, 11  

par(mfrow=c(3,4))
for(i in c(2:10, 12:13)) {
  dat<-filter(ternprod_long, site==i, !is.na(chickprod))
  mod2<-with(dat, nls((1/chickprod) ~ betastar + alphastar*(1/HerringSSB), 
                      start=list(alphastar = 1.0 , betastar = 1.0)))
  
  with(dat, plot(HerringSSB, chickprod, xlim=c(0,max(HerringSSB)), ylim=c(0,max(chickprod, na.rm=T)), main=paste("BH test site ", i)))
  
  alphaBH<-1/coef(mod2)[1]
  betaBH<-1/coef(mod2)[1]*coef(mod2)[2]
  with(dat, lines(seq(0,max(HerringSSB),max(HerringSSB)/100), alphaBH*(seq(0,max(HerringSSB),max(HerringSSB)/100))/(1+(betaBH*seq(0,max(HerringSSB),max(HerringSSB)/100))), col="blue"))
  
}

#Herring age 1 rec 
par(mfrow=c(3,4))
for(i in c(2:10, 12:13)) {
  dat<-filter(ternprod_long, site==i, !is.na(chickprod))
  mod2<-with(dat, nls((1/chickprod) ~ betastar + alphastar*(1/HerringAge1Rec), 
                      start=list(alphastar = 1.0 , betastar = 1.0)))
  
  with(dat, plot(HerringAge1Rec, chickprod, xlim=c(0,max(HerringAge1Rec)), ylim=c(0,max(chickprod, na.rm=T)), main=paste("BH test site ", i)))
  
  alphaBH<-1/coef(mod2)[1]
  betaBH<-1/coef(mod2)[1]*coef(mod2)[2]
  with(dat, lines(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100), alphaBH*(seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))/(1+(betaBH*seq(0,max(HerringAge1Rec),max(HerringAge1Rec)/100))), col="blue"))
  
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(p2, p3)
