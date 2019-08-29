rm(list=ls())

# Tern data
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
#library(minpack.lm)

#setwd("~/Documents/0_Data/MSE/HerringMSE")
setwd("~/Data/Projects/MSE/HerringMSE")
#ternprod<-read.csv("TernProdHerring.csv")
#ternpopproddiet<-read.csv("TernPopProdDiet.csv")
ternpopproddiet<-read.csv("GOMTernPopsProdDiet_Corrected.csv")

names(ternpopproddiet)
#reconcile names
names(ternpopproddiet)[5]<-"productivity"

#get rid of trailing white spaces
ternpopproddiet$majority.diet <- sub(" +$", "", ternpopproddiet$majority.diet)

#ternprod_long <- ternprod %>%
#  gather(spp_site, chickprod, Common.1:Arctic.13) %>%
#  separate(spp_site, into=c("species", "site"), sep="\\.")

herring <- read.csv("herringassessout2015.csv")
herring <- herring[,1:6]

herringMENH <- read.csv("ME-NHtrawl_herring.csv")

ternpopproddiet_herring <- left_join(ternpopproddiet, herring)
ternpopproddiet_herring <- left_join(ternpopproddiet_herring, herringMENH)


###########################################################################################
#  Exloring the data, just plotting by colony
###########################################################################################

ternpops <- ggplot(ternpopproddiet, aes(x=yr, y=breedpairs, colour=species)) +
  geom_point() +
  ggtitle("Tern breeding pairs GOM")

ternpops + facet_wrap(~colony)
ggsave("ternpops.png")

ternprod <- ggplot(ternpopproddiet, aes(x=yr, y=productivity, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  ggtitle("Tern productivity GOM")

ternprod + facet_wrap(~colony)
ggsave("ternprod.png")

terndiet <- ggplot(ternpopproddiet, aes(x=yr, y=prop.herring, colour=species)) +
  geom_point() +
  ggtitle("Herring proportion in tern diet GOM")

terndiet + facet_wrap(~colony)
ggsave("terndiet.png")

ternproddiet <- ggplot(ternpopproddiet, aes(x=prop.herring, y=productivity, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Herring proportion in diet and tern productivity GOM")

ternproddiet + facet_wrap(~colony) #+ theme(text=element_text(size=16))

ggsave("ternproddiet.png")

# what is productivity by major diet item?
#with(ternpopproddiet, boxplot(productivity ~ majority.diet, las=3))
#abline(h=1.0, lty=3)
#title("productivity by majority diet item")

majorterndiet <- ternpopproddiet %>% 
  filter(!majority.diet %in% c("", "hake and sandlance",
                             "hake and butterfish",
                             "herring and hake",
                             "herring present"))

dietitemall <- ggplot(ternpopproddiet, aes(x=majority.diet, y=productivity, colour=species)) +
  geom_hline(yintercept = 1) + 
  geom_boxplot()

dietitemall + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))                   

dietitem <- ggplot(majorterndiet, aes(x=majority.diet, y=productivity, colour=species)) +
  geom_hline(yintercept = 1) + 
  geom_boxplot()

dietitem + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))                   

dietitem + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_jitter(alpha=0.3) +
  geom_boxplot(alpha=0)

#-------------------------------------------------------------------------------
# is the relationship between herring in diet and tern productivity significant?
# common terns only since the model is mainly for them
# correlation test (positive only), 2 of 13 colonies with p<0.05

colonydat <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$prop.herring, .$productivity,alternative="greater", method="spearman", na.action=na.rm)))

#colonydat
#Source: local data frame [13 x 6]
#Groups: colony [13]
#
#colony                 estimate statistic    p.value
#<fctr>                   <dbl>     <dbl>      <dbl>
#1  Eastern Egg Rock  0.27840489 959.72149 0.11730349
#2           Jenny I -0.09278361 743.09285 0.63373696
#3      Machias Seal  0.62870323  81.68529 0.01913475
#4    Matinicus Rock  0.44858432 200.71531 0.06208615
#5           Metinic -0.34505495 612.00000 0.88952018
#6           Monomoy  0.01801875  54.99095 0.48470769
#7      Outer green   0.06611595 339.93379 0.41504051
#8       Petit Manan  0.25889590 604.74095 0.15783610
#9            Pond I  0.23235294 522.00000 0.19258403
#10           Seal I  0.38727858 223.03060 0.09553945
#11           Ship I  0.30952381  58.00000 0.23090278
#12         Stratton  0.51719444 175.74122 0.03514742
#13 White and Seavey -0.05361932 590.02682 0.57526155

filter(colonydat, p.value<0.05)
#Source: local data frame [2 x 6]
#Groups: colony [2]
#
#colony  estimate statistic    p.value                          method
#<fctr>     <dbl>     <dbl>      <dbl>                          <fctr>
#1 Machias Seal 0.6287032  81.68529 0.01913475 Spearman's rank correlation rho
#2     Stratton 0.5171944 175.74122 0.03514742 Spearman's rank correlation rho
# ... with 1 more variables: alternative <fctr>


#linear models are non-significant slopes for all colonies using just common terns

colonydat <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerDiet = lm(productivity ~ prop.herring, data = .))
colonycoef <- tidy(colonydat, fitHerDiet)
filter(colonycoef, p.value<0.05)

#Source: local data frame [11 x 6]
#Groups: colony [11]
#
#colony        term  estimate std.error statistic      p.value
#<fctr>       <chr>     <dbl>     <dbl>     <dbl>        <dbl>
#  1  Eastern Egg Rock (Intercept) 0.9295303 0.1234340  7.530587 5.734943e-07
#2           Jenny I (Intercept) 1.5943009 0.3916569  4.070657 1.145995e-03
#3    Matinicus Rock (Intercept) 0.8072639 0.1101296  7.330127 1.485102e-05
#4           Metinic (Intercept) 1.6770953 0.3141998  5.337671 1.771248e-04
#5           Monomoy (Intercept) 1.3370583 0.3008719  4.443945 6.740113e-03
#6      Outer green  (Intercept) 1.3739975 0.2637232  5.209998 2.899567e-04
#7       Petit Manan (Intercept) 0.7291428 0.2138888  3.408980 3.885520e-03
#8            Pond I (Intercept) 0.9114894 0.2538450  3.590732 2.952280e-03
#9            Seal I (Intercept) 0.7564574 0.1127727  6.707804 3.340944e-05
#10         Stratton (Intercept) 1.0008011 0.2157381  4.638963 7.178557e-04
#11 White and Seavey (Intercept) 1.3550055 0.1815945  7.461709 4.753354e-06

# the below groups both species, exaggerating trends where common terns have higher diet props than arctic
# combined with generally better productivity for common terns, so don't use this

colonydat <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerDiet = lm(productivity ~ prop.herring, data = .))
colonycoef <- tidy(colonydat, fitHerDiet)
filter(colonycoef, p.value<0.05)

#Source: local data frame [14 x 6]
#Groups: colony [12]
#
#colony         term   estimate  std.error statistic      p.value
#<fctr>        <chr>      <dbl>      <dbl>     <dbl>        <dbl>
#  1  Eastern Egg Rock  (Intercept)  0.7692069 0.07754518  9.919468 2.762048e-11
#2           Jenny I  (Intercept)  1.5943009 0.39165686  4.070657 1.145995e-03
#3      Machias Seal  (Intercept)  0.3024373 0.09718273  3.112048 5.079933e-03
#4      Machias Seal prop.herring  0.5515016 0.19556477  2.820046 9.970689e-03
#5    Matinicus Rock  (Intercept)  0.7201252 0.06537636 11.015071 7.043355e-12
#6           Metinic  (Intercept)  1.4497208 0.20337784  7.128214 1.152015e-07
#7           Metinic prop.herring -1.1746318 0.45797032 -2.564864 1.619727e-02
#8           Monomoy  (Intercept)  1.3370583 0.30087190  4.443945 6.740113e-03
#9      Outer green   (Intercept)  1.3739975 0.26372322  5.209998 2.899567e-04
#10      Petit Manan  (Intercept)  0.6034167 0.13243781  4.556227 7.181420e-05
#11           Pond I  (Intercept)  0.9114894 0.25384501  3.590732 2.952280e-03
#12           Seal I  (Intercept)  0.7986442 0.05457768 14.633165 6.350738e-15
#13         Stratton  (Intercept)  1.0008011 0.21573813  4.638963 7.178557e-04
#14 White and Seavey  (Intercept)  1.3550055 0.18159453  7.461709 4.753354e-06

#---------------------------------------------------------------------------------
# more plotting: is tern diet related to herring population?


terndietherring <- ggplot(ternpopproddiet_herring, aes(x=HerringAge1Rec, y=prop.herring, colour=species)) +
  geom_point() +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Age 1 Herring and herring in tern diet GOM")

terndietherring + facet_wrap(~colony) #+ theme(text=element_text(size=16))
ggsave("terndietherringrec.png")

terndietherring2 <- ggplot(ternpopproddiet_herring, aes(x=HerringTotalB, y=prop.herring, colour=species)) +
  geom_point() +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Herring total B and herring in tern diet GOM")

terndietherring2 + facet_wrap(~colony) #+ theme(text=element_text(size=16))
ggsave("terndietherringtotB.png")

terndietherring3 <- ggplot(ternpopproddiet_herring, aes(x=HerringSSB, y=prop.herring, colour=species)) +
  geom_point() +
  #stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Herring SSB and herring in tern diet GOM")

terndietherring3 + facet_wrap(~colony) #+ theme(text=element_text(size=16))
ggsave("terndietherringSSB.png")

terndietherring4 <- ggplot(ternpopproddiet_herring, aes(x=Nmean, y=prop.herring, colour=species)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("ME NH trawl herring N and herring in tern diet GOM")

terndietherring4 + facet_wrap(~colony) + theme(text=element_text(size=16))

terndietherring5 <- ggplot(ternpopproddiet_herring, aes(x=Wtmean, y=prop.herring, colour=species)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("ME NH trawl herring Wt and herring in tern diet GOM")

terndietherring5 + facet_wrap(~colony) + theme(text=element_text(size=16))

#-----------------------------------------------------------------------------------
# are any of the relationships between herring and tern diet proportion significant?

colonydietcor <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$HerringAge1Rec,.$prop.herring, alternative="greater", method="spearman", na.action=na.rm))) 
filter(colonydietcor, p.value<0.05)

#Source: local data frame [4 x 6]
#Groups: colony [4]
#
#colony  estimate statistic      p.value
#<fctr>     <dbl>     <dbl>        <dbl>
#1 Matinicus Rock 0.5012716 142.63631 0.0484339003
#2   Outer green  0.6713287  94.00000 0.0102209392
#3         Seal I 0.8521338  42.28973 0.0002159835
#4       Stratton 0.5314685 134.00000 0.0396514697
# ... with 2 more variables: method <fctr>, alternative <fctr>

colonydietcor <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$HerringTotalB,.$prop.herring, alternative="greater", method="spearman", na.action=na.rm))) 
filter(colonydietcor, p.value<0.05)

#Source: local data frame [5 x 6]
#Groups: colony [5]
#
#colony  estimate statistic     p.value                          method
#<fctr>     <dbl>     <dbl>       <dbl>                          <fctr>
#1        Jenny I 0.4750000  294.0000 0.037935041 Spearman's rank correlation rho
#2   Machias Seal 0.8454545   34.0000 0.001038844 Spearman's rank correlation rho
#3 Matinicus Rock 0.5339633  133.2865 0.036873435 Spearman's rank correlation rho
#4   Outer green  0.6013986  114.0000 0.021403434 Spearman's rank correlation rho
#5         Seal I 0.5563518  126.8834 0.030150013 Spearman's rank correlation rho
## ... with 1 more variables: alternative <fctr>

colonydietcor <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$HerringSSB,.$prop.herring, alternative="greater", method="spearman", na.action=na.rm))) 
filter(colonydietcor, p.value<0.05)

#Source: local data frame [4 x 6]
#Groups: colony [4]
#
#colony  estimate statistic     p.value
#<fctr>     <dbl>     <dbl>       <dbl>
# 1 Eastern Egg Rock 0.4705884  603.5292 0.021003985
#2          Jenny I 0.6178571  214.0000 0.008155134
#3     Machias Seal 0.6454545   78.0000 0.018500523
#4           Ship I 0.9642857    2.0000 0.001388889
# ... with 2 more variables: method <fctr>, alternative <fctr>

colonydietcor <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$Nmean,.$prop.herring, alternative="greater", method="spearman", na.action=na.rm))) 
filter(colonydietcor, p.value<0.05)

#Source: local data frame [0 x 6]
#Groups: colony [0]


colonydietcor <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(tidy(cor.test(.$Wtmean,.$prop.herring, alternative="greater", method="spearman", na.action=na.rm)))
filter(colonydietcor, p.value<0.05)

#Source: local data frame [0 x 6]
#Groups: colony [0]


colonydat1 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerRec = lm(prop.herring ~ HerringAge1Rec, data = .))
colony1coef <- tidy(colonydat1, fitHerRec)
filter(colony1coef, p.value<0.05)

#Source: local data frame [10 x 6]
#Groups: colony [10]
#
#colony           term     estimate    std.error statistic      p.value
#(fctr)          (chr)        (dbl)        (dbl)     (dbl)        (dbl)
#1  Eastern Egg Rock    (Intercept) 1.828337e-01 5.979646e-02  3.057600 7.122529e-03
#2           Jenny I    (Intercept) 3.193621e-01 6.492480e-02  4.918953 2.804797e-04
#3    Matinicus Rock HerringAge1Rec 4.428688e-09 1.476315e-09  2.999825 1.334764e-02
#4           Metinic    (Intercept) 5.012099e-01 9.782816e-02  5.123371 2.518470e-04
#5      Outer green     (Intercept) 1.829509e-01 6.808271e-02  2.687186 2.280960e-02
#6       Petit Manan    (Intercept) 4.694850e-01 9.686972e-02  4.846560 2.589248e-04
#7            Pond I    (Intercept) 1.721888e-01 6.248205e-02  2.755813 1.546740e-02
#8            Seal I HerringAge1Rec 1.186358e-08 1.283970e-09  9.239759 3.264872e-06
#9            Ship I HerringAge1Rec 1.322569e-08 4.975440e-09  2.658196 4.498128e-02
#10         Stratton HerringAge1Rec 6.239316e-09 2.520065e-09  2.475855 3.277463e-02

colonydat1 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerRec = lm(prop.herring ~ HerringAge1Rec, data = .))
colony1coef <- tidy(colonydat1, fitHerRec)
filter(colony1coef, p.value<0.05)

#Source: local data frame [9 x 6]
#Groups: colony [9]
#
#colony           term     estimate    std.error statistic      p.value
#(fctr)          (chr)        (dbl)        (dbl)     (dbl)        (dbl)
#1 Eastern Egg Rock    (Intercept) 1.097248e-01 4.119321e-02  2.663664 1.231282e-02
#2          Jenny I    (Intercept) 3.193621e-01 6.492480e-02  4.918953 2.804797e-04
#3          Metinic    (Intercept) 4.867950e-01 6.295199e-02  7.732798 2.572001e-08
#4     Outer green     (Intercept) 1.829509e-01 6.808271e-02  2.687186 2.280960e-02
#5      Petit Manan    (Intercept) 4.159017e-01 6.544624e-02  6.354860 5.193903e-07
#6           Pond I    (Intercept) 1.721888e-01 6.248205e-02  2.755813 1.546740e-02
#7           Seal I HerringAge1Rec 8.388527e-09 1.225501e-09  6.844976 2.361422e-07
#8           Ship I HerringAge1Rec 1.322569e-08 4.975440e-09  2.658196 4.498128e-02
#9         Stratton HerringAge1Rec 6.239316e-09 2.520065e-09  2.475855 3.277463e-02

colonydat2 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerB = lm(prop.herring ~ HerringTotalB, data = .))
colony2coef <- tidy(colonydat2, fitHerB)
filter(colony2coef, p.value<0.05)

#Source: local data frame [9 x 6]
#Groups: colony [8]
#
#colony          term      estimate    std.error statistic      p.value
#(fctr)         (chr)         (dbl)        (dbl)     (dbl)        (dbl)
#1          Jenny I HerringTotalB  2.252319e-07 8.780522e-08  2.565131 0.0235096001
#2     Machias Seal   (Intercept) -1.049394e+00 3.693783e-01 -2.840974 0.0193723502
#3     Machias Seal HerringTotalB  1.016394e-06 2.405940e-07  4.224522 0.0022245389
#4          Metinic   (Intercept)  8.550026e-01 2.234442e-01  3.826470 0.0024109747
#5     Outer green  HerringTotalB  2.284862e-07 8.727454e-08  2.618017 0.0256849725
#6      Petit Manan   (Intercept)  9.426054e-01 2.207835e-01  4.269364 0.0007783701
#7           Seal I HerringTotalB  2.821979e-07 1.085193e-07  2.600440 0.0264714902
#8           Ship I HerringTotalB  5.459466e-07 1.867170e-07  2.923926 0.0328635733
#9 White and Seavey HerringTotalB  2.458012e-07 1.021727e-07  2.405742 0.0331694282

colonydat2 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerB = lm(prop.herring ~ HerringTotalB, data = .))
colony2coef <- tidy(colonydat2, fitHerB)
filter(colony2coef, p.value<0.05)

#Source: local data frame [9 x 6]
#Groups: colony [8]
#
#colony          term      estimate    std.error statistic      p.value
#(fctr)         (chr)         (dbl)        (dbl)     (dbl)        (dbl)
#1          Jenny I HerringTotalB  2.252319e-07 8.780522e-08  2.565131 2.350960e-02
#2     Machias Seal HerringTotalB  5.835127e-07 2.104765e-07  2.772342 1.141503e-02
#3          Metinic   (Intercept)  7.952334e-01 1.508982e-01  5.270000 1.475415e-05
#4          Metinic HerringTotalB -2.650043e-07 9.728860e-08 -2.723899 1.117338e-02
#5     Outer green  HerringTotalB  2.284862e-07 8.727454e-08  2.618017 2.568497e-02
#6      Petit Manan   (Intercept)  7.947924e-01 1.591657e-01  4.993490 2.372882e-05
#7           Seal I HerringTotalB  1.869721e-07 6.442728e-08  2.902064 7.294513e-03
#8           Ship I HerringTotalB  5.459466e-07 1.867170e-07  2.923926 3.286357e-02
#9 White and Seavey HerringTotalB  2.458012e-07 1.021727e-07  2.405742 3.316943e-02

colonydat3 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerSSB = lm(prop.herring ~ HerringSSB, data = .))
colony3coef <- tidy(colonydat3, fitHerSSB)
filter(colony3coef, p.value<0.05)

#Source: local data frame [6 x 6]
#Groups: colony [5]
#
#colony        term      estimate    std.error statistic      p.value
#(fctr)       (chr)         (dbl)        (dbl)     (dbl)        (dbl)
#1 Eastern Egg Rock  HerringSSB  3.471541e-07 1.590096e-07  2.183227 4.332927e-02
#2          Jenny I  HerringSSB  5.943246e-07 1.558033e-07  3.814583 2.147045e-03
#3          Metinic (Intercept)  5.252152e-01 1.905080e-01  2.756920 1.737806e-02
#4     Outer green   HerringSSB  5.561164e-07 1.896708e-07  2.932008 1.498916e-02
#5      Petit Manan (Intercept)  9.107092e-01 1.650417e-01  5.518055 7.574612e-05
#6      Petit Manan  HerringSSB -6.860021e-07 2.784406e-07 -2.463728 2.731241e-02

colonydat3 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerSSB = lm(prop.herring ~ HerringSSB, data = .))
colony3coef <- tidy(colonydat3, fitHerSSB)
filter(colony3coef, p.value<0.05)

#Source: local data frame [5 x 6]
#Groups: colony [4]
#
#colony        term      estimate    std.error statistic      p.value
#(fctr)       (chr)         (dbl)        (dbl)     (dbl)        (dbl)
#1      Jenny I  HerringSSB  5.943246e-07 1.558033e-07  3.814583 2.147045e-03
#2      Metinic (Intercept)  4.773676e-01 1.245674e-01  3.832202 6.887200e-04
#3 Outer green   HerringSSB  5.561164e-07 1.896708e-07  2.932008 1.498916e-02
#4  Petit Manan (Intercept)  7.863331e-01 1.212669e-01  6.484317 3.631744e-07
#5  Petit Manan  HerringSSB -5.406096e-07 2.045885e-07 -2.642424 1.295363e-02

#---------------------------------------------------------------------------------------------
# more plotting, is the herring population related to tern productivity? warning! causation 
# may not be evident given the above


ternprodherring <- ggplot(ternpopproddiet_herring, aes(x=HerringAge1Rec, y=productivity, colour=species)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Age 1 Herring and tern productivity GOM")

ternprodherring + facet_wrap(~colony) + theme(text=element_text(size=16))

ternprodherring2 <- ggplot(ternpopproddiet_herring, aes(x=HerringTotalB, y=productivity, colour=species)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Herring total B and tern productivity GOM")

ternprodherring2 + facet_wrap(~colony) + theme(text=element_text(size=16))

ternprodherring3 <- ggplot(ternpopproddiet_herring, aes(x=HerringSSB, y=productivity, colour=species)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  #stat_smooth(method = "loess", formula = y ~ x, size = 1) +
  ggtitle("Est Herring SSB and tern productivity GOM")

ternprodherring3 + facet_wrap(~colony) + theme(text=element_text(size=16))

#---------------------------------------------------------------------------------------------
# are any of the relationships between herring and tern productivity significant?

# LIKELY SPURIOUS given that only Stratton had a sig. positive relationship
# between herring in diet and productivity!

colonydat1 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerRec = lm(productivity ~ HerringAge1Rec, data = .))
colony1coef <- tidy(colonydat1, fitHerRec)
filter(colony1coef, p.value<0.05)

#Source: local data frame [11 x 6]
#Groups: colony [11]
#
#colony        term  estimate std.error statistic      p.value
#<fctr>       <chr>     <dbl>     <dbl>     <dbl>        <dbl>
#  1  Eastern Egg Rock (Intercept) 0.9311534 0.1285918  7.241157 1.379104e-06
#2           Jenny I (Intercept) 1.1756832 0.2061274  5.703671 9.770739e-06
#3    Matinicus Rock (Intercept) 0.7850854 0.1454009  5.399455 1.602247e-04
#4           Metinic (Intercept) 1.0190105 0.2918110  3.492022 3.974390e-03
#5           Monomoy (Intercept) 1.4108492 0.1471015  9.590992 8.631382e-08
#6      Outer green  (Intercept) 1.3427314 0.2133379  6.293919 8.978800e-05
#7       Petit Manan (Intercept) 0.9195563 0.1122817  8.189725 8.542256e-09
#8            Pond I (Intercept) 0.9665425 0.2743972  3.522421 3.078593e-03
#9            Seal I (Intercept) 0.8719750 0.1070947  8.142092 1.289047e-07
#10         Stratton (Intercept) 0.9273018 0.2486285  3.729669 1.824412e-03
#11 White and Seavey (Intercept) 1.3574677 0.1787971  7.592224 1.085852e-06

colonydat1 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerRec = lm(productivity ~ HerringAge1Rec, data = .))
colony1coef <- tidy(colonydat1, fitHerRec)
filter(colony1coef, p.value<0.05)

#Source: local data frame [12 x 6]
#Groups: colony [12]
#
#colony        term  estimate  std.error statistic      p.value
#<fctr>       <chr>     <dbl>      <dbl>     <dbl>        <dbl>
#  1  Eastern Egg Rock (Intercept) 0.8443736 0.09637870  8.760998 3.984155e-10
#2           Jenny I (Intercept) 1.1756832 0.20612745  5.703671 9.770739e-06
#3      Machias Seal (Intercept) 0.2922311 0.10200413  2.864895 6.759861e-03
#4    Matinicus Rock (Intercept) 0.7199000 0.08965519  8.029653 4.567455e-09
#5           Metinic (Intercept) 0.8980569 0.19139446  4.692178 6.433060e-05
#6           Monomoy (Intercept) 1.4108492 0.14710149  9.590992 8.631382e-08
#7      Outer green  (Intercept) 1.3427314 0.21333790  6.293919 8.978800e-05
#8       Petit Manan (Intercept) 0.8283966 0.09295657  8.911652 2.060891e-11
#9            Pond I (Intercept) 0.8764488 0.24040167  3.645768 1.510530e-03
#10           Seal I (Intercept) 0.8536076 0.06661479 12.814086 5.602823e-15
#11         Stratton (Intercept) 0.6754455 0.19082029  3.539694 1.328851e-03
#12 White and Seavey (Intercept) 1.1976760 0.15959322  7.504554 2.850337e-08


colonydat2 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerB = lm(productivity ~ HerringTotalB, data = .))
colony2coef <- tidy(colonydat2, fitHerB)
filter(colony2coef, p.value<0.05)

#Source: local data frame [7 x 6]
#Groups: colony [7]
#
#colony          term     estimate    std.error statistic
#<fctr>         <chr>        <dbl>        <dbl>     <dbl>
#  1 Eastern Egg Rock   (Intercept) 7.484166e-01 3.193034e-01  2.343904
#2          Monomoy   (Intercept) 8.683906e-01 3.740497e-01  2.321592
#3     Outer green    (Intercept) 1.055014e+00 4.643734e-01  2.271910
#4      Petit Manan   (Intercept) 1.115116e+00 1.882850e-01  5.922489
#5           Seal I   (Intercept) 1.119680e+00 2.758297e-01  4.059316
#6           Ship I HerringTotalB 9.222767e-07 3.589396e-07  2.569448
#7 White and Seavey   (Intercept) 1.188591e+00 4.568467e-01  2.601728

colonydat2 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerB = lm(productivity ~ HerringTotalB, data = .))
colony2coef <- tidy(colonydat2, fitHerB)
filter(colony2coef, p.value<0.05)

#Source: local data frame [8 x 6]
#Groups: colony [8]
#
#colony          term     estimate    std.error statistic
#<fctr>         <chr>        <dbl>        <dbl>     <dbl>
#  1 Eastern Egg Rock   (Intercept) 8.328106e-01 2.408276e-01  3.458119
#2   Matinicus Rock   (Intercept) 4.866580e-01 2.087795e-01  2.330967
#3          Metinic   (Intercept) 1.058238e+00 4.547483e-01  2.327085
#4          Monomoy   (Intercept) 8.683906e-01 3.740497e-01  2.321592
#5     Outer green    (Intercept) 1.055014e+00 4.643734e-01  2.271910
#6      Petit Manan   (Intercept) 1.126677e+00 1.716550e-01  6.563616
#7           Seal I   (Intercept) 1.055051e+00 1.674511e-01  6.300651
#8           Ship I HerringTotalB 9.222767e-07 3.589396e-07  2.569448

colonydat3 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  filter(species == "Common") %>%
  do(fitHerSSB = lm(productivity ~ HerringSSB, data = .))
colony3coef <- tidy(colonydat3, fitHerSSB)
filter(colony3coef, p.value<0.05)

#Source: local data frame [13 x 6]
#Groups: colony [10]
#
#colony        term      estimate    std.error statistic
#<fctr>       <chr>         <dbl>        <dbl>     <dbl>
#  1  Eastern Egg Rock (Intercept)  4.893150e-01 2.062242e-01  2.372734
#2  Eastern Egg Rock  HerringSSB  7.624986e-07 3.395357e-07  2.245710
#3           Jenny I (Intercept)  1.332724e+00 4.268335e-01  3.122352
#4    Matinicus Rock  HerringSSB  1.001856e-06 3.617730e-07  2.769294
#5           Monomoy (Intercept)  6.236351e-01 2.456170e-01  2.539055
#6           Monomoy  HerringSSB  1.157736e-06 4.203991e-07  2.753896
#7      Outer green  (Intercept)  1.628220e+00 3.702485e-01  4.397643
#8       Petit Manan (Intercept)  9.911370e-01 1.727882e-01  5.736138
#9            Pond I (Intercept)  1.428785e+00 4.777275e-01  2.990796
#10           Seal I (Intercept)  7.820303e-01 2.149920e-01  3.637485
#11           Ship I (Intercept) -8.114037e-01 3.448419e-01 -2.352973
#12           Ship I  HerringSSB  2.272879e-06 5.642173e-07  4.028375
#13 White and Seavey (Intercept)  1.395731e+00 3.271488e-01  4.266350

colonydat3 <- ternpopproddiet_herring %>%
  group_by(colony) %>%
  do(fitHerSSB = lm(productivity ~ HerringSSB, data = .))
colony3coef <- tidy(colonydat3, fitHerSSB)
filter(colony3coef, p.value<0.05)

#Source: local data frame [14 x 6]
#Groups: colony [11]
#
#colony        term      estimate    std.error statistic
#<fctr>       <chr>         <dbl>        <dbl>     <dbl>
#  1  Eastern Egg Rock (Intercept)  5.177596e-01 1.744492e-01  2.967968
#2           Jenny I (Intercept)  1.332724e+00 4.268335e-01  3.122352
#3      Machias Seal  HerringSSB  6.479592e-07 2.891545e-07  2.240876
#4    Matinicus Rock (Intercept)  3.452314e-01 1.417609e-01  2.435308
#5    Matinicus Rock  HerringSSB  7.124122e-07 2.419406e-07  2.944574
#6           Monomoy (Intercept)  6.236351e-01 2.456170e-01  2.539055
#7           Monomoy  HerringSSB  1.157736e-06 4.203991e-07  2.753896
#8      Outer green  (Intercept)  1.628220e+00 3.702485e-01  4.397643
#9       Petit Manan (Intercept)  9.703511e-01 1.539284e-01  6.303913
#10           Pond I (Intercept)  9.183981e-01 4.015240e-01  2.287281
#11           Seal I (Intercept)  7.424305e-01 1.306713e-01  5.681665
#12           Ship I (Intercept) -8.114037e-01 3.448419e-01 -2.352973
#13           Ship I  HerringSSB  2.272879e-06 5.642173e-07  4.028375
#14 White and Seavey (Intercept)  8.339439e-01 2.818465e-01  2.958859

##################################################################################################
# total tern pop dataset to try to derive overall stock recruit relationship 
##################################################################################################

terntotalpop <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

write.csv(terntotalpop, "terntotalpop_corrected.csv")
terntotalpop_old <- read.csv("terntotalpop.csv")

# Monomoy Island is a third of the population but eats mainly sandlance and is south of GOM proper--leave out?

terntotpop_noMon <- ternpopproddiet %>%
  group_by(yr, species) %>%
  filter(colony != "Monomoy") %>%
  filter(yr<2016) %>%
  select(breedpairs, nfledged) %>%
  filter(!is.na(nfledged)) %>%
  summarise(totpairs = sum(breedpairs, na.rm=T), totfledge = sum(nfledged, na.rm=T))

write.csv(terntotpop_noMon, "terntotpop_noMon_corrected.csv")
terntotpop_noMon_old <- read.csv("terntotpop_noMon.csv")

ternpops2 <- ggplot(terntotalpop, aes(x=yr, y=totpairs, colour=species))+
  geom_point() +
  ggtitle("Tern total breeding pairs and estimated fledglings GOM")

ternpops2 +   geom_line(aes(x=yr, y=totfledge, colour=species))

ternpops2 + geom_line(aes(x=yr, y=totpairs, colour=species), data=terntotalpop_old)

ternpops2_noMon <- ggplot(terntotpop_noMon, aes(x=yr, y=totpairs, colour=species))+
  geom_point() +
  ggtitle("Tern total breeding pairs and estimated fledglings GOM no Monomoy")

ternpops2_noMon +   geom_line(aes(x=yr, y=totfledge, colour=species))

ternpops2_noMon + geom_line(aes(x=yr, y=totpairs, colour=species), data=terntotpop_noMon_old)

terntotalpop19982015 <- terntotalpop %>%
  filter(yr>1997)

ternpoptrend <- ggplot(terntotalpop19982015, aes(x=yr, y=totpairs, colour=species))

ternpoptrend +   geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) 

# plot with actual data and simulated data from pop model
# WARNING only works with PredN from MSE_TernLink.R in memory
simTern <- data.frame(yr=seq(1998,2015,1), PredN=PredN[6:23])

ternpoptrend +   geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  geom_point(data=simTern, aes(x=yr, y=PredN, colour="Simulated"))


commontotalpop19982015 <- terntotalpop19982015 %>%
  filter(species=="Common")

commonpoptrend <- lm(totpairs ~ yr, data=commontotalpop19982015)

commonpoptrend

#Call:
#  lm(formula = totpairs ~ yr, data = commontotalpop19982015)
#
#Coefficients:
#  (Intercept)           yr  
#-692682.1        353.2  


with(commontotalpop19982015, mean(totpairs))
#[1] 16047.28

commonpoptrend$coefficients[2]/with(commontotalpop19982015, mean(totpairs))

#0.02201101  
  
terntotalpop19982015_noMon <- terntotpop_noMon %>%
  filter(yr>1997)

ternpoptrend_noMon <- ggplot(terntotalpop19982015_noMon, aes(x=yr, y=totpairs, colour=species))

ternpoptrend_noMon +   geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1)

commontotalpop19982015_noMon <- terntotalpop19982015_noMon %>%
  filter(species=="Common")

commonpoptrend_noMon <- lm(totpairs ~ yr, data=commontotalpop19982015_noMon)

commonpoptrend_noMon

#Call:
#  lm(formula = totpairs ~ yr, data = commontotalpop19982015_noMon)
#
#Coefficients:
#  (Intercept)           yr  
#-502026.5        254.6  

with(commontotalpop19982015_noMon, mean(totpairs))
#[1] 8815.167
commonpoptrend_noMon$coefficients[2]/with(commontotalpop19982015_noMon, mean(totpairs))

#0.02888129  

# Some recruitment functions and assumptions

# Rec(=totfledge*.1) = (4*steep*Rmax*totpairs)/(totpairsmax*(1-steep)+totpairs(5*steep-1))
# pop dat from BNA, New England 45K pairs in 1930s = totpairsmax
# surv rate fledge to age 4 (nesting @3) 0.06-0.1 or 0.07-.13
# mean prod 1-2 fledglings per nesting pair so call Rmax 45K*2*.1 = 9K
# initial est of steepness 0.2

BHsteep <- function(Sdat, steep, Rmax, Smax){
  4*steep*Rmax*Sdat/
    (Smax*(1-steep)+Sdat*(5*steep-1))
}

#maybe this can have depensation? or better to fit separate functions for common and arctic?
BHsteepDep <- function(Sdat, steep, Rmax, Smax, Dep){
  4*steep*Rmax*(Sdat^Dep)/
    (Smax*(1-steep)+(Sdat^Dep)*(5*steep-1))
}

ternrec_noMon <- ggplot(terntotpop_noMon, aes(x=totpairs, y=totfledge*.1, colour=species)) +
  geom_point() +
  scale_y_continuous(limits = c(0,5000)) +
  scale_x_continuous(limits = c(0,45000)) +
  ggtitle("Tern stock and recruitment assuming 10% fledgling->adult survival")

# playing with parameter values without monomoy island colony
ternrec_noMon + 
  stat_function(fun=BHsteep, args=list(steep=0.201, Rmax=5000, Smax=45000), col="black") +
  stat_function(fun=BHsteep, args=list(steep=0.3, Rmax=5000, Smax=45000), col="blue") +
  stat_function(fun=BHsteep, args=list(steep=0.4, Rmax=5000, Smax=45000), col="red") +
  stat_function(fun=BHsteep, args=list(steep=0.5, Rmax=5000, Smax=45000), col="green")

ternrec_noMon + 
  stat_function(fun=BHsteep, args=list(steep=0.201, Rmax=5000, Smax=45000), col="black") +
  stat_function(fun=BHsteep, args=list(steep=0.201, Rmax=2500, Smax=45000), col="blue") +
  stat_function(fun=BHsteep, args=list(steep=0.3, Rmax=2500, Smax=45000), col="red") +
  stat_function(fun=BHsteep, args=list(steep=0.4, Rmax=2500, Smax=45000), col="green")

ternrec_noMon + 
  stat_function(fun=BHsteep, args=list(steep=0.201, Rmax=5000, Smax=45000), col="black") +
  stat_function(fun=BHsteep, args=list(steep=0.201, Rmax=3000, Smax=45000), col="blue") +
  stat_function(fun=BHsteepDep, args=list(steep=0.401, Rmax=5000, Smax=45000, Dep=1.2), col="red") +
  stat_function(fun=BHsteepDep, args=list(steep=0.401, Rmax=3000, Smax=45000,  Dep=1.2), col="green")

steeplist<-seq(0.201, 0.991, 0.05)
Rmaxlist<-seq(2000, 10000, 500)
Smaxlist<-seq(10000, 50000, 5000)
Deplist<-seq(0.5,1.5,0.1)

#try simple linear and loess fits compared to guess at BHsteep pars

ternrec_noMon + 
  stat_function(fun=BHsteepDep, args=list(steep=0.201, Rmax=5000, Smax=45000, Dep=1), col="black") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1)

ternrec_noMon + 
  stat_function(fun=BHsteepDep, args=list(steep=0.201, Rmax=5000, Smax=45000, Dep=1), col="black") +
  stat_smooth(method = "loess", formula = y ~ x, size = 1)

#maybe recruitment fitting for common terns would be better with monomoy, more x range

ternrec <- ggplot(terntotalpop, aes(x=totpairs, y=totfledge*.1, colour=species)) +
  geom_point() +
  scale_y_continuous(limits = c(0,5000)) +
  scale_x_continuous(limits = c(0,45000)) +
  ggtitle("Tern stock and recruitment assuming 10% fledgling->adult survival")

ternrec + 
  stat_function(fun=BHsteep, args=list(steep=0.27, Rmax=4500, Smax=45000), col="black") +
  stat_smooth(method = "lm", formula = y ~ x, size = 1)

ternrec + 
  stat_function(fun=BHsteepDep, args=list(steep=0.25, Rmax=4500, Smax=45000, Dep=1), col="black") +
  stat_smooth(method = "loess", formula = y ~ x, size = 1)



  stat_function(fun=BHsteepDep, args=list(steep=steeplist, Rmax=5000, Smax=45000, Dep=1)) 
  

#BHsteep(terntotpop_noMon$totpairs, steep=0.2, Rmax=9000, Smax=45000)

with(terntotpop_noMon, plot(totpairs, totfledge*.1))#, xlim=c(0,45000), ylim=c(0,10000)))
lines(seq(0, 45000, by=100), BHsteep(seq(0, 45000, by=100), steep=0.201, Rmax=5000, Smax=45000))

#with(terntotpop_noMon, points(totpairs, BHsteep(totpairs, steep=0.25, Rmax=2500, Smax=30000), pch="x"))

terntotpop_com <- terntotalpop %>%
  filter(species=="Common") %>%

#these arent working
ternBH <- nls(totfledge*.1 ~ BHsteep(totpairs, steep, Rmax, Smax), terntotpop_com,
              start=list(steep=0.25, Rmax=4500, Smax=45000))

ternBH <- nls(totfledge*.1 ~ BHsteep(totpairs, steep, Rmax, Smax), commontotalpop19982015,
              start=list(steep=0.25, Rmax=4500, Smax=45000))

##### Tern SR fit to full dataset (possible missing values <1998)
#try a different parameterization with only 2 pars
BHfun <- function(x,a,b){
  (x/a)/(1+(1/a*b)*x)
}

alldatBH <- terntotalpop %>%
  group_by(species) %>%
  filter(!is.na(totfledge)) %>%
  do(mod2 = nls((1/(totfledge*0.1)) ~ betastar + alphastar*(1/totpairs), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

spcoefBH <- tidy(alldatBH, mod2)

spcoefBH
#Source: local data frame [4 x 6]
#Groups: species [2]
#
#species      term      estimate    std.error  statistic      p.value
#<fctr>     <chr>         <dbl>        <dbl>      <dbl>        <dbl>
#1  Arctic alphastar  7.9851619025 1.9727515128  4.0477282 7.551535e-04
#2  Arctic  betastar  0.0024994673 0.0007634703  3.2738240 4.218001e-03
#3  Common alphastar  9.7860199801 0.6094627512 16.0567975 5.727951e-16
#4  Common  betastar -0.0001570145 0.0004843356 -0.3241853 7.481244e-01

spcoefBH.1 <- spcoefBH %>%
  group_by(species) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

spcoefBH.1
#Source: local data frame [2 x 3]
#Groups: species [2]
#
#species alphastar      betastar
#*  <fctr>     <dbl>         <dbl>
# 1  Arctic  7.985162  0.0024994673
# 2  Common  9.786020 -0.0001570145

pBHspdat <- left_join(terntotalpop, spcoefBH.1) %>%
  mutate(BHfit = BHfun(totpairs, alphastar, betastar))


ternrec + 
  #geom_smooth(aes(x=totpairs, y=BHfit, col=species), data=pBHspdat, fullrange=T)+
  stat_function(fun=BHfun, args=list(a= 9.786020,b=-0.0001570145), col="black")+
  stat_function(fun=BHfun, args=list(a= 9.786020,b=0), col="green")+
  stat_function(fun=BHfun, args=list(a= 7.985162,b=0.0024994673), col="red")
  
##### Tern SR fit to 1998-2015 dataset, all islands reporting

alldatBH_9815 <- terntotalpop19982015 %>%
  group_by(species) %>%
  filter(!is.na(totfledge)) %>%
  do(mod2 = nls((1/(totfledge*0.1)) ~ betastar + alphastar*(1/totpairs), 
                start=list(alphastar = 1.0 , betastar = 1.0), data=.))

spcoefBH_9815 <- tidy(alldatBH_9815, mod2)

#Source: local data frame [4 x 6]
#Groups: species [2]
#
#species      term     estimate    std.error statistic     p.value
#<fctr>     <chr>        <dbl>        <dbl>     <dbl>       <dbl>
# 1  Arctic alphastar 1.092357e+01 4.1193347934  2.651780 0.017405871
# 2  Arctic  betastar 1.470393e-03 0.0012620068  1.165123 0.261046526
# 3  Common alphastar 5.445551e+00 1.8628439439  2.923246 0.009948929
# 4  Common  betastar 2.266453e-04 0.0001225097  1.850019 0.082867156

spcoefBH_9815.1 <- spcoefBH_9815 %>%
  group_by(species) %>%
  select(term, estimate) %>%
  spread(key=term, value=estimate) 

#Source: local data frame [2 x 3]
#Groups: species [2]
#
#species alphastar     betastar
#*  <fctr>     <dbl>        <dbl>
# 1  Arctic 10.923570 0.0014703930
# 2  Common  5.445551 0.0002266453

spcoefBH_9815.2 <- spcoefBH_9815.1 %>%
  group_by(species) %>%
  mutate(alpha = 1/alphastar) %>%
  mutate(beta = 1/alphastar*betastar) %>%
  mutate(steep = 0.2/(1-(1-alphastar)*0.8)) %>%
  mutate(Smax = (1-alphastar)/betastar) %>%
  mutate(Rmax = 1/betastar) %>%
  #mutate(Rmax2 = alpha/beta) # testing; same as above
  mutate(steep2 = 0.2*((alphastar/betastar)+Smax)/((alphastar/betastar) + 0.2*(alphastar/betastar)))



pBHspdat2 <- left_join(terntotalpop, spcoefBH_9815.1) %>%
  mutate(BHfit = BHfun(totpairs, alphastar, betastar))

ternrec + 
  #geom_smooth(aes(x=totpairs, y=BHfit, col=species), data=pBHspdat2, fullrange=T)
  stat_function(fun=BHfun, args=list(a= 5.445551,b=0.0002266453), col="green")+
  stat_function(fun=BHfun, args=list(a= 10.923570,b=0.0014703930), col="red")+
  stat_function(fun=BHfun, args=list(a= 10.923570,b=0), col="black")

# compare all fits to full vs complete shorter dataset
ternrec + 
  stat_function(fun=BHfun, args=list(a= 9.786020,b=-0.0001570145), col="cyan3")+
  stat_function(fun=BHfun, args=list(a= 9.786020,b=0), col="cyan3", lty=3)+
  stat_function(fun=BHfun, args=list(a= 7.985162,b=0.0024994673), col="red")+
  stat_function(fun=BHfun, args=list(a= 5.445551,b=0.0002266453), col="darkgreen")+
  stat_function(fun=BHfun, args=list(a= 10.923570,b=0.0014703930), col="darkred")+
  stat_function(fun=BHfun, args=list(a= 10.923570,b=0), col="darkred", lty=3)

#compare fit pars to converted steepness parameterization
BH<-function(x,alpha,beta){
  (alpha*x)/(1+beta*x)
}

ternrec +
  stat_function(fun=BHfun, args=list(a= 5.445551,b=0.0002266453), col="green", lwd=2)+
  #stat_function(fun=BH, args=list(alpha=0.18363616, beta=4.162027e-05), col="black") +
  stat_function(fun=BHfun, args=list(a= 9.786020,b=-0.0001570145), col="cyan3")+
  stat_function(fun=BHfun, args=list(a= 9.786020,b=0), col="cyan3", lty=3)+
  stat_function(fun=BHsteep, args=list(steep=0.41, Rmax=2900, Smax=45000), col="black") +
  stat_function(fun=BHsteep, args=list(steep=0.26, Rmax=4500, Smax=45000), col="black")

######################################################################
# herring pop plots
######################################################################
with(herring, plot(yr, HerringTotalB, ylim=c(0,max(HerringAge1Rec/10, na.rm=T)), type="l"))
with(herring, lines(yr, HerringSSB, col="blue"))
with(herring, points(yr, HerringAge1Rec/10))

ternpop_herringpop <- left_join(terntotalpop, herring)
ternpopnomon_herringpop <- left_join(terntotpop_noMon, herring)
  
popsplot <- ggplot(ternpop_herringpop, aes(x=yr, y=totfledge, colour=species)) +
  geom_point() +
  geom_line(aes(x=yr, y=HerringTotalB/100), color="black") +
  geom_point(aes(x=yr, y=HerringAge1Rec/1000), color="black") +
  ggtitle("Tern fledgelings and assessed herring population")

popsplot + theme(text=element_text(size=16))

popsplot <- ggplot(ternpopnomon_herringpop, aes(x=yr, y=totfledge, colour=species)) +
  geom_point() +
  geom_line(aes(x=yr, y=HerringTotalB/100), color="black") +
  geom_point(aes(x=yr, y=HerringAge1Rec/1000), color="black") +
  ggtitle("Tern fledgelings (no Monomoy) and assessed herring population")

popsplot + theme(text=element_text(size=16))

totprodherr <- ggplot(ternpop_herringpop, aes(x=HerringTotalB, y=totfledge/totpairs, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  ggtitle("Tern production and assessed herring tot B")

totprodherr + stat_smooth(method = "lm", formula = y ~ x, size = 1) 

lmtot <- lm(totfledge/totpairs ~ HerringTotalB, species=="Common", data=ternpop_herringpop)

totprodnomonherr <- ggplot(ternpopnomon_herringpop, aes(x=HerringTotalB, y=totfledge/totpairs, colour=species)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  ggtitle("Tern production (no Monomoy) and assessed herring tot B")

totprodnomonherr + stat_smooth(method = "lm", formula = y ~ x, size = 1) 

HerrAlpha <- 1.09
HerringThresh <- 400000

recmult<-(HerrAlpha*(ternpopnomon_herringpop$HerringTotalB/HerringThresh))/((HerrAlpha-1)+(ternpopnomon_herringpop$HerringTotalB/HerringThresh))

totprodnomonherr + geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5)

lmtotnomon <- lm(totfledge/totpairs ~ HerringTotalB, species=="Common", data=ternpopnomon_herringpop)

totprodnomonherr + geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5) +
  geom_line(aes(x=HerringTotalB, y=9.439e-01 + 1.210e-07*HerringTotalB), color="black", lty=3)

totprodnomonherr + 
  scale_y_continuous(limits = c(0,1.8)) +
  scale_x_continuous(limits = c(0,2500000)) +
  geom_line(aes(x=HerringTotalB, y=recmult), lwd=1.5) +
  geom_line(aes(x=HerringTotalB, y=9.439e-01 + 1.210e-07*HerringTotalB), color="black", lty=3) +
  geom_line(aes(x=seq(0,2500000,length.out=51), y=(HerrAlpha*(seq(0,2500000,length.out=51)/HerringThresh))/((HerrAlpha-1)+(seq(0,2500000,length.out=51)/HerringThresh))))

comprodherr <- ternpop_herringpop %>%
  group_by(species) %>%
  filter(!is.na(totfledge)) %>%
  do(mod2 = nls((totfledge) ~ mult*((predalpha*HerringTotalB)/(predbeta+HerringTotalB)), 
                start=list(mult = 1.0), data=.))

comprodpar <- tidy(comprodherr, mod2)

comprodpar



##########################################################################################################
# Basic exploratory plotting below here before I had the full dataset--ignore
##########################################################################################################



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


