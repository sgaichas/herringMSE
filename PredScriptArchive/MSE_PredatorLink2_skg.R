rm(list=ls())
#Some functions and code to simulate predators depending on prey abundance

#######Some Functions#########
###Function to calculate unfished "SSBR", and abundance or biomass, or equilibrium characteristics given exploitation rate for delay-difference dynamics with Beverton-Holt SR.
#function(unfishedspr,desired unfished N or B, steepness,units of numbs or bio?, annual natural mortality, annual fishery exploitation rate, ford walford parameter intercept, ford walford slope param, do unfished or exploited equilib?)
unfishspr<-function(predunfishspr=NULL,targzero=NULL,steep=NULL,numbs_weightb=NULL,annualnatmort=NULL,annualexploit=NULL,FWalpha=NULL,FWrho=NULL,unfishflag=NULL) {
  #see delay-diff chapter of hilborn and walters and common pubs on steepness (Mangel et al 2010)
  predrzero<-1/predunfishspr*targzero  #unfished recruitment=f(unfishedSSBR,unfished N or B specified by user
  predalpha<-(4*steep*predrzero)/(5*steep-1) #alpha of BH SR #1/beta? alpha in Mangel et al 2010 is (targzero/predrzero)*((1-steep)/(4*steep)) #
  predbeta<-((targzero/predrzero)*((1-steep)/(4*steep)))/((5*steep-1)/(4*steep*predrzero)) #beta of BH SR #alpha/beta? beta ibid (5*steep-1)/(4*steep*predrzero)
  #predalpha<- (targzero/predrzero)*((1-steep)/(4*steep)) #direct from Mangel et al 2010
  #predbeta<-(5*steep-1)/(4*steep*predrzero) #direct from Mangel et al 2010
  if(numbs_weightb==1){ #solving in units of numbers or biomass
    kappa<-1-(1-annualnatmort)*(1-annualexploit) #growth-survival constant (pg 339 Hilborn and Walters)
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium N
    Ceq<-eqvalue*annualexploit #equilibrium catch
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Ceq  #equilibrium spawning stock
  } else {
    kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/FWalpha  #growth-survival constant (pg 339 Hilborn and Walters)
    #kappa<-((1-(1+FWrho)*(1-annualnatmort)*(1-annualexploit))+(FWrho*(1-annualnatmort)^2*(1-annualexploit)^2))/(Recwt-((Prerecwt*FWrho)*(1-annualnatmort)*(1-annualexploit)))
    eqvalue<-(predalpha/kappa)-(predbeta/(1-annualexploit)) #equilibrium Biomass
    Yeq<-eqvalue*annualexploit #equilibrium yield
    Req<-kappa*eqvalue #equilibrium recruitment
    SSeq<-eqvalue-Yeq #equilibrium spawning stock
    Neq<-Req/(1-(1-annualnatmort)*(1-annualexploit)) #equilibrium N
    Weq<-eqvalue/Neq #equilibrium average weight
  }  
  if(unfishflag==0){ #if doing unfished then just return unfished spr to optimize function; required for equlibrium calcs
  return((targzero-eqvalue)^2)
  } else { #if not doing unfished (exploit>0) then return various equilibrium values
   if(numbs_weightb==1) { #units of numbers or biomass?
     return(data.frame("Neq"=eqvalue,"Ceq"=Ceq,"Req"=Req,"SSeq"=SSeq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }else {
     return(data.frame("Beq"=eqvalue,"Yeq"=Yeq,"Req"=Req,"SSeq"=SSeq,"Neq"=Neq,"Weq"=Weq,"predrzero"=predrzero,"predalpha"=predalpha,"predbeta"=predbeta))
   }
  }
}
#####end function unfishspr#####

####User Inputs###########

#require(plyr)
require(dplyr)

setwd("~/Data/Projects/MSE/HerringMSE")
#setwd("~/Documents/0_Data/MSE/HerringMSE")

preyB<-read.table("Unadj1111TotBioSimYear.txt", header=T)  #will need to look for different scenario numbers in title and folder, splice out and label outputs here
preysim<-read.table("Unadj1111NAASimYear.txt", header=T)
#names(preysim)<-c("a1","a2","a3","a4","a5","a6", "a7","a8","simN")
#preyNsim<-transmute(preysim, preyN=a1+a2+a3+a4+a5+a6+a7+a8, simN=Sim)
preyNsim<-transmute(preysim, preyN=Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8, Sim=Sim)

#nyears<-28 #number of years
by_simN<-group_by(preyNsim, Sim)
Nyears<-summarise(by_simN, n=n())

preychar<-read.table("Unadj1111SimCharAA.txt", header=T)
preychar<-cbind(preychar, Age=rep(1:8,max(preychar$Sim)))

preychar2<-read.table("Unadj1111SimCharAA_OldWt.txt", header=T)
preychar2<-cbind(preychar2, Age=rep(1:8,max(preychar2$Sim)))

#preyAvgWt<-#for each Sim, for ages >2, multiply N at age by Fselectivity and Weight, then sum it. divide by sum (Fselectivity*N at age)
BFTforage<-select(preysim, Age3:Sim)

BFTforage<- cbind(BFTforage, Yr=rep(1:Nyears$n, max(preychar$Sim)))

BFTforage <- reshape(BFTforage, 
             varying = c("Age3", "Age4", "Age5", "Age6", "Age7", "Age8"), 
             v.names = "N",
             timevar = "Age", 
             times = c(3,4,5,6,7,8), 
             direction = "long")

BFTforagechar <- left_join(BFTforage, preychar)

BFTforagewtage <- BFTforagechar %>%
                     group_by(Sim,Yr) %>%
                     mutate(preyAgeAvgWt=N*Weight*Fselectivity, preyAgeselN=N*Fselectivity) %>%
                     summarise(preyAvgWt = sum(preyAgeAvgWt)/sum(preyAgeselN))

  
# Get vector of pars for 4 pred types (groundfish, tuna, seabird, marine mammal)
# two herring dependence types each, "sensitive" and "robust"

#testing WBFT pars (see HerringMSEpredspars.xlsx)

numbs_weight  <-  2  #2 #1=numbers, else weight.  Will predator be tracked just in numbers or numbers and weight?  Weight allows for growth effects.
Predzero      <-  6.69E+04#5.90E+05 # #5.60E+09  #50000 #predator unfished total N or biomass
Predsteep     <-  1.0#0.51 # #0.45  #0.8 #predator steepness
PredannualA   <-  0.14  #0.6 #predator annual natural mortality rate (0-1); serves as max if time varying
Predexploit   <-  0.079  #0.3 #predator annual exploitation rate (0-1)
FWalpha       <-  0.020605  #0.243 #ford-walford plot intercept parameter if doing biomass units tons
FWrho         <-  0.9675  #1.45 #ford-walford plot slope parameter if doing biomass units; serves as max if time varying

Recwt         <-  0.041  #age 6 not 1 0.00278  #avg wt at recruitment age, t (wt at PredRecdelay)
Prerecwt      <-  0.027  #age 5 not 0 0.00047  #avg wt at rec age - 1, t

PredN         <-  c() #for predator abundance vector
PredN[1]      <-  111864  #7000 #initial abundance in year 1
PredB         <-  c() #for predator biomass vector if doing in those units
PredB[1]      <-  27966  #6000 #initial biomass in year 1 
PredRecdelay  <-  1  #1 #delay in years for when recruits in year y are added to population

preypredrec<-1 #strength of effect of prey N on predator recruitment (like Plaganyi and Butterworth) >=1; 1=no effect
preyprednatm<-0 #strength of effect of prey N on predator annual natural mort (0-1) 0=no effect
preypredgrow<-1.1 #strengh of effect of prey N on predator growth, >=1, 1=no effect

BFT<-TRUE

doplotsbysim<-TRUE


#####End user inputs########

#Do predator unfished equilibrium calculation and get unfished spr
predunfishsprfxn<-optimize(f=unfishspr,interval=c(0,50),targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=0,FWalpha=FWalpha,FWrho=FWrho,unfishflag=0,tol=0.0000001)
predunfishspr<-predunfishsprfxn$minimum
#send unfished spr through function again with an exploitation rate to get equlibrium conditions
PredEquilib<-unfishspr(predunfishspr=predunfishspr,targzero=Predzero,steep=Predsteep,numbs_weightb=numbs_weight,annualnatmort=PredannualA,annualexploit=Predexploit,FWalpha=FWalpha,FWrho=FWrho,unfishflag=1)
predalpha<-PredEquilib$predalpha #predator BH SR parm
predbeta<-PredEquilib$predbeta #predator BH SR parm
Req<-PredEquilib$Req
Weq<-PredEquilib$Weq

if(BFT){
  #BFT test only Porch and Lauretta 2016 PLoS ONE
  #predunfishspr<-0.7
  #predalpha<- #9.18E+05 #predator BH SR parm eq 7 tau = 1.2
  #predbeta<- #1.48E+05 #predator BH SR parm eq 7 tau = 1.2
  SSB_MSY<-13226  #63102
  BFTGrowThresh<-0.15
}
  

## just a quick loop to generate prey numbers for testing.  Could be any prey t.s. This justmade up.
#preyNzero<-11000
#preyN<-c()
#preyN[1]<-1000 #n in year 1
#for(y in 1:(nyears-1)){
#   if((y+1)<=(nyears/2)) {
#     preyN[y+1]<-preyN[y]+(preyN[y]*0.2) #increase N by 20% each year for half of years
#   } else {
#     preyN[y+1]<-preyN[y]-(preyN[y]*0.2) #decline by 20% each year for last half of years
#   }
#   if(preyN[y+1]<0) {preyN[y+1]==0} #make sure we don't go negative
#}
####end generic preyN loop

# add loop for multiple preds (inside or outside?)

# add loop for multiple sims per pred
for (h in 1:max(preyNsim$Sim)){
  preyN<-preyNsim$preyN[preyNsim$Sim==h]
  nyears<-Nyears$n[h]
  ## WARNING HACK need actual unfished herring B
  #preyNzero<-max(preyN)
  #sum unfished N across ages for this sim in preychar
  preyNzero<-sum(preychar$UnfishedNAA[preychar$Sim==h])
  preyTotB<-preyB$TotalBio[preyB$Sim==h]
  preyAvgWt<-BFTforagewtage$preyAvgWt[BFTforagewtage$Sim==h]

  ###predator dynamics###
  #stuff I need
  Recmult<-(preypredrec*(preyN/preyNzero))/((preypredrec-1)+(preyN/preyNzero)) #fraction of expected predator BH recruitment dependent on prey N
  PredAnnualNatM<-PredannualA*exp(-(preyN/preyNzero)*preyprednatm) #needed if annual nat mort is time varying
  TotalSurv<-(1-PredAnnualNatM)*(1-Predexploit) #total annual predator survival
  catch<-c()
  Spawn<-c() #spawners in numbers or biomass depending
  yield<-c()
  if(BFT){ #base growth on prey avereage weight; use generalized logistic with lower bound on growth rate (97% of 102% of FWrho to start, trying to center on FWrho)
    AnnualAlpha<-(0.9*FWalpha) + ((1.1*FWalpha) - (0.9*FWalpha))/(1+exp((1-preypredgrow)*(100*(preyAvgWt-BFTGrowThresh)/BFTGrowThresh))) #alpha changes with herring avg wt, not slope
    } else { #base growth on abundance of prey
    AnnualGrowParm<-FWrho*((preypredgrow*(preyN/preyNzero))/((preypredgrow-1)+(preyN/preyNzero))) #if grow time varies then this is annual FW slope
    }
  Predrec<-c()
  Predrec[1:PredRecdelay]<-Req #set recruitment in initial years at equlibrium to account for delay
  ###year loop for predator dynamics
  for(y in 1:(nyears-1)){
    if(numbs_weight==1){
      catch[y]<-PredN[y]*Predexploit
      Spawn[y]<-PredN[y]-catch[y]
      } else {
      yield[y]<-PredB[y]*Predexploit
      Spawn[y]<-PredB[y]-yield[y]
     }
   Predrec[y+PredRecdelay]<-Recmult[y]*((predalpha*Spawn[y])/(predbeta+Spawn[y])) #SR
   #Predrec[y+PredRecdelay]<-Recmult[y]*(Spawn[y])/(predalpha+(predbeta*Spawn[y])) #SR from Mangel
   PredN[y+1]<-PredN[y]*TotalSurv[y]+Predrec[y+1]
    if(numbs_weight!=1){ #only do biomass calcs if requested in those units
      if(BFT){
        PredB[y+1]<-TotalSurv[y]*(AnnualAlpha[y]*PredN[y]+FWrho*PredB[y])+AnnualAlpha[y]*Predrec[y+1]
      } else {
        PredB[y+1]<-TotalSurv[y]*(FWalpha*PredN[y]+AnnualGrowParm[y]*PredB[y])+FWalpha*Predrec[y+1]
      }
    }  
  } #end y loop for predators
  
  if(h==1){
    predOut<-data.frame(Sim=h, Yr=c(1:nyears), PredB=PredB, PredN=PredN, PredRec=Predrec, PredAvWt=PredB/PredN, PredB_status=PredB/SSB_MSY, PredAvWt_status=(PredB/PredN)/Weq )
  } else {
    predOutn<-data.frame(Sim=h, Yr=c(1:nyears), PredB=PredB, PredN=PredN, PredRec=Predrec, PredAvWt=PredB/PredN, PredB_status=PredB/SSB_MSY, PredAvWt_status=(PredB/PredN)/Weq )
    predOut<-bind_rows(predOut, predOutn)
  }
  
  
  if(doplotsbysim){
    par(mfrow=c(3,2))
    par(mar=c(4,4,2,2)+0.1)
    par(oma=c(2,2,2,0))
    plot(PredN,type='l',col="black",xlab="Year",ylab="Predator Abundance",lwd=2, ylim=c(0, max(PredN)))
    plot(preyN,Recmult,col="black",xlab="Prey Abundance",ylab="Recruitment Fraction",lwd=2)
    title(paste("preypredrec = ", preypredrec), line=-2)
    plot(preyN,PredAnnualNatM,col="black",xlab="Prey Abundance",ylab="Annual Natural Mortality Rate",lwd=2)
    title(paste("preyprednatm = ", preyprednatm), line=-2)
    if(BFT){
      plot(preyAvgWt,AnnualAlpha,col="black",xlab="Prey Avg Wt",ylab="Growth Intercept",lwd=2)
      title(paste("preypredgrow = ", preypredgrow), line=-2)
    } else {
      plot(preyN,AnnualGrowParm,col="black",xlab="Prey Abundance",ylab="Growth Rate",lwd=2)
      title(paste("preypredgrow = ", preypredgrow), line=-2)
    }
    plot(preyN, type="l", xlab="Year", ylab="Prey Abundance", lwd=2)
    
    if(numbs_weight != 1){ #only do if biomass units requested
    plot(PredB,type='l',col="black",xlab="Year",ylab="Predator Biomass",lwd=2, ylim=c(0, max(PredB))) 
    abline(h=SSB_MSY)
    par(new=T)
    plot(PredB/PredN, type="l", col="blue", axes=F, xlab=NA, ylab=NA, ylim=c(0, max(PredB/PredN)))
    axis(side = 4)
    mtext(side = 4, line = 3, 'Pop avg wt')
    lines(PredB/PredN, type='l',col="blue")
    abline(h=Weq, col="blue", lty=3)
    }
    mtext(paste("Sim ", h, sep=""), outer=T, side=3)
  } # end doplots
} # end loop over prey sims

write.table(predOut, "test.txt") #put control rule id in filename

par(mfrow=c(3,1))

plot(predOut$Yr[predOut$Sim==1], predOut$PredB[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredB", xlab="")
for (j in 1:max(predOut$Sim)){
  lines(predOut$Yr[predOut$Sim==j], predOut$PredB[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
}

plot(predOut$Yr[predOut$Sim==1], predOut$PredN[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredN", xlab="")
for (j in 1:max(predOut$Sim)){
  lines(predOut$Yr[predOut$Sim==j], predOut$PredN[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
}

plot(predOut$Yr[predOut$Sim==1], predOut$PredAvWt_status[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredAvWt", xlab="")
for (j in 1:max(predOut$Sim)){
  lines(predOut$Yr[predOut$Sim==j], predOut$PredAvWt_status[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
}
abline(h=1.0, col="blue", lwd=3)
mtext(paste("Control Rule ", cr, sep=""), outer=T, side=3)



library(ggplot2)

ggplot(predOut, aes(x=Yr, y=PredB, group=Sim)) +
  geom_line()
ggplot(predOut, aes(x=Yr, y=PredN, group=Sim)) +
  geom_line()
ggplot(predOut, aes(x=Yr, y=PredAvWt_status, group=Sim)) +
  geom_line()


