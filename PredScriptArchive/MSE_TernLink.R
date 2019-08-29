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

#### Common Tern Gulf of Maine pars (see HerringMSEpredspars.xlsx for derivation

numbs_weight  <-  1  #2 #1=numbers, else weight.  Will predator be tracked just in numbers or numbers and weight?  Weight allows for growth effects.
Predzero      <-  45000 #5.90E+05 # #5.60E+09  #50000 #predator unfished total N or biomass
Predsteep     <-  0.26 #0.25  #0.8 #predator steepness
PredannualA   <-  0.1  #0.6 #predator annual natural mortality rate (0-1); serves as max if time varying
Predexploit   <-  0.00  #0.3 #predator annual exploitation rate (0-1)
FWalpha       <-  0.00015#1.5 #  #0.243 #ford-walford plot intercept parameter if doing biomass units tons
FWrho         <-  0.0  #1.45 #ford-walford plot slope parameter if doing biomass units; serves as max if time varying

#Recwt         <-  0.041  #age 6 not 1 0.00278  #avg wt at recruitment age, t (wt at PredRecdelay)
#Prerecwt      <-  0.027  #age 5 not 0 0.00047  #avg wt at rec age - 1, t

PredN         <-  c() #for predator abundance vector
PredN[1]      <-  3000  #7000 #initial abundance in year 1
PredB         <-  c() #for predator biomass vector if doing in those units
PredB[1]      <-  1.5  #6000 #initial biomass in year 1 
PredRecdelay  <-  4  #1 #delay in years for when recruits in year y are added to population

preypredrec<-1.09 #strength of effect of prey N on predator recruitment (like Plaganyi and Butterworth) >=1; 1=no effect
preyprednatm<-0 #strength of effect of prey N on predator annual natural mort (0-1) 0=no effect
preypredgrow<-1 #strengh of effect of prey N on predator growth, >=1, 1=no effect

BFT<-F
COTE<-T

if(COTE){
  #pars below for alternative steepness=0.41, same B0=45000
  #predunfishspr<-15
  #predalpha<- 4686
  #predbeta<- 25286
  COTEprodThresh<-400000
}

doplotsbysim<-F
dosummaryplots<-T


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


###################################################################

#### Read in Herring "data" for each OM and Control Rule ###########

require(dplyr)

#setwd("~/Data/Projects/MSE/HerringMSE")
#setwd("~/Documents/0_Data/MSE/HerringMSE")

#Set up directory to save output and graphics
#ResultsDirec<-c('//net/home5/jderoba/MSE/HiM_LowSteep_NoAssBias_OldWt/')
PredDirec<-c('~/Data/Projects/MSE/HerringMSE/Terns')

#OMDirec<-c('//net/home5/jderoba/MSE/HiM_LowSteep_NoAssBias_RecWt/')   #Change for different OM combinations HOW TO ACCESS FROM MY COMPUTER?
OMDirec<-c('//net/mse/')
#OMtype<-c('HiM_LowSteep_NoAssBias_OldWt') #done first
OMlist<-c('LoM_HiSteep_NoAssBias_OldWt','HiM_LowSteep_NoAssBias_OldWt','LoM_HiSteep_NoAssBias_RecWt', 'HiM_LowSteep_NoAssBias_RecWt')
CRtype<-c("BB", "BB3yr", "BB5yr", "BB3yrPerc", "CC", "CCC")

for(OMtype in OMlist){

dir.create(PredDirecOM<-file.path(PredDirec,paste(OMtype)))

for(crt in CRtype){
  dir.create(PredDirecCR<-file.path(PredDirecOM,paste(crt)))
  PreyDirecCR<-paste(OMDirec,OMtype,"/",crt,"/", sep="")
  
  #get control rule number from filenames
  library(stringr)
  k<-list.files(path=PreyDirecCR,pattern="\\SimCharAA.txt$")
  # prepare regular expression
  regexp <- "[[:digit:]]+"
  
  CRnum<-str_extract(k, regexp)
    
  for(crnum in CRnum){

     filename1<-paste("Unadj", crnum, "NAASimYear.txt", sep="")
     filename2<-paste("Unadj", crnum, "SimCharAA.txt", sep="")
     filename3<-paste("Unadj", crnum, "TotBioSimYear.txt", sep="")
    
     NAAfile<-paste(PreyDirecCR,"/",filename1, sep="")
     charfile<-paste(PreyDirecCR,"/",filename2, sep="")
     totBfile<-paste(PreyDirecCR,"/",filename3, sep="")

     preyB<-read.table(totBfile, header=T)  
     preysim<-read.table(NAAfile, header=T)
     preyNsim<-transmute(preysim, preyN=Age1+Age2+Age3+Age4+Age5+Age6+Age7+Age8, Sim=Sim)
     
     by_simN<-group_by(preyNsim, Sim)
     Nyears<-summarise(by_simN, n=n())
     
     preychar<-read.table(charfile, header=T)
     preychar<-cbind(preychar, Age=rep(1:8,max(preychar$Sim)))
     
     #preyAvgWt<-#for each Sim, for ages <2, keep for tern forage
     #Ternforage<-select(preysim, Age1, Sim)
     
     #BASE ON TOTAL B
     Ternforage<-cbind(preyB, Yr=rep(1:Nyears$n, max(preychar$Sim)))
     
     #Ternforage<- cbind(Ternforage, Yr=rep(1:Nyears$n, max(preychar$Sim)))
     
     
     # add loop for multiple sims per pred
     for (h in 1:max(preyNsim$Sim)){
       #preyN for terns is just age 1
       #preyN<-Ternforage$Age1[Ternforage$Sim==h]
       #preyN for terns is totalB
       preyN<-Ternforage$TotalBio[Ternforage$Sim==h]
       #preyN<-preyNsim$preyN[preyNsim$Sim==h] #this would be if all ages were forage
       nyears<-Nyears$n[h]
       #preyNzero for terns is just age1 unfished
       #preyNzero<-preychar$UnfishedNAA[preychar$Age==1 & preychar$Sim==h]
       #preyNzero for terns is the threshold total B where prod drops, FIXED based on GOM data
       preyNzero<-COTEprodThresh

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
       
       #construct and append to dataframe for all sims
       if(numbs_weight==1){
         if(h==1){
           #predOut<-data.frame(Sim=h, Yr=c(1:nyears), PredN=PredN, PredRec=c(Predrec[5:153], 0), PredProd=c(Predrec[5:153]*10/PredN[1:149], 0))
           predOut<-data.frame(Sim=h, Yr=c(1:nyears), PredN=PredN, PredRec=Predrec[4:153], PredProd=c(0,Predrec[5:153]*10/PredN[1:149]))
         } else {
           #predOutn<-data.frame(Sim=h, Yr=c(1:nyears), PredN=PredN, PredRec=c(Predrec[5:153], 0), PredProd=c(Predrec[5:153]*10/PredN[1:149], 0))
           predOutn<-data.frame(Sim=h, Yr=c(1:nyears), PredN=PredN, PredRec=Predrec[4:153], PredProd=c(0, Predrec[5:153]*10/PredN[1:149]))
           predOut<-bind_rows(predOut, predOutn)
         }
         
       } else {
         if(h==1){
           predOut<-data.frame(Sim=h, Yr=c(1:nyears), PredB=PredB, PredN=PredN, PredRec=Predrec, PredAvWt=PredB/PredN, PredB_status=PredB/SSB_MSY, PredAvWt_status=(PredB/PredN)/Weq )
         } else {
           predOutn<-data.frame(Sim=h, Yr=c(1:nyears), PredB=PredB, PredN=PredN, PredRec=Predrec, PredAvWt=PredB/PredN, PredB_status=PredB/SSB_MSY, PredAvWt_status=(PredB/PredN)/Weq )
           predOut<-bind_rows(predOut, predOutn)
         }
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
         if(numbs_weight == 1){ #only do for numbers 
           plot(Predrec[5:153]*10/PredN[1:149], type="l", xlab="Year", ylab="Predator Productivity", lwd=2, ylim=c(0,2))
           abline(h=1, col="blue", lty=3)
         }
         
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
     
     write.table(predOut, paste(PredDirecCR,"/",crnum,"predBNAvWtStatus.txt", sep="")) #put control rule id in filename
                 
     if(dosummaryplots) {
       par(mfrow=c(3,1))
       if(numbs_weight != 1){ #only do if biomass units requested
         plot(predOut$Yr[predOut$Sim==1], predOut$PredB[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredB", xlab="")
         for (j in 1:max(predOut$Sim)){
           lines(predOut$Yr[predOut$Sim==j], predOut$PredB[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
         }
       }
       plot(predOut$Yr[predOut$Sim==1], predOut$PredN[predOut$Sim==1], ylim=c(0,max(predOut$PredN)),type="l", col = rgb(0, 0, 0, 0.3), ylab="PredN", xlab="")
       for (j in 1:max(predOut$Sim)){
         lines(predOut$Yr[predOut$Sim==j], predOut$PredN[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
       }
       if(numbs_weight != 1){ #only do if biomass units requested  
         plot(predOut$Yr[predOut$Sim==1], predOut$PredAvWt_status[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredAvWt", xlab="")
         for (j in 1:max(predOut$Sim)){
           lines(predOut$Yr[predOut$Sim==j], predOut$PredAvWt_status[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
         }
         abline(h=1.0, col="blue", lwd=3)
         mtext(paste("Control Rule ", crnum, sep=""), outer=T, side=3)
       }
       if(numbs_weight == 1){ #only do if numbers units requested
         plot(predOut$Yr[predOut$Sim==1], predOut$PredProd[predOut$Sim==1], type="l", col = rgb(0, 0, 0, 0.3), ylab="PredProd", xlab="")
         for (j in 1:max(predOut$Sim)){
           lines(predOut$Yr[predOut$Sim==j], predOut$PredProd[predOut$Sim==j], type="l", col = rgb(0, 0, 0, 0.3))
         }
         abline(h=1.0, col="blue", lwd=3)
         mtext(paste("Control Rule ", crnum, sep=""), outer=T, side=3)
       }
     }
     
  }#end crnum loop over control rule variants 
  
}#end CRtype loop over control rule type folders

}#end OMtype list over operating model folders
