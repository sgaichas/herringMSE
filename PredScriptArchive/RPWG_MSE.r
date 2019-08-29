rm(list=ls(all=TRUE))  #Remove all objects that currently exist in memory
#CodeDirect<-"C:\\Herring\\MSE ABC Control Rule\\MSE\\"  #location of functions.r and userinputs.r
#CodeDirect<-"C:\\MSE ABC Control Rule\\MSE\\"
CodeDirect<-c('//net/home5/jderoba/MSE/')
source(paste(CodeDirect,"UserInputs.r",sep="")) #source user inputs
source(paste(CodeDirect,"SimsFunction.R",sep="")) #source the simulation loop code

#########################Go get user inputs##################
inputs<-userinputs()
###############################

########harvest Control rule stuff
CtrlRule<-inputs$CtrlRule
FracBmsyThreshHi<-inputs$FracBmsyThreshHi # For control rule; Threshold to switch from Fmsy as target F to linear decline to zero
FracBmsyThreshLo<-inputs$FracBmsyThreshLo #For control rule; Level of SSB where target F set to 0
FracFtarg<-inputs$FracFtarg #fraction of Fmsy that serves as max target F in control rule
useminmax<-inputs$useminmax #use min and max quota levels for biobased CR?
minquota<-inputs$minquota #min quota if min max invoked
maxquota<-inputs$maxquota #max quota if min max invoked
usequotavar<-inputs$usequotavar #control interannual variation in quota? 0=no; >0 yes
quotavar<-inputs$quotavar #fraction that quota can vary year to year if using 'usequotavar'
Thresh<-inputs$Thresh #fraction of Bmsy where target yield set to 0; for proportional threshold
Prop<-inputs$Prop #fraction of fully selected Fmsy to be harvested when population above Thresh (converted to exploitation rate in sim code)
Propc<-inputs$Propc #fraction of MSY for constant catch control rule
useFcap<-inputs$useFcap #use conditional constant catch with max F?  0=no >0 is yes
Fcapprop<-inputs$Fcapprop #fraction of fully selected Fmsy to use for max F in conditional constant catch
OverfishedSSB<-inputs$OverfishedSSB  #overfished definition for %collapse calc.  Not used in control rule; only as measure of collapse
quotablock<-inputs$quotablock #quota will be set to same value for this number of years (control rule only applied this frequently)
################################

if(CtrlRule==1){ #biomass based
for(f in 1:(length(FracFtarg))) { #max f for biomass based ctrl rule loop
for(u in 1:(length(FracBmsyThreshHi))){ #upper threshold loop
  for(l in 1:u){ #lower threshold loop
  c<-ifelse(l<=length(Propc),l,length(Propc)) #c is only for constant catch but function expects something and so this isn't important; just filler
  sims(u=u,l=l,f=f,c=c,CtrlRule=CtrlRule,FracBmsyThreshHi=FracBmsyThreshHi[u],FracBmsyThreshLo=FracBmsyThreshLo[l],FracFtarg=FracFtarg[f],useminmax=useminmax,minquota=minquota,maxquota=maxquota,usequotavar=usequotavar,quotavar=quotavar,Thresh=Thresh,Prop=Prop,Propc=Propc[c],useFcap=useFcap,Fcapprop=Fcapprop,OverfishedSSB=OverfishedSSB,quotablock=quotablock,CodeDirect=CodeDirect)
  } #end l loop
} #end u loop
} #end f loop
} else if(CtrlRule==3){ #constant catch
  for(c in 1:(length(Propc))) {
  u<-ifelse(c<=length(FracBmsyThreshHi),c,length(FracBmsyThreshHi))
  l<-ifelse(c<=length(FracBmsyThreshLo),c,length(FracBmsyThreshLo))
  f<-ifelse(c<=length(FracBmsyThreshLo),c,length(FracFtarg))
  sims(u=u,l=l,f=f,c=c,CtrlRule=CtrlRule,FracBmsyThreshHi=FracBmsyThreshHi[u],FracBmsyThreshLo=FracBmsyThreshLo[l],FracFtarg=FracFtarg[f],useminmax=useminmax,minquota=minquota,maxquota=maxquota,usequotavar=usequotavar,quotavar=quotavar,Thresh=Thresh,Prop=Prop,Propc=Propc[c],useFcap=useFcap,Fcapprop=Fcapprop,OverfishedSSB=OverfishedSSB,quotablock=quotablock,CodeDirect=CodeDirect)
  } # end c loop
} else { print ("Poop: unknown control rule") }

print("Simulations Complete")
  
