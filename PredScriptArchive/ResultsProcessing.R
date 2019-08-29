rm(list=ls(all=TRUE))  #Remove all objects that currently exist in memory
#CodeDirect<-"F:\\Herring\\MSE ABC Control Rule\\MSE\\"  #location of functions.r and userinputs.r
CodeDirect<-c('//net/home5/jderoba/MSE/')
source(paste(CodeDirect,"UserInputs.r",sep="")) #source user inputs

wantyears<-50 #number of years to take median over for each sim for results

#########################Go get user inputs##################
inputs<-userinputs()
###############################

########harvest Control rule stuff
nages<-inputs$nages
nyears<-inputs$nyears
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
nsims<-inputs$nsims
################################
quartilefun<-function(data=NA,quart=NA){ #function to return desired quartile from desired data
 quartile<-quantile(x=data,probs=quart,na.rm=TRUE)
 return(quartile)
}
ResFun<-function(u=NA,l=NA,f=NA,c=NA,CtrlRule=NA,FracBmsyThreshHi=NA,FracBmsyThreshLo=NA,FracFtarg=NA,useminmax=NA,minquota=NA,maxquota=NA,usequotavar=NA,quotavar=NA,Thresh=NA,Prop=NA,Propc=NA,useFcap=NA,Fcapprop=NA,quotablock=NA,nsims=NA,wantyears=NA){
  ctrlrulenum<-paste(u,l,f,c,sep="") #results recorded this way and so referenced this way here
  ####SSB
  SSBSimsMed<-matrix(data=NA,nrow=nsims,ncol=7) #holds median among 'wantyears' for each sim
  SSBSimYear<-read.table(paste(paste(readfile,ctrlrulenum,sep=""),"SSBSimYear.txt",sep=""),header=T)  #read SSB by year and sim
  SSBSimYear$SSBrelSSBmsy<-SSBSimYear$SSB/SSBSimYear$Bmsy #SSB over SSBmsy
  SSBSimYear$SSBpropSSBmsy<-ifelse(SSBSimYear$SSBrelSSBmsy<1,1,0) #if B<Bmsy then 1, else 0.
  SSBSimYear$SSBprophalfSSBmsy<-ifelse(SSBSimYear$SSBrelSSBmsy<0.5,1,0) #if B<0.5Bmsy  then 1, else 0.
  SSBSimYear$SSBrelSSBzero<-SSBSimYear$SSB/SSBSimYear$Bzero #SSB over SSBzero
  SSBSimYear$SSBpropthreeSSBzero<-ifelse(SSBSimYear$SSBrelSSBzero<0.3,1,0) #if B<0.3Bzero then 1, else 0.
  SSBSimYear$SSBprop75SSBzero<-ifelse(SSBSimYear$SSBrelSSBzero<0.75,1,0) #if B<0.75Bzero then 1, else 0.
  ####End SSB
  ####F
  TrueFSimsMed<-matrix(data=NA,nrow=nsims,ncol=1) #holds median among 'wantyears' for each sim
  TrueFSimYear<-read.table(paste(paste(readfile,ctrlrulenum,sep=""),"TrueFSimYear.txt",sep=""),header=T)  #read fully sel F by year and sim
  TrueFSimYear$FrelFmsy<-TrueFSimYear$TrueFullF/TrueFSimYear$Fmsy #F relative to Fmsy  
  TrueFSimYear$FpropFmsy<-ifelse(TrueFSimYear$FrelFmsy>1,1,0) #if F>Fmsy then 1, else 0
  ####End F
  ####Yield and NatDeaths
  TotYieldSimMed<-matrix(data=NA,nrow=nsims,ncol=6) #holds median among 'wantyears' for each sim
  TotYieldSimYear<-read.table(paste(paste(readfile,ctrlrulenum,sep=""),"TotYieldSimYear.txt",sep=""),header=T)  #read Yield by year and sim
  TotYieldSimYear$YieldrelMSY<-TotYieldSimYear$Yield/TotYieldSimYear$MSY #yield relative to MSY
  TotYieldSimYear$Closure<-ifelse(TotYieldSimYear$TargetQuota==0,1,0) #if fishery closed (i.e. targquota=0) then 1, else 0.
  TotYieldSimYear$YieldrelNatDeath<-TotYieldSimYear$Yield/TotYieldSimYear$NatDeaths #yield relative to natural deaths
  ####End Yield
  ####NAA
  NAASimMed<-matrix(data=NA,nrow=nsims,ncol=(nages+4))
  NAASimYear<-read.table(paste(paste(readfile,ctrlrulenum,sep=""),"NAASimYear.txt",sep=""),header=T)  #read NAA by year and sim
  NAASimYear$TotN<-rowSums(NAASimYear)-NAASimYear$Sim #total abundance each sim and year
  NAASimYear$PropOne<-NAASimYear$Age1/NAASimYear$TotN #proportion age one
  NAASimYear$PropTwo<-NAASimYear$Age2/NAASimYear$TotN #proportion age two
  NAASimYear$PropOneTwo<-NAASimYear$PropOne+NAASimYear$PropTwo #proportion age one and two
  ####End NAA
  ####Total biomass for Surplus Production
  SPSimsMed<-matrix(data=NA,nrow=nsims,ncol=1) #holds median among 'wantyears' for each sim
  TotBioSimYear<-read.table(paste(paste(readfile,ctrlrulenum,sep=""),"TotBioSimYear.txt",sep=""),header=T)  #read TotBio by year and sim
  ####End Surplus Production
  for(s in 1:nsims){ 
    ####SSB
    SSBSimYearSub<-SSBSimYear[SSBSimYear$Sim==s,c("SSB","SSBrelSSBmsy","SSBrelSSBzero","SSBpropSSBmsy","SSBprophalfSSBmsy","SSBpropthreeSSBzero","SSBprop75SSBzero")] #extract just results for 1 sim
    SSBSimsMed[s,1]<-median(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSB"],na.rm=T) #median among year for each sim
    SSBSimsMed[s,2]<-median(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBrelSSBmsy"],na.rm=T) #median among year for each sim 
    SSBSimsMed[s,3]<-median(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBrelSSBzero"],na.rm=T) #median among year for each sim
    SSBSimsMed[s,4]<-sum(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBpropSSBmsy"],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where B<Bmsy
    SSBSimsMed[s,5]<-sum(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBprophalfSSBmsy"],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where B<Bmsy
    SSBSimsMed[s,6]<-sum(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBpropthreeSSBzero"],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where B<0.3Bzero
    SSBSimsMed[s,7]<-sum(SSBSimYearSub[c((nrow(SSBSimYearSub)-wantyears):nrow(SSBSimYearSub)),"SSBprop75SSBzero"],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where B<0.75Bzero
    ####End SSB
    ####F
    TrueFSimYearSub<-TrueFSimYear[TrueFSimYear$Sim==s,"FpropFmsy"] #extract just results for 1 sim
    TrueFSimsMed[s,1]<-sum(TrueFSimYearSub[c((length(TrueFSimYearSub)-(wantyears+1)):(length(TrueFSimYearSub)-1))],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where F>Fmsy
    ####End F
    ####Yield and Nat Deaths
    TotYieldSimYearSub<-TotYieldSimYear[TotYieldSimYear$Sim==s,c("Yield","NatDeaths","YieldrelMSY","Closure","YieldrelNatDeath")] #extract just results for 1 sim
    TotYieldSimYearSub[1:(nrow(TotYieldSimYearSub)-2),"Yvar"]<-(TotYieldSimYearSub[2:(nrow(TotYieldSimYearSub)-1),"Yield"]-TotYieldSimYearSub[1:(nrow(TotYieldSimYearSub)-2),"Yield"])^2
    TotYieldSimMed[s,1]<-median(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+1)):(nrow(TotYieldSimYearSub)-1)),"Yield"],na.rm=T) #median among year for each sim
    TotYieldSimMed[s,2]<-median(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+1)):(nrow(TotYieldSimYearSub)-1)),"NatDeaths"],na.rm=T) #median among year for each sim
    TotYieldSimMed[s,3]<-median(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+1)):(nrow(TotYieldSimYearSub)-1)),"YieldrelMSY"],na.rm=T) #median among year for each sim
    TotYieldSimMed[s,4]<-median(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+1)):(nrow(TotYieldSimYearSub)-1)),"YieldrelNatDeath"],na.rm=T) #median among year for each sim
    TotYieldSimMed[s,5]<-sum(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+1)):(nrow(TotYieldSimYearSub)-1)),"Closure"],na.rm=T)/(wantyears+1) #number of last wantyears of each sim where fishery closed
    TotYieldSimMed[s,6]<-sqrt(mean(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+2)):(nrow(TotYieldSimYearSub)-2)),"Yvar"],na.rm=T))/mean(TotYieldSimYearSub[c((nrow(TotYieldSimYearSub)-(wantyears+2)):(nrow(TotYieldSimYearSub)-2)),"Yield"],na.rm=T)
    ####End Yield
    ####NAA
    NAASimYearSub<-NAASimYear[NAASimYear$Sim==s,c(paste("Age",seq(1,nages,1),sep=""),"TotN","PropOne","PropTwo","PropOneTwo")] #extract just results for 1 sim
    for(a in 1:nages){
      NAASimMed[s,a]<-median(NAASimYearSub[c((nrow(NAASimYearSub)-wantyears):nrow(NAASimYearSub)),paste("Age",a,sep="")],na.rm=T) #median among year for each sim      
    }
    NAASimMed[s,(nages+1)]<-median(NAASimYearSub[c((nrow(NAASimYearSub)-wantyears):nrow(NAASimYearSub)),"TotN"],na.rm=T) #median among year for each sim 
    NAASimMed[s,(nages+2)]<-median(NAASimYearSub[c((nrow(NAASimYearSub)-wantyears):nrow(NAASimYearSub)),"PropOne"],na.rm=T) #median among year for each sim      
    NAASimMed[s,(nages+3)]<-median(NAASimYearSub[c((nrow(NAASimYearSub)-wantyears):nrow(NAASimYearSub)),"PropTwo"],na.rm=T) #median among year for each sim      
    NAASimMed[s,(nages+4)]<-median(NAASimYearSub[c((nrow(NAASimYearSub)-wantyears):nrow(NAASimYearSub)),"PropOneTwo"],na.rm=T) #median among year for each sim      
    ####End NAA
    ####Surplus Production
    TotBioSimYearSub<-TotBioSimYear[SSBSimYear$Sim==s,"TotalBio"] #extract just results for 1 sim
    TotBioSimYearSub<-data.frame(cbind(TotBioSimYearSub,TotYieldSimYearSub$Yield))
    names(TotBioSimYearSub)<-c("TotalBio","Yield")
    TotBioSimYearSub[2:nyears,"SP"]<-TotBioSimYearSub[2:nyears,"TotalBio"]-TotBioSimYearSub[1:(nyears-1),"TotalBio"]+TotBioSimYearSub[1:(nyears-1),"Yield"]
    SPSimsMed[s,1]<-median(TotBioSimYearSub[c((nrow(TotBioSimYearSub)-wantyears):nrow(TotBioSimYearSub)),"SP"],na.rm=T) #median among year for each sim    
    ####End Surp Prod
  } #end sims loop
  
  colnames(SSBSimsMed)<-c("MedSSB","MedSSBrelSSBmsy","MedSSBrelSSBzero","PropSSBrelSSBmsy","PropSSBrelhalfSSBmsy","PropSSBrel3SSBzero","PropSSBrel75SSBzero")
  colnames(TrueFSimsMed)<-c("PropFrelFmsy")
  colnames(TotYieldSimMed)<-c("Yield","NatDeaths","YieldrelMSY","YieldrelNatDeath","PropClosure","Yvar")
  colnames(NAASimMed)<-c(paste("MedAge",seq(1,nages,1),sep=""),"TotN","PropOne","PropTwo","PropOneTwo")
  colnames(SPSimsMed)<-c("MedSP")
  ###medians###
  for(a in 1:nages){ #loop to calculate median among sims of NAge for each age and then gets added to result table
    if(a==1) { NAAMedRes<-median(NAASimMed[,paste("MedAge",a,sep="")])
    } else { 
      NAAMedRes<-data.frame(NAAMedRes,median(NAASimMed[,paste("MedAge",a,sep="")]))      
    } #end else
  } #end for
  Res<-data.frame(median(SSBSimsMed[,"MedSSB"],na.rm=T),median(SSBSimsMed[,"MedSSBrelSSBmsy"],na.rm=T),median(SSBSimsMed[,"MedSSBrelSSBzero"],na.rm=T),median(SSBSimsMed[,"PropSSBrelSSBmsy"],na.rm=T),median(SSBSimsMed[,"PropSSBrelhalfSSBmsy"],na.rm=T),median(SSBSimsMed[,"PropSSBrel3SSBzero"],na.rm=T),median(SSBSimsMed[,"PropSSBrel75SSBzero"],na.rm=T),median(TrueFSimsMed[,"PropFrelFmsy"],na.rm=T),median(TotYieldSimMed[,"Yield"],na.rm=T),median(TotYieldSimMed[,"Yvar"],na.rm=T),median(TotYieldSimMed[,"NatDeaths"],na.rm=T),median(TotYieldSimMed[,"YieldrelMSY"],na.rm=T),median(TotYieldSimMed[,"YieldrelNatDeath"],na.rm=T),median(TotYieldSimMed[,"PropClosure"],na.rm=T),NAAMedRes,median(NAASimMed[,"PropOne"],na.rm=T),median(NAASimMed[,"PropTwo"],na.rm=T),median(NAASimMed[,"PropOneTwo"],na.rm=T),median(SPSimsMed[,"MedSP"],na.rm=T),CtrlRule,FracBmsyThreshHi,FracBmsyThreshLo,FracFtarg,useminmax,minquota,maxquota,usequotavar,quotavar,Thresh,Prop,Propc,useFcap,Fcapprop,quotablock) #SSB resuls with CR characteristics
  names(Res)<-c("MedianSSB","MedSSBrelSSBmsy","MedSSBrelSSBzero","PropSSBrelSSBmsy","PropSSBrelhalfSSBmsy","PropSSBrel3SSBzero","PropSSBrel75SSBzero","PropFrelFmsy","Yield","Yvar","NatDeaths","YieldrelMSY","YieldrelNatDeath","PropClosure",paste("NAge",seq(1,nages,1),sep=""),"PropAgeOne","PropAgeTwo","PropAgeOneTwo","SurpProd","CtrlRule","FracBmsyThreshHi","FracBmsyThreshLo","FracFtarg","useminmax","minquota","maxquota","usequotavar","quotavar","Thresh","Prop","Propc","useFcap","Fcapprop","quotablock")
  if(CtrlRule==1 & f==1 & u==1 & l==1 || CtrlRule==3 & c==1){ #write the summary results
    write.table(Res,paste(readfile,"SummaryResults.txt",sep=""),append=FALSE,row.names=FALSE,col.names=TRUE,quote=TRUE)   
  } else {
    write.table(Res,paste(readfile,"SummaryResults.txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE)   
  }  
  ###25th quartiles
  for(a in 1:nages){ #loop to calculate median among sims of NAge for each age and then gets added to result table
    if(a==1) { NAA25Res<-quartilefun(data=NAASimMed[,paste("MedAge",a,sep="")],quart=0.25)
    } else { 
      NAA25Res<-data.frame(NAA25Res,quartilefun(data=NAASimMed[,paste("MedAge",a,sep="")],quart=0.25))      
    } #end else
  } #end for
  Res25<-data.frame(quartilefun(data=SSBSimsMed[,"MedSSB"],quart=0.25),quartilefun(data=SSBSimsMed[,"MedSSBrelSSBmsy"],quart=0.25),quartilefun(data=SSBSimsMed[,"MedSSBrelSSBzero"],quart=0.25),quartilefun(data=SSBSimsMed[,"PropSSBrelSSBmsy"],quart=0.25),quartilefun(data=SSBSimsMed[,"PropSSBrelhalfSSBmsy"],quart=0.25),quartilefun(data=SSBSimsMed[,"PropSSBrel3SSBzero"],quart=0.25),quartilefun(data=SSBSimsMed[,"PropSSBrel75SSBzero"],quart=0.25),quartilefun(data=TrueFSimsMed[,"PropFrelFmsy"],quart=0.25),quartilefun(data=TotYieldSimMed[,"Yield"],quart=0.25),quartilefun(data=TotYieldSimMed[,"Yvar"],quart=0.25),quartilefun(data=TotYieldSimMed[,"NatDeaths"],quart=0.25),quartilefun(data=TotYieldSimMed[,"YieldrelMSY"],quart=0.25),quartilefun(data=TotYieldSimMed[,"YieldrelNatDeath"],quart=0.25),quartilefun(data=TotYieldSimMed[,"PropClosure"],quart=0.25),NAA25Res,quartilefun(data=NAASimMed[,"PropOne"],quart=0.25),quartilefun(data=NAASimMed[,"PropTwo"],quart=0.25),quartilefun(data=NAASimMed[,"PropOneTwo"],quart=0.25),quartilefun(data=SPSimsMed[,"MedSP"],quart=0.25),CtrlRule,FracBmsyThreshHi,FracBmsyThreshLo,FracFtarg,useminmax,minquota,maxquota,usequotavar,quotavar,Thresh,Prop,Propc,useFcap,Fcapprop,quotablock) #SSB resuls with CR characteristics
  names(Res25)<-c("MedianSSB","MedSSBrelSSBmsy","MedSSBrelSSBzero","PropSSBrelSSBmsy","PropSSBrelhalfSSBmsy","PropSSBrel3SSBzero","PropSSBrel75SSBzero","PropFrelFmsy","Yield","Yvar","NatDeaths","YieldrelMSY","YieldrelNatDeath","PropClosure",paste("NAge",seq(1,nages,1),sep=""),"PropAgeOne","PropAgeTwo","PropAgeOneTwo","SurpProd","CtrlRule","FracBmsyThreshHi","FracBmsyThreshLo","FracFtarg","useminmax","minquota","maxquota","usequotavar","quotavar","Thresh","Prop","Propc","useFcap","Fcapprop","quotablock")  
  if(CtrlRule==1 & f==1 & u==1 & l==1 || CtrlRule==3 & c==1){ #write the summary results
    write.table(Res25,paste(readfile,"SummaryResults25.txt",sep=""),append=FALSE,row.names=FALSE,col.names=TRUE,quote=TRUE)   
  } else {
    write.table(Res25,paste(readfile,"SummaryResults25.txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE)   
  }  
  ###75th quartiles
  for(a in 1:nages){ #loop to calculate median among sims of NAge for each age and then gets added to result table
    if(a==1) { NAA75Res<-quartilefun(data=NAASimMed[,paste("MedAge",a,sep="")],quart=0.75)
    } else { 
      NAA75Res<-data.frame(NAA75Res,quartilefun(data=NAASimMed[,paste("MedAge",a,sep="")],quart=0.75))      
    } #end else
  } #end for
  Res75<-data.frame(quartilefun(data=SSBSimsMed[,"MedSSB"],quart=0.75),quartilefun(data=SSBSimsMed[,"MedSSBrelSSBmsy"],quart=0.75),quartilefun(data=SSBSimsMed[,"MedSSBrelSSBzero"],quart=0.75),quartilefun(data=SSBSimsMed[,"PropSSBrelSSBmsy"],quart=0.75),quartilefun(data=SSBSimsMed[,"PropSSBrelhalfSSBmsy"],quart=0.75),quartilefun(data=SSBSimsMed[,"PropSSBrel3SSBzero"],quart=0.75),quartilefun(data=SSBSimsMed[,"PropSSBrel75SSBzero"],quart=0.75),quartilefun(data=TrueFSimsMed[,"PropFrelFmsy"],quart=0.75),quartilefun(data=TotYieldSimMed[,"Yield"],quart=0.75),quartilefun(data=TotYieldSimMed[,"Yvar"],quart=0.75),quartilefun(data=TotYieldSimMed[,"NatDeaths"],quart=0.75),quartilefun(data=TotYieldSimMed[,"YieldrelMSY"],quart=0.75),quartilefun(data=TotYieldSimMed[,"YieldrelNatDeath"],quart=0.75),quartilefun(data=TotYieldSimMed[,"PropClosure"],quart=0.75),NAA75Res,quartilefun(data=NAASimMed[,"PropOne"],quart=0.75),quartilefun(data=NAASimMed[,"PropTwo"],quart=0.75),quartilefun(data=NAASimMed[,"PropOneTwo"],quart=0.75),quartilefun(data=SPSimsMed[,"MedSP"],quart=0.75),CtrlRule,FracBmsyThreshHi,FracBmsyThreshLo,FracFtarg,useminmax,minquota,maxquota,usequotavar,quotavar,Thresh,Prop,Propc,useFcap,Fcapprop,quotablock) #SSB resuls with CR characteristics
  names(Res75)<-c("MedianSSB","MedSSBrelSSBmsy","MedSSBrelSSBzero","PropSSBrelSSBmsy","PropSSBrelhalfSSBmsy","PropSSBrel3SSBzero","PropSSBrel75SSBzero","PropFrelFmsy","Yield","Yvar","NatDeaths","YieldrelMSY","YieldrelNatDeath","PropClosure",paste("NAge",seq(1,nages,1),sep=""),"PropAgeOne","PropAgeTwo","PropAgeOneTwo","SurpProd","CtrlRule","FracBmsyThreshHi","FracBmsyThreshLo","FracFtarg","useminmax","minquota","maxquota","usequotavar","quotavar","Thresh","Prop","Propc","useFcap","Fcapprop","quotablock")  
  if(CtrlRule==1 & f==1 & u==1 & l==1 || CtrlRule==3 & c==1){ #write the summary results
    write.table(Res75,paste(readfile,"SummaryResults75.txt",sep=""),append=FALSE,row.names=FALSE,col.names=TRUE,quote=TRUE)   
  } else {
    write.table(Res75,paste(readfile,"SummaryResults75.txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE)   
  }  
} #end ResFun


if(CtrlRule==1){ #biomass based
  for(f in 1:(length(FracFtarg))) { #max f for biomass based ctrl rule loop
    for(u in 1:(length(FracBmsyThreshHi))){ #upper threshold loop
      for(l in 1:u){ #lower threshold loop
        c<-ifelse(l<=length(Propc),l,length(Propc)) #c is only for constant catch but function expects something and so this isn't important; just filler
        ResFun(u=u,l=l,f=f,c=c,CtrlRule=CtrlRule,FracBmsyThreshHi=FracBmsyThreshHi[u],FracBmsyThreshLo=FracBmsyThreshLo[l],FracFtarg=FracFtarg[f],useminmax=useminmax,minquota=minquota,maxquota=maxquota,usequotavar=usequotavar,quotavar=quotavar,Thresh=Thresh,Prop=Prop,Propc=Propc[c],useFcap=useFcap,Fcapprop=Fcapprop,quotablock=quotablock,nsims=nsims,wantyears=(wantyears-1))
      } #end l loop
    } #end u loop
  } #end f loop
} else if(CtrlRule==3){ #constant catch
  for(c in 1:(length(Propc))) {
    u<-ifelse(c<=length(FracBmsyThreshHi),c,length(FracBmsyThreshHi))
    l<-ifelse(c<=length(FracBmsyThreshLo),c,length(FracBmsyThreshLo))
    f<-ifelse(c<=length(FracBmsyThreshLo),c,length(FracFtarg))  
    ResFun(u=u,l=l,f=f,c=c,CtrlRule=CtrlRule,FracBmsyThreshHi=FracBmsyThreshHi[u],FracBmsyThreshLo=FracBmsyThreshLo[l],FracFtarg=FracFtarg[f],useminmax=useminmax,minquota=minquota,maxquota=maxquota,usequotavar=usequotavar,quotavar=quotavar,Thresh=Thresh,Prop=Prop,Propc=Propc[c],useFcap=useFcap,Fcapprop=Fcapprop,quotablock=quotablock,nsims=nsims,wantyears=(wantyears-1))
  } # end c loop
} else { print ("Poop: unknown control rule") }

if(FALSE){
  resdirectplay<-"C:\\Herring\\MSE ABC Control Rule\\MSE\\"  #location of functions.r and userinputs.r
  resname.a<-"HiM_LowSteep_NoAssBias_RecWt\\BB3yr\\"
  resname.b<-"HiM_LowSteep_NoAssBias_RecWt\\BB3yrPerc\\"
  resname.c<-"HiM_LowSteep_NoAssBias_RecWt\\CCC\\"
  resfile<-"UnadjSummaryResults.txt"
  res.a<-read.table(paste(paste(resdirectplay,resname.a,sep=""),resfile,sep=""),header=T)
  res.b<-read.table(paste(paste(resdirectplay,resname.b,sep=""),resfile,sep=""),header=T)
  res.c<-read.table(paste(paste(resdirectplay,resname.c,sep=""),resfile,sep=""),header=T)
  #head(res.a)
  
  #install.packages('rggobi')
  library(rggobi)
  #res.a2<-res.a[res.a$YieldrelMSY>0.85,]
  #ggobi(res.a2)
  
  ggobi(res.a)
  a<-ggobi_get()$res.a
  glyph_colour(a)<-ifelse(a$Yvar<0.5 & a$YieldrelMSY>0.9 & a$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(a)<-ifelse(a$Yvar<0.5 & a$YieldrelMSY>0.9 & a$MedSSBrelSSBzero>0.5,6,2)
  
  ggobi(res.b)
  b<-ggobi_get()$res.b
  glyph_colour(b)<-ifelse(b$Yvar<0.5 & b$YieldrelMSY>0.8 & b$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(b)<-ifelse(b$Yvar<0.5 & b$YieldrelMSY>0.8 & b$MedSSBrelSSBzero>0.5,6,2)
  
  ggobi(res.c)
  c<-ggobi_get()$res.c
  glyph_colour(c)<-ifelse(c$Yvar<0.5 & c$YieldrelMSY>0.9 & c$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(c)<-ifelse(c$Yvar<0.5 & c$YieldrelMSY>0.9 & c$MedSSBrelSSBzero>0.5,6,2)
  
  comp<-ggobi()  
  comp["BB3"]<-res.a
  comp$BB3P<-res.b
  glyph_colour(comp$BB3)<-ifelse(comp$BB3$Yvar<0.5 & comp$BB3$YieldrelMSY>0.85 & comp$BB3$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(comp$BB3)<-ifelse(comp$BB3$Yvar<0.5 & comp$BB3$YieldrelMSY>0.85 & comp$BB3$MedSSBrelSSBzero>0.5,3,1)
  glyph_colour(comp$BB3P)<-ifelse(comp$BB3P$Yvar<0.5 & comp$BB3P$YieldrelMSY>0.85 & comp$BB3P$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(comp$BB3P)<-ifelse(comp$BB3P$Yvar<0.5 & comp$BB3P$YieldrelMSY>0.85 & comp$BB3P$MedSSBrelSSBzero>0.5,3,1)
  
}