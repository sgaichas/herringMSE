  rm(list=ls(all=TRUE))
  #direct<-"C:\\Herring\\MSE ABC Control Rule\\MSE\\"  #location of functions.r and userinputs.r
  direct<-"F:\\Herring\\MSE ABC Control Rule\\MSE\\" 
  #resdirectpred<-"F:\\Herring\\MSE ABC Control Rule\\MSE\\TunaTest\\TunaResultsSummaries\\ResultsSummaries\\"
  resfile<-"UnadjSummaryResults.txt" #the summary results file name (medians of metrics)
  
  #OMnames<-c("HiM_LowSteep_NoAssBias_OldWt","HiM_LowSteep_NoAssBias_RecWt","LoM_HiSteep_NoAssBias_OldWt","LoM_HiSteep_NoAssBias_RecWt","HiM_LowSteep_AssBias_OldWt","HiM_LowSteep_AssBias_RecWt","LoM_HiSteep_AssBias_OldWt","LoM_HiSteep_AssBias_RecWt") #names of operating model file names
  OMnames<-c("HiM_LowSteep_NoAssBias_OldWt","HiM_LowSteep_NoAssBias_RecWt","LoM_HiSteep_NoAssBias_OldWt","LoM_HiSteep_NoAssBias_RecWt") #names of operating model file names
  #OMnames<-c("HiM_LowSteep_AssBias_OldWt","HiM_LowSteep_AssBias_RecWt","LoM_HiSteep_AssBias_OldWt","LoM_HiSteep_AssBias_RecWt") #names of operating model file names 
  assign("OMnames", OMnames, envir = .GlobalEnv)
  ctrlrulenames<-c("BB","BB3yr","BB5yr","BB3yrPerc","CC","CCC") #names of ctrl rule file names
  assign("ctrlrulenames", ctrlrulenames, envir = .GlobalEnv)
  OMwant<-2 #position in OMnames of desired OM for which detailed CR plots will be produced
  assign("OMwant", OMwant, envir = .GlobalEnv)
  
  #function to go to a location and read MSE results for given OM and CR
  getres<-function(locale=NULL){
    res<-read.table(locale,header=T)
    return(res)
  }
 
  source(paste(direct,"geteconfxn.R",sep=""))  ##get econ results
  source(paste(direct,"getpredatorfxn.R",sep=""))  ##get econ results
  
  yaxiswant<-c("YieldrelMSY","YieldrelMSY","YieldrelMSY","PropSSBrelhalfSSBmsy","Yield","YieldrelMSY", "AvPropYrs_okBstatusgf","AvPropYrs_goodProd_Targplus", "AvPropYrs_goodAvWtstatus" ) #for tradeoff plots; desired y-axis; must be same size as xaxiswant
  xaxiswant<-c("MedSSBrelSSBzero","Yvar","PropSSBrelhalfSSBmsy","Yvar","MedianSSB","PropSSBrel3SSBzero","MedSSBrelSSBzero","MedSSBrelSSBzero","MedSSBrelSSBzero") #as with yaxiswant, but for xaxis.
  yaxisnames<-c("Yield/MSY","Yield/MSY","Yield/MSY","Prob Overfished","Yield","Yield/MSY","Prob GF >0.5Bmsy","Prob Tern Prod >1","Prob Wt > Avg")
  xaxisnames<-c("SSB/Unfished","Variation in Yield","Prob Overfished","Variation in Yield","SSB","Prob SSB < 30%Unfished","SSB/Unfished","SSB/Unfished","SSB/Unfished") 
  yaxisval<-c(0.8,0.8,0.8,0.0,100000,0.8,1.0,0.9,1.0)
  xaxisval<-c(0.60,0.35,0.0,0.35,300000,0.0,0.60,0.60,0.60)
  
  #yaxiswant<-c("Yield","MedPredNtern","MedPredN_Curstatus") #for tradeoff plots; desired y-axis; must be same size as xaxiswant
  #xaxiswant<-c("MedianSSB","MedianSSB","MedSSBrelSSBzero") #as with yaxiswant, but for xaxis.
  #yaxisnames<-c("Yield","MedPredNTern","Tern CurStatus")
  #xaxisnames<-c("MedianSSB","MedianSSB","SSB/Unfished") 
  
  #function to make boxplot
  boxplotfun<-function(results=NULL,reswanted=NULL,marwant=NULL,aba=NULL,abb=NULL,mainname=NULL,ylims=NULL,cexaxis=NULL) {
  par(mar=marwant) #set desired margins
  boxplot(data=results,results[,names(results) %in% reswanted],las=2.2,main=mainname,ylim=ylims,cex.axis=cexaxis)
  abline(a=aba,b=abb,lty=2)
  }
  
  
  #loop to read in results that occur in given combo of OM and CR and place in single dataframe
  for(OM in 1:length(OMnames)) {
  for(cr in 1:length(ctrlrulenames)) {
    name<-paste(paste(paste(direct,OMnames[OM],sep=""),ctrlrulenames[cr],sep="\\"),resfile,sep="\\")
    tempres<-getres(locale=name)
    tempres$OM<-OMnames[OM]
    tempres$CR<-ctrlrulenames[cr]
    if(OM==1 & cr==1) {
      allres<-tempres
    } else {
      allres<-rbind(allres,tempres)
    }
  }
  }
                              
  #merge econ and predator results with all results
    #library(dplyr)
    #anti_join(poopb,poop)
  allres<-merge(allres,allecon)
  allres<-merge(allres,allpredres)
  #saveRDS(allres,file=paste0(direct,"allres.rds"))
  assign("allres", allres, envir = .GlobalEnv)
#######################################################
  #plot results across all OMs and CRs and see what varies
  windows(height=10,width=10)  
  boxplotfun(results=allres,reswanted=c("Yvar","YieldrelNatDeath","YieldrelMSY","MedSSBrelSSBzero","MedPredAvWt_status"),marwant=c(10.5,4.1,.2,2.1),aba=1,abb=0,mainname="",ylims=c(0,2),cexaxis=1.2)
  #savePlot(filename=paste(direct,paste("Plots","AllRes_RelResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllRes_RelResults_AssBias.png",sep="\\"),sep=""),type="png") 
  windows(height=10,width=10)  
  boxplotfun(results=allres,reswanted=c("PropSSBrelhalfSSBmsy","PropFrelFmsy","PropClosure","PropAgeOneTwo"),marwant=c(12.5,4.1,.2,2.1),aba=0,abb=0,mainname="",ylims=c(0,1),cexaxis=1.2)
  #savePlot(filename=paste(direct,paste("Plots","AllRes_PropResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllRes_PropResults_AssBias.png",sep="\\"),sep=""),type="png") 
  windows(height=10,width=10)  
  boxplotfun(results=allres,reswanted=c("MedianSSB","Yield","NatDeaths"),marwant=c(8,6,.2,2.1),aba=0,abb=0,mainname="",ylims=c(0,1500000),cexaxis=1.2)
  #savePlot(filename=paste(direct,paste("Plots","AllRes_AbsResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllRes_AbsResults_AssBias.png",sep="\\"),sep=""),type="png") 
  
  graphics.off()
  
  #plot results by OM among all CRs
  #layout.show(graph)
  for(OM. in 1:length(OMnames)){
    allres.sub<-allres[allres$OM==OMnames[OM.],]
    if(OM. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,4,1),1,4,byrow=F),respect=F)   #creates graphics "matrix"
    } else {
      dev.set(which=devices[1])
    }
      boxplotfun(results=allres.sub,reswanted=c("Yvar","YieldrelNatDeath","YieldrelMSY","MedSSBrelSSBzero","MedPredAvWt_status"),marwant=c(12,4.1,2,2.1),aba=1,abb=0,mainname=OMnames[OM.],ylims=c(0,2),cexaxis=1.4)
    
    if(OM. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,4,1),1,4,byrow=F),respect=F)   #creates graphics "matrix"
    } else {
      dev.set(which=devices[2])
    }
      boxplotfun(results=allres.sub,reswanted=c("PropSSBrelhalfSSBmsy","PropFrelFmsy","PropClosure","PropAgeOneTwo"),marwant=c(13.5,4.1,2,2.1),aba=0,abb=0,mainname=OMnames[OM.],ylims=c(0,1),cexaxis=1.4)

    if(OM. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,4,1),1,4,byrow=F),respect=F)   #creates graphics "matrix"
      devices<-dev.list()
    } else {      
      dev.set(which=devices[3])
    }
      boxplotfun(results=allres.sub,reswanted=c("MedianSSB","Yield","NatDeaths"),marwant=c(10,6,2,2.1),aba=0,abb=0,mainname=OMnames[OM.],ylims=c(0,1500000),cexaxis=1.4)
  }
  dev.set(which=devices[1])
  #savePlot(filename=paste(direct,paste("Plots","AllResByOM_RelResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllResByOM_RelResults_AssBias.png",sep="\\"),sep=""),type="png") 
  dev.set(which=devices[2])
  #savePlot(filename=paste(direct,paste("Plots","AllResByOM_PropResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllResByOM_PropResults_AssBias.png",sep="\\"),sep=""),type="png") 
  dev.set(which=devices[3])
  #savePlot(filename=paste(direct,paste("Plots","AllResByOM_AbsResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllResByOM_AbsResults_AssBias.png",sep="\\"),sep=""),type="png") 
  
  
graphics.off()  
  
  #plot results by CR among all OMs
  #layout.show(graph)
  for(CR. in 1:length(ctrlrulenames)){
    allres.sub<-allres[allres$CR==ctrlrulenames[CR.],]
    if(CR. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,6,1),1,6,byrow=F),respect=F)   #creates graphics "matrix"
    } else {
      dev.set(which=devices[1])
    }
    boxplotfun(results=allres.sub,reswanted=c("Yvar","YieldrelNatDeath","YieldrelMSY","MedSSBrelSSBzero","MedPredAvWt_status"),marwant=c(12,4.1,2,2.1),aba=1,abb=0,mainname=ctrlrulenames[CR.],ylims=c(0,2),cexaxis=1.4)
    
    if(CR. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,6,1),1,6,byrow=F),respect=F)   #creates graphics "matrix"
    } else {
      dev.set(which=devices[2])
    }
    boxplotfun(results=allres.sub,reswanted=c("PropSSBrelhalfSSBmsy","PropFrelFmsy","PropClosure","PropAgeOneTwo"),marwant=c(13.5,4.1,2,2.1),aba=0,abb=0,mainname=ctrlrulenames[CR.],ylims=c(0,1),cexaxis=1.4)
    
    if(CR. == 1) { 
      windows(height=8,width=15)
      graph<-layout(matrix(seq(1,6,1),1,6,byrow=F),respect=F)   #creates graphics "matrix"
      devices<-dev.list()
    } else {      
      dev.set(which=devices[3])
    }
    boxplotfun(results=allres.sub,reswanted=c("MedianSSB","Yield","NatDeaths"),marwant=c(10,6,2,2.1),aba=0,abb=0,mainname=ctrlrulenames[CR.],ylims=c(0,1500000),cexaxis=1.4)
  }
  dev.set(which=devices[1])
  #savePlot(filename=paste(direct,paste("Plots","AllResByCR_RelResults_NoAssBias.png",sep="\\"),sep=""),type="png")
  savePlot(filename=paste(direct,paste("Plots","AllResByCR_RelResults_AssBias.png",sep="\\"),sep=""),type="png") 
  dev.set(which=devices[2])
  #savePlot(filename=paste(direct,paste("Plots","AllResByCR_PropResults_NoAssBias.png",sep="\\"),sep=""),type="png") 
  savePlot(filename=paste(direct,paste("Plots","AllResByCR_PropResults_AssBias.png",sep="\\"),sep=""),type="png") 
  dev.set(which=devices[3])
  #savePlot(filename=paste(direct,paste("Plots","AllResByCR_AbsResults_NoAssBias.png",sep="\\"),sep=""),type="png")
  savePlot(filename=paste(direct,paste("Plots","AllResByCR_AbsResults_AssBias.png",sep="\\"),sep=""),type="png") 
  
  graphics.off()  
###########################################################  
  
 #tradeoff plot functions
  tradeoffplot<-function(xdat=NULL,ydat=NULL,typeb=NULL,ylimb=NULL,xlimb=NULL,pchb=NULL,xlabb=NULL,ylabb=NULL,cex.labb=NULL,cex.axisb=NULL,colb=NULL){
    plot(xdat,ydat,type=typeb,ylim=ylimb,xlim=xlimb,pch=pchb,xlab=xlabb,ylab=ylabb,cex.lab=cex.labb,cex.axis=cex.axisb,col=colb)
  }
  tradeofflines<-function(xdat=NULL,ydat=NULL,typeb=NULL,pchb=NULL,cexb=NULL,colb=NULL){
    lines(xdat,ydat,type=typeb,pch=pchb,col=colb,cex=cexb)
  }

#function to plot specified CR for each OM to view differences among tradeoffs caused by OM  
CR.OM<-function(xwantb=NULL,ywantb=NULL,xnameb=NULL,ynameb=NULL,ctrlwant=NULL,Hiwant=NULL,Lowant=NULL,Fwant=NULL,Propcwant=NULL){
  #windows(height=10,width=10)
  for(OM. in 1:length(OMnames)){
   if(ctrlwant<5) { subres<-allres[allres$FracBmsyThreshHi==Hiwant & allres$FracBmsyThreshLo==Lowant & allres$FracFtarg==Fwant & allres$OM==OMnames[OM.] & allres$CR==ctrlrulenames[ctrlwant],] 
                    #print(OMnames[OM.])
                    #print(subres)
                    
   } else {
     subres<-allres[allres$Propc==Propcwant & allres$OM==OMnames[OM.] & allres$CR==ctrlrulenames[ctrlwant],] 
   }   
   if(OM.==1){   
   tradeoffplot(xdat=as.numeric(subres[xwantb]),ydat=as.numeric(subres[ywantb]),typeb="p",ylimb=c(0,110000),xlimb=c(0,500000),pchb=19,xlabb=xnameb,ylabb=ynameb,cex.labb=1.5,cex.axisb=1.3,colb="black")
   legend
   } else {
     col<-OM.
   tradeofflines(xdat=as.numeric(subres[xwantb]),ydat=as.numeric(subres[ywantb]),typeb="p",pchb=19,cexb=1,colb=col) 
   }
  }
  legend("topright",legend=OMnames,pch=19,col=seq(1:length(OMnames)))
}


    #########Detailed plots for just one OM
TradeOffPlots<-function(xwantb=NULL,ywantb=NULL,xnameb=NULL,ynameb=NULL,wantcolor=NULL,MSYwant=NULL,SSBwant=NULL){  
  Statusquo<-allres[allres$FracBmsyThreshHi==0.5 & allres$FracBmsyThreshLo==0.0 & allres$FracFtarg==0.9 & allres$OM==OMnames[OMwant] & allres$CR==ctrlrulenames[2],]
  lessthans<-c("Yvar","PropSSBrelhalfSSBmsy","PropSSBrel3SSBzero")
  if(ywantb=="MedPredNtern"){ylimb<-c(0,50000)}else{ylimb<-c(0,1)} #max(allres$Yvar)
  if(ywantb=="Yield"){ylimb<-c(0,140000)}else{ylimb<-c(0,1)}
  #ylimb<-NULL
  #ylimb<-c(0,1)
  
  if(xwantb=="nr" | xwantb=="gr" | xwantb=="MedianSSB" | xwantb=="MedSSBrelSSBmsy"){xlimb<-NULL}else{xlimb<-c(0,1)}
  for(CR. in 1:(length(ctrlrulenames)-1)){
    if(CR. < (length(ctrlrulenames)-1)) {      
      allres.sub<-allres[allres$OM==OMnames[OMwant] & allres$CR==ctrlrulenames[CR.],]      
      if(wantcolor==FALSE){
       colors<-"black"
       
      } else { colors<-ifelse(allres.sub$YieldrelMSY >= max(MSYwant) & allres.sub$MedSSBrelSSBzero >= SSBwant,"green1",ifelse(allres.sub$YieldrelMSY >= min(MSYwant) & allres.sub$YieldrelMSY < max(MSYwant) & allres.sub$MedSSBrelSSBzero >=SSBwant,"skyblue","black"))
        
      }
      
      if(CR.==1){
      tradeoffplot(xdat=as.matrix(allres.sub[xwantb]),ydat=as.matrix(allres.sub[ywantb]),typeb="p",ylimb=ylimb,xlimb=xlimb,pchb=19,xlabb=xnameb,ylabb=ynameb,cex.labb=1.5,cex.axisb=1.3,colb=colors)
      } else {tradeoffplot(xdat=as.matrix(allres.sub[xwantb]),ydat=as.matrix(allres.sub[ywantb]),typeb="p",ylimb=ylimb,xlimb=xlimb,pchb=19,xlabb=xnameb,ylabb=" ",cex.labb=1.5,cex.axisb=1.3,colb=colors) }
      title(main=ctrlrulenames[CR.])
      tradeofflines(xdat=as.numeric(Statusquo[xwantb]),ydat=as.numeric(Statusquo[ywantb]),typeb="p",pchb=17,cexb=1.3,colb="red")
      if(ywantb=="Yvar" | ywantb=="YieldrelMSY"){
      axis(side=4,labels=seq(0,100000,20000),at=(seq(0,100000,20000)/100000),cex.axis=1.3)
      } #put absolute scale to variation in yield or yield/MSY   
    } else {
      allres.sub<-allres[allres$OM==OMnames[OMwant] & allres$CR==ctrlrulenames[CR.],]
      if(wantcolor==FALSE){
        colors<-"black"
        
      } else { colors<-ifelse(allres.sub$YieldrelMSY >= max(MSYwant) & allres.sub$MedSSBrelSSBzero >= SSBwant,"green1",ifelse(allres.sub$YieldrelMSY >= min(MSYwant) & allres.sub$YieldrelMSY < max(MSYwant) & allres.sub$MedSSBrelSSBzero >=SSBwant,"skyblue","black"))
        
      }
      
      tradeoffplot(xdat=as.matrix(allres.sub[xwantb]),ydat=as.matrix(allres.sub[ywantb]),typeb="p",ylimb=ylimb,xlimb=xlimb,pchb=19,xlabb=xnameb,ylabb=" ",cex.labb=1.5,cex.axisb=1.3,colb=colors)
      title(main=paste(ctrlrulenames[CR.],ctrlrulenames[CR.+1],sep=" and "))
      allres.sub<-allres[allres$OM==OMnames[OMwant] & allres$CR==ctrlrulenames[CR.+1],]
      if(wantcolor==FALSE){
        colors<-"black"
        
      } else { colors<-ifelse(allres.sub$YieldrelMSY >= max(MSYwant) & allres.sub$MedSSBrelSSBzero >= SSBwant,"green1",ifelse(allres.sub$YieldrelMSY >= min(MSYwant) & allres.sub$YieldrelMSY < max(MSYwant) & allres.sub$MedSSBrelSSBzero >=SSBwant,"skyblue","black"))
        
      }
      
      tradeofflines(xdat=as.matrix(allres.sub[xwantb]),ydat=as.matrix(allres.sub[ywantb]),typeb="p",pchb=23,cexb=1,colb=colors)
      tradeofflines(xdat=as.numeric(Statusquo[xwantb]),ydat=as.numeric(Statusquo[ywantb]),typeb="p",pchb=17,cexb=1.3,colb="red")
    }
  }
}
####################################################
  
  
  ##loop calling to function to make plots of a given HCR for each OM 
  ctrlwant=4 #which control desired; number is position in ctrlrulenames
  Hiwant=0 #If a BB control rule, then what hi threshold as fraction of SSBmsy
  Lowant=0 #If a BB control rule, then what lo threshold as fraction of SSBmsy
  Fwant=0.8 #If a BB control rule, what is max target F as frac of Fmsy
  Propcwant=0.8 #If a CC/CCC control rule then what fraction of MSY  
  windows(height=8,width=15)
  graph<-layout(matrix(seq(1,length(xaxiswant),1),1,length(xaxiswant),byrow=F),respect=F)   #creates graphics "matrix"
  for(i in 1:length(xaxiswant)) {    
    CR.OM(xwantb=xaxiswant[i],ywantb=yaxiswant[i],xnameb=xaxisnames[i],ynameb=yaxisnames[i],ctrlwant=ctrlwant,Hiwant=Hiwant,Lowant=Lowant,Fwant=Fwant,Propcwant=Propcwant)
    CR.OMplotname<-paste(paste(ctrlrulenames[ctrlwant],paste("Hi",paste(Hiwant,paste("Lo",paste(Lowant,paste("Fcap",paste(Fwant,paste("CC",Propcwant,sep=""),sep="_"),sep=""),sep="_"),sep=""),sep="_"),sep=""),sep="_"),".png",sep="")
    #savePlot(filename=paste(direct,paste("Plots",CR.OMplotname,sep="\\"),sep=""),type="png") 
}
  graphics.off()
  
  #loop to make desired tradeoff plots for OMwant  
  for(i in 1:length(xaxiswant)) {
    windows(height=8,width=15)
    graph<-layout(matrix(seq(1,5,1),1,5,byrow=F),respect=F)   #creates graphics "matrix"
    #layout.show(graph)
    TradeOffPlots(xwantb=xaxiswant[i],ywantb=yaxiswant[i],xnameb=xaxisnames[i],ynameb=yaxisnames[i],wantcolor=TRUE,MSYwant=c(0.7,0.8),SSBwant=0.6)
    Tradeplotname<-paste(paste(xaxiswant[i],paste(yaxiswant[i],"_TRADE_FollowY_SSB",sep=""),sep="VS"),OMnames[OMwant],sep="_")
    savePlot(filename=paste(direct,paste("Plots",Tradeplotname,sep="\\"),sep=""),type="png")
  }
  graphics.off()
################################################################

################################################################
###graph CR shape and desired CR tradeoff plots
  CRplotfun<-function(FracHi=NULL,FracLo=NULL,maxFtarg=NULL,newplot=NULL,c=NULL){
    Ftarg<-c()
    SSB<-seq(0,1.5,0.05)
    for(i in 1:(length(SSB))){
    if(SSB[i]>=FracHi){Ftarg[i]<-maxFtarg}else{Ftarg[i]<-(maxFtarg)*((SSB[i]-FracLo)/(FracHi-FracLo))} #set target fully selected F based on a control rule (See Katsukawa 2004 Fisheries Science)
    if (Ftarg[i]<0) Ftarg[i]<-0  #when SSBstatus is less than SSBthreshlo then Ftarg is negative based on control rule.  This fixes that problem.
    }
    SSB<-data.frame(SSB,Ftarg)
    if(!is.null(dev.list())) {par(new=newplot) }
    plot(SSB,type='l',xlab="SSB/SSBmsy",ylab="F/Fmsy",lwd=2,cex.axis=1.3,cex.lab=1.4,ylim=c(0,1),lty=c)
  }
###Desired CR characteristics  
#ForCRplot<-allres[allres$FracBmsyThreshHi==0.5 & allres$FracBmsyThreshLo==0.0 & allres$FracFtarg==0.9 & allres$OM==OMnames[2] & allres$CR==ctrlrulenames[2],]  
  ForCRplot<-allres[allres$MedPredNtern<28000 & allres$OM==OMnames[4] & allres$CR==ctrlrulenames[4],]   
    #allres[allres$MedSSBrelSSBzero>0.5 & allres$YieldrelMSY>0.8 & allres$Yvar>0.88 & allres$OM==OMnames[2] & allres$CR==ctrlrulenames[1],]   
FracHis<-ForCRplot$FracBmsyThreshHi
FracLos<-ForCRplot$FracBmsyThreshLo
maxFs<-ForCRplot$FracFtarg
  
  windows(height=10,width=10)
for(c in 1:(length(maxFs))){
CRplotfun(FracHi=FracHis[c],FracLo=FracLos[c],maxFtarg=maxFs[c],c=c,newplot=TRUE)
}
  #for plotting and saving status quo CR shape
  #lines(x=c(0.0,0.5),y=c(0.9,0.9),lty=2,lwd=2)
  #lines(x=c(0.5,0.5),y=c(0.0,0.9),lty=2,lwd=2)
  #savePlot(filename=paste(direct,paste("Plots","StatusQuoCRshape.png",sep="\\"),sep=""),type="png")
graphics.off()
  
  windows(height=8,width=15)
  graph<-layout(matrix(seq(1,length(xaxiswant),1),1,length(xaxiswant),byrow=F),respect=F)   #creates graphics "matrix"
for(i in 1:length(xaxiswant)) {  
tradeoffplot(xdat=as.matrix(ForCRplot[xaxiswant[i]]),ydat=as.matrix(ForCRplot[yaxiswant[i]]),typeb="p",ylimb=c(0,max(ForCRplot[yaxiswant[i]])),xlimb=c(0,1),pchb=19,xlabb=xaxisnames[i],ylabb=yaxisnames[i],cex.labb=1.5,cex.axisb=1.3,colb="black")
}
graphics.off()  
###################################################################  
if(FALSE)  {
  resname.b<-"HiM_LowSteep_AssBias_RecWt\\CCC\\"
  #resname.c<-"HiM_LowSteep_NoAssBias_RecWt\\CCC\\"
  
  resfile25<-"UnadjSummaryResults25.txt"
  resfile75<-"UnadjSummaryResults75.txt"
  resfilepred<-"HiM_LowSteep_NoAssBias_RecWt_BB3yrPerc_predSum.txt"
  
  res.pred<-read.table(paste(resdirectpred,resfilepred,sep=""),header=T)
  names(res.pred)[names(res.pred) %in% c("minmaxquota.","LimQuotaVar.","AnnQuotaVarFrac","ConCatchPropMSY","CondtlCatch.","MaxF_fracFmsy")]<-c("useminmax","usequotavar","quotavar","Propc","useFcap","Fcapprop")
  
  res.a<-read.table(paste(paste(resdirectplay,resname.a,sep=""),resfile,sep=""),header=T)
  res.a25<-read.table(paste(paste(resdirectplay,resname.a,sep=""),resfile25,sep=""),header=T)
  res.a75<-read.table(paste(paste(resdirectplay,resname.a,sep=""),resfile75,sep=""),header=T)
  res.a$Percentile<-'Median' 
  res.a25$Percentile<-"TwoFive"
  res.a75$Percentile<-"SevFive"
  #res.a<-rbind(res.a,res.a25,res.a75)  
  #head(res.a)
  
  res.poo<-merge(res.a,res.pred)
  res.poo<-res.poo[,order(colnames(res.poo),decreasing=T)]
  
  res.b<-read.table(paste(paste(resdirectplay,resname.b,sep=""),resfile,sep=""),header=T)
  res.b25<-read.table(paste(paste(resdirectplay,resname.b,sep=""),resfile25,sep=""),header=T)
  res.b75<-read.table(paste(paste(resdirectplay,resname.b,sep=""),resfile75,sep=""),header=T)
  res.b$Percentile<-'Median' 
  res.b25$Percentile<-"TwoFive"
  res.b75$Percentile<-"SevFive"
  #res.b<-rbind(res.b,res.b25,res.b75)  
  
  res.a<-res.a[res.a$FracBmsyThreshHi %in% c(0,0.5,0.7,1) & res.a$FracBmsyThreshLo %in% c(0,0.5,0.7,1) & res.a$FracFtarg %in% c(0.2,0.5,0.7),]
  res.b<-res.b[res.b$FracBmsyThreshHi %in% c(0,0.5,0.7,1) & res.b$FracBmsyThreshLo %in% c(0,0.5,0.7,1) & res.b$FracFtarg %in% c(0.2,0.5,0.7),]  
  
  par(mar=c(9.2,4.1,.2,2.1))
  windows(height=10,width=10)  
  graph<-layout(matrix(seq(1,6,1),2,3,byrow=F),respect=F)   #creates graphics "matrix" 
  #layout.show(graph)
  
  
  
  #install.packages('rggobi')
  library(rggobi)
  #res.a2<-res.a[res.a$YieldrelMSY>0.85,]
  #ggobi(res.a2)
  
  ggobi(res.poo)
  
  ggobi(res.a)
  a<-ggobi_get()$res.a
  glyph_colour(a)<-ifelse(a$Percentile=="Median",3,1)
  
  glyph_colour(a)<-ifelse(a$Yvar<0.5 & a$YieldrelMSY>0.9 & a$MedSSBrelSSBzero>0.5,2,1)
  glyph_size(a)<-ifelse(a$Yvar<0.5 & a$YieldrelMSY>0.9 & a$MedSSBrelSSBzero>0.5,6,2)
  
  ggobi(res.b)
  b<-ggobi_get()$res.b
  glyph_colour(b)<-ifelse(b$Percentile=="Median",3,1)
  
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
  
  
  #a<-ggobi(res.a)
  #b<-ggobi(res.b)
  
}