#25th percentile
resfile25<-"UnadjSummaryResults25.txt" #the summary results file name (medians of metrics)

for(OM in 1:length(OMnames)) {
  for(cr in 1:length(ctrlrulenames)) {
    name<-paste(paste(paste(direct,OMnames[OM],sep=""),ctrlrulenames[cr],sep="\\"),resfile25,sep="\\")
    tempres<-getres(locale=name)
    tempres$OM<-OMnames[OM]
    tempres$CR<-ctrlrulenames[cr]
    if(OM==1 & cr==1) {
      allres25<-tempres
    } else {
      allres25<-rbind(allres25,tempres)
    }
  }
}

names(allres25)[1:26]<-paste(names(allres25[1:26]),"25",sep="_")

#75th percentile
resfile75<-"UnadjSummaryResults75.txt" #the summary results file name (medians of metrics)

for(OM in 1:length(OMnames)) {
  for(cr in 1:length(ctrlrulenames)) {
    name<-paste(paste(paste(direct,OMnames[OM],sep=""),ctrlrulenames[cr],sep="\\"),resfile75,sep="\\")
    tempres<-getres(locale=name)
    tempres$OM<-OMnames[OM]
    tempres$CR<-ctrlrulenames[cr]
    if(OM==1 & cr==1) {
      allres75<-tempres
    } else {
      allres75<-rbind(allres75,tempres)
    }
  }
}

names(allres75)[1:26]<-paste(names(allres75[1:26]),"75",sep="_")
