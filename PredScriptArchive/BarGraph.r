folder<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\"
folder1<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\BioBased_A\\"
folder2<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\BioBased_B\\"
folder3<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\BioBased_C\\"
folder4<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\CF_75Fmsy\\"
folder5<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\PropThresh_A\\"
folder6<-"C:\\NEFMC ABC WG\\MSE\\RetroManuCode\\ResultsWithTimeM_BiasLinkM\\TrueLenfest\\"

BioARes<-read.table(paste(folder1,"UnadjResults.txt",sep=""),header=T)
BioBRes<-read.table(paste(folder2,"UnadjResults.txt",sep=""),header=T)
BioCRes<-read.table(paste(folder3,"UnadjResults.txt",sep=""),header=T)
CF75Res<-read.table(paste(folder4,"UnadjResults.txt",sep=""),header=T)
PropT_ARes<-read.table(paste(folder5,"UnadjResults.txt",sep=""),header=T)
TrueLenRes<-read.table(paste(folder6,"UnadjResults.txt",sep=""),header=T)

SSBYYvarmed<-rbind(BioARes['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioBRes['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioCRes['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],CF75Res['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],
PropT_ARes['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],TrueLenRes['Median',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')]) ##median results
rownames(SSBYYvarmed)<-c("BioBasedA","BioBasedB","BioBasedC","75%Fmsy","PropThreshA","Lenfest")
SSBYYvar25<-rbind(BioARes['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioBRes['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioCRes['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],CF75Res['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],
PropT_ARes['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],TrueLenRes['25',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')])
rownames(SSBYYvar25)<-c("BioBasedA","BioBasedB","BioBasedC","75%Fmsy","PropThreshA","Lenfest")
SSBYYvar75<-rbind(BioARes['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioBRes['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],BioCRes['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],CF75Res['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],
PropT_ARes['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')],TrueLenRes['75',c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')])
rownames(SSBYYvar75)<-c("BioBasedA","BioBasedB","BioBasedC","75%Fmsy","PropThreshA","Lenfest")

windows(width=18,height=13) #create big plot window
ResPlot<-barplot(as.matrix(SSBYYvarmed[,]),beside=TRUE,names.arg=c("SSB/UnfishedSSB","SSB/SSBmsy","Years<0.4UnfishedSSB","Y/MSY","IAV","Years Closed"),lwd=2,cex.axis=1.2,cex.names=1.2,border=FALSE,ylab="Relative Results (%)",cex.lab=1.2,ylim=c(0,max(SSBYYvar75[,c('SSBRelBzero','SSBRelBmsy','YieldRelMSY','Yvar','FreqFisheryClosed')])),
   col=c("deepskyblue","darkgoldenrod2","yellow","darkgray","green","darkmagenta"),legend = rownames(SSBYYvarmed))   #make a barplot ,angle=c(30,60),density=c(4,4)
abline(h=0,lwd=2)
abline(a=25,b=0,lwd=2,col="gray",lty=2)
abline(a=50,b=0,lwd=2,col="gray",lty=2)
abline(a=75,b=0,lwd=2,col="gray",lty=2)
abline(a=100,b=0,lwd=2,col="gray",lty=2)
segments(ResPlot,as.matrix(SSBYYvar25[,c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')]),ResPlot,as.matrix(SSBYYvar75[,c('SSBRelBzero','SSBRelBmsy','X.OverfishedSSB','YieldRelMSY','Yvar','FreqFisheryClosed')]),lwd=2)  #error bars based on interquartiles

savePlot(paste(folder,"Results_barplot.tiff",sep=""),type="tiff")
graphics.off()



if(FALSE) #does multiple plots in same window; good for different scales in plots
{
dev.new(width=16,height=13) #create big plot window
resplot<-layout(matrix(c(1,2),1,2,byrow=F),respect=F)  #setup plot window matrix for plotting SSB and yield in one, and Yvar in the next
#SSB and yield plot
SSBYPlot<-barplot(as.matrix(SSBYYvarmed[,c('SSB','Yield')]),beside=TRUE,lwd=2,cex.axis=1.2,cex.names=1.2,border=FALSE,ylab="Metric Tons",cex.lab=1.2,ylim=c(0,max(SSBYYvar75[,c('SSB','Yield')])),
   col=c("deepskyblue","darkgoldenrod2","yellow"),legend = rownames(SSBYYvarmed))   #make a barplot ,angle=c(30,60),density=c(4,4)
abline(h=0,lwd=2)
segments(SSBYPlot,as.matrix(SSBYYvar25[,c('SSB','Yield')]),SSBYPlot,as.matrix(SSBYYvar75[,c('SSB','Yield')]),lwd=2)  #error bars based on interquartiles
#Yvar plot
IAVdat<-as.matrix(SSBYYvarmed[,c('Yvar')])   #for unknown R reasons, next two lines required for R to auto make xaxis label as it does above without these 2 lines
colnames(IAVdat)<-c("IAV")    #see comment 1 line up.
YvarPlot<-barplot(IAVdat,beside=TRUE,lwd=2,cex.axis=1.2,cex.names=1.2,border=FALSE,cex.lab=1.2,ylab="",ylim=c(0,max(SSBYYvar75[,c('Yvar')])),
   col=c("deepskyblue","darkgoldenrod2","yellow"))   #make a barplot ,angle=c(30,60),density=c(4,4)
abline(h=0,lwd=2)
segments(YvarPlot,as.matrix(SSBYYvar25[,c('Yvar')]),YvarPlot,as.matrix(SSBYYvar75[,c('Yvar')]),lwd=2)  #error bars based on interquartiles

savePlot(paste(folder,"SSB_Y_Yvar_barplot.tiff",sep=""),type="tiff")
} #end if false


