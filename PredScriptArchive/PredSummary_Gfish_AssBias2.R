rm(list=ls())
# summarize predator model results in same format as herring results
# S. Gaichas, herring MSE August 2016

# be in the correct predator OM and CRtype folder 
# find the equivalent prey OM and CRtype folder on network drive
# check that control rule numbers in predator folder match those in SimCharacteristics file in prey folder
# make a file with filename control rule numbers and sim characteristics?

require(dplyr)

wantyears<-50 #number of years to take median over for each sim for results

PredDirec<-c('~/Data/Projects/MSE/HerringMSE/Groundfish')
#PredDirec<-c('~/Documents/0_Data/MSE/HerringMSE/Tuna')

OMDirec<-c('//net/mse/')
#OMtype<-c('HiM_LowSteep_NoAssBias_OldWt') #done first
#OMlist<-c('HiM_LowSteep_NoAssBias_OldWt','LoM_HiSteep_NoAssBias_OldWt', 'LoM_HiSteep_NoAssBias_RecWt', 'HiM_LowSteep_NoAssBias_RecWt')
OMlist<-c('HiM_LowSteep_AssBias_OldWt','LoM_HiSteep_AssBias_OldWt')
#OMlist<-c('LoM_HiSteep_AssBias_RecWt', 'HiM_LowSteep_AssBias_RecWt')
CRtype<-c("BB", "BB3yr", "BB5yr", "BB3yrPerc", "CC", "CCC")

for(OMtype in OMlist){
  
  PredDirecOM<-file.path(PredDirec,paste(OMtype))
  
  for(crt in CRtype){
    PredDirecCR<-file.path(PredDirecOM,paste(crt))
    PreyDirecCR<-paste(OMDirec,OMtype,"/",crt,"/", sep="")
    
    #get control rule number from filenames
    library(stringr)
    k<-list.files(path=PreyDirecCR,pattern="\\SimCharacteristics.txt$")
    # prepare regular expression
    regexp <- "[[:digit:]]+"
    
    CRnum<-str_extract(k, regexp)
   
    #for each CRnum, 
    for(crnum in CRnum){
      #get median of last "wantyears" years in each of the nsims
      predOut <-read.table(paste(PredDirecCR,"/",crnum,"predBNAvWtStatus.txt", sep=""))
      simMed <- predOut %>%
        group_by(Sim) %>%
        filter(Yr>(max(Yr)-wantyears)) %>%
        summarize(medPredB = median(PredB, na.rm=T),
                  medPredN = median(PredN, na.rm=T),
                  medPredRec = median(PredRec, na.rm=T),
                  medPredAvWt = median(PredAvWt, na.rm=T),
                  medPredB_status = median(PredB_status, na.rm=T),
                  nPredB_okstat = sum(PredB_status>=0.5, na.rm=T),
                  medPredAvWt_status = median(PredAvWt_status, na.rm=T),
                  nPredAvWt_goodstat = sum(PredAvWt_status>=1, na.rm=T),
                  NyrsForMed = wantyears-sum(is.na(PredB)))
      
      #then take the median and other quants of the nsim medians to get pred metric medians
      crMed <- data.frame(CRnum = crnum, 
                          Q25PredB = quantile(simMed$medPredB, probs=0.25, na.rm=T),
                          MedPredB = median(simMed$medPredB, na.rm=T),
                          Q75PredB = quantile(simMed$medPredB, probs=0.75, na.rm=T),
                          Q25PredN = quantile(simMed$medPredN, probs=0.25, na.rm=T),
                          MedPredN = median(simMed$medPredN, na.rm=T),
                          Q75PredN = quantile(simMed$medPredN, probs=0.75, na.rm=T),
                          Q25PredRec = quantile(simMed$medPredRec, probs=0.25, na.rm=T),
                          MedPredRec = median(simMed$medPredRec, na.rm=T),
                          Q75PredRec = quantile(simMed$medPredRec, probs=0.75, na.rm=T),
                          Q25PredAvWt = quantile(simMed$medPredAvWt, probs=0.25, na.rm=T),
                          MedPredAvWt = median(simMed$medPredAvWt, na.rm=T),
                          Q75PredAvWt = quantile(simMed$medPredAvWt, probs=0.75, na.rm=T),
                          Q25PredB_status = quantile(simMed$medPredB_status, probs=0.25, na.rm=T),
                          MedPredB_status = median(simMed$medPredB_status, na.rm=T),
                          Q75PredB_status = quantile(simMed$medPredB_status, probs=0.75, na.rm=T),
                          MinYrs_okBStatus = min(simMed$nPredB_okstat, na.rm=T),
                          AvPropYrs_okBstatus = mean(simMed$nPredB_okstat, na.rm=T)/wantyears,
                          Q25PredAvWt_status = quantile(simMed$medPredAvWt_status, probs=0.25, na.rm=T),
                          MedPredAvWt_status = median(simMed$medPredAvWt_status, na.rm=T),
                          Q75PredAvWt_status = quantile(simMed$medPredAvWt_status, probs=0.75, na.rm=T),
                          MinYrs_goodAvWtstatus = min(simMed$nPredAvWt_goodstat, na.rm=T),
                          AvPropYrs_goodAvWtstatus = mean(simMed$nPredAvWt_goodstat, na.rm=T)/wantyears,
                          NsimsForQuants = dim(simMed)[1]-sum(is.na(simMed$medPredB)))
        
      #then get CRnum characteristics
      crChar <- read.table(paste(PreyDirecCR,"/","Unadj", crnum, "SimCharacteristics.txt", sep=""), header=T)
      crChar <- distinct(select(crChar, -sim)) #should return 1 row
      crChar <- select(crChar, CtrlRule:TermBias) #keep only control rule
      
      #align with predator median results--this is one line in the pred results summary file
      crRes <- cbind(crMed, crChar)
      
      #append lines for every CRnum in the folder
      if(crnum == CRnum[1]){
        Res <- crRes
      } else {
        Resn <- crRes
        Res <- bind_rows(Res, Resn)
      }
      
    } #end crnum loop over control rule files
    
    write.table(Res, paste(PredDirecCR,"/",OMtype,"_",crt,"_","predSum.txt", sep=""), row.names = F)
    
  } #end crt loop of control rule type folers
} #end om loop over om folders