# summarize predator model results in same format as herring results
# S. Gaichas, herring MSE August 2016

# be in the correct predator OM and CRtype folder 
# find the equivalent prey OM and CRtype folder on network drive
# check that control rule numbers in predator folder match those in SimCharacteristics file in prey folder
# make a file with filename control rule numbers and sim characteristics?

require(dplyr)

wantyears<-50 #number of years to take median over for each sim for results

PredDirec<-c('~/Data/Projects/MSE/HerringMSE/Terns')
#PredDirec<-c('~/Documents/0_Data/MSE/HerringMSE/Terns')

OMDirec<-c('//net/mse/')
#OMtype<-c('HiM_LowSteep_NoAssBias_OldWt') #done first
OMlist<-c('HiM_LowSteep_NoAssBias_OldWt','LoM_HiSteep_NoAssBias_OldWt', 'LoM_HiSteep_NoAssBias_RecWt', 'HiM_LowSteep_NoAssBias_RecWt')
#OMlist<-c('HiM_LowSteep_AssBias_OldWt','LoM_HiSteep_AssBias_OldWt', 'LoM_HiSteep_AssBias_RecWt', 'HiM_LowSteep_AssBias_RecWt')
CRtype<-c("BB", "BB3yr", "BB5yr", "BB3yrPerc", "CC", "CCC")

TernCurrentPop <- 16000 #average pop in nesting pairs including Monomoy 1998-2015
TernK <- 45000  #carrying capacity in nesting pairs assumed based on historical "New England" pop BNA
TernProdTarg <- 1.0 #stated productivity target of F&W, results in eq pop
TernProdThresh <- 0.8 #suggested lower prod threshold at stakeholder workshop

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
      predOut <-read.table(paste(PredDirecCR,"/",crnum,"predBNAvWtStatus.txt", sep=""))
      #calculate status for all sims
      predOut <- predOut %>%
        mutate(PredNstatCur = PredN/TernCurrentPop,
               PredNstatK = PredN/TernK,
               PredProdstatTarg = PredProd/TernProdTarg,
               PredProdstatThresh = PredProd/TernProdThresh)
      
      #get median of last "wantyears" years in each of the nsims
      simMed <- predOut %>%
        group_by(Sim) %>%
        filter(Yr>(max(Yr)-wantyears)) %>%
        summarize(medPredN = median(PredN, na.rm=T),
                  medPredRec = median(PredRec, na.rm=T),
                  medPredProd = median(PredProd, na.rm=T),
                  medPredN_Curstatus = median(PredNstatCur, na.rm=T),
                  nPredN_CurPlus = sum(PredNstatCur>=1, na.rm=T),
                  medPredN_Kstatus = median(PredNstatK, na.rm=T),
                  medPredProd_Targ = median(PredProdstatTarg, na.rm=T),
                  nPredProd_Targplus = sum(PredProdstatTarg>=1, na.rm=T),
                  medPredProd_Thresh = median(PredProdstatThresh, na.rm=T),
                  nPredProd_Threshplus = sum(PredProdstatThresh>=1, na.rm=T),
                  NyrsForMed = wantyears-sum(is.na(PredN)))
      
      #then take the median and other quants of the nsim medians to get pred metric medians
      crMed <- data.frame(CRnum = crnum, 
                          Q25PredN = quantile(simMed$medPredN, probs=0.25, na.rm=T),
                          MedPredN = median(simMed$medPredN, na.rm=T),
                          Q75PredN = quantile(simMed$medPredN, probs=0.75, na.rm=T),
                          Q25PredRec = quantile(simMed$medPredRec, probs=0.25, na.rm=T),
                          MedPredRec = median(simMed$medPredRec, na.rm=T),
                          Q75PredRec = quantile(simMed$medPredRec, probs=0.75, na.rm=T),
                          Q25PredProd = quantile(simMed$medPredProd, probs=0.25, na.rm=T),
                          MedPredProd = median(simMed$medPredProd, na.rm=T),
                          Q75PredProd = quantile(simMed$medPredProd, probs=0.75, na.rm=T),
                          Q25PredN_Curstatus = quantile(simMed$medPredN_Curstatus, probs=0.25, na.rm=T),
                          MedPredN_Curstatus = median(simMed$medPredN_Curstatus, na.rm=T),
                          Q75PredN_Curstatus = quantile(simMed$medPredN_Curstatus, probs=0.75, na.rm=T),
                          MinYrs_goodCurPlus = min(simMed$nPredN_CurPlus, na.rm=T),
                          AvPropYrs_goodCurPlus = mean(simMed$nPredN_CurPlus, na.rm=T)/wantyears,
                          Q25PredProd_Targ = quantile(simMed$medPredProd_Targ, probs=0.25, na.rm=T),
                          MedPredProd_Targ = median(simMed$medPredProd_Targ, na.rm=T),
                          Q75PredProd_Targ = quantile(simMed$medPredProd_Targ, probs=0.75, na.rm=T),
                          MinYrs_goodProd_Targplus = min(simMed$nPredProd_Targplus, na.rm=T),
                          AvPropYrs_goodProd_Targplus = mean(simMed$nPredProd_Targplus, na.rm=T)/wantyears,
                          Q25PredProd_Thresh = quantile(simMed$medPredProd_Thresh, probs=0.25, na.rm=T),
                          MedPredProd_Thresh = median(simMed$medPredProd_Thresh, na.rm=T),
                          Q75PredProd_Thresh = quantile(simMed$medPredProd_Thresh, probs=0.75, na.rm=T),
                          MinYrs_goodProd_Threshplus = min(simMed$nPredProd_Threshplus, na.rm=T),
                          AvPropYrs_goodProd_Threshplus = mean(simMed$nPredProd_Threshplus, na.rm=T)/wantyears,
                          NsimsForQuants = dim(simMed)[1]-sum(is.na(simMed$medPredN)))
        
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