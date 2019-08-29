# read in all summarized predator results, add om and cr columms
# save as rds files for further visualization in VisPredResults.R
# November 2016

rm(list=ls())

library(tools)
library(dplyr)
library(ggplot2)

PredDirec<-c('~/Data/Projects/MSE/HerringMSE')

#combine extra metric summaries for each predator; Dec 15 2016

Pred <- c("Tuna", "Terns", "Groundfish")
ResultsDir <- c("ExtraResultsSummaries_Tuna", "ExtraResultsSummaries_Terns", "ExtraResultsSummaries_Groundfish")

for (i in 1:length(Pred)){
  PredResultsDir <- paste(PredDirec,Pred[i],ResultsDir[i], sep="/")
  all <- list.files(path=PredResultsDir)
  
  ldf = lapply(all, function(x) {
    
    dat = read.table(file.path(PredResultsDir,x), header=TRUE)
    
    colvals = strsplit(x, "_")
    
    # Add columns for om and cr types 
    omcr <- data.frame(M=rep(colvals[[1]][1], dim(dat)[1]), 
                       Steep=rep(colvals[[1]][2], dim(dat)[1]), 
                       Bias=rep(colvals[[1]][3], dim(dat)[1]), 
                       Wt=rep(colvals[[1]][4], dim(dat)[1]), 
                       CRtype=rep(colvals[[1]][5], dim(dat)[1]))
    
    
    #merge them
    dat <- bind_cols(dat, omcr)
    
    return(dat)
  })
  
  df = bind_rows(ldf)
  
  predall<-df
  saveRDS(predall, paste0(Pred[i],"extra.rds"))
  
}


#old code below originally used to combine

#############################################################################

#Pred <- c("Tuna", "Terns")
#ResultsDir <- c("ResultsSummaries_all_Tuna", "ResultsSummaries_Terns")
Pred <- c("Groundfish")
ResultsDir <- c("ResultsSummaries_Groundfish")

PredResultsDir <- paste(PredDirec,Pred,ResultsDir, sep="/")

allgf <-list.files(path=PredResultsDir)

allgfcols <- strsplit(allgf, "_")

#dataset <- ldply(file.path(PredResultsDir[1], alltuna), read.table, header=TRUE)

listtxt <- allgf
#colvals <- alltunacols

ldf = lapply(listtxt, function(x) {
  
  dat = read.table(file.path(PredResultsDir[1],x), header=TRUE)
  
  colvals = strsplit(x, "_")
  
  # Add columns for om and cr types 
  omcr <- data.frame(M=rep(colvals[[1]][1], dim(dat)[1]), 
                     Steep=rep(colvals[[1]][2], dim(dat)[1]), 
                     Bias=rep(colvals[[1]][3], dim(dat)[1]), 
                     Wt=rep(colvals[[1]][4], dim(dat)[1]), 
                     CRtype=rep(colvals[[1]][5], dim(dat)[1]))
  
  
  #merge them
  dat <- bind_cols(dat, omcr)
  
  return(dat)
})


df = bind_rows(ldf)

gfishall<-df
saveRDS(gfishall, "gfishall.rds")


PredResultsDir <- paste(PredDirec,Pred,ResultsDir, sep="/")

alltuna <-list.files(path=PredResultsDir[1])

alltunacols <- strsplit(alltuna, "_")

#dataset <- ldply(file.path(PredResultsDir[1], alltuna), read.table, header=TRUE)

listtxt <- alltuna
#colvals <- alltunacols

ldf = lapply(listtxt, function(x) {
  
  dat = read.table(file.path(PredResultsDir[1],x), header=TRUE)
  
  colvals = strsplit(x, "_")
  
  # Add columns for om and cr types 
  omcr <- data.frame(M=rep(colvals[[1]][1], dim(dat)[1]), 
                     Steep=rep(colvals[[1]][2], dim(dat)[1]), 
                     Bias=rep(colvals[[1]][3], dim(dat)[1]), 
                     Wt=rep(colvals[[1]][4], dim(dat)[1]), 
                     CRtype=rep(colvals[[1]][5], dim(dat)[1]))
  
  
  #merge them
  dat <- bind_cols(dat, omcr)
  
  return(dat)
})


df = bind_rows(ldf)

tunaall<-df
saveRDS(tunaall, "tunaall.rds")


allterns <-list.files(path=PredResultsDir[2])

allterncols <- strsplit(allterns, "_")

#dataset <- ldply(file.path(PredResultsDir[1], alltuna), read.table, header=TRUE)

listtxt <- allterns
#colvals <- alltunacols

ldf = lapply(listtxt, function(x) {
  
  dat = read.table(file.path(PredResultsDir[2],x), header=TRUE)
  
  colvals = strsplit(x, "_")
  
  # Add columns for om and cr types 
  omcr <- data.frame(M=rep(colvals[[1]][1], dim(dat)[1]), 
                     Steep=rep(colvals[[1]][2], dim(dat)[1]), 
                     Bias=rep(colvals[[1]][3], dim(dat)[1]), 
                     Wt=rep(colvals[[1]][4], dim(dat)[1]), 
                     CRtype=rep(colvals[[1]][5], dim(dat)[1]))
  
  
  #merge them
  dat <- bind_cols(dat, omcr)
  
  return(dat)
})


df = bind_rows(ldf)

ternsall<-df
saveRDS(ternsall, "ternsall.rds")

