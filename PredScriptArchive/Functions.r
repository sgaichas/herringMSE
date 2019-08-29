############################# a function to make a vector out of a constant and a matrix out of a vector
fill<-function(x,n,m){
  if(length(x)==1) x<-rep(x,m) #make a vector out of a constant
  if(is.matrix(x)==FALSE) return(matrix(rep(x,n),ncol=m,nrow=n,byrow=TRUE)) #make a matrix out of a vector
  #A little ticky when you are running the unfished algorith, and need to get 1000 years out
  #of a time varying vector - this will hold the last year's values constant for the remaining years
  if(dim(x)[1]<n) {
    return(rbind(x,matrix(rep(x[dim(x)[1],],(n-dim(x)[1])),ncol=m,nrow=(n-dim(x)[1]),byrow=TRUE)))
  }
  if(dim(x)[1]==n & dim(x)[2]==m) return(x)
}
#############################

#########################Given M, Wt, Mat, steepness, and B0; returns R0, alpha and beta of BevHolt or Ricker.
get.SR.pars<-function(M=NULL,Wt=NULL,Mat=NULL,steep=NULL
          ,B0=NULL,type=NULL){

     dim.var<-dim(Mat)[2]
     if(is.null(dim.var)) dim.var<-length(Mat) #Mat can be an array or vector...
     #if(length(M)==1) M<-fill(M,n=dim.var,m=1)

     nvec<-c()
     nvec[1]<-1 #set the population to one so you don't have to rescale SSBR
     for(i in 2:(dim.var-1)){nvec[i]<-nvec[i-1]*exp(-M[i-1])}  #unfished abundance vector
     nvec[dim.var]<-nvec[(dim.var-1)]*exp(-M[(dim.var-1)])*(1/(1-exp(-M[dim.var]))) #plus group!
     SSBR0<-sum(Mat*Wt*nvec)  #unfished SSBR

     if(type==1) { #Beverton-Holt stock recruit
      R0<-(1/SSBR0)*B0
      Reca<-(B0/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
      Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta
      return(data.frame("R0"=R0,"alpha"=Reca,"beta"=Recb,"SSBR0"=SSBR0))
     } else if(type==2){ #Ricker
      beta<-(log(steep)-log(0.2))/(0.8*B0)
      alpha<-(exp(beta*B0))/SSBR0
      return(data.frame("R0"=NULL,"alpha"=alpha,"beta"=beta,"SSBR0"=SSBR0)    )
     }
}
###################################

################YPR function used to calculate MSY parms, and Fx% levels
max.ypr<-function(F,M,Wt,Mat,selectivity.F,B0,type,med.recr,steep,SSBR0,B.MSYflag){
  #function to be maximized for MSY calculation
  nvec<-c()
  nvec[1]<-1 #set the population to one so you don't have to rescale SSBR
  Faa<-F*selectivity.F
  Z<-M+Faa #total mortality
  for(i in 2:(length(Mat)-1)){nvec[i]<-nvec[i-1]*exp(-Z[i-1])} #ages up to oldest
  nvec[length(Mat)]<-nvec[(length(Mat)-1)]*exp(-Z[length(Mat)-1])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  SSBR<-sum(Mat*Wt*nvec)

  if(type==1) { #Beverton-Holt stock recruit
      R0<-(1/SSBR0)*B0
      Reca<-(B0/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
      Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta
      SSBhat<-(SSBR-Reca)/Recb  #equlibrium SSB
      Rhat<-SSBhat/SSBR        #equilibrium recruitment
  } else if(type==2){ #Ricker
      beta<-(log(steep)-log(0.2))/(0.8*B0)
      alpha<-(exp(beta*B0))/SSBR0
      SSBhat<-(log(alpha)-log(1/SSBR))/(beta)  #equlibrium SSB
      Rhat<-alpha*SSBhat*exp(-beta*SSBhat)     #equilibrium recruitment
  } else (Rhat<-med.recr) #Random recruitment - just use the median for the year being evaluated
  yield<-sum((Faa/Z)*nvec*Wt*(1-exp(-1*Z)))
  if(B.MSYflag) {         #see MSY function for need for flag
    return(SSBhat)
  } else return(yield*Rhat)
}
###################################################################################

########################Function to return MSY, FMSY.  Calls to max.ypr
MSY.find<-function(M=NULL,Wt=NULL,Mat=NULL,steep=NULL,selectivity.F=NULL
			    ,B0=NULL,type=NULL,med.recr=NULL,SSBR0=NULL,B.MSYflag=NULL)
#need B.MSYflag because it would screw up optimize function to return something from max.ypr other than MSY (i.e., yield*Rhat)
{
     MSY<-optimize(f=max.ypr,interval=c(0,10),M,Wt,Mat,selectivity.F,B0,type,med.recr,steep,SSBR0,B.MSYflag,maximum=TRUE)
     return(data.frame("MSY"=MSY$objective,"F.MSY"=MSY$maximum))
}
###################################################################################

######################Function to find F and NAA from a given level of unfished SSB (not SSBR); similar to max.ypr
FSSB0<-function(F,M,Wt,Mat,selectivity.F,B0,type,med.recr,steep,SSBR0,SSBfrac,EqNAAflag){
  #function to be minimized to find level of F that will produce a specified percent of unfished SSB (not SSBR)
  nvec<-c()
  nvec[1]<-1 #set the population to one so you don't have to rescale SSBR
  Faa<-F*selectivity.F
  Z<-M+Faa #total mortality
  for(i in 2:(length(Mat)-1)){nvec[i]<-nvec[i-1]*exp(-Z[i-1])} #ages up to oldest
  nvec[length(Mat)]<-nvec[(length(Mat)-1)]*exp(-Z[length(Mat)-1])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  SSBR<-sum(Mat*Wt*nvec)

  if(type==1) { #Beverton-Holt stock recruit
      R0<-(1/SSBR0)*B0
      Reca<-(B0/R0)*((1-steep)/(4*steep)) #Bev-Holt alpha parm (See Mangel et al 2010 Fish and Fisheries)
      Recb<-(5*steep-1)/(4*steep*R0)  #Bev-Holt Beta
      SSBhat<-(SSBR-Reca)/Recb  #equlibrium SSB
      Rhat<-SSBhat/SSBR        #equilibrium recruitment
      NAA<-c()
      NAA[1]<-Rhat      #number at age-1 equals equilibrium recruitment
      for(i in 2:(length(Mat)-1)){NAA[i]<-NAA[i-1]*exp(-Z[i-1])} #ages up to oldest
      NAA[length(Mat)]<-NAA[(length(Mat)-1)]*exp(-Z[length(Mat)-1])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  } else if(type==2){ #Ricker
      beta<-(log(steep)-log(0.2))/(0.8*B0)
      alpha<-(exp(beta*B0))/SSBR0
      SSBhat<-(log(alpha)-log(1/SSBR))/(beta)  #equlibrium SSB
      Rhat<-alpha*SSBhat*exp(-beta*SSBhat)     #equilibrium recruitment
      NAA<-c()
      NAA[1]<-Rhat      #number at age-1 equals equilibrium recruitment
      for(i in 2:(length(Mat)-1)){NAA[i]<-NAA[i-1]*exp(-Z[i-1])} #ages up to oldest
      NAA[length(Mat)]<-NAA[(length(Mat)-1)]*exp(-Z[length(Mat)-1])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  } else {
      Rhat<-med.recr #Random recruitment - just use the median for the year being evaluated
      NAA<-c()
      NAA[1]<-Rhat      #number at age-1 equals equilibrium recruitment
      for(i in 2:(length(Mat)-1)){NAA[i]<-NAA[i-1]*exp(-Z[i-1])} #ages up to oldest
      NAA[length(Mat)]<-NAA[(length(Mat)-1)]*exp(-Z[length(Mat)-1])*(1/(1-exp(-Z[length(Mat)]))) #plus group!
  }
  
  ypr<-sum((Faa/Z)*nvec*Wt*(1-exp(-1*Z)))
  EqYield<-ypr*Rhat
  if(EqNAAflag) {
    return(c(data.frame("NAA"=NAA,"ZAA"=Z),"EqYield"=EqYield))
  } else return((SSBhat-(SSBfrac*B0))*(SSBhat-(SSBfrac*B0)))
}

#######################Function to return the NAA, ZAA, and else for a specified level of F_unfishedSSB (not SSBR0); calls to FSSB0
PercB0<-function(M=NULL,Wt=NULL,Mat=NULL,steep=NULL,selectivity.F=NULL,B0=NULL,type=NULL,med.recr=NULL,SSBR0=NULL,SSBfrac=NULL,EqNAAflag=NULL)
 {
  PercB0<-optimize(f=FSSB0,interval=c(0,10),M,Wt,Mat,selectivity.F,B0,type,med.recr,steep,SSBR0,SSBfrac,EqNAAflag,maximum=FALSE)
  return(data.frame("PercB0objfxn"=PercB0$objective,"F.SSBfrac"=PercB0$minimum))
 }
#########################################################

###########################Stock-recruitment options
stock.recruit<-function(SSB1=NULL,Reca=NULL,Recb=NULL,type=NULL)
{
      if(type==1) { #Beverton-Holt stock recruit
      return((SSB1/(Reca+Recb*SSB1)) )
     } else if(type==2){ #Ricker
      return(alpha*SSB1*exp(-beta*SSB1))
     }  else (print(paste("ERROR - unrecognized stock recruit model, type=", type,"must be 1, 2",sep=" ")))
}
################################################################# end SR

#########################################Newton Raphson###############
Solve.F<-function(Fnew=NULL,Fsel=NULL,M=NULL,Wt=NULL,NAA=NULL,quotatarg=NULL)
{
  for (n in 1:7) #Newton Raphson loop
  {
   delta=Fnew*0.0001
   fofFa=sum(((Fnew*Fsel)/(Fnew*Fsel+M))*NAA*Wt*(1-exp(-1*(Fnew*Fsel+M))))
   fprimeFhi=sum((((Fnew+delta)*Fsel)/((Fnew+delta)*Fsel+M))*NAA*Wt*(1-exp(-1*((Fnew+delta)*Fsel+M))))
   fprimeFlo=sum((((Fnew-delta)*Fsel)/((Fnew-delta)*Fsel+M))*NAA*Wt*(1-exp(-1*((Fnew-delta)*Fsel+M))))

   fofF=fofFa-quotatarg
   fprimeF=(fprimeFhi-fprimeFlo)/(2.0*delta)
   Fnew=Fnew-(fofF/fprimeF)
   if (Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF))
   if (Fnew>3) Fnew=3
  } #end newton raphson for loop
  return(Fnew)
}
###################################end newton raphson