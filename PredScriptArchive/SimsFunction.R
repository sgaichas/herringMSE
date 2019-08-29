sims<-function(u=NA,l=NA,f=NA,c=NA,CtrlRule=NA,FracBmsyThreshHi=NA,FracBmsyThreshLo=NA,FracFtarg=NA,useminmax=NA,minquota=NA,maxquota=NA,usequotavar=NA,quotavar=NA,Thresh=NA,Prop=NA,Propc=NA,useFcap=NA,Fcapprop=NA,OverfishedSSB=NA,quotablock=NA,CodeDirect=NA)  {
  
  CodeDirect<-CodeDirect  #location of functions.r and userinputs.r
  source(paste(CodeDirect,"Functions.r",sep="")) #source required functions
  source(paste(CodeDirect,"UserInputs.r",sep="")) #source user inputs
  
  #########################Go get user inputs##################
  inputs<-userinputs()
  nages<-inputs$nages
  nyears<-inputs$nyears
  nsims<-inputs$nsims
  startseed<-inputs$startseed  #For specifying random number seed.  This set up ensures differences aren't due to random number generation
  ###############################
  
  #####Recruitment stuff########
  RecType<-inputs$RecType #Can only do Bev Holt 12/11/14; BevHolt=1; Ricker =2
  ############################
  
  #######Inputs if doing retrospective sims; If no retro then check over quick, but don't worry about it########
  AdjSwitch<-inputs$AdjSwitch  #####Switch that defines whether target fullly selected F is based on Mohn's Rho adjusted or unadjusted SSB.  0 is undadjusted, else adjusts for Mohn's Rho
  peels<-inputs$peels    #number of peels for mohn's rho calculation used by managers
  SSBSimYearDiffSwitch<-inputs$SSBSimYearDiffSwitch #Switch that turns on or off (0 is off; else is on) the calculation of the %difference in SSB between adj and unaj for each year and simulation.
  ###Turn SSBSimYearDiffSwitch on if either the adjusted or unadjusted sims of a given pair has already been run.  For example,
  ###Run Unadjusted sims with retro above the truth and switch off.  Then run adjusted sims with switch on.  The code using the switch requires results from both adj and unaj sims.
  ####################################end getting user inputs################
  
  ###########setup standard storage stuff#############
  TrueF<-matrix(NA,nrow=nyears,ncol=1) #F actually applied to the stock after assessment and implementation error
  TotBio<-matrix(NA,nrow=nyears,ncol=1) #total biomass in each year
  SSB<-matrix(NA,nrow=nyears,ncol=1) #SSB in each year
  SSBAA<-matrix(NA,nrow=nyears,ncol=nages) #SSB at age in each year
  SSBretro<-matrix(NA,nrow=nyears,ncol=nyears) #assessed SSB in each year
  Nretro<-matrix(NA,nyears,nyears)  #matrix of assessed abundance in each year with retrospective
  NretroTermAA<-matrix(NA,nyears,nages) #matrix of terminal year's assessed abundance at age with retrospective.  Necessary for quota setting (see below quotatarg).
  FAA<-matrix(NA,nrow=nyears,ncol=nages)  #F at age in each year
  exploitAA<-matrix(NA,nrow=nyears,ncol=nages)  #F at age in each year
  quotatarg<-matrix(NA,nrow=nyears,ncol=1) #quota target set by managers each year
  ZAA<-matrix(NA,nrow=nyears,ncol=nages) #total mortality at age in each year
  NAA<-matrix(NA,nrow=nyears,ncol=nages)  #numbers at age in each year
  TotN<-matrix(NA,nrow=nyears,ncol=1) #sum of true abundance at age in each year
  TotCatch<-matrix(NA,nrow=nyears,ncol=1) #vector of total catch each year
  TotYield<-matrix(NA,nrow=nyears,ncol=1) #vector of yield each year (CAA * Wt)
  TotNatDeath<-matrix(NA,nrow=nyears,ncol=1) #vector of natural deaths, much like yield
  YieldVar<-matrix(NA,nrow=nyears,ncol=1) #squared difference between yield this year and yield last year for interannual variation in yield calc (A'mar et al 2009 CJFAS)
  RealTargFDiff<-matrix(NA,nrow=nyears,ncol=1) #Percent difference between the target F and the true/realized F
  NerrRan<-matrix(NA,nrow=nyears,ncol=1) #year specific autocorrelated assessment errors
  NerrRanNorm<-matrix(NA,nrow=nyears,ncol=1) #year specific random variables for assessment errors
  Recerr<-matrix(NA,nrow=nyears,ncol=1) #year specific normal R.V. for recruitment process errors
  Recerrauto<-matrix(NA,nrow=nyears,ncol=1) #year specific, autocorr, recruitment process errors
  TrueM<-matrix(NA,nrow=nyears,ncol=nages) #year and age specific M; filled using "separability" assumption below
  Merrauto<-matrix(NA,nrow=nyears,ncol=1) #year specific rand walk process errors around mean M
  ############################################
  
  ##########setup results storage stuff###############
  OverfishedSSBSimYear<-matrix(NA,nsims,nyears) #1 if given year in a given sim is overfished, else 0.  Used for %overfished (see below).
  ClosureSimYear<-matrix(NA,nsims,nyears) #1 if given year in a given sim has closed fishery, else 0.  Used for %fishery closure (see below).
  ###########################################
  
  #############setup retro sim storage stuff that's only relevant if running with retros
  QuotaAdjDiff<-matrix(NA,nrow=nyears,ncol=1) #Percent difference in the quota with and without adjusting for Mohn's Rho in each year
  MohnsRho<-matrix(NA,nyears,1) #Mohn's rho as estimated by assessment biologists using x number of peels specified above.
  MohnsRhoStep<-matrix(NA,peels,1) #See calculation of Mohn's Rho below.  Mohn's rho values for each tip for each of the peels.  Avg of these equals the Mohn's Rho used by managers.
  MeanQuotaAdjDiffSims<-matrix(NA,nsims,1) #mean percent difference between the quota with and without adjusting for Mohn's Rho.
  QuotaDiffByYearMed<-matrix(NA,nyears,2) #median % difference among sims in quota between adjusting and not adjusting for each year
  QuotaDiffByYear25<-matrix(NA,nyears,2) #25th % difference among sims in quota between adjusting and not adjusting for each year
  QuotaDiffByYear75<-matrix(NA,nyears,2) #75th % difference among sims in quota between adjusting and not adjusting for each year
  SSBSimYearDiffMed<-matrix(NA,nyears,2) #median among sims of %diff in SSB between adj and unadj sims for each year
  SSBSimYearDiff25<-matrix(NA,nyears,2) #25th among sims of %diff in SSB between adj and unadj sims for each year
  SSBSimYearDiff75<-matrix(NA,nyears,2) #75th among sims of %diff in SSB between adj and unadj sims for each year
  #############
  
  #############Simulation loop
  for (s in 1:nsims)
  {
    if(s==1) {seed=startseed} else {seed=seed+100}  #Changes the random number seed among simulations; startseed specified in user inputs.
    set.seed(seed) #Sets random number seed; This set up ensures differences aren't due to random number generation
    
    ##########################################Read in user inputs each sim because some or all can be specified as random variables.
    inputs<-userinputs() #run user inputs because some inputs are random variable that must change each sim
    M<-inputs$M    #fully selected natural mortality;
    Msel<-inputs$Msel #how relative M changes with age, like selectivity
    Mrho<-inputs$Mrho #degree of auto corr for M AR1 process error
    Msd<-inputs$Msd #std dev of M process error
    Mat<-inputs$Mat #maturity at age
    Wt<-inputs$Wt  #wt at age
    Fsel<-inputs$Fsel #fishery selectivity at age
    SSBfrac<-inputs$SSBfrac #proportion of unfished SSB (not SSBR) at which to find equilibrium NAA for starting each sim in year 1; desired level of depletion in year 1.
    steepness<-inputs$steepness #Maybe link this to M,
    Bzero<-inputs$Bzero  #unfished spawning stock biomass
    Recsd<-inputs$Recsd #standard deviation of lognormal recruitment process error each simulation
    Recrho<-inputs$Recrho #autocorrelation of recruitment process errors
    NerrRho<-inputs$NerrRho #level of autocorrelation in assessment error (terminal year estimates)
    Nerrsd<-inputs$Nerrsd #standard deviation of assessment error (terminal year estimates)
    Impsd<-inputs$Impsd #standard deviation of implementation error
    MohnsRhoInit<-inputs$MohnsRhoInit #Mohn's Rho values for year one to number of peels.  These years don't have the adequate assessed values to calculate a Mohn's Rho and so must be fixed. Perhaps should be changed to the mean Mohn's Rho suggested by whatever the mohnsrhofan values are defined below.
    mohnsrhotermA<-inputs$mohnsrhoterm #"mohn's rho" for terminal year estimates of each "assessment"
    mohnsrhoterm<-mohnsrhotermA #done so that bias (mohnsrhoterm) can be linked to error in M
    mohnsrhofan<-inputs$mohnsrhofan
    NerrFansd<-inputs$NerrFansd
    ############################################################
    
    ###########################Call function to calculate Rzero and parms of recruitment relationship, given various other life history characteristics
    Rzero<-get.SR.pars(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,B0=Bzero,type=RecType)$R0     #unfished recruitment
    Reca<-get.SR.pars(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,B0=Bzero,type=RecType)$alpha   #recruitment parameter
    Recb<-get.SR.pars(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,B0=Bzero,type=RecType)$beta    #recruitment parameter
    SSBR0<-get.SR.pars(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,B0=Bzero,type=RecType)$SSBR0   #unfished SSBR
    ###########################
    
    ###########################MSY reference points based on S-R and other parms
    Fmsy<-MSY.find(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,SSBR0=SSBR0,B.MSYflag=F)$F.MSY
    MSY<-MSY.find(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,SSBR0=SSBR0,B.MSYflag=F)$MSY
    Bmsy<-max.ypr(F=Fmsy,M=M*Msel,Wt=Wt,Mat=Mat,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,steep=steepness,SSBR0=SSBR0,B.MSYflag=T)
    SSBthresh<-Bmsy*FracBmsyThreshHi # For BB control rule; Threshold to switch from Fmsy as target F to linear decline to zero
    SSBthreshlo<-Bmsy*FracBmsyThreshLo #For BB control rule; Level of SSB where target F set to 0
    NAAmsy<-FSSB0(F=Fmsy,M=M*Msel,Wt=Wt,Mat=Mat,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,steep=steepness,SSBR0=SSBR0,SSBfrac=(Bmsy/Bzero),EqNAAflag=T)$NAA #NAA at MSY
    NAAzero<-FSSB0(F=0,M=M*Msel,Wt=Wt,Mat=Mat,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,steep=steepness,SSBR0=SSBR0,SSBfrac=(Bmsy/Bzero),EqNAAflag=T)$NAA #NAA at unfished
    SSBthreshPT<-Bmsy*Thresh #for prop. threshold control rule;   level of SSB below which quota=0
    BcutPT<-(NAAmsy*Wt)*Thresh #for prop. threshold; the cutoff at age in total biomass; needed because fishing is on total bio and not just SSB.  So cut defined by SSB, but applied in total bio
    Fcap<-Fcapprop*Fmsy
    ##########################
    
    ########################### level of F that produces the specified level of SSB depletion and the equilibrium NAA at that F; mostly for characterisitics in year 1
    F.SSBfrac<-PercB0(M=M*Msel,Wt=Wt,Mat=Mat,steep=steepness,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,SSBR0=SSBR0,SSBfrac=SSBfrac,EqNAAflag=F)$F.SSBfrac
    EqNAA<-FSSB0(F=F.SSBfrac,M=M*Msel,Wt=Wt,Mat=Mat,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,steep=steepness,SSBR0=SSBR0,SSBfrac=SSBfrac,EqNAAflag=T)$NAA
    EqYield<-FSSB0(F=F.SSBfrac,M=M*Msel,Wt=Wt,Mat=Mat,selectivity.F=Fsel,B0=Bzero,type=RecType,med.recr=0,steep=steepness,SSBR0=SSBR0,SSBfrac=SSBfrac,EqNAAflag=T)$EqYield
    ############################
    
    #Sim Counter for keeping track of number of sims when running.
    if(s==1) { write.table(s,paste(readfile2,"SimCount.txt",sep=""),append=FALSE,row.names=FALSE,col.names=FALSE,quote=FALSE) } else { write.table(s,paste(readfile2,"SimCount.txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE) }
    
    #Fill years one to start each sim from the same place
    NAA[1,]<-EqNAA  #defines NAA in year one
    TotN[1,1]=sum(NAA[1,])
    SSB[1,1]=sum(NAA[1,]*Mat*Wt)
    TotBio[1,1]=sum(NAA[1,]*Wt)
    Nretro[1,1]=sum(NAA[1,])*(mohnsrhoterm+1)
    NretroTermAA[1,]=NAA[1,]*(mohnsrhoterm+1)
    SSBretro[1,1]=sum(NAA[1,]*Mat*Wt)*(mohnsrhoterm+1)
    quotatarg[1,1]=EqYield
    TrueM[1,]<-M*Msel
    
    #loop to fill in Mohn's Rho values for year one to number of peels.  These years don't have the adequate assessed values to calculate a Mohn's Rho and so must be fixed here.
    for (q in 1:peels)
    {
      MohnsRho[q,1]=MohnsRhoInit
    }
    
    quotacounter<-1 #for multiyear quota setting; see control rules below
    
    #############################Start year loop
    for (i in 2:nyears) #year loop; Coded such that reading i-1 as "this year" makes more sense
    {
      Imperr=rnorm(1,0,Impsd) #Implementation error for each year (see below)
      NerrFan=rnorm(1,0,NerrFansd) #Assessment error for estimates other than terminal year (see below)
      NerrRanNorm[i,1]=rnorm(1,0,Nerrsd) #Random normal variable for use in assessment error below
      Recerr[i,1]<-rnorm(1,0,Recsd) #year specfic normal R.V. for autocorr errors below
      Merr<-rnorm(1,0,Msd)
      if (i==2) {
        NerrRan[i,1]=NerrRanNorm[i,1]*((Nerrsd)/(sqrt(1-(NerrRho*NerrRho)))) #First year of sim comes from stationary distribution.
        Recerrauto[i,1]<-Recerr[i,1]*((Recsd)/(sqrt(1-(Recrho*Recrho))))  #First year of sim comes from stationary distribution.
        Merrauto[i,1]<-Merr*((Msd)/(sqrt(1-(Mrho*Mrho))))  #First year of sim comes from stationary distribution.
      } else {
        NerrRan[i,1]=NerrRho*NerrRan[i-1,1]+(sqrt(1-(NerrRho*NerrRho))*NerrRanNorm[i,1]) #for year specific, autocorrelated, assessment error. 
        Recerrauto[i,1]<-Recrho*Recerrauto[i-1,1]+(sqrt(1-(Recrho*Recrho))*Recerr[i,1]) #year specific autorcorrelated recruitment process errors
        Merrauto[i,1]<-Mrho*Merrauto[i-1,1]+(sqrt(1-(Mrho*Mrho))*Merr) #year specific autorcorrelated M process errors
      } #for year specific, autocorrelated, assessment error.  First year of sim comes from stationary distribution.
      
      TrueM[i,]<-(M*exp(Merrauto[i,1]-((Msd*Msd)/2)))*Msel #specify M each year
      mohnsrhoterm<-mohnsrhotermA+((M-max(TrueM[i,]))/max(TrueM[i,])) #use TrueM[i] because applied to N[i] farther below.  Want bias in same year.
      
      ifelse (AdjSwitch==0,SSBstatus<-SSBretro[i-1,i-1],SSBstatus<-SSBretro[i-1,i-1]/(MohnsRho[i-1,1]+1))  #Based on switch value, uses either SSB unadjusted or Mohn'sRho adjust to set Target F in control rule
      
      #Control rules:
      if (CtrlRule==1) {      #biomass based
        ifelse (SSBstatus>=SSBthresh,Ftarg<-Fmsy*FracFtarg,Ftarg<-(Fmsy*FracFtarg)*((SSBstatus-SSBthreshlo)/(SSBthresh-SSBthreshlo))) #set target fully selected F based on a control rule (See Katsukawa 2004 Fisheries Science)
        if (Ftarg<0) Ftarg<-0  #when SSBstatus is less than SSBthreshlo then Ftarg is negative based on control rule.  This fixes that problem.
        if (TotN[i-1,1]==0) Ftarg<-0 #Assumes that no fishing occurs after extinction.  Without this statement the Newton Raphson iterations explode.
        FAA[i-1,]=Ftarg*Fsel  #Determine target F at age as product of fully selected F and selectivity at age
        ####Set a quota. Note that quota in year i is based on the assessment in year i-1.  This creates a lag in the assessment and the quota.  The lag could be altered if desired.  Quota in biomass set by managers based on target F and terminal year assessed abundance; either adjusted or unadjusted for Mohn's Rho depending on AdjSwitch.
        if(quotacounter>(quotablock-1)) { #used for multi year quotas
          ifelse (AdjSwitch==0, quotatarg[i,1]<-sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*NretroTermAA[i-1,]*Wt*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))),quotatarg[i,1]<-sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*((NretroTermAA[i-1,]*Wt)/(MohnsRho[i-1,1]+1))*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))))
          quotacounter<-1
        } else {
          quotatarg[i,1]<-quotatarg[i-1,1] #if quotacounter < # of multi year quota desired set equal to prior year
          quotacounter<-quotacounter+1 #increment counter
        }
        #if desired, set min and max levels of quota
        if(useminmax>0){
          if(quotatarg[i-1,1]<minquota) { quotatarg[i-1,1]<-minquota }
          if(quotatarg[i-1,1]>maxquota) { quotatarg[i-1,1]<-maxquota }
        }  
        if(usequotavar>0){ #control degree of interannual variation in quota?
          if(quotatarg[i-1,1]>0) { #needed because if quota =0 then it would stay 0 forever
          if(quotatarg[i,1]>(quotatarg[i-1,1]+(quotatarg[i-1,1]*quotavar))) { #if quota next year > than some amount of quota this year then
             quotatarg[i,1]<-(quotatarg[i-1,1]+(quotatarg[i-1,1]*quotavar)) #quota next year = some fraction of last quota
          }
          if(quotatarg[i,1]<(quotatarg[i-1,1]-(quotatarg[i-1,1]*quotavar))) { #if quota next year < than some amount of quota this year then
            quotatarg[i,1]<-(quotatarg[i-1,1]-(quotatarg[i-1,1]*quotavar)) #quota next year = some fraction of last quota
          }
          } #end if quota > 0 
        } #end if for controlling interannual varition in quota
        ####End quota setting.
        ###Newton Raphson iterations to determine the F that would result in the quota being removed from the true population; needed due to lag between assess. yr and quotatarg yr (see above)
        #if (Ftarg==0) { Fnew<-0}
        if (quotatarg[i-1,1]==0) { Fnew<-0 }
        else {
          Fnew<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i-1,1])
        }
        ###End Newton Raphson
        ############################################
      } else if (CtrlRule==2) {  #proportional threshold
        if (SSBstatus>SSBthreshPT) {
          exploitAA[i-1,]=((Fmsy*Prop*Fsel)/((Fmsy*Prop*Fsel)+TrueM[i-1,]))*(1-exp(-((Fmsy*Prop*Fsel)+TrueM[i-1,]))) #convert instantaneous F to annual rate
          ####Set a quota. Note that quota in year i is based on the assessment in year i-1.  This creates a lag in the assessment and the quota.  The lag could be altered if desired.  Quota in biomass set by managers based on target F and terminal year assessed abundance; either adjusted or unadjusted for Mohn's Rho depending on AdjSwitch.
          ifelse (AdjSwitch==0, quotatarg[i,1]<-sum(exploitAA[i-1,]*((NretroTermAA[i-1,]*Wt)-BcutPT)),quotatarg[i,1]<-sum(exploitAA[i-1,]*(((NretroTermAA[i-1,]*Wt)/(MohnsRho[i-1,1]+1))-BcutPT)))
          ####End quota setting.
          ###Newton Raphson iterations to determine the F that would result in the quota being removed from the entire true population; accounting for assessment lag; it's Biomass based calculations in reverse
          if (quotatarg[i-1,1]==0) { Ftarg<-0}
          else {
            Ftarg<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i,1])
          }
          FAA[i-1,]=Ftarg*Fsel
          ###End Newton Raphson
          
          ###Newton Raphson iterations to determine the F that would result in the quota being removed from the entire true population; needed due to lag between assess. yr and quotatarg yr (see above)
          if (quotatarg[i-1,1]==0) { Fnew<-0}
          else {
            Fnew<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i-1,1])
          }
          ###End Newton Raphson
        } else {
          exploitAA[i-1,]=0
          ####Set a quota. Note that quota in year i is based on the assessment in year i-1.  This creates a lag in the assessment and the quota.  The lag could be altered if desired.  Quota in biomass set by managers based on target F and terminal year assessed abundance; either adjusted or unadjusted for Mohn's Rho depending on AdjSwitch.
          ifelse (AdjSwitch==0, quotatarg[i,1]<-sum(exploitAA[i-1,]*((NretroTermAA[i-1,]*Wt)-BcutPT)),quotatarg[i,1]<-sum(exploitAA[i-1,]*(((NretroTermAA[i-1,]*Wt)/(MohnsRho[i-1,1]+1))-BcutPT)))
          ####End quota setting.
          Fnew<-0
          Ftarg<-0
          FAA[i-1,]=Ftarg*Fsel
        }
      } else if (CtrlRule==3) { #constant catch
        quotatarg[i,1]<-Propc*MSY #set quota; same each year
        #newton raphsons to find F values
        if (quotatarg[i-1,1]==0) { Ftarg<-0}
        else {
          Ftarg<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NretroTermAA[i-1,],quotatarg=quotatarg[i,1])
        }
        ###End Newton Raphson
        if(useFcap>0){
        if(Ftarg>Fcap){  ###conditional constant catch; is target F > than desired max F?
          Ftarg<-Fcap #if yes, fish at Fcap, much like any F based control rule
          if (TotN[i-1,1]==0) Ftarg<-0 #Assumes that no fishing occurs after extinction.  Without this statement the Newton Raphson iterations explode.          
          FAA[i-1,]=Ftarg*Fsel
          ifelse (AdjSwitch==0, quotatarg[i,1]<-sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*NretroTermAA[i-1,]*Wt*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))),quotatarg[i,1]<-sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*((NretroTermAA[i-1,]*Wt)/(MohnsRho[i-1,1]+1))*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))))
          ###Newton Raphson iterations to determine the F that would result in the quota being removed from the true population; needed due to lag between assess. yr and quotatarg yr (see above)
          if (Ftarg==0) { Fnew<-0}
          else {
            Fnew<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i-1,1])
          }
          ###End Newton Raphson
        } else {        
        ###Newton Raphson iterations to determine the F that would result in the quota being removed from the entire true population; needed due to lag between assess. yr and quotatarg yr (see above)
        if (quotatarg[i-1,1]==0) { Fnew<-0}
        else {
          Fnew<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i-1,1])
        }
        } #end Fcap if else
        } else {
          ###Newton Raphson iterations to determine the F that would result in the quota being removed from the entire true population; needed due to lag between assess. yr and quotatarg yr (see above)
          if (quotatarg[i-1,1]==0) { Fnew<-0}
          else {
            Fnew<-Solve.F(Fnew=0.1,Fsel=Fsel,M=TrueM[i-1,],Wt=Wt,NAA=NAA[i-1,],quotatarg=quotatarg[i-1,1])
          }  
        } #end if useFcap if else
      } else { print ("Poop: unknown control rule") }
      ##End Control Rules
      
      if(SSB[i-1,1]<OverfishedSSB) {OverfishedSSBSimYear[s,i-1]=1} else {OverfishedSSBSimYear[s,i-1]=0} #1 if year i in sim s is overfished, else 0.  Summed later for % overfished.
      if(quotatarg[i-1,1]==0) {ClosureSimYear[s,i-1]=1} else {ClosureSimYear[s,i-1]=0} #1 if year i in sim s has closed fishery, else 0.  Summed later for %Closure.
      
      #Line below: Percent difference in quota between adjusting for retro and not adjusting.  Equals: ((Qadj-Qunadj)/Qunadj)*100
      if (TotN[i-1,1]>0) {QuotaAdjDiff[i-1,1]=((sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*((NretroTermAA[i-1,]*Wt)/(MohnsRho[i-1,1]+1))*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,]))))-sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*NretroTermAA[i-1,]*Wt*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))))/sum((FAA[i-1,]/(FAA[i-1,]+TrueM[i-1,]))*NretroTermAA[i-1,]*Wt*(1-exp(-1*(FAA[i-1,]+TrueM[i-1,])))))*100 }
      else { QuotaAdjDiff[i-1,1]=0 }
      
      TrueF[i-1,1]=Fnew*exp(Imperr-((Impsd*Impsd)/2))   #Application of Implementation error
      if (Ftarg>0) {RealTargFDiff[i-1,1]=((TrueF[i-1,1]-Ftarg)/Ftarg)*100} else {RealTargFDiff[i-1,1]=0} #percent difference between True and Target F in each year
      TotCatch[i-1,1]=sum(((TrueF[i-1,1]*Fsel)/(TrueF[i-1,1]*Fsel+TrueM[i-1,]))*NAA[i-1,]*(1-exp(-1*((TrueF[i-1,1]*Fsel)+TrueM[i-1,]))))  #Total Catch (millions of fish) each year
      TotYield[i-1,1]=sum(((TrueF[i-1,1]*Fsel)/(TrueF[i-1,1]*Fsel+TrueM[i-1,]))*NAA[i-1,]*Wt*(1-exp(-1*((TrueF[i-1,1]*Fsel)+TrueM[i-1,]))))  #Total Yield each year
      TotNatDeath[i-1,1]=sum(((TrueM[i-1,])/(TrueF[i-1,1]*Fsel+TrueM[i-1,]))*NAA[i-1,]*Wt*(1-exp(-1*((TrueF[i-1,1]*Fsel)+TrueM[i-1,]))))  #Total Natural Deaths each year
      if (i>3) YieldVar[i-1,1]=(TotYield[i-1,1]-TotYield[i-2,1])^2 #squared difference between yield this year and yield last year; Year must >3 to accomodate i-2 term.
      ZAA[i-1,]=(TrueF[i-1,1]*Fsel)+TrueM[i-1,]   #calculate total mortality at age (Z) as sum of F at age and M
      for (j in 1:nages)  #age loop
      {
        if (j==1) { NAA[i,j]<-stock.recruit(SSB1=SSB[i-1,1],Reca=Reca,Recb=Recb,type=RecType)*exp(Recerrauto[i,1]-((Recsd*Recsd)/2)) }  #Recruitment at age 1 with bias corrected error
        else { if (j<nages) { NAA[i,j]=NAA[i-1,j-1]*exp(-ZAA[i-1,j-1]) }  #Fill abundance at age matrix on the diagonals
               else { if (j==nages) { NAA[i,j]=(NAA[i-1,j-1]*exp(-ZAA[i-1,j-1]))+(NAA[i-1,nages]*exp(-ZAA[i-1,nages])) }}}  #Plus Group
        NAA[i,j]=round(NAA[i,j],digits=6) #Rounds N at age values to six digits to accomodate units of millions.  So, when abundance is less than 1, it is set equal to zero.
        SSBAA[i,j]=NAA[i,j]*Mat[j]*Wt[j] #SSB at age in each year
      } #end age loop
      
      Nretro[i,i]=(sum(NAA[i,])*(mohnsrhoterm+1))*exp(NerrRan[i,1]-((Nerrsd*Nerrsd)/2)) #Terminal year estimates for each "assessment" with bias corrected error
      NretroTermAA[i,]=NAA[i,]*(mohnsrhoterm+1)*exp(NerrRan[i,1]-((Nerrsd*Nerrsd)/2))  #Terminal year estimates at age for each assessment.  Needed for quota setting (see quotatarg).  These calculations summed among ages equals applying mohn's rho to the total true N.  i.e., sum of products equals product with sum.
      SSBretro[i,i]=(sum(NAA[i,]*Mat*Wt)*(mohnsrhoterm+1))*exp(NerrRan[i,1]-((Nerrsd*Nerrsd)/2))  #Terminal year estimates of SSB for each "assessment"
      for (k in 1:(i-1)) #loop for retrospective "fan".  Applies to all years of "assessment" except terminal year
      {
        Nretro[i-k,i]=((Nretro[i-k,i-1]/(mohnsrhofan[k]+1)))*exp(NerrFan-((NerrFansd*NerrFansd)/2))
        SSBretro[i-k,i]=((SSBretro[i-k,i-1]/(mohnsrhofan[k]+1)))*exp(NerrFan-((NerrFansd*NerrFansd)/2))
      }
      
      TotBio[i,1]=sum(NAA[i,]*Wt) #True total biomass each year
      SSB[i,1]=sum(SSBAA[i,]) #Calculate true total SSB in each year (sum across ages)
      if (SSB[i,1]<1) SSB[i,1]=0 #Once SSB is less 1 kmt (units of SSB are already in 000's of mt), SSB is set to zero and the stock collapses.
      TotN[i,1]=sum(NAA[i,])  #Total true abundance in each year (sum across ages)
      
      #Loop for calculating Mohn's Rho as assessment scientists would do using assessed SSB.
      if (i>peels) {
        for (p in 1:peels)
        {
          MohnsRhoStep[p,1]=(SSBretro[i-p,i-p]-SSBretro[i-p,i])/SSBretro[i-p,i]
        }
        MohnsRho[i,1]=mean(MohnsRhoStep)
      }
      ##End loop for calculating Mohn's Rho as assessment scientists would.
      
    } #end year loop
    
    #If-Else below writes years specific results for each simulation into a table (each table nsims x nyears big).  Coded hear to cut down on run time.  Doing in loops above slowed sims.
    if(s==1)
    {
      ctrlrulenum<-paste(u,l,f,c,sep="")
      TotBiob<-data.frame(TotBio,s)
      write.table(TotBiob,paste(readfile,paste(ctrlrulenum,"TotBioSimYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("TotalBio","Sim"),quote=TRUE) #Total biomass in each year of each sim
      SSBb<-data.frame(SSB,s,Bzero,Bmsy)
      write.table(SSBb,paste(readfile,paste(ctrlrulenum,"SSBSimYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("SSB","Sim","Bzero","Bmsy"),quote=TRUE) #SSB in each year of each sim
      TotYieldb<-data.frame(TotYield,TotNatDeath,quotatarg,s,MSY)
      write.table(TotYieldb,paste(readfile,paste(ctrlrulenum,"TotYieldSimYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("Yield","NatDeaths","TargetQuota","Sim","MSY"),quote=TRUE) #TotYield in each year of each sim
      QuotaAdjDiffb<-data.frame(QuotaAdjDiff,s)
      write.table(QuotaAdjDiffb,paste(readfile,paste(ctrlrulenum,"QuotaDiffByYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("QuotaDiffw_wo_rhoadj","Sim"),quote=TRUE) #Percent change in quota with adjustment for each year of each sim
      RealTargFDiffb<-data.frame(RealTargFDiff,s)
      write.table(RealTargFDiffb,paste(readfile,paste(ctrlrulenum,"RealTargFDiffByYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("RealTargFdiff","Sim"),quote=TRUE) #Percent difference between real and target F for each year of each sim
      write.table(array(c(s,nsims,nyears,nages,M,Mrho,Msd,SSBfrac,RecType,steepness,Bzero,Recsd,Recrho,CtrlRule,FracBmsyThreshHi,FracBmsyThreshLo,FracFtarg,useminmax,minquota,maxquota,usequotavar,quotavar,Thresh,Prop,Propc,useFcap,Fcapprop,quotablock,NerrRho,Nerrsd,Impsd,AdjSwitch,mohnsrhoterm),dim=(c(1,33))),paste(readfile,paste(ctrlrulenum,"SimCharacteristics.txt",sep=""),sep=""),row.names=FALSE,col.names=c("sim","nsims","nyear","nages","M","Mrho","Msd","SSBfrac","RecType","steep","Bzero","RecSD","RecAutoCorr","CtrlRule","FracBmsyThreshHi","FracBmsyThreshLo","FracFtarg","minmaxquota?","minquota","maxquota","LimQuotaVar?","AnnQuotaVarFrac","Thresh","Prop","ConCatchPropMSY","CondtlCatch?","MaxF_fracFmsy","quotablock","AsserrRho","AsserrSD","ImperrSD","AdjSwitch","TermBias"),quote=FALSE)     #sim characteristics
      SimCharAA<-data.frame(rep(s,nages),Msel,Mat,Wt,Fsel,NAAzero)
      write.table(SimCharAA,paste(readfile,paste(ctrlrulenum,"SimCharAA.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("Sim","Mselectivity","Maturity","Weight","Fselectivity","UnfishedNAA"),quote=TRUE) #NAA in each year of each sim
      NAAb<-data.frame(NAA,s)
      write.table(NAAb,paste(readfile,paste(ctrlrulenum,"NAASimYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c(paste("Age",seq(1,nages,1),sep=""),"Sim"),quote=TRUE) #NAA in each year of each sim
      TrueFb<-data.frame(TrueF,s,Fmsy)
      write.table(TrueFb,paste(readfile,paste(ctrlrulenum,"TrueFSimYear.txt",sep=""),sep=""),append=FALSE,row.names=FALSE,col.names=c("TrueFullF","Sim","Fmsy"),quote=TRUE) #TrueF in each year of each sim      
    } else {
      ctrlrulenum<-paste(u,l,f,c,sep="")
      TotBiob<-data.frame(TotBio,s)
      write.table(TotBiob,paste(readfile,paste(ctrlrulenum,"TotBioSimYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #Total biomass in each year of each sim
      SSBb<-data.frame(SSB,s,Bzero,Bmsy)
      write.table(SSBb,paste(readfile,paste(ctrlrulenum,"SSBSimYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #SSB in each year of each sim
      TotYieldb<-data.frame(TotYield,TotNatDeath,quotatarg,s,MSY)
      write.table(TotYieldb,paste(readfile,paste(ctrlrulenum,"TotYieldSimYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #TotYield in each year of each sim
      QuotaAdjDiffb<-data.frame(QuotaAdjDiff,s)
      write.table(QuotaAdjDiffb,paste(readfile,paste(ctrlrulenum,"QuotaDiffByYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #Percent change in quota with adjustment for each year of each sim
      RealTargFDiffb<-data.frame(RealTargFDiff,s)
      write.table(RealTargFDiffb,paste(readfile,paste(ctrlrulenum,"RealTargFDiffByYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #Percent difference between real and target F for each year of each sim
      write.table(array(c(s,nsims,nyears,nages,M,Mrho,Msd,SSBfrac,RecType,steepness,Bzero,Recsd,Recrho,CtrlRule,FracBmsyThreshHi,FracBmsyThreshLo,FracFtarg,useminmax,minquota,maxquota,usequotavar,quotavar,Thresh,Prop,Propc,useFcap,Fcapprop,quotablock,NerrRho,Nerrsd,Impsd,AdjSwitch,mohnsrhoterm),dim=(c(1,33))),paste(readfile,paste(ctrlrulenum,"SimCharacteristics.txt",sep=""),sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)   #sim characteristics
      SimCharAA<-data.frame(rep(s,nages),Msel,Mat,Wt,Fsel,NAAzero)
      write.table(SimCharAA,paste(readfile,paste(ctrlrulenum,"SimCharAA.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #NAA in each year of each sim
      NAAb<-data.frame(NAA,s)
      write.table(NAAb,paste(readfile,paste(ctrlrulenum,"NAASimYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #NAA in each year of each sim
      TrueFb<-data.frame(TrueF,s,Fmsy)
      write.table(TrueFb,paste(readfile,paste(ctrlrulenum,"TrueFSimYear.txt",sep=""),sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=TRUE) #TrueF in each year of each sim
    }
    
  } #######End Sims loop
  
  write.table(OverfishedSSBSimYear,paste(readfile,paste(ctrlrulenum,"OverfishedSSBSimYear.txt",sep=""),sep=""),col.names=c(paste("Year",seq(1,nyears,1),sep="-"))) #write table for whether SSB overfished in each year and sim (rows are sims; cols are years)
  write.table(ClosureSimYear,paste(readfile,paste(ctrlrulenum,"ClosureSimYear.txt",sep=""),sep=""),col.names=c(paste("Year",seq(1,nyears,1),sep="-"))) #write table for whether fishery is close in each year and sim (rows are sims; cols are years)
  
  #If the switch is turned on, this if else statement records the %difference in SSB for each year and simulation between adj and unadj sims
  if (SSBSimYearDiffSwitch!=0)
  { SSBSimYearAdj<-read.table(paste(readfile2,"AdjSSBSimYear.txt",sep=""),header=F)
    SSBSimYearUnadj<-read.table(paste(readfile2,"UnadjSSBSimYear.txt",sep=""),header=F)
    SSBSimYearDiff=((SSBSimYearAdj-SSBSimYearUnadj)/SSBSimYearUnadj)*100
  } ###End if
  
  if(FALSE){ ##not needed for herring MSE
    #loop records the median percent difference of various metrics between adjusting or not adjusting for the retro for each year among all sims
    QuotaDiffByYear<-read.table(paste(readfile,"QuotaDiffByYear.txt",sep=""),header=F)
    RealTargFDiffByYear<-read.table(paste(readfile,"RealTargFDiffByYear.txt",sep=""),header=F)
    for (i in 2:nyears)
    {
      QuotaDiffByYearMed[i,1]=median(QuotaDiffByYear[,i])
      QuotaDiffByYearMed[i,2]=i
      RealTargFDiffByYearMed[i,1]=median(RealTargFDiffByYear[,i])
      RealTargFDiffByYearMed[i,2]=i
      RealTargFDiffByYear25[i,1]=quantile(RealTargFDiffByYear[,i],probs=0.25,na.rm=TRUE)
      RealTargFDiffByYear25[i,2]=i
      QuotaDiffByYear25[i,1]=quantile(QuotaDiffByYear[,i],probs=0.25,na.rm=TRUE)
      QuotaDiffByYear25[i,2]=i
      RealTargFDiffByYear75[i,1]=quantile(RealTargFDiffByYear[,i],probs=0.75,na.rm=TRUE)
      RealTargFDiffByYear75[i,2]=i
      QuotaDiffByYear75[i,1]=quantile(QuotaDiffByYear[,i],probs=0.75,na.rm=TRUE)
      QuotaDiffByYear75[i,2]=i
      OverfishedSSBByYear[i,1]=(sum(OverfishedSSBSimYear[,i])/nsims)*100
      OverfishedSSBByYear[i,2]=i
      if(SSBSimYearDiffSwitch!=0) {SSBSimYearDiffMed[i,1]=median(SSBSimYearDiff[,i])
                                   SSBSimYearDiffMed[i,2]=i
                                   SSBSimYearDiff25[i,1]=quantile(SSBSimYearDiff[,i],probs=0.25,na.rm=TRUE)
                                   SSBSimYearDiff25[i,2]=i
                                   SSBSimYearDiff75[i,1]=quantile(SSBSimYearDiff[,i],probs=0.75,na.rm=TRUE)
                                   SSBSimYearDiff75[i,2]=i
      } #end if
    }     #end loop for recording median percent
  } #close if false
} #end sims function