#script for visualizing summary herring MSE results August 2016

preysum<-read.table("UnadjSummaryResults.txt", header=T)
names(preysum)
#[1] "MedianSSB"            "MedSSBrelSSBmsy"      "MedSSBrelSSBzero"     "PropSSBrelSSBmsy"     "PropSSBrelhalfSSBmsy" "PropSSBrel3SSBzero"   "PropSSBrel75SSBzero" 
#[8] "PropFrelFmsy"         "Yield"                "Yvar"                 "NatDeaths"            "YieldrelMSY"          "YieldrelNatDeath"     "PropClosure"         
#[15] "NAge1"                "NAge2"                "NAge3"                "NAge4"                "NAge5"                "NAge6"                "NAge7"               
#[22] "NAge8"                "PropAgeOne"           "PropAgeTwo"           "PropAgeOneTwo"        "SurpProd"             "CtrlRule"             "FracBmsyThreshHi"    
#[29] "FracBmsyThreshLo"     "FracFtarg"            "useminmax"            "minquota"             "maxquota"             "usequotavar"          "quotavar"            
#[36] "Thresh"               "Prop"                 "Propc"                "useFcap"              "Fcapprop"             "quotablock"          

preysum25<-read.table("UnadjSummaryResults25.txt", header=T)
names(preysum25)


preysum75<-read.table("UnadjSummaryResults75.txt", header=T)
names(preysum75)

plot.new()
plot(preysum$MedianSSB)
plot(preysum$MedSSBrelSSBmsy)
plot(preysum$MedSSBrelSSBzero)

