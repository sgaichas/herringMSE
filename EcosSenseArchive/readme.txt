updates made to GOM model for herring MSE
November-December 2016
recalculated all marine mammal diets based on Laurel Smith et al.'s 2015 paper
see Data/Projects/MSE/MammalDat/Mammals_forDietUpdate.xls and GOM_EMAXdietmatrixupdate.csv
directly modified the header file for C code \
new file is GOM_mammaldietupdate.h

changing diets of course knocked balance off for small pelagics_other and small pelagics_squids
changed biomass of both as noted in Mammals_forDietUpdate.xls
changes PB for other pelagics too also as noted in above

calculated new EEs in Ecopath using altered model
entered new B, PB, and EE directly in header file

HUGE CHANGE required to get flatlines in dynamic model (although EEs are balanced)
set detritus fate for all things going into discards to 0 except fishery is set to 1
increase B of detritus to 100
all detritus fate for all things going to detritus is 1 aside from discards and detritus themselves
discards now way out of balance (17)
so nobody eats discards aside from seabirds
adjusted seabird diet so that discards are <1%, benthivores 2% (hake) and juvenile fish 10% (hake again)


ran all new scenarios with this file.
                                                                    
