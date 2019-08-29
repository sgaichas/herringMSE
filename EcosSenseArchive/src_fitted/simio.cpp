// -----------------------------------------------------------------------------

int rate_output(struct SimRun *v){

double ssq;

int i, link, sp;
//float *fit_vector;
char outFile[80];					   
FILE *dmpptr;
const char* w = "w";
int rlookup[NUM_GROUPS+1][NUM_GROUPS+1];

           ssq = v->fit_fittot, v->fit_forcetot, v->fit_diettot;

           sprintf(outFile,"%s_rates.csv",OUTFOLDER);
           dmpptr=fopen(outFile,w);
           
           fprintf(dmpptr,"Seed,System,Species,Year,ssqCheck,");
           fprintf(dmpptr,"PB,QB,F,M2,Mzero,Mtot,Force,B,");
           for (link=1; link<=NUM_LIVING; link++){
               //if (v->RRR.rpar_PreyFrom[link] == sp){
                   fprintf(dmpptr,"Mby_%s,",path_species[link]);
               //}
           }
           fprintf(dmpptr,"\n");    
           memset(rlookup,0,(NUM_GROUPS+1)*(NUM_GROUPS+1)*sizeof(int));
           for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
              rlookup[v->RRR.rpar_PreyTo[link]][v->RRR.rpar_PreyFrom[link]]=link;
           }
             int yy;
             for (yy=0; yy<=fit_Years; yy++){
             for (sp=1; sp<=NUM_LIVING; sp++){
                 //fprintf(dmpptr[sp],"Seed,System,Species,Year,ssqOrig,ssqCheck,");
                   fprintf(dmpptr,"%d,%d,%s,%d,%g,",v->thread_ID,0,path_species[sp],
                                                         yy+fit_StartYear,ssq);
                 //fprintf(dmpptr[sp],"F,Mtot,Mzero,QB,PB");
                   fprintf(dmpptr,"%g,%g,%g,%g,%g,%g,%g,%g,",
                                                  v->out_PB[sp][yy],
                                                  v->out_QB[sp][yy],
                                                  v->out_TotF[sp][yy],
                                                  v->out_M2[sp][yy],
                                                  v->RRR.rpar_MzeroMort[sp],
                                                  v->out_MM[sp][yy],
                                                  v->fit_bioanom[sp][yy]/v->out_BB[sp][yy],
                                                  v->out_BB[sp][yy]);


                 for (link=1; link<=NUM_LIVING; link++){
                     //if (v->RRR.rpar_PreyFrom[link] == sp){
                         fprintf(dmpptr,"%g,",v->out_LinksM[rlookup[link][sp]][yy]);
                     //}
                 }

                 fprintf(dmpptr,"\n");
              }
              }
           
           
           fclose(dmpptr);
                
}
// -----------------------------------------------------------------------------
#define HEXOUT(ff) fp1=(unsigned int *)&(ff);sprintf(hex,"%X%X",fp1[0],fp1[4]);\
fprintf(fptr,"%s,%d,%d,%s,\n",mID,i,j,hex);
#define NAMEIT(ss) sprintf(mID,"%s",ss);
//----------------------------------------------------------------------
void giant_dump(struct SimRun *v, char* fbase){

    const char* w = "w";
	  FILE *fptr;   
    char outFile[80];
    char hex[80];
    char mID[80];
    int i,j;
    unsigned int *fp1;
		    
		         sprintf(outFile,"%s%s_giant.csv",OUTFOLDER,fbase);
             fptr=fopen(outFile,w);

//struct SimRun{
//
//    struct randvec rng;   // The random number state holder for this run
//    int    thread_ID;     // An ID number for convenience
fprintf(fptr,"thread_ID,%d,\n",v->thread_ID);
//    FILE  *dump;          // A File to dump results to
//    struct RatePar RRR;   // The RATE PARAMETERS for the run
//    
//    //float fit_EST[MAX_COLS+1][MAX_YEARS+1];
//    float *fit_vector;
//    float **fit_EST;
//    float **fit_diet_EST; //[NUM_PREDPREYLINKS+1][MAX_YEARS+1];
//    //float fit_SSQ[MAX_COLS+1];
//    float **fit_SSQ;
//    float *fit_diet_SSQ;
//    float **fit_bioanom;
//
//    float  *fit_vec_control;
//    float **fit_ssq_control;
//    //float fit_bioanom[NUM_GROUPS+1][MAX_YEARS+1];
//    //float 
//    //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
//    float **out_BB;
//    float **out_CC;
//    float **out_RR;
//    float **out_EE;
//    float **out_SSB;
//    //float out_CC[NUM_GROUPS+1][MAX_YEARS+1];
//    float **out_MM; //float out_MM[NUM_GROUPS+1][MAX_YEARS+1];
//    float **out_M2; //float out_MM[NUM_GROUPS+1][MAX_YEARS+1];
//    float **out_PB;  //out_PB[NUM_GROUPS+1];
//    float **out_QB;  //out_QB[NUM_GROUPS+1];
//    float **out_LinksM; //out_LinksM[NUM_PREDPREYLINKS+1];
//    float **out_TotF;   //out_TotF[NUM_GROUPS+1];
//    float out_PP;
//
//    float **force_bymort; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
//
//    float  fit_SSQraw[NUM_GROUPS+MAX_COLS+1];
//    float  fit_SSQfit[NUM_GROUPS+MAX_COLS+1];
//    float  fit_SSQwt[NUM_GROUPS+MAX_COLS+1];
//    float  fit_FORCEraw[NUM_GROUPS+1];    
//
//    float state_BB[NUM_GROUPS+1];
NAMEIT("state_BB")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->state_BB[j])}

//    float state_Ftime[NUM_GROUPS+1];
//    //float NageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
//    //float WageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];  
//    float **NageS, **WageS;
//    float state_NN[NUM_GROUPS+1];
//    float stanzaPred[NUM_GROUPS+1];
//    float stanzaBasePred[NUM_GROUPS+1];
//    
//		float vBM[MAX_SPLIT+1];
//    float recruits[MAX_SPLIT+1];
//		float RzeroS[MAX_SPLIT+1];
//    //float WWa[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
//    float **WWa;
//    float stanzaGGJuv[MAX_SPLIT+1];
//    float stanzaGGAdu[MAX_SPLIT+1];
//    float baseEggsStanza[MAX_SPLIT+1];
//    float EggsStanza[MAX_SPLIT+1];
//    float SpawnBio[MAX_SPLIT+1];
//    float baseSpawnBio[MAX_SPLIT+1];
//    float Rbase[MAX_SPLIT+1];
//    
//    //float baseSpawnEnergy[MAX_SPLIT+1];
//    
////  hardcoding only two stanzas for now; these dimension them by months
//    int firstMoJuv[MAX_SPLIT+1], lastMoJuv[MAX_SPLIT+1]; 
//    int firstMoAdu[MAX_SPLIT+1], lastMoAdu[MAX_SPLIT+1];
//    //float SplitAlpha[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
//    float **SplitAlpha;
//    float RscaleSplit[MAX_SPLIT+1];  //not used right now
//
//// OUTPUTS
//   //float out_BB[NUM_GROUPS+1][MAX_YEARS*STEPS_PER_YEAR+1];
//   //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
//   //float out_nn[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
//   //float out_ww[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
//   //float out_CC[NUM_GROUPS+1][MAX_YEARS+1];
//   
//   float fish_Effort[NUM_GROUPS+1];
//   float TerminalF[NUM_GROUPS+1];
//
//// FOR RECORDING SENSE RESULTS
//   int DISCARD_YEAR;
//	 int discarded[NUM_GROUPS+1];
//	 //float month_noise[NUM_GROUPS + 1];
//
//// Variables prefaced by deriv_ are derivative parts (e.g. current consumption)
//   float deriv_TotGain[NUM_GROUPS + 1];
NAMEIT("deriv_TotGain")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_TotGain[j])}
//   float deriv_TotLoss[NUM_GROUPS + 1];
NAMEIT("deriv_TotLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_TotLoss[j])}
//   float deriv_LossPropToB[NUM_GROUPS + 1];
NAMEIT("deriv_LossPropToB")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_LossPropToB[j])}
//   float deriv_LossPropToQ[NUM_GROUPS + 1];
NAMEIT("deriv_LossPropToQ")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_LossPropToQ[j])}
//   float deriv_ConsMat[NUM_GROUPS + 1][NUM_GROUPS + 1];
//   float deriv_DerivT[NUM_GROUPS + 1];
NAMEIT("deriv_DerivT")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_DerivT[j])}
//   float deriv_dyt[NUM_GROUPS + 1];
NAMEIT("deriv_dyt")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_dyt[j])}
//   float deriv_biomeq[NUM_GROUPS + 1];
NAMEIT("deriv_biomeq")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_biomeq[j])}
//
//   float deriv_FoodGain[NUM_GROUPS + 1];
NAMEIT("deriv_FoodGain")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_FoodGain[j])}
//   float deriv_DetritalGain[NUM_GROUPS + 1];
NAMEIT("deriv_DetritalGain")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_DetritalGain[j])}
//   float deriv_FoodLoss[NUM_GROUPS + 1];
NAMEIT("deriv_FoodLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_FoodLoss[j])}
//   float deriv_UnAssimLoss[NUM_GROUPS + 1];
NAMEIT("deriv_UnAssimLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_UnAssimLoss[j])}
//   float deriv_ActiveRespLoss[NUM_GROUPS + 1];
NAMEIT("deriv_ActiveRespLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_ActiveRespLoss[j])}
//   //float deriv_PassiveRespLoss[NUM_GROUPS + 1];
//NAMEIT("deriv_PassiveRespLoss")
//   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_PassiveRespLoss[j])}
//   float deriv_MzeroLoss[NUM_GROUPS + 1];
NAMEIT("deriv_MzeroLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_MzeroLoss[j])}
//   float deriv_DetritalLoss[NUM_GROUPS + 1];
NAMEIT("deriv_DetritalLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_DetritalLoss[j])}
//   float deriv_TotDetOut[NUM_GROUPS + 1];
NAMEIT("deriv_TotDetOut")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_TotDetOut[j])}   
//   float deriv_FishingLoss[NUM_GROUPS + 1];
NAMEIT("deriv_FishingLoss")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_FishingLoss[j])}
//   float deriv_FishingGain[NUM_GROUPS + 1];
NAMEIT("deriv_FishingGain")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_FishingGain[j])}
//   float deriv_FishingThru[NUM_GROUPS + 1];
NAMEIT("deriv_FishingThru")
   i=0; for (j=0; j<=NUM_GROUPS; j++){HEXOUT(v->deriv_FishingThru[j])}
//}; 

//    v->out_BB      = matrix(0,NUM_GROUPS,0,fit_Years);
//    v->out_CC      = matrix(0,NUM_GROUPS,0,fit_Years);
//    v->out_RR      = matrix(0,juv_N,0,fit_Years);
//    v->out_EE      = matrix(0,juv_N,0,fit_Years);
//    v->out_SSB     = matrix(0,juv_N,0,fit_Years);
//
//    v->out_MM      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
//    v->out_M2      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
//    v->out_PB      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
//    v->out_QB      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
//    v->out_LinksM  = matrix(0,NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
//    v->out_TotF    = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
//
//    v->fit_vector  = vector(0,NDIM);
//    v->fit_EST     = matrix(0,fitN,0,fit_Years);
//    v->fit_diet_EST= matrix(0,NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
//    v->fit_diet_SSQ= vector(0,NUM_PREDPREYLINKS);
//    v->fit_SSQ     = matrix(0,4,0,NUM_GROUPS);
//    v->fit_bioanom = matrix(0,NUM_GROUPS,0,fit_Years);
//
//    v->fit_vec_control = vector(0,NDIM);
//    v->fit_ssq_control = matrix(0,4,0,NUM_GROUPS);
//
//    v->NageS       = matrix(0,juv_N,0,MAX_MONTHS_STANZA);
//    v->WageS       = matrix(0,juv_N,0,MAX_MONTHS_STANZA);
//    v->WWa         = matrix(0,juv_N,0,MAX_MONTHS_STANZA);
//    v->SplitAlpha  = matrix(0,juv_N,0,MAX_MONTHS_STANZA);
//    
//    v->force_bymort = matrix(0,NUM_GROUPS,0,MAX_YEARS*STEPS_PER_YEAR);
//
//     for (sp=0; sp<=NUM_GROUPS; sp++){
//		     for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++){
//				  //force_byprey[sp][t]=1.0;
//				  v->force_bymort[sp][t]=force_bymort[sp][t];
//				  //force_byrecs[sp][t]=1.0;
//				 }
//		 }

fclose(fptr);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int output_run(struct SimRun *v, char* fbase){
    int t,sp,i,s, ageMo,gr;
    int pred,prey;
    float est,obs,sd,raw,craw;
    const char* w = "w";
	  FILE *fptr;   
    char outFile[80];
    char tmpstr[80];
    
		// Save biomass to a file
		         //sprintf(outFile,"%s/%s_bio.csv",OUTFOLDER,fbase);
		         sprintf(outFile,"%s%s_bio.csv",OUTFOLDER,fbase);
             fptr=fopen(outFile,w);
             fprintf(fptr,"Year,");
             for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                 if(v->discarded[sp]){fprintf(fptr,"dd%s,",path_species[sp]);}
                 else             {fprintf(fptr,"%s,",path_species[sp]);}
						 }   fprintf(fptr,"\n");
             //for (t=0; t<=fit_Years*STEPS_PER_YEAR; t+=12){
             for (t=0; t<=fit_Years; t++){
                 fprintf(fptr,"%d,",t + fit_StartYear);
                 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                     fprintf(fptr,"%g,",v->out_BB[sp][t]);}   fprintf(fptr,"\n");
						 }
						 fclose(fptr);
/*
		// Save catch to a file
		         //sprintf(outFile,"%s/%s_bio.csv",OUTFOLDER,fbase);
		         sprintf(outFile,"%s%s_cat.csv",OUTFOLDER,fbase);
             fptr=fopen(outFile,w);
             fprintf(fptr,"Year,");
             for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                 if(v->discarded[sp]){fprintf(fptr,"dd%s,",path_species[sp]);}
                 else             {fprintf(fptr,"%s,",path_species[sp]);}
						 }   fprintf(fptr,"\n");
             //for (t=0; t<=fit_Years*STEPS_PER_YEAR; t+=12){
             for (t=0; t<=fit_Years; t++){
                 fprintf(fptr,"%d,",t + fit_StartYear);
                 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                     fprintf(fptr,"%g,",v->out_CC[sp][t]);}   fprintf(fptr,"\n");
						 }
						 fclose(fptr);


		// Save mort to a file
		         //sprintf(outFile,"%s/%s_bio.csv",OUTFOLDER,fbase);
		         sprintf(outFile,"%s%s_mort.csv",OUTFOLDER,fbase);
             fptr=fopen(outFile,w);
             fprintf(fptr,"Year,");
             for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                 if(v->discarded[sp]){fprintf(fptr,"dd%s,",path_species[sp]);}
                 else             {fprintf(fptr,"%s,",path_species[sp]);}
						 }   fprintf(fptr,"\n");
             //for (t=0; t<=fit_Years*STEPS_PER_YEAR; t+=12){
             for (t=0; t<=fit_Years; t++){
                 fprintf(fptr,"%d,",t + fit_StartYear);
                 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                     fprintf(fptr,"%g,",v->out_MM[sp][t]);}   fprintf(fptr,"\n");
						 }
						 fclose(fptr);
   						     
   // Save juvenile/adults
          sprintf(outFile,"%s%s_rec.csv",OUTFOLDER,fbase);
          fptr=fopen(outFile,w);
          for (i=1; i<=juv_N; i++){
              fprintf(fptr,"%dSB,%dSE,%dEggs,%dNrec,%dWjuv,%dWadu,%dBjuv,%dBadu,",i,i,i,i,i,i,i,i);
              fprintf(fptr,"%dNjuv,%dNadu,%dFoodJuv,%dNFoodAdu,%dGGjuv,%dGGadu,%dPredJuv,%dPredAdu,",i,i,i,i,i,i,i,i);
             }
          fprintf(fptr,"\n");          
          for (t=0; t<=fit_Years*STEPS_PER_YEAR; t++){
					   for (i=1; i<=juv_N; i++){
					     fprintf(fptr,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,",
					      v->SpawnBio[i][t],
                out_SpawnEnergy[i][t],
                out_eggs[i][t],
                out_Nrec[i][t],
                out_WavgJuv[i][t],
                out_WavgAdu[i][t],
								out_BB[juv_JuvNum[i]][t],
								out_BB[juv_AduNum[i]][t],
								out_BB[juv_JuvNum[i]][t]/out_WavgJuv[i][t],
								out_BB[juv_AduNum[i]][t]/out_WavgAdu[i][t],
								out_foodgain[juv_JuvNum[i]][t],
								out_foodgain[juv_AduNum[i]][t],
                out_GGJuv[i][t],
                out_GGAdu[i][t],
                out_stanzaPredJuv[i][t],
                out_stanzaPredAdu[i][t]
								);
             }
             fprintf(fptr,"\n");
          }
          fclose(fptr);
          
    */
      // Save Juvenile/Recruits as minimal number
      /*
          sprintf(outFile,"%s%s_rec.csv",OUTFOLDER,fbase);
          fptr=fopen(outFile,w);
          for (i=1; i<=juv_N; i++){
              fprintf(fptr,"%dSB,%dSE,%dEggs,%dNrec,%dWjuv,%dWadu,%dNjuv,%dNadu,",i,i,i,i,i,i,i,i);
              //fprintf(fptr,"%dNjuv,%dNadu,%dFoodJuv,%dNFoodAdu,%dGGjuv,%dGGadu,%dPredJuv,%dPredAdu,",i,i,i,i,i,i,i,i);
             }
          fprintf(fptr,"\n");          
          for (t=11; t<=fit_Years*STEPS_PER_YEAR; t+=12){
					   for (i=1; i<=juv_N; i++){
					     fprintf(fptr,"%g,%g,%g,%g,%g,%g,%g,%g,",
					      out_SpawnBio[i][t],
                out_SpawnEnergy[i][t],
                out_eggs[i][t],
                out_Nrec[i][t],
                out_WavgJuv[i][t],
                out_WavgAdu[i][t],
								out_BB[juv_JuvNum[i]][t]/out_WavgJuv[i][t],
								out_BB[juv_AduNum[i]][t]/out_WavgAdu[i][t]
								);
             }
             fprintf(fptr,"\n");
          }
          fclose(fptr);          
          
    
		// Save age structure to a file (only for one species numbered in Adams-Basforth)
		         sprintf(outFile,"%s%s_age.csv",OUTFOLDER,fbase);
             fptr=fopen(outFile,w);
             for (t=0; t<=fit_Years*STEPS_PER_YEAR; t++)
                 fprintf(fptr,"%dN,%dW,",t,t);
             fprintf(fptr,"\n");
             //hard coded for POLLOCK species 3 in the stanza_in file
             i=3;
             for (ageMo=firstMoJuv[i]; ageMo<=lastMoAdu[i]; ageMo++){
                  for (t=0; t<=fit_Years*STEPS_PER_YEAR; t++){
						         fprintf(fptr,"%g,%g,",out_nn[ageMo][t],out_ww[ageMo][t]);
								     }
                  fprintf(fptr,"\n");
             }
						 fclose(fptr);   

    */
    // Output diet fitting
       if (dietN>0){
          sprintf(outFile,"%s%s_fit_diets.csv",OUTFOLDER,fbase);
				  fptr=fopen(outFile,w);
          fprintf(fptr,"N,Series,Type,Pred,Prey,Mult,Output,");
          fprintf(fptr,"RawFit,Wt,LastFit,");
        for (t=0; t<=fit_lastyear; t++){
            fprintf(fptr,"y%d,",t+fit_StartYear);
            }
        fprintf(fptr,"\n");     
				//cout << dietN << " yata!" << endl;     
        for (s=1; s<=dietN; s++){
            pred = int(fit_diet_lookup[s]/1000); 
            prey = fit_diet_lookup[s] - pred*1000; 
        
            if (prey<=NUM_LIVING){
               sprintf(tmpstr,"%s",path_species[prey]);
            }
            else{
               sprintf(tmpstr,"other");
            }    
            //cout << dietN << " ya!" << endl;
            // OBSERVED
               fprintf(fptr,"%d,%s,%d,%s,%s,%g,OBSERVED,",
                      s, fit_diet_names[s].c_str(),fit_diet_type[s],
                      path_species[pred],tmpstr,
                      0.0);
               fprintf (fptr,"%g,%g,%g,",
                       0.0,0.0,v->fit_diet_SSQ[s]
               );
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",fit_diet_OBS[s][t]);
               }
               fprintf(fptr,"\n");               
            
            // ESTIMATED
               fprintf(fptr,"%d,%s,%d,%s,%s,%g,ESTIMATED,",
                      s, fit_diet_names[s].c_str(),fit_diet_type[s],
                      path_species[pred],tmpstr,
                      0.0);
               fprintf (fptr,"%g,%g,%g,",
                       0.0,0.0,v->fit_diet_SSQ[s]
               );
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",v->fit_diet_EST[s][t]);
               }
               fprintf(fptr,"\n");            
        }

       fclose(fptr);   
       }
    // Output fitting data format #2
        sprintf(outFile,"%s%s_fit_new.csv",OUTFOLDER,fbase);
				fptr=fopen(outFile,w);
        fprintf(fptr,"N,Series,Type,ApplyTo,Mult,Output,");
        fprintf(fptr,"RawFit,Wt,LastFit,");
        for (t=0; t<=fit_lastyear; t++){
            fprintf(fptr,"y%d,",t+fit_StartYear);
            }
        fprintf(fptr,"\n");
        
        
        for (s=1; s<=NUM_LIVING+NUM_DEAD; s++){
        			 //est    = v->out_BB[s][t];
				       //cest   = v->out_CC[s][t];
        
               fprintf(fptr,"%d,%s,%d,%s,%g,%s,",
                      s,path_species[s],10,path_species[s],0.0,"ESTBIO");
               fprintf (fptr,"%g,%g,%g,",0.0,0.0,0.0);
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",v->out_BB[s][t]);
               }
               fprintf(fptr,"\n");				

               fprintf(fptr,"%d,%s,%d,%s,%g,%s,",
                      s,path_species[s],11,path_species[s],0.0,"ESTCAT");
               fprintf (fptr,"%g,%g,%g,",0.0,0.0,0.0);
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",v->out_CC[s][t]);
               }
               fprintf(fptr,"\n");
				}
				
        for (s=1; s<=fitN; s++){
            // OBSERVED
               fprintf(fptr,"%d,%s,%d,%s,%g,OBSERVED,",
                      s, fit_names[s].c_str(),fit_type[s],
                      path_species[fit_ApplyTo[s]],
                      fit_Norm[s]);
               fprintf (fptr,"%g,%g,%g,",
                       v->fit_SSQraw[s],v->fit_SSQwt[s],v->fit_SSQfit[s]
               );
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",fit_OBS[s][t]);
               }
               fprintf(fptr,"\n");               
            
            // ESTIMATED
               //fprintf(fptr,"%d,%s,%d,%s,%g,ESTIMATED,",
               //       s, fit_names[s].c_str(),fit_type[s],
               //       path_species[fit_ApplyTo[s]],
               //       0.0);
               //fprintf (fptr,"%g,%g,%g,",
               //        v->fit_SSQraw[s],v->fit_SSQwt[s],v->fit_SSQfit[s]
               //);
               //for (t=0; t<=fit_lastyear; t++){
               //    fprintf(fptr,"%g,",v->fit_EST[s][t]);
               //}
               //fprintf(fptr,"\n");
							 
						// SD	             
               fprintf(fptr,"%d,%s,%d,%s,%g,SD,",
                      s, fit_names[s].c_str(),fit_type[s],
                      path_species[fit_ApplyTo[s]],
                      0.0);
               fprintf (fptr,"%g,%g,%g,",
                       v->fit_SSQraw[s],v->fit_SSQwt[s],v->fit_SSQfit[s]
               );
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",fit_SD[s][t]);
               }
               fprintf(fptr,"\n"); 

        }

        // FORCED
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
               fprintf(fptr,"%d,FORCED,%d,%s,%g,FORCED,",
                      1000+sp, -1,
                      path_species[sp],
                      0.0);
               fprintf (fptr,"%g,%g,%g,",
                       v->fit_FORCEraw[sp],
                       v->fit_ssq_control[3][sp],
                       v->fit_SSQ[3][sp]
               );
               for (t=0; t<=fit_lastyear; t++){
                   fprintf(fptr,"%g,",v->fit_bioanom[sp][t]);
               }
               fprintf(fptr,"\n");
        }
           
        fclose(fptr);          

    // Output fitting data
    int link;
    int rlookup[NUM_GROUPS+1][NUM_GROUPS+1];
    float forced, cest, cobs,cpat,cattot,eggs,recs,ssb;
    if (!v->DISCARD_YEAR){
        sprintf(outFile,"%s%s_fit.csv",OUTFOLDER,fbase);
				fptr=fopen(outFile,w);
		    
		    fprintf (fptr,"Num,Species,Guild,Year,Est,Path,Obs,Raw,ForcedBio,ForcedBioRaw,SD,");
				fprintf (fptr,"Cest,CatPath,CatObs,CatRaw,ForcedCat,ForcedCatRaw,ForceAnom,Eggs,RecBio,SSB,");

        fprintf(fptr,"PB,QB,F,M2,Mzero,Mtot,Force,B,");
           for (link=1; link<=NUM_LIVING; link++){
                   fprintf(fptr,"Mby_%s,",path_species[link]);
           }

        fprintf (fptr,"CleanCatch,");
        for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
            fprintf(fptr,"%s,",path_species[gr]);
        }

        //fprintf(dmpptr,"\n"); 
        fprintf(fptr,"\n");

           memset(rlookup,0,(NUM_GROUPS+1)*(NUM_GROUPS+1)*sizeof(int));
           for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
              rlookup[v->RRR.rpar_PreyTo[link]][v->RRR.rpar_PreyFrom[link]]=link;
           }
        
		    int Bshown[NUM_GROUPS+1];
		    int Cshown[NUM_GROUPS+1];
		    int Rshown[NUM_GROUPS+1];
		 // Now Time Series         
        memset(Bshown,0,(NUM_GROUPS+1)*sizeof(int));   
        memset(Cshown,0,(NUM_GROUPS+1)*sizeof(int)); 
        memset(Rshown,0,(NUM_GROUPS+1)*sizeof(int)); 
		    for (s=1; s<=fitN; s++){
		        if (fit_type[s] == 0){
                Bshown[fit_ApplyTo[s]] = s;
            }
		        if (fit_type[s] == -100){
                Bshown[fit_ApplyTo[s]] = s;
            }
		        if (fit_type[s] == 6){
                Cshown[fit_ApplyTo[s]] = s;
            }            
        }
        for (i=1; i<=juv_N; i++){
             //cout << i << " " << juv_AduNum[i] << endl;
             Rshown[juv_AduNum[i]]=i;
        }
        


        float pat;

        float gearcatch[NUM_GROUPS+1][NUM_GROUPS+1];
        float caught;
        int links;
		    for (t=0; t<=fit_lastyear; t++){
        // Now apply effort to determine catch
           memset(gearcatch,0,(NUM_GROUPS+1)*(NUM_GROUPS+1)*sizeof(float));
					 //for (sp=0; sp<=NUM_GROUPS; sp++){
           //  for (gr=0; gr<=NUM_GROUPS; gr++){
           //      gearcatch[gr][sp]=0.0;
           //  }
           //}
           
           for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){
				       sp = v->RRR.rpar_FishingFrom[links];
				       gr   = v->RRR.rpar_FishingThrough[links];
				       //dest = v->RRR.rpar_FishingTo[links];
				       gearcatch[gr][sp] += v->RRR.rpar_FishingQ[links] * v->out_Effort[gr][t] * v->out_BB[sp][t]; 
		       }	


           for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
		        //for (t=0; t<=150; t++){
		          //fprintf(fptr,"%d,", t+fit_StartYear);
				      //est = v->out_BB[bioApplyTo[s]][t*STEPS_PER_YEAR+MEASURE_MONTH];
				      est    = v->out_BB[sp][t];
				      cest   = v->out_CC[sp][t];
				      forced =v->fit_bioanom[sp][t];
              				      
				      if (Bshown[sp]>0){ obs = fit_OBS[Bshown[sp]][t];
				                         raw = fit_RAW[Bshown[sp]][t];
                                 sd  = fit_SD[Bshown[sp]][t];}
				      else{              obs = 0.0;
							                   raw = 0.0; 
                                 sd=1.0;}
                                 
				      if (Cshown[sp]>0){ cobs = fit_OBS[Cshown[sp]][t];
							                   craw = fit_RAW[Cshown[sp]][t];
							}
				      else{              cobs = 0.0;
							                   craw = 0.0;
							}
                            
              if  ( ((t+fit_StartYear)>=PATH_YEAR_START) && ((t+fit_StartYear)<=PATH_YEAR_END)){
                        pat = path_BB[sp];
                        cpat = 0.0;
                        for (gr=0; gr<=NUM_GEARS; gr++){cpat+=path_Catch[gr][sp] + path_Discards[gr][sp];}                            
              }				     
              else{     pat  =0.0; cpat = 0.0;}
              
              if (Rshown[sp]>0) {eggs = v->out_EE[Rshown[sp]][t];  
                                 recs = v->out_RR[Rshown[sp]][t];
                                 ssb  = v->out_SSB[Rshown[sp]][t];
                                 }
              else{              eggs = 0.0;  recs = 0.0; ssb=0.0; }
              
              fprintf(fptr,"%d,%s,%s,%i,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,",
                     sp,path_species[sp],guildlist[PathNames[sp]].c_str(),
										 t+fit_StartYear,
                     est*AREA_MULT,pat*AREA_MULT,obs*AREA_MULT,raw,
										 FORCED_BIO[sp][t],FORCED_BIO_RAW[sp][t],sd,
										 cest*AREA_MULT,cpat*AREA_MULT,cobs*AREA_MULT,craw,
										 FORCED_CATCH[sp][t],FORCED_CATCH_RAW[sp][t],
										 forced,eggs,recs,ssb);

              // KYA Inserted lines to output production rates here 12/7/09
                   fprintf(fptr,"%g,%g,%g,%g,%g,%g,%g,%g,",
                                                  v->out_PB[sp][t],
                                                  v->out_QB[sp][t],
                                                  v->out_TotF[sp][t],
                                                  v->out_M2[sp][t],
                                                  v->RRR.rpar_MzeroMort[sp],
                                                  v->out_MM[sp][t],
                                                  v->fit_bioanom[sp][t]/v->out_BB[sp][t],
                                                  v->out_BB[sp][t]);


                   for (link=1; link<=NUM_LIVING; link++){
                             fprintf(fptr,"%g,",v->out_LinksM[rlookup[link][sp]][t]);
                   }
										 	 
              fprintf (fptr,"%g,",v->out_Effort[sp][t] * v->out_BB[sp][t]*AREA_MULT);
              for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
                  fprintf(fptr,"%g,",gearcatch[gr][sp]*AREA_MULT);
              }
              
              fprintf(fptr,"\n");       
				      //cout << fit_OBS[s][t] << " ";
				      //if (obs>0){fprintf(fptr,"%s,%s,%g,%i,%g,%g,%g,%g,\n",
				      //           fit_names[s].c_str(),path_species[fit_ApplyTo[s]],
				      //           v->fit_SSQ[s],t+fit_StartYear,
              //           obs,sd,est,forced);}
				      //else      {fprintf(fptr,"%s,%s,%g,%i,,,%g,%g,\n",
				      //           fit_names[s].c_str(),path_species[fit_ApplyTo[s]],
				      //           v->fit_SSQ[s],t+fit_StartYear,
              //           est,forced);
              //}
						}  
 						//cout << endl;
						//fprintf(fptr,"\n");
				}
					/*
				float cattot;
				for (s=1; s<=NUM_LIVING+NUM_DEAD; s++){
            if (!Bshown[s]){
                for (t=0; t<=fit_lastyear; t++){
                    est = bioMat[s][t];
                    if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
												  ((t+fit_StartYear)<=PATH_YEAR_END)){
                        obs = path_BB[s];    
                        }
                    else {obs=0.0;}
                    sd  = 1;                                     
				      if (obs>0){fprintf(fptr,"%s,%s,%g,%i,%g,%g,%g,\n",
				                 "PathBio",path_species[s],
				                 0.0,t+fit_StartYear,
                         obs,sd,est);}
				      else      {fprintf(fptr,"%s,%s,%g,%i,,,%g,\n",
				                 "PathBio",path_species[s],
				                 0.0,t+fit_StartYear,
                         est);                    
                
                }
            }
          }
          if (!Cshown[s]){
             cattot=0.0;
             for (gr=0; gr<=NUM_GEARS; gr++){
                  cattot+=path_Catch[gr][s] + path_Discards[gr][s];
             }
             if (cattot>0.0){
                for (t=0; t<=fit_lastyear; t++){
                    est = catMat[s][t];
                    if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
												  ((t+fit_StartYear)<=PATH_YEAR_END)){
                        obs = cattot;    
                        }
                    else {obs=0.0;}
                    sd  = 1;
				      if (obs>0){fprintf(fptr,"%s,%s,%g,%i,%g,%g,%g,\n",
				                 "PathCat",path_species[s],
				                 0.0,t+fit_StartYear,
                         obs,sd,est);}
				      else      {fprintf(fptr,"%s,%s,%g,%i,,,%g,\n",
				                 "PathCat",path_species[s],
				                 0.0,t+fit_StartYear,
                         est);                    
                
                }
            }             
             
             }  
          }
        
        
        }
				*/
        		    
/*
     // First line: titles
				fprintf(fptr,"Year,");
        for (s=1; s<=fitN; s++){
            //fprintf(fptr,"OBS%s,EST%s,sd,",fit_names[s],path_species[fit_ApplyTo[s]]);
            //fprintf(fptr,"OBS%s,EST%s,sd,",path_species[fit_ApplyTo[s]],path_species[fit_ApplyTo[s]]);
            fprintf(fptr,"OBS%s,EST%s,sd,",fit_names[s].c_str(),path_species[fit_ApplyTo[s]]);
        }   	
				fprintf(fptr,"\n");
				
		 // Second Line:  SSQs
        fprintf(fptr,"SSQ,");
        for (s=1; s<=fitN; s++){
              fprintf(fptr,"v%g,v%g,sd,",v->fit_SSQ[s],v->fit_SSQ[s]);
        }			
				fprintf(fptr,"\n");
    
		 // Now Time Series            
		    for (t=0; t<=fit_lastyear; t++){
		        fprintf(fptr,"%d,", t+fit_StartYear);
		        for (s=1; s<=fitN; s++){
				      //est = v->out_BB[bioApplyTo[s]][t*STEPS_PER_YEAR+MEASURE_MONTH];
				      est = v->fit_EST[s][t];
				      obs = fit_OBS[s][t];
				      sd  = fit_SD[s][t];
				      //cout << fit_OBS[s][t] << " ";
				      if (obs>0){fprintf(fptr,"%g,%g,%g,",obs,est,sd);}
				      else      {fprintf(fptr,",%g,,",est);}
						}  
 						//cout << endl;
						fprintf(fptr,"\n");
				}
*/				
				
				fclose(fptr);
		}
}
// -----------------------------------------------------------------------------
/*
int huge_dump(struct SimRun *v)
{

    struct randvec rng;   // The random number state holder for this run
    int    thread_ID;     // An ID number for convenience
    FILE  *dump;          // A File to dump results to
    struct RatePar RRR;   // The RATE PARAMETERS for the run
    
    //float fit_EST[MAX_COLS+1][MAX_YEARS+1];
    float *fit_vector;
    float **fit_EST;
    float **fit_diet_EST; //[NUM_PREDPREYLINKS+1][MAX_YEARS+1];
    //float fit_SSQ[MAX_COLS+1];
    double **fit_SSQ;
    float *fit_diet_SSQ;
    double **fit_bioanom;

    float  *fit_vec_control;
    float **fit_ssq_control;
    //float fit_bioanom[NUM_GROUPS+1][MAX_YEARS+1];
    //float 
    //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
    float **out_BB;
    float **out_CC;
    float **out_RR;
    float **out_EE;
    float **out_SSB;
    //float out_CC[NUM_GROUPS+1][MAX_YEARS+1];
    float **out_MM; //float out_MM[NUM_GROUPS+1][MAX_YEARS+1];
    float **out_M2; //float out_MM[NUM_GROUPS+1][MAX_YEARS+1];
    float **out_PB;  //out_PB[NUM_GROUPS+1];
    float **out_QB;  //out_QB[NUM_GROUPS+1];
    float **out_LinksM; //out_LinksM[NUM_PREDPREYLINKS+1];
    float **out_TotF;   //out_TotF[NUM_GROUPS+1];
    float out_PP;

    float **force_bymort; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];

    double  fit_SSQraw[NUM_GROUPS+MAX_COLS+1];
    double  fit_SSQfit[NUM_GROUPS+MAX_COLS+1];
    double  fit_SSQwt[NUM_GROUPS+MAX_COLS+1];
    double  fit_FORCEraw[NUM_GROUPS+1];    

    float state_BB[NUM_GROUPS+1];
    float state_Ftime[NUM_GROUPS+1];
    //float NageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    //float WageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];  
    float **NageS, **WageS;
    float state_NN[NUM_GROUPS+1];
    float stanzaPred[NUM_GROUPS+1];
    float stanzaBasePred[NUM_GROUPS+1];
    
		float vBM[MAX_SPLIT+1];
    float recruits[MAX_SPLIT+1];
		float RzeroS[MAX_SPLIT+1];
    //float WWa[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    float **WWa;
    float stanzaGGJuv[MAX_SPLIT+1];
    float stanzaGGAdu[MAX_SPLIT+1];
    float baseEggsStanza[MAX_SPLIT+1];
    float EggsStanza[MAX_SPLIT+1];
    float SpawnBio[MAX_SPLIT+1];
    float baseSpawnBio[MAX_SPLIT+1];
    float Rbase[MAX_SPLIT+1];
    
    //float baseSpawnEnergy[MAX_SPLIT+1];
    
//  hardcoding only two stanzas for now; these dimension them by months
    int firstMoJuv[MAX_SPLIT+1], lastMoJuv[MAX_SPLIT+1]; 
    int firstMoAdu[MAX_SPLIT+1], lastMoAdu[MAX_SPLIT+1];
    //float SplitAlpha[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    float **SplitAlpha;
    float RscaleSplit[MAX_SPLIT+1];  //not used right now

// OUTPUTS
   //float out_BB[NUM_GROUPS+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
   //float out_nn[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_ww[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_CC[NUM_GROUPS+1][MAX_YEARS+1];
   
   float fish_Effort[NUM_GROUPS+1];
   float TerminalF[NUM_GROUPS+1];

// FOR RECORDING SENSE RESULTS
   int DISCARD_YEAR;
	 int discarded[NUM_GROUPS+1];
	 //float month_noise[NUM_GROUPS + 1];

// Variables prefaced by deriv_ are derivative parts (e.g. current consumption)
   float deriv_TotGain[NUM_GROUPS + 1];
   float deriv_TotLoss[NUM_GROUPS + 1];
   float deriv_LossPropToB[NUM_GROUPS + 1];
   float deriv_LossPropToQ[NUM_GROUPS + 1];
   float deriv_ConsMat[NUM_GROUPS + 1][NUM_GROUPS + 1];
   float deriv_DerivT[NUM_GROUPS + 1];
   float deriv_dyt[NUM_GROUPS + 1];
   float deriv_biomeq[NUM_GROUPS + 1];

   float deriv_FoodGain[NUM_GROUPS + 1];
   float deriv_DetritalGain[NUM_GROUPS + 1];
   float deriv_FoodLoss[NUM_GROUPS + 1];
   float deriv_UnAssimLoss[NUM_GROUPS + 1];
   float deriv_ActiveRespLoss[NUM_GROUPS + 1];
   //float deriv_PassiveRespLoss[NUM_GROUPS + 1];
   float deriv_MzeroLoss[NUM_GROUPS + 1];
   float deriv_DetritalLoss[NUM_GROUPS + 1];
   float deriv_TotDetOut[NUM_GROUPS + 1];
   
   float deriv_FishingLoss[NUM_GROUPS + 1];
   float deriv_FishingGain[NUM_GROUPS + 1];
   float deriv_FishingThru[NUM_GROUPS + 1];

return 0;
}
*/
//-----------------------------------------------------------------------------
int read_climate(){

int cols, rows, iout,sp,t;
int series_index[MAX_COLS+1], series_type[MAX_COLS+1];
string sLine, ins;
double flout;

  // FORMAT OF CSV INPUT FILE MUST BE 
	// Col 1 should be timestep (numbering arbitrary)
	// Row 1:  Series Name
	// Row 2:  Group Number (0 for prim prod)
  // Row 3:  Series Type (1 = prey forcing)
  //                     (2 = mort forcing)
  //                     (3 = rec forcing)
	// ALL LINES after that:  Climate value by MONTH
	
  // Set Prey Forcing to 1.0
     for (sp=0; sp<=NUM_GROUPS; sp++){
		     for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++){
				  force_byprey[sp][t]=1.0;
				  force_bymort[sp][t]=1.0;
				  force_byrecs[sp][t]=1.0;
				  force_bysearch[sp][t]=1.0;
				 }
		 }
		 
  // Open this file and declare parser that breaks lines into fields
     ifstream infile(CLIMATE_FILE);
     if (!infile){cout << "No climate loaded, unable to read " 
		                   << CLIMATE_FILE << endl;  return 1;
		 }else {cout << "Climate loaded from " << CLIMATE_FILE << endl;}
     CSVParser parser;

  // read the first line of the file in and pass it to the parser
  // then count the number of columns in the top line
     getline(infile,sLine);  parser << sLine;
     cols=0;  ins="in";
     while (ins != ""){
		     parser >> ins;
		     if (ins == "") continue;
		     cols++;
		 }
	// read the second line and get the group indices
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){
		     parser >> series_index[sp];
		 }

	// read the third line and get the forcing type
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){
		     parser >> series_type[sp];
		 }
		 
  // read the rest of the lines and put in the forcing    
     t=0;
     while (!infile.eof()) {
           t++;
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;    
           parser >> ins; // Get the first column from parser and ignore
           for (sp=1; sp<cols; sp++){
		            parser >> flout;
		            if (t<=MAX_YEARS*STEPS_PER_YEAR){
                switch (series_type[sp]){
								    case 1:
								         force_byprey[series_index[sp]][t] += flout;
								         break;
								    case 2:
								         force_bymort[series_index[sp]][t] += flout;
								         break;								    
								    case 3:
								         force_byrecs[series_index[sp]][t] += flout;
								         break;
										case 4:
										     force_bysearch[series_index[sp]][t] += flout;
												 break;	
								    default:
								         break;
								}
								}
		       }
     }

  // Close the file
     infile.close();
}

// -----------------------------------------------------------------------------

int read_juveniles(){

int cols,iout,i,t;
string sLine, ins;
double flout;

     ifstream infile(JUV_FILE); 
     if (!infile){cout << "No juveniles loaded, unable to read " 
		                   << JUV_FILE << endl;  return 1;
		 }else {cout << "Juveniles loaded from " << JUV_FILE << endl;}
     CSVParser parser;

     getline(infile,sLine);
		 parser << sLine;
     cols=0;  ins="in";
     while (ins != ""){
         parser >> ins;
		     if (ins == "") continue;
		     cols++;
     }
     // REC_CHANGE
     if (cols!=14){cout << "Juvenile file should have 14 columns" << endl;
		             return 1;}

		 juv_N = 0;
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;
					 juv_N++;
					 parser >> ins;    juv_Names[juv_N] = ins;
           parser >> iout;   juv_JuvNum[juv_N] = (int)iout;
           parser >> iout;   juv_AduNum[juv_N] = (int)iout; 
		  		 parser >> iout;   juv_RecAge[juv_N] = (int)iout;
		  		 parser >> iout;   juv_RecMonth[juv_N] = (int)iout;
           parser >> flout;  juv_VonBK[juv_N] = (float)flout;  
	 				 parser >> flout;  juv_aduEqAgeZ[juv_N] = (float)flout;
	 				 parser >> flout;  juv_juvEqAgeZ[juv_N] = (float)flout;
            parser >> flout; juv_VonBD[juv_N] =  (float)flout;
           parser >> flout;  Wmat50[juv_N] = (float)flout;
           parser >> flout;  Wmat001[juv_N] = (float)flout;
           parser >> flout;  Amat50[juv_N] = (float)flout * 12.;
           parser >> flout;  Amat001[juv_N] = (float)flout * 12.;           
           parser >> flout;  juv_RecPower[juv_N] = (float)flout;
                   // REC_CHANGE
           //WmatWinf[i] = juv_Wmat[i]/juv_Winf[i]; 
           WmatSpread[juv_N] = -(Wmat001[juv_N]- Wmat50[juv_N])/MIN_REC_FACTOR;
           AmatSpread[juv_N] = -(Amat001[juv_N]- Amat50[juv_N])/MIN_REC_FACTOR;
        // END REC_CHANGE
           //cout << juv_Names[juv_N] << " " << juv_JuvNum[juv_N] << " " <<  juv_AduNum[juv_N]
           //<<  " " << juv_RecAge[juv_N] <<  " " << juv_VonBK[juv_N] <<  " " << juv_RecPower[juv_N]
           //<<  " " << juv_BAB[juv_N] <<  " " << juv_Wmat[juv_N] <<  " " << juv_Winf[juv_N] << endl;
		 }
		 // END REC_CHANGE		            
     cout << juv_N << endl;
     infile.close();
     
     for (i=1; i<=juv_N; i++){
         for (t=0; t<=MAX_YEARS; t++){
             FORCED_REC[i][t]=0.0;
         }
     }
     return 0;
}

// -----------------------------------------------------------------------------
int old_read_fitting(){

int cols, rows, iout,sp,t,c,iin, gr,j,i;
string sLine, ins;
double flout; 
float cattot;
float fit_renormI[NUM_LIVING+MAX_COLS+1];
float fit_renormB[NUM_LIVING+MAX_COLS+1];
int   fit_renormN[NUM_LIVING+MAX_COLS+1];

float FORCED_BIO_renormI[NUM_GROUPS+1];
int   FORCED_BIO_renormN[NUM_GROUPS+1];
float FORCED_CAT_renormI[NUM_GROUPS+1];
int   FORCED_CAT_renormN[NUM_GROUPS+1];

float FORCED_REC_renormI[NUM_GROUPS+1];
int   FORCED_REC_renormN[NUM_GROUPS+1];
int marked[NUM_GROUPS+1];

   string in_series_names[NUM_LIVING+MAX_COLS+1];
   int    in_series_index[NUM_LIVING+MAX_COLS+1];
   int    in_series_type[NUM_LIVING+MAX_COLS+1];
   int    in_pedigree[NUM_LIVING+MAX_COLS+1];
  // FORMAT OF CSV INPUT FILE MUST BE 
	// Col 1 should be timestep (numbering arbitrary)
	// Row 1:  Series Name
	// Row 2:  Group Number (0 for prim prod)
  // Row 3:  Series Type 

	// ALL LINES after that:  FITTING VALUE BY YEARS
	
	
  // Set Prey Forcing to 1.0
     for (sp=0; sp<=NUM_GROUPS; sp++){
         //bio_renorm[sp]    = -1.;
         marked[sp]=0;
         FORCE_PEDIGREE[sp] = 0;
         FORCE_BIO_FLAG[sp]=0;
         fit_bioStart[sp]  = 0.0;
         FORCED_TARGET[sp] = 0;
         FORCED_FTARGET[sp] = 0;
         FORCED_BIO_renormI[sp]=0.0;
         FORCED_BIO_renormN[sp]=0;
         FORCED_CAT_renormI[sp]=0.0;
         FORCED_CAT_renormN[sp]=0;
         FORCED_REC_renormI[sp]=0.0;
         FORCED_REC_renormN[sp]=0;
		     for (t=0; t<=MAX_YEARS; t++){
		          if (sp<=NUM_LIVING+NUM_DEAD){FORCED_CATCH[sp][t]=0.0; 
		                                       FORCED_CATCH_RAW[sp][t]=0.0; 
							                             FORCED_FRATE[sp][t]=0.0; 
							}
		          else                        {FORCED_CATCH[sp][t]=1.0; 
		                                       FORCED_CATCH_RAW[sp][t]=1.0; 
							                             FORCED_FRATE[sp][t]=1.0;
							}


              FORCED_EFFORT[sp][t] = -1.0; 
		          FORCED_BIO[sp][t]=0.0;
		          FORCED_BIO_RAW[sp][t]=0.0;
		          if (sp<=juv_N){
		            FORCED_REC[sp][t]=-1.0;
		          }
				 }
		 }
		 //bioN=0; catN=0;
		 fitN=0;
     for (c=0; c<=MAX_COLS; c++){
		     for (t=0; t<=MAX_YEARS; t++){
		          fit_renormI[c] = 0.0;
		          fit_renormB[c] = 0.0;
		          fit_renormN[c] = 0;
              fit_OBS[c][t] = 0.0;
              fit_RAW[c][t] = -90.0;
              fit_SD[c][t]  = 1.0;		 
		 //	      bioOBS[c][t]=0.0;  bioSD[c][t]=1.0;  bioEST[c][t]=0.0;
		 //		      catOBS[c][t]=0.0;  catSD[c][t]=1.0;  catEST[c][t]=0.0;
         }
      }
    
  // Open this file and declare parser that breaks lines into fields
     ifstream infile(FITTING_FILE); 
     if (!infile){cout << "No fitting loaded, unable to read " 
		                   << FITTING_FILE << endl;  return 1;
		 }else {cout << "Fitting loaded from " << FITTING_FILE << endl;}
     CSVParser parser;

  // read the first line of the file in and pass it to the parser
  // then count the number of columns in the top line
     getline(infile,sLine);
		 parser << sLine;
     cols=0;  ins="in";
     while (ins != ""){
         parser >> ins;
		     in_series_names[cols] = ins;
		     if (ins == "") continue;
		     cols++;
		 }
	// read the second line and get the group indices
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){parser >> in_series_index[sp];}

	// read the next line and get the pedigree
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){parser >> in_pedigree[sp];}

	// read the third line and get the forcing type
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){parser >> in_series_type[sp];}	 

  // read the rest of the lines and put in the forcing    
     t=0;
     fit_lastyear=0;
     //fit_max=0;
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;    
           parser >> iin; // Get the first column from parser and ignore
           if (t==0){fit_StartYear=iin; }
           fitN=0;
           // KYA 4/8/08 WHY HERE fit year???  Move down below flout read
           //if (flout>0.){fit_lastyear=t;}
           for (c=1; c<cols; c++){
		            parser >> flout;
		            if (flout>0.){fit_lastyear=t;}
                switch (in_series_type[c]){

                       //Case 0:  return_fittype = "Relative Biomass" 'Biomass relative (CPUE)
                       //Case -100:  Absolute Biomass with no normalization
                       case -100:
                       case 0:
                          fitN++;
                          fit_names[fitN]     = in_series_names[c];
                          fit_type[fitN]      = in_series_type[c];
													fit_ApplyTo[fitN]   = in_series_index[c];
													fit_ApplyFrom[fitN] = 0;
													marked[fit_ApplyTo[fitN]]=1;
													//if (flout>0.){fit_lastyear=t;}
													fit_OBS[fitN][t]  = flout;
													fit_RAW[fitN][t]  = flout;
													if(in_series_type[c]!=-100){
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              fit_renormI[fitN]  += flout;
                              fit_renormB[fitN]  = path_BB[fit_ApplyTo[fitN]];
                              fit_renormN[fitN] ++;
                              } 
                          }
                          if ((fit_bioStart[fit_ApplyTo[fitN]] <= 0. ) && (flout > 0.))
													    {
													      fit_bioStart[fit_ApplyTo[fitN]] = flout; 
															}

                       break;
                       case 6:
                          fitN++;
                          fit_names[fitN]     = in_series_names[c];
                          fit_type[fitN]      = in_series_type[c];
													fit_ApplyTo[fitN]   = in_series_index[c];
													fit_ApplyFrom[fitN] = 0;
                          fit_OBS[fitN][t]    = flout;
													fit_RAW[fitN][t]    = flout;  

                          //if (flout>0.){fit_lastyear=t;}
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              fit_renormI[fitN]  += flout;
                              cattot=0.0;
                              for (gr=0; gr<=NUM_GEARS; gr++){
                                   cattot+=path_Catch[gr][fit_ApplyTo[fitN]] + path_Discards[gr][fit_ApplyTo[fitN]];
                              }
                              fit_renormB[fitN]  = cattot; //path_BB[fit_ApplyTo[fitN]];
                              fit_renormN[fitN] ++;
                              }                           
                          
                                                
                       break;

                       //Case 1:  return_fittype = "Absolute Biomass" 'Biomass absolute
                       //Case 2:  return_fittype = "unk" 'Time Forcing, in climate instead
                       case 3:  //return_fittype = "Forced Gear effort" 'Effort data by gear type       
                            if (in_series_index[c] > NUM_LIVING+NUM_DEAD){
											         FORCED_EFFORT[in_series_index[c]][t] =  flout;
											      }
											 case 4:  //return_fittype = "F by Species" 'F by biomass pool
                          if (in_series_index[c] <= NUM_LIVING+NUM_DEAD){
											       FORCED_FRATE[in_series_index[c]][t] +=  flout;
											       }
											    else{ if (in_series_index[c] <= NUM_GROUPS){
											              FORCED_FRATE[in_series_index[c]][t] = flout;
											          }
											          else {
														         gr = (int)(in_series_index[c]/1000);
														         sp  = in_series_index[c] - gr * 1000;
														         FORCED_FRATE[gr][t] = flout;
														         FORCED_FTARGET[gr] = sp;
														    }     
											    }
											 break;                       
                       
                       //Case 5:  return_fittype = "Z for fitting" 'Z by pool for fitting
                       //Case 6:  return_fittype = "Absolute Catch" 'Catches by pool for fitting

                       
                       //Case 7:  return_fittype = "unk" 'Mean Body weight from Martell's brain
                  //case 101: //return_fittype = "Recruitment" 'Kerim's B/Bstart for no good reason

                  //Case 102: return_fittype  = "Eggs"
                  case -102:
                  case 102:
                          j=0;
                          for (i=1; i<=juv_N; i++){
                              if ((juv_JuvNum[i]==in_series_index[c]) || (juv_AduNum[i]==in_series_index[c])){
                                 j=i;
                              }
                          }
                          if (flout>=0.0){
                             FORCED_REC[j][t]=flout;
                             }
                          if(in_series_type[c]!=-102){
													   if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													         ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													         (flout>0.)
													       )
													       {
                                 FORCED_REC_renormI[j]  += flout;
                                 FORCED_REC_renormN[j] ++;
                                 }
                             }                                       
                       break;                
                // Case 105: R/S anomaly
                
                
                // Case -1: return_fittype = "Forced Biomass" 'Pool biomass forcing
                // Case -101: non-normalized Forced Biomass                
                    case -101:
                    case -1:
                       FORCE_PEDIGREE[in_series_index[c]] = in_pedigree[c];
                       FORCED_BIO[in_series_index[c]][t] = flout;
                       FORCED_BIO_RAW[in_series_index[c]][t] = flout;
                       FORCE_BIO_FLAG[in_series_index[c]]=1;
                       marked[in_series_index[c]]=1;
                       if(in_series_type[c]!=-101){
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              FORCED_BIO_renormI[in_series_index[c]]  += flout;
                              FORCED_BIO_renormN[in_series_index[c]] ++;
                              }          
                       }
											 else{
											    FORCED_BIO[in_series_index[c]][t] /= AREA_MULT;
											 }           

                       break;
                
								// Case 103: return_fittype = "Gear Targeted Catch"
                    case 103:
                    case -103:
                       //if (in_series_index[c] <= NUM_LIVING+NUM_DEAD){
											 //   FORCED_CATCH[in_series_index[c]][t] +=  flout;
											 //   }
											 if (in_series_index[c] <= NUM_GROUPS){
											     gr = in_series_index[c];
											     sp = in_series_index[c];
											     //FORCED_CATCH[gr][t] +=  flout;
											     //FORCED_TARGET[gr] = sp;
											 }											 
											 else{ //if (in_series_index[c] <= NUM_GROUPS){
											       //    FORCED_CATCH[in_series_index[c]][t] = flout;
											       //}
											       //else {
														      gr = (int)(in_series_index[c]/1000);
														      sp  = in_series_index[c] - gr * 1000;
														      //FORCED_CATCH[gr][t] = flout;
														      //FORCED_TARGET[gr] = sp;
														 //}      
											 }
											 FORCED_CATCH[gr][t] = flout;
											 FORCED_TARGET[gr]   = sp;
											 FORCED_CATCH_RAW[FORCED_TARGET[gr]][t] = flout;
											 if(in_series_type[c]!=-103){
											 if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              FORCED_CAT_renormI[gr]  += flout;
                              FORCED_CAT_renormN[gr] ++;
                              }  
											 }
											 else{
											       FORCED_CATCH[in_series_index[c]][t] /= AREA_MULT;
											 }
											 

											 break;
                       //Case Else: return_fittype = "unk"
								    default:
								         break;
								}
		       }
		 //if (fit_max<fitN){fit_max=fitN;}
		 t++;	
     //cout << endl;       
     }
     fit_Years=t;
     MAX_THRESHOLD(fit_Years,MAX_YEARS);
     cout << "fit_Years: " << fit_Years << " Fit LastYear: " << fit_lastyear << endl;
     //fitN=fit_max;
  // Close the file
     infile.close();
     
  // Add Path fitting to the list!
     for (c=1; c<=NUM_LIVING; c++){
         fitN++;
         fit_names[fitN]     = "Path";
         fit_type[fitN]      = 50;
				 fit_ApplyTo[fitN]   = c;
				 fit_ApplyFrom[fitN] = 0;
				 for (t=PATH_YEAR_START; t<=PATH_YEAR_END; t++){
				     fit_OBS[fitN][t-fit_StartYear]  = path_BB[c];
         }
         fit_renormI[fitN]  = 1.0;
         fit_renormB[fitN]  = 1.0;
         fit_renormN[fitN]  = 0;
    }

     for (c=1; c<=fitN; c++){
          if (fit_renormN[c]){
              fit_renormI[c] /= (float)fit_renormN[c];
              if (fit_type[c]==0){ 
               if (fit_bioStart[fit_ApplyTo[c]] > 0){
                  //cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
                  //cout << " " << fit_renormI[c];
                  fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]];
                  fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
               }
              }
    		      for (t=0; t<=fit_lastyear; t++){
    		          fit_OBS[c][t] *= fit_renormB[c]/fit_renormI[c];
    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  fit_Norm[c] = fit_renormB[c]/fit_renormI[c];
    		  }
    		  else{
              fit_Norm[c]=1.0;
          }     		  
     }
     
     for (c=1; c<=juv_N; c++){
          if (FORCED_REC_renormN[c]){
              FORCED_REC_renormI[c] /= (float)FORCED_REC_renormN[c];
              //if (fit_bioStart[fit_ApplyTo[c]] > 0){
              //    cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
              //    cout << " " << fit_renormI[c];
              //    fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
              //    cout << " " << fit_bioStart[fit_ApplyTo[c]];
              //    fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
              //    cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
              //}
    		      for (t=0; t<=fit_lastyear; t++){
    		          FORCED_REC[c][t] *= 1.0/FORCED_REC_renormI[c];
    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  }
     } 
     

    
    for (c=1; c<=NUM_GROUPS; c++){
         cattot=0.0;
         // Case of non-gear applied catch
         if ((c<=NUM_LIVING+NUM_DEAD) && (FORCED_TARGET[c]!=0)){
            for (gr=0; gr<=NUM_GEARS; gr++){
                cattot+=path_Catch[gr][FORCED_TARGET[c]] + path_Discards[gr][FORCED_TARGET[c]];
                }
         }              
         // Case of gear applied catch
         if ((c>NUM_LIVING+NUM_DEAD) && (FORCED_TARGET[c]!=0)){
             gr = c - NUM_LIVING - NUM_DEAD;
             cattot = path_Catch[gr][FORCED_TARGET[c]] + path_Discards[gr][FORCED_TARGET[c]];
         }
         // cout << c << " c " << endl;
         if (FORCED_CAT_renormN[c]){
            FORCED_CAT_renormI[c] /= (float)FORCED_CAT_renormN[c];
    		    //cout << c << " " << path_species[c] << " " << FORCED_CAT_renormN[c] 
            //     << " " << cattot << " " << FORCED_CAT_renormI[c] <<  endl;
    		    for (t=0; t<=fit_lastyear; t++){
    		         FORCED_CATCH[c][t] *= cattot/FORCED_CAT_renormI[c];


    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		    }         
         }
     }
     
     for (c=1; c<=NUM_GROUPS; c++){
          if (FORCED_BIO_renormN[c]){
              FORCED_BIO_renormI[c] /= (float)FORCED_BIO_renormN[c];
              //if (fit_bioStart[fit_ApplyTo[c]] > 0){
                  //cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
                  //cout << " " << fit_renormI[c];
                  //fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]];
                  //fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
                  //cout << c << " " << path_BB[c] << " " << FORCED_BIO_renormI[c] 
                  //<< " " << path_BB[c]/FORCED_BIO_renormI[c] << endl;
              //}
    		      for (t=0; t<=fit_lastyear; t++){
    		          FORCED_BIO[c][t] *= path_BB[c]/FORCED_BIO_renormI[c];
    		          
    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  }
     }
     
 // THIS LINE IS USED to fix uforced groups to not change
    // for (c=1; c<=NUM_GROUPS; c++){
    //     if(!marked[c]){
    //		      for (t=0; t<=fit_lastyear; t++){
    //		          FORCED_BIO[c][t] = path_BB[c];
    //		      }         
    //    }
    // }
     //cout << fit_lastyear << endl;         
               
}

//------------------------------------------------------------------------------
int read_diets(){

int cols, rows, iout,sp,t,c,iin, gr,j,i,pred,prey;
string sLine, ins;
double flout; 
char tmpstr[80];

   string in_series_names[NUM_LIVING+MAX_COLS+1];
   int    in_series_index[NUM_LIVING+MAX_COLS+1];
   int    in_series_type[NUM_LIVING+MAX_COLS+1];
  // FORMAT OF CSV INPUT FILE MUST BE 
	// Col 1 should be timestep (numbering arbitrary)
	// Row 1:  Series Name
	// Row 2:  Group Number (0 for prim prod)
  // Row 3:  Series Type 

	// ALL LINES after that:  FITTING VALUE BY YEARS
		
		 dietN = NUM_LIVING;
		 
		 for (c=1; c<=NUM_LIVING; c++){
         for (t=0; t<=MAX_YEARS; t++){
		       fit_diet_OBS[c][t] = 1.0;
		     }
		     sprintf(tmpstr,"%s_diet_other",path_species[c]);
         fit_diet_names[c].assign(tmpstr);
         fit_diet_type[c]        = 40;
				 fit_diet_lookup[c]      = c*1000 + 999;     
		}
     for (c=NUM_LIVING+1; c<=NUM_LIVING+NUM_PREDPREYLINKS; c++){
		     for (t=0; t<=MAX_YEARS; t++){
       		 fit_diet_OBS[c][t]  = 0.0;
          }
      }
    
  // Open this file and declare parser that breaks lines into fields
     ifstream infile(DIET_FILE); 
     if (!infile){cout << "No diet loaded, unable to read " 
		                   << DIET_FILE << endl;  return 1;
		 }else {cout << "Diets loaded from " << DIET_FILE << endl;}
     CSVParser parser;

  // read the first line of the file in and pass it to the parser
  // then count the number of columns in the top line
     getline(infile,sLine);
		 parser << sLine;
     cols=0;  ins="in";
     while (ins != ""){
         parser >> ins;
		     in_series_names[cols] = ins;
		     if (ins == "") continue;
		     cols++;
		 }
	// read the second line and get the group indices
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){parser >> in_series_index[sp];}
	// read the third line and get the forcing type
     getline(infile,sLine);  parser << sLine;    
     parser >> ins; // Get the first column from parser and ignore
     for (sp=1; sp<cols; sp++){parser >> in_series_type[sp];}	 
  // read the rest of the lines and put in the forcing    
     t=0;
     //fit_lastyear=0;
     //fit_max=0;
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;    
           parser >> iin; // Get the first column from parser and ignore
           //if (t==0){fit_StartYear=iin; }
           dietN=NUM_LIVING;
           //if (flout>0.){fit_lastyear=t;}
           for (c=1; c<cols; c++){
		            parser >> flout;
                switch (in_series_type[c]){
                       case 40:
												  pred = int(in_series_index[c]/1000); 
                          prey = in_series_index[c] - pred*1000;
                          
                          if (prey<=NUM_LIVING){
                             dietN++;
                             fit_diet_names[dietN]     = in_series_names[c];
                             fit_diet_type[dietN]      = in_series_type[c];
													   fit_diet_lookup[dietN]  = in_series_index[c];
													   
													   fit_diet_OBS[dietN][t]  = flout;
													   fit_diet_OBS[pred][t]  -=flout;
														 if (fit_diet_OBS[pred][t]<-0.0000001){
														     cout << "Error: Diet of pred " << pred
														          << "in year " << t << " more than 1.0"
																			<< "prey was " << prey << " flout " << flout 
																			<< "other " << fit_diet_OBS[pred][t]  
																			<< endl;
																 exit(1);
																 }
                          }
													//      int dietN;
                          //      int fit_diet_lookup[NUM_GROUPS+1][NUM_GROUPS+1];
                          //      float fit_diet_obs[NUM_PREDPREYLINKS][MAX_YEARS+1];
                          //      float fit_diet_est[NUM_PREDPREYLINKS][MAX_YEARS+1];
                       break;
                       
                       default:
                       break;
                }
           }
    t++;
    }
}
//------------------------------------------------------------------------------

int new_read_fitting(){

int cols, rows, iout,sp,t,c,iin, gr,j,i,cmatch;
string sLine, ins;

string Series, Source;
int TypeCode, ApplyTo;
double DefaultCV, Weight;
double Renorm;
double flout; 
float cattot;
float fit_renormI[NUM_LIVING+MAX_COLS+1];
float fit_renormB[NUM_LIVING+MAX_COLS+1];
int   fit_renormN[NUM_LIVING+MAX_COLS+1];

float FORCED_BIO_renormI[NUM_GROUPS+1];
int   FORCED_BIO_renormN[NUM_GROUPS+1];
float FORCED_CAT_renormI[NUM_GROUPS+1];
int   FORCED_CAT_renormN[NUM_GROUPS+1];

float FORCED_REC_renormI[NUM_GROUPS+1];
int   FORCED_REC_renormN[NUM_GROUPS+1];
int marked[NUM_GROUPS+1];

   string in_series_names[NUM_LIVING+MAX_COLS+1];
   int    in_series_index[NUM_LIVING+MAX_COLS+1];
   int    in_series_type[NUM_LIVING+MAX_COLS+1];
   int    in_pedigree[NUM_LIVING+MAX_COLS+1];
  // FORMAT OF CSV INPUT FILE MUST BE 
	// Col 1 should be timestep (numbering arbitrary)
	// Row 1:  Series Name
	// Row 2:  Group Number (0 for prim prod)
  // Row 3:  Series Type 

	// ALL LINES after that:  FITTING VALUE BY YEARS
	
	
  // Set Prey Forcing to 1.0
     for (sp=0; sp<=NUM_GROUPS; sp++){
         //bio_renorm[sp]    = -1.;
         marked[sp]=0;
         FORCE_PEDIGREE[sp] = 0;
         FORCE_BIO_FLAG[sp]=0;
         fit_bioStart[sp]  = 0.0;
         FORCED_TARGET[sp] = 0;
         FORCED_FTARGET[sp] = 0;
         FORCED_BIO_renormI[sp]=0.0;
         FORCED_BIO_renormN[sp]=0;
         FORCED_CAT_renormI[sp]=0.0;
         FORCED_CAT_renormN[sp]=0;
         FORCED_REC_renormI[sp]=0.0;
         FORCED_REC_renormN[sp]=0;
		     for (t=0; t<=MAX_YEARS; t++){
		          if (sp<=NUM_LIVING+NUM_DEAD){FORCED_CATCH[sp][t]=0.0; 
		                                       FORCED_CATCH_RAW[sp][t]=0.0; 
							                             FORCED_FRATE[sp][t]=0.0; 
							}
		          else                        {FORCED_CATCH[sp][t]=1.0; 
		                                       FORCED_CATCH_RAW[sp][t]=1.0; 
							                             FORCED_FRATE[sp][t]=1.0;
							}


              FORCED_EFFORT[sp][t] = -1.0; 
		          FORCED_BIO[sp][t]=0.0;
		          FORCED_BIO_RAW[sp][t]=0.0;
		          if (sp<=juv_N){
		            FORCED_REC[sp][t]=-1.0;
		          }
				 }
		 }
		 //bioN=0; catN=0;
		 fitN=0;
     for (c=0; c<=MAX_COLS; c++){
		     for (t=0; t<=MAX_YEARS; t++){
		          fit_renormI[c] = 0.0;
		          fit_renormB[c] = 0.0;
		          fit_renormN[c] = 0;
              fit_OBS[c][t] = 0.0;
              fit_RAW[c][t] = -90.0;
              fit_SD[c][t]  = 1.0;		 
		 //	      bioOBS[c][t]=0.0;  bioSD[c][t]=1.0;  bioEST[c][t]=0.0;
		 //		      catOBS[c][t]=0.0;  catSD[c][t]=1.0;  catEST[c][t]=0.0;
         }
      }
    
  // Open this file and declare parser that breaks lines into fields
     ifstream infile(FITTING_FILE); 
     if (!infile){cout << "No fitting loaded, unable to read " 
		                   << FITTING_FILE << endl;  return 1;
		 }else {cout << "Fitting loaded from " << FITTING_FILE << endl;}
     CSVParser parser;

  // read the first line of the file in and pass it to the parser
  // then count the number of columns in the top line
     getline(infile,sLine);
		 parser << sLine;
	// Trim the first set of columns
	   for (i=1; i<=7; i++){parser >> ins;}
     fit_lastyear=0;  ins="in";
     while (ins != ""){
         parser >> ins;
		     if (ins == "") continue;
		     iin=atoi(ins.c_str());
		     if (fit_lastyear==0){fit_StartYear=iin;}
		     fit_lastyear++;
		 }

		 fitN = 0;		 
		 while (!infile.eof()) {
		      getline(infile,sLine);  parser << sLine;
					parser >> Series;
					parser >> Source;
          parser >> TypeCode;
					parser >> ApplyTo;
          parser >> DefaultCV;
					parser >> Weight;
          parser >> Renorm;
          
          switch (TypeCode){
					
					   case 0:
					   case 100:
                   if (DefaultCV>0){
									   fitN++;
                     fit_names[fitN]           = Series;
                     fit_type[fitN]            = TypeCode;
									   fit_ApplyTo[fitN]         = ApplyTo;
									   fit_ApplyFrom[fitN]       = 0;
									   marked[fit_ApplyTo[fitN]] = 1;
                     for (t=0; t<fit_lastyear; t++){
		                      parser >> flout;										 
													fit_OBS[fitN][t]  = flout;
													fit_RAW[fitN][t]  = flout;
													fit_SD[fitN][t]   = flout * DefaultCV;
													if(TypeCode != -100){
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              fit_renormI[fitN]  += flout;
                              fit_renormB[fitN]  = path_BB[fit_ApplyTo[fitN]];
                              fit_renormN[fitN] ++;
                              } 
                          }
                          if ((fit_bioStart[fit_ApplyTo[fitN]] <= 0. ) && (flout > 0.))
													    {
													      fit_bioStart[fit_ApplyTo[fitN]] = flout; 
															}
										 }													 					      
					         }
					         else{
					           // Find the first mean time series that this can apply to!
					           cmatch=0;
					           for (c=1; c<=fitN; c++){
										      if ((ApplyTo == fit_ApplyTo[c]) && (TypeCode == fit_type[c])){
													   cmatch=c;
													}
										 }
                     for (t=0; t<fit_lastyear; t++){
		                      parser >> flout;
													fit_SD[cmatch][t] = flout;
										 }	
									 }
					
					       break;
					   
                       case 6:
                          fitN++;
                          fit_names[fitN]     = Series;
                          fit_type[fitN]      = TypeCode;
													fit_ApplyTo[fitN]   = ApplyTo;
													fit_ApplyFrom[fitN] = 0;
																							
													for (t=0; t<fit_lastyear; t++){
		                         parser >> flout;
                             fit_OBS[fitN][t]    = flout;
													   fit_RAW[fitN][t]    = flout;  
                             fit_SD[fitN][t]     = flout * DefaultCV;
                          //if (flout>0.){fit_lastyear=t;}
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              fit_renormI[fitN]  += flout;
                              cattot=0.0;
                              for (gr=0; gr<=NUM_GEARS; gr++){
                                   cattot+=path_Catch[gr][fit_ApplyTo[fitN]] + path_Discards[gr][fit_ApplyTo[fitN]];
                              }
                              fit_renormB[fitN]  = cattot; //path_BB[fit_ApplyTo[fitN]];
                              fit_renormN[fitN] ++;
                              }                           
                          } 
                                                
                       break;

                       //Case 1:  return_fittype = "Absolute Biomass" 'Biomass absolute
                       //Case 2:  return_fittype = "unk" 'Time Forcing, in climate instead
                       case 3:  //return_fittype = "Forced Gear effort" 'Effort data by gear type       
                            if (ApplyTo > NUM_LIVING+NUM_DEAD){                            
                            	for (t=0; t<fit_lastyear; t++){
		                             parser >> flout;
											           FORCED_EFFORT[ApplyTo][t] =  flout;
											         }
											      }
											 case 4:  //return_fittype = "F by Species" 'F by biomass pool
                          if (ApplyTo <= NUM_LIVING+NUM_DEAD){
                            	for (t=0; t<fit_lastyear; t++){
		                             parser >> flout;                          
											           FORCED_FRATE[ApplyTo][t] +=  flout * Renorm;
											        }
											       }
											    else{ if (ApplyTo <= NUM_GROUPS){
                            	   for (t=0; t<fit_lastyear; t++){
		                                parser >> flout;         											           
											              FORCED_FRATE[ApplyTo][t] = flout;
											            }
											          }
											          else {
														         gr = (int)(ApplyTo/1000);
														         sp  = ApplyTo - gr * 1000;
                            	       for (t=0; t<fit_lastyear; t++){
		                                     parser >> flout; 
														             FORCED_FRATE[gr][t] = flout;
														         }
														         FORCED_FTARGET[gr] = sp;
														    }     
											    }
											 break;                       
                       
                       //Case 5:  return_fittype = "Z for fitting" 'Z by pool for fitting
                       //Case 6:  return_fittype = "Absolute Catch" 'Catches by pool for fitting

                       
                       //Case 7:  return_fittype = "unk" 'Mean Body weight from Martell's brain
                  //case 101: //return_fittype = "Recruitment" 'Kerim's B/Bstart for no good reason

                  //Case 102: return_fittype  = "Eggs"
                  case -102:
                  case 102:
                          j=0;
                          for (i=1; i<=juv_N; i++){
                              if ((juv_JuvNum[i]==ApplyTo) || (juv_AduNum[i]==ApplyTo)){
                                 j=i;
                              }
                          }
                          
                          for (t=0; t<fit_lastyear; t++){
		                          parser >> flout; 
                              if (flout>=0.0){FORCED_REC[j][t]=flout;}
                              
                             
                          if(in_series_type[c]!=-102){
													   if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													         ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													         (flout>0.)
													       )
													       {
                                 FORCED_REC_renormI[j]  += flout;
                                 FORCED_REC_renormN[j] ++;
                                 }
                               }
														 }                                       
                       break;                
                // Case 105: R/S anomaly
                
                
                // Case -1: return_fittype = "Forced Biomass" 'Pool biomass forcing
                // Case -101: non-normalized Forced Biomass                
                    case -101:
                    case -1:
                     for (t=0; t<fit_lastyear; t++){
		                   parser >> flout;                     
                    
                       FORCE_PEDIGREE[ApplyTo] = int(Weight);
                       FORCED_BIO[ApplyTo][t] = flout;
                       FORCED_BIO_RAW[ApplyTo][t] = flout;
                       FORCE_BIO_FLAG[ApplyTo]=1;
                       marked[ApplyTo]=1;
                       if(in_series_type[c]!=-101){
													if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              FORCED_BIO_renormI[ApplyTo]  += flout;
                              FORCED_BIO_renormN[ApplyTo] ++;
                              }          
                       }
											 else{
											    FORCED_BIO[ApplyTo][t] /= AREA_MULT;
											 }           
                       }
                       break;
                
								// Case 103: return_fittype = "Gear Targeted Catch"
                    case 103:
                    case -103:
                       //if (ApplyTo <= NUM_LIVING+NUM_DEAD){
											 //   FORCED_CATCH[ApplyTo][t] +=  flout;
											 //   }
											 if (ApplyTo <= NUM_GROUPS){
											     gr = ApplyTo;
											     sp = ApplyTo;
											     //FORCED_CATCH[gr][t] +=  flout;
											     //FORCED_TARGET[gr] = sp;
											 }											 
											 else{ //if (ApplyTo <= NUM_GROUPS){
											       //    FORCED_CATCH[ApplyTo][t] = flout;
											       //}
											       //else {
														      gr = (int)(ApplyTo/1000);
														      sp  = ApplyTo - gr * 1000;
														      //FORCED_CATCH[gr][t] = flout;
														      //FORCED_TARGET[gr] = sp;
														 //}      
											 }
											 for (t=0; t<fit_lastyear; t++){
		                   parser >> flout; 
											 
											 FORCED_CATCH[gr][t] = flout;
											 FORCED_TARGET[gr]   = sp;
											 FORCED_CATCH_RAW[FORCED_TARGET[gr]][t] = flout;
											 if(in_series_type[c]!=-103){
											 if  ( ((t+fit_StartYear)>=PATH_YEAR_START) &&
													      ((t+fit_StartYear)<=PATH_YEAR_END)   &&
													      (flout>0.)
													    )
													    {
                              FORCED_CAT_renormI[gr]  += flout;
                              FORCED_CAT_renormN[gr] ++;
                              }  
											 }
											 else{
											       FORCED_CATCH[ApplyTo][t] /= AREA_MULT;
											 }
											 }

											 break;
                       //Case Else: return_fittype = "unk"
								    default:
								         break;
								}
		       }
		 //if (fit_max<fitN){fit_max=fitN;}
     fit_lastyear--;
     
     fit_Years = fit_lastyear + 1;
     MAX_THRESHOLD(fit_Years,MAX_YEARS);
     cout << "fit_Years: " << fit_Years+1 << " Fit LastYear: " << fit_lastyear+1 << endl;
     //fitN=fit_max;
  // Close the file
     infile.close();
     
  // Add Path fitting to the list!
     for (c=1; c<=NUM_LIVING; c++){
         fitN++;
         fit_names[fitN]     = "Path";
         fit_type[fitN]      = 50;
				 fit_ApplyTo[fitN]   = c;
				 fit_ApplyFrom[fitN] = 0;
				 for (t=PATH_YEAR_START; t<=PATH_YEAR_END; t++){
				     fit_OBS[fitN][t-fit_StartYear]  = path_BB[c];
         }
         fit_renormI[fitN]  = 1.0;
         fit_renormB[fitN]  = 1.0;
         fit_renormN[fitN]  = 0;
    }

     for (c=1; c<=fitN; c++){
          if (fit_renormN[c]){
              fit_renormI[c] /= (float)fit_renormN[c];
              if (fit_type[c]==0){ 
               if (fit_bioStart[fit_ApplyTo[c]] > 0){
                  //cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
                  //cout << " " << fit_renormI[c];
                  fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]];
                  fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
               }
              }
    		      for (t=0; t<=fit_lastyear; t++){
    		          fit_OBS[c][t] *= fit_renormB[c]/fit_renormI[c];
    		          fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  fit_Norm[c] = fit_renormB[c]/fit_renormI[c];
    		  }
    		  else{
              fit_Norm[c]=1.0;
          }     		  
     }
     
     for (c=1; c<=juv_N; c++){
          if (FORCED_REC_renormN[c]){
              FORCED_REC_renormI[c] /= (float)FORCED_REC_renormN[c];
              //if (fit_bioStart[fit_ApplyTo[c]] > 0){
              //    cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
              //    cout << " " << fit_renormI[c];
              //    fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
              //    cout << " " << fit_bioStart[fit_ApplyTo[c]];
              //    fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
              //    cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
              //}
    		      for (t=0; t<=fit_lastyear; t++){
    		          FORCED_REC[c][t] *= 1.0/FORCED_REC_renormI[c];
    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  }
     } 
     

    
    for (c=1; c<=NUM_GROUPS; c++){
         cattot=0.0;
         // Case of non-gear applied catch
         if ((c<=NUM_LIVING+NUM_DEAD) && (FORCED_TARGET[c]!=0)){
            for (gr=0; gr<=NUM_GEARS; gr++){
                cattot+=path_Catch[gr][FORCED_TARGET[c]] + path_Discards[gr][FORCED_TARGET[c]];
                }
         }              
         // Case of gear applied catch
         if ((c>NUM_LIVING+NUM_DEAD) && (FORCED_TARGET[c]!=0)){
             gr = c - NUM_LIVING - NUM_DEAD;
             cattot = path_Catch[gr][FORCED_TARGET[c]] + path_Discards[gr][FORCED_TARGET[c]];
         }
         // cout << c << " c " << endl;
         if (FORCED_CAT_renormN[c]){
            FORCED_CAT_renormI[c] /= (float)FORCED_CAT_renormN[c];
    		    //cout << c << " " << path_species[c] << " " << FORCED_CAT_renormN[c] 
            //     << " " << cattot << " " << FORCED_CAT_renormI[c] <<  endl;
    		    for (t=0; t<=fit_lastyear; t++){
    		         FORCED_CATCH[c][t] *= cattot/FORCED_CAT_renormI[c];


    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		    }         
         }
     }
     
     for (c=1; c<=NUM_GROUPS; c++){
          if (FORCED_BIO_renormN[c]){
              FORCED_BIO_renormI[c] /= (float)FORCED_BIO_renormN[c];
              //if (fit_bioStart[fit_ApplyTo[c]] > 0){
                  //cout << fit_ApplyTo[c] << " " << fit_bioStart[fit_ApplyTo[c]];
                  //cout << " " << fit_renormI[c];
                  //fit_bioStart[fit_ApplyTo[c]] /= fit_renormI[c];
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]];
                  //fit_bioStart[fit_ApplyTo[c]] = log(fit_bioStart[fit_ApplyTo[c]]);
                  //cout << " " << fit_bioStart[fit_ApplyTo[c]] << endl;
                  //cout << c << " " << path_BB[c] << " " << FORCED_BIO_renormI[c] 
                  //<< " " << path_BB[c]/FORCED_BIO_renormI[c] << endl;
              //}
    		      for (t=0; t<=fit_lastyear; t++){
    		          FORCED_BIO[c][t] *= path_BB[c]/FORCED_BIO_renormI[c];
    		          
    		          //fit_SD[c][t]  *= fit_renormB[c]/fit_renormI[c];
    		      }
    		  }
     }
     
 // THIS LINE IS USED to fix uforced groups to not change
    // for (c=1; c<=NUM_GROUPS; c++){
    //     if(!marked[c]){
    //		      for (t=0; t<=fit_lastyear; t++){
    //		          FORCED_BIO[c][t] = path_BB[c];
    //		      }         
    //    }
    // }
     //cout << fit_lastyear << endl;         
               
}


// -----------------------------------------------------------------------------


int read_guilds(){

int sp, cols, inint;
string sLine, ins, species, guild, instr;

       for (sp=1; sp<=NUM_GROUPS; sp++){
            PathNames[sp].assign(path_species[sp]);             
       }

  // Open this file and declare parser that breaks lines into fields
     ifstream infile("guild_lookup.csv"); 
     if (!infile){cout << "WARNING: unable to read guild_lookup.csv, " 
		                   << "not outputting guild info." << endl;  return 1;
		 }else {cout << "Loading guilds from guild_lookup.csv"  << endl;}
     CSVParser parser;		   
		 
     getline(infile,sLine);
		 parser << sLine;
     cols=0;  ins="in";
     while (ins != ""){
         parser >> ins;
		     if (ins == "") continue;
		     cols++;
		 }		   
     if (cols<2){cout << "Not enough columns in guild_lookup.csv" << endl;
                 return 1;}
     
		 NGuilds=0;            
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;    
		       parser >> species;
		       parser >> guild;   guildlist[species]   = guild;
		       parser >> instr;   trophiclist[species] = instr;
		       parser >> inint;   GuildOrder[inint]    = guild;
		       parser >> instr;   dattype[species]     = instr;
		       parser >> inint;   spcolor[species]     = inint;
		       parser >> instr;   // dummy

		       if (!guildcount[guild]){cout << "Making Guild " << guild << endl; 
					                         //GuildOrder[NGuilds]=guild;
					                         NGuilds++;
																	 guildnum[guild] = NGuilds;
					                         //sprintf(webFile,"./%s/Guild%.3d_index.html",web,guildnum[guild]);
																	 //filelist[guild]=fopen(webFile,w);
																	 //html_header(filelist[guild],guild.c_str());
																	 }
		       guildcount[guild]++;
		 }
  // Close the file
     infile.close();

}

//------------------------------------------------------------------------------
int commParse (int argc, char **argv)
  {

	int i, j, k;
  int error;
	float tmpflt;
  error=0;
  strcpy(OUTFOLDER,"./out");
  strcpy(INFOLDER,"./in");
  strcpy(CLIMATE_FILE,"no_climate.csv");
  strcpy(JUV_FILE,"no_juvs.csv");
  strcpy(FITTING_FILE,"no_fits.csv");
  strcpy(DIET_FILE,"no_diets.csv");
  strcpy(FIT_VECTOR_FILE, "no_start_fit.csv");
  strcpy(FIT_CONTROL_FILE, "no_start_control.csv");
  
	for (i = 1; i < argc; i++) { 
	 if (argv[i][0] != '-'){cout << "Option " << argv[i] << " unknown." << endl; error=1;}
	 else {
	  //cout << argv[i][1] << endl;
		switch (argv[i][1]) {
								
			// -M[char] Run Mode Switch
			case 'M':
			  RUNMODE = argv[i][2];
			  switch (RUNMODE){
			       case 'M':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][2] << " Must specify submode for mode runmode M"; error=1;}
						   else{
                   SUBMODE = argv[i][3];
               }
               cout << "setting run mode " << RUNMODE << " submode " << SUBMODE << endl;
							 break;
						 default:     
			          cout << "setting run mode " << RUNMODE << " with no submode" << endl;
			          if (RUNMODE=='l'){
								   if(strlen(&argv[i][3])>0){sscanf(&argv[i][3],"%d",&LONG_CYCLES);}
				           cout << "LONG_CYCLES set to " << LONG_CYCLES << endl;
									 }
			          break;
			  }
			  break;
			  
			// -S[n] run n random systems.  
			// Default is 0, which runs the base system.				
			case 'S':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " systems to run error."; error=1;}
			  else {sscanf(&argv[i][2],"%d",&SYSTEMS_TO_RUN);
				      cout << "Systems to Run set to " << SYSTEMS_TO_RUN << endl;
				}
			  break;

      // -B[n] burn time is n years.  Default 0 (no systems discarded).
			// During burn time, if a species collapses the system is discarded
			case 'B':
			  //cout << "f this" << endl;
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " burn time error."; error=1;}
			  else {sscanf(&argv[i][2],"%d",&BURN_TIME);
				      cout << "Burn time set to " << BURN_TIME << endl;
				}
			  break;

			case 'A':
			  //cout << "f this" << endl;
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " area multiplier error."; error=1;}
			  else {sscanf(&argv[i][2],"%f",&AREA_MULT);
				      cout << "Area Multiplier set to " << AREA_MULT << endl;
				}
			  break;


      case 'N':
               if   ((strlen(&argv[i][2])>79) || (strlen(&argv[i][2])<1)) 
						        {cout << argv[i][2] << " Noise error."; error=1;}
						   else {sscanf(&argv[i][2],"%f",&NOISE_RANGE);
						        cout << "NOISE_RANGE set to " << NOISE_RANGE << endl;
						        }       
        break;
      // -V[s] -V[h] Turn on Scramble [s] or handling time [h]
      case 'V':
         switch (argv[i][2]){
            case 'x':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][3] << " Vul error."; error=1;}
						   else {sscanf(&argv[i][3],"%f",&MSCRAMBLE);
						        cout << "MSCRAMBLE set to " << MSCRAMBLE << endl;
						        }
						   break;
						case 'd':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][3] << " Vul error."; error=1;}
						   else {sscanf(&argv[i][3],"%f",&MHANDLE);
						        cout << "MHANDLE set to " << MHANDLE << endl;
						        }
						   break;            
				    case 's':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][3] << " Vul error."; error=1;}
						   else {sscanf(&argv[i][3],"%f",&ScrambleSelfWt);
						        cout << "ScrambleSelfWt set to " << ScrambleSelfWt << endl;
						        }
						   break;
				    case 'h':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][3] << " Vul error."; error=1;}
						   else {sscanf(&argv[i][3],"%f",&HandleSelfWt);
						        cout << "HandleSelfWt set to " << HandleSelfWt << endl;
						        }
						   break;				    
				    case 'v':
               if   ((strlen(&argv[i][3])>79) || (strlen(&argv[i][3])<1)) 
						        {cout << argv[i][3] << " Vul error."; error=1;}
						   else {sscanf(&argv[i][3],"%f",&VULVAR);
						        cout << "Vulvar set to " << VULVAR << endl;
						        }
						   break;	
				    default:
				       break;
				    }
				 break;
                
      // -O[folder] Save output files to [folder]
			case 'O':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " outfile too long."; error=1;}
			  else {strcpy(OUTFOLDER,&argv[i][2]);
				      cout << "writing to folder " << OUTFOLDER << endl;
				}
			  break;
      // -I[folder] Read input from [folder]
			case 'I':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " outfile too long."; error=1;}
			  else {strcpy(INFOLDER,&argv[i][2]);
				      cout << "Reading from folder " << INFOLDER << endl;
				}
			  break;


			//case 'P':

      // -C[file.csv] read climate from [file.csv], must be comma-delimited
			case 'C':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " climate file too long."; error=1;}
			  else {strcpy(CLIMATE_FILE,&argv[i][2]);
				      cout << "Climate from " << CLIMATE_FILE << endl;
				}
			  break;

      // -J[file.csv] read juvenile parameters from [file.csv], must be comma-delimited
			case 'J':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " juv file too long."; error=1;}
			  else {strcpy(JUV_FILE,&argv[i][2]);
				      cout << "Juveniles from " << JUV_FILE << endl;
				}
			  break;

      // -F[file.csv] read fitting parameters from [file.csv], must be comma-delimited
      case 'F':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " fitting file too long."; error=1;}
			  else {strcpy(FITTING_FILE,&argv[i][2]);
				      cout << "Fitting from " << FITTING_FILE << endl;
				}
      	break;

      // -D[file.csv] read diet fitting parameters from [file.csv], must be comma-delimited
      case 'D':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " diet fitting file too long."; error=1;}
			  else {strcpy(DIET_FILE,&argv[i][2]);
				      cout << "Diet fitting from " << DIET_FILE << endl;
				}
      	break;

      // -P[file.csv] read start fitting vectors from [file.csv], must be comma-delimited
      case 'P':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " start fit too long."; error=1;}
			  else {strcpy(FIT_VECTOR_FILE,&argv[i][2]);
				      cout << "Fitting starts from " << FIT_VECTOR_FILE << endl;
				}
      	break;

      // -K[file.csv] read fitting control vectors from [file.csv], must be comma-delimited
      case 'K':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " start fit control too long."; error=1;}
			  else {strcpy(FIT_CONTROL_FILE,&argv[i][2]);
				      cout << "Fit controls from " << FIT_CONTROL_FILE << endl;
				}
      	break;

      // -R[n] Random seed as an integer.  Default is 0.
			// range is 0-511, numbers outside this range are modded to get values inside range.
			// Number is lookup in RandomSeeds.h, and each refers to a 256bit random # generated
			// from atmospheric noise (http://Random.org)
			case 'R':
			  if   (strlen(&argv[i][2])>79){cout << argv[i][2] << " rand seed error."; error=1;}
			  else {sscanf(&argv[i][2],"%lu",&rseed);	//SKG changed to %lu from %d to fix compiler warning June 2012
				      cout << "Random Seed set to " << rseed << endl;
				}
			  break;

			case '\0':
				break;

			/* Illegal options */
			
			default:
				cout << "Option " << argv[i][1] << " unknown." << endl; 
				error=1;
				break;
		  }
	 }
	}
	return error;
}						


