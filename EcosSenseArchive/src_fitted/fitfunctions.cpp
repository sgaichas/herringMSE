//
//#define NUM_VARS   5
//#define NUM_ROWS   NUM_LIVING
//#define NDIM       NUM_VARS*NUM_ROWS     // Number of parameters in maximization vector
#define GTOL      0.01                   // Min %change of a parameter at final solution step

#define EPS       3.0e-8
#define TOLX      (4*EPS)
#define STPMX     1.0          //past runs at 0.1, recent at 0.5
#define ALF       1.0e-4
#define BIGSSQ    1000000.0
#define DERIV_STEP 0.1
//#define GTOL      1.0e-4
#define ONEHALFLOGTWOPI 0.918938533204672

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
                  
int ITMAX   =  200; 
double **dSdX;
int TESTFLAG = 0;
FILE *tptr;
int BigIt =0;
int rr[NDIM+1];

const char *VarNames[NUM_VARS] = {"PB","M0",
                            "PyHand","PdHand",
                            "PyVul","PdVul",
                            "PySwitch","PdSwitch"};

const char *FitNames[4]       = {"BioFit","CatFit","PathFith","ForceFit"};
// -----------------------------------------------------------------------------
void load_control(struct SimRun *v)
{
int sp,i;
string sout, sLine;
//char STEP_FILE[80];
double flout;

     //sprintf(STEP_FILE,"%s_fitVectorsStart.csv",OUTFOLDER);
     
     for (i=0; i<=NDIM; i++){v->fit_vec_control[i]=1.0;}
     for (sp=0; sp<=NUM_GROUPS; sp++){
          for (i=0; i<=4; i++){
              v->fit_ssq_control[i][sp] = 1.0;
          }
     }
     ifstream infile(FIT_CONTROL_FILE); 
     if (!infile){cout << "No step parameters loaded, unable to read " 
		                   << FIT_CONTROL_FILE << endl;  return;
		 }else {cout << "Fit controller loaded from " << FIT_CONTROL_FILE << endl;}
     CSVParser parser;

     getline(infile,sLine);
		 parser << sLine;

     sp=0;
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;
		       parser >> sout;  parser >> sout;
					 sp++;
					 for (i=0; i<NUM_VARS; i++){
               parser >> flout; 
               v->fit_vec_control[i*NUM_LIVING+sp] = flout;
               //cout << i*NUM_LIVING+sp << flout << endl;
           }
           for (i=0; i<=4; i++){
               parser >> flout;
               v->fit_ssq_control[i][sp] = flout;
           }

		 }
		 // END REC_CHANGE		            
     infile.close();

}
// -----------------------------------------------------------------------------
void bin_outstep(int st, double *p, struct SimRun *v){
char outFile[80]; 
FILE *binptr;

  sprintf(outFile,"%s_fitVectorsStep%d.bin",OUTFOLDER,st);
  binptr=fopen(outFile,"wb");

     //for (sp=1; sp<=NUM_LIVING; sp++){
     //    fprintf(fptr,"%d,%s,",sp,path_species[sp]);
     //    for (i=0; i<NUM_VARS; i++){fprintf(fptr,"%g,",p[i*NUM_LIVING+sp]);}
     //    for (i=0; i<=4; i++){fprintf (fptr,"%g,", v->fit_SSQ[i][sp]);}
     //    fprintf (fptr,"\n");
     //}
     //fclose(fptr);

     fwrite( p , 1, sizeof(double)*(NDIM+1), binptr);
		 fclose(binptr);
}
// -----------------------------------------------------------------------------
void bin_instep(struct SimRun *v){
FILE *binptr;

  load_control(v);
  binptr=fopen(FIT_VECTOR_FILE,"rb");
  if (binptr==NULL){cout << "No Bin Vector file loaded: " << FIT_VECTOR_FILE << endl;return;}
	else {cout << "Binary fit parameters loaded from " << FIT_VECTOR_FILE << endl;}
  fread(v->fit_vector,1, sizeof(double)*(NDIM+1), binptr);
  fclose(binptr);

}

// ----------------------------------------------------------------------------
void load_instep(struct SimRun *v)
{
int sp,i;
string sout, sLine;
//char STEP_FILE[80];
double flout;

     //sprintf(STEP_FILE,"%s_fitVectorsStart.csv",OUTFOLDER);
     load_control(v);

     for (i=0; i<=NDIM; i++){v->fit_vector[i]=0.0;}
     ifstream infile(FIT_VECTOR_FILE); 
     if (!infile){cout << "No step parameters loaded, unable to read " 
		                   << FIT_VECTOR_FILE << endl;  return;
		 }else {cout << "Fit parameters loaded from " << FIT_VECTOR_FILE << endl;}
     CSVParser parser;

     getline(infile,sLine);
		 parser << sLine;

     sp=0;
     while (!infile.eof()) {
           getline(infile,sLine); if (sLine == "") continue; 
		       parser << sLine;
		       parser >> sout;  parser >> sout;
					 sp++;
					 for (i=0; i<NUM_VARS; i++){
               parser >> flout; 
               v->fit_vector[i*NUM_LIVING+sp] = flout;
           }

		 }
		 // END REC_CHANGE		            
     infile.close();
     
}

// -----------------------------------------------------------------------------

void outstep(int st, double *p, struct SimRun *v){
int sp,i,c;
char outFile[80]; 
FILE *fptr;

  // SAVE diagnostics and final p values in fitVectors file 
     sprintf(outFile,"%s_fitVectorsStep%d.csv",OUTFOLDER,st);
     fptr=fopen(outFile,"w");
     fprintf(fptr,"Num,Species,");
     for (i=0; i<NUM_VARS; i++){fprintf (fptr,"p%d,",i);}
     for (i=0; i<=4; i++){fprintf (fptr,"SSQ%d,",i);}
     fprintf (fptr,"\n");
     
     for (sp=1; sp<=NUM_LIVING; sp++){
         fprintf(fptr,"%d,%s,",sp,path_species[sp]);
         for (i=0; i<NUM_VARS; i++){fprintf(fptr,"%.36f,",p[i*NUM_LIVING+sp]);}
         for (i=0; i<=4; i++){fprintf (fptr,"%.36f,", v->fit_SSQ[i][sp]);}
         fprintf (fptr,"\n");
     }
     fclose(fptr);
}

//----------------------------------------------------------------------
  // SAVE diagnostics and final p values in fitVectors file
void outderiv(struct SimRun *v, FILE *fptr, double *p, double *dp, char* rlabel){ 
  int i,sp,c;
  double endTot;
  
      for (i=0; i<NUM_VARS; i++){      
        for (sp=1; sp<=NUM_ROWS; sp++){
          fprintf(fptr,"%s,%s,%d,%f,%f,%f,",
					        rlabel,path_species[sp],i,p[sp+NUM_ROWS*i],
									                        v->fit_vec_control[sp+NUM_ROWS*i],
																					dp[sp+NUM_ROWS*i]);
          endTot=0.0;
					for (c=1; c<=fitN+NUM_LIVING+NUM_DEAD; c++){
              if((c<=fitN) || (FORCE_BIO_FLAG[c-fitN])){
               fprintf(fptr,"%f,",dSdX[sp+NUM_ROWS*i][c]);
               endTot+=dSdX[sp+NUM_ROWS*i][c];
              }
          }
          fprintf(fptr,"%g,\n",endTot);
        }
      }
     
}
//------------------------------------------------------------------------------

void apply_vector_to_rates(struct SimRun *v, double q[])
{
int i, sp, links, prey, pred;
double PreyX[NUM_LIVING+NUM_DEAD+1],      PredX[NUM_LIVING+NUM_DEAD+1];
double PredHandle[NUM_LIVING+NUM_DEAD+1], PreyHandle[NUM_LIVING+NUM_DEAD+1];
double PredExp[NUM_LIVING+NUM_DEAD+1],    PreyExp[NUM_LIVING+NUM_DEAD+1];
double x[NDIM+1];     
     
     for (i=0; i<=NDIM; i++){
         x[i]=q[i];
         if (x[i]>40.0){x[i]=40.0;}
     }
     
     for (sp=1; sp<=NUM_LIVING; sp++){
      //oldvals[i] = rpar_ActiveRespFrac[i];
         i=NUM_LIVING*0 + sp;
         //cout << i << ","<< sp << ","<< v->RRR.rpar_ActiveRespFrac[sp];
         if (path_QB[sp] > EPSILON){
            	v->RRR.rpar_ActiveRespFrac[sp]  = 1.0 - ( path_PB[sp] * 
					                                (( 2.0 * exp(x[i]) /(1+exp(x[i])))) 
					                                /path_QB[sp])
							                          - v->RRR.rpar_UnassimRespFrac[sp]; 
         MIN_THRESHOLD(v->RRR.rpar_ActiveRespFrac[sp],0.01);
         }
         
         //i=NUM_LIVING*1 + sp;
         //v->RRR.rpar_B_BaseRef[sp] *= exp(x[i]);
         
        i=NUM_LIVING*1 + sp;
        //v->RRR.rpar_MzeroMort[sp] *= (1.0 + x[i]);
        v->RRR.rpar_MzeroMort[sp] *= (( 2.0 * exp(x[i]) /(1+exp(x[i]))));        
//        MIN_THRESHOLD(v->RRR.rpar_MzeroMort[sp] , 0.0);
//        cout << i << " " << x[i] << " " <<  (( 2.0 * exp(x[i]) /(1+exp(x[i])))) << endl;

// Make relationship between x and actual vulnerabilties a logistic
// curve that is near-linear in a range that for VTOL = 0.1, is near
// half the limits
#define VTOL 0.1                  // sets linearity range
#define VLIM 10.0                 // sets min and max
#define LAMBDA 0.219722457733622  // -log(VTOL/(1-VTOL))/VLIM
#define VFUNC(r) -VLIM + 2.0*VLIM*exp(LAMBDA*(r))/(1.0 + exp(LAMBDA*(r)));   

        i=NUM_LIVING*2 + sp;
        PreyHandle[sp] = VFUNC(x[i]);  

        i=NUM_LIVING*3 + sp;
        PredHandle[sp] = VFUNC(x[i]); 

        i=NUM_LIVING*4 + sp;
        PreyX[sp] = VFUNC(x[i]);
         
        i=NUM_LIVING*5 + sp;
        PredX[sp] = VFUNC(x[i]);   

        i=NUM_LIVING*6 + sp;
        PreyExp[sp] = VFUNC(x[i]);
        
        i=NUM_LIVING*7 + sp;
        PredExp[sp] = VFUNC(x[i]);
        //v->RRR.rpar_HandleSwitch[sp] = 0.5 + (( 2.0 * exp(x[i]) /(1+exp(x[i])))); 

     }
     PreyX[0] = 0;
     PredX[0] = 0;
     PreyHandle[0] =0;
     PredHandle[0] =0;
     PredExp[0]=0;
     PreyExp[0]=0;
    for (sp=NUM_LIVING+1; sp<=NUM_LIVING+NUM_DEAD; sp++){
     PreyX[sp] = 0;
     PredX[sp] = 0;
     PreyHandle[sp] =0;
     PredHandle[sp] =0;
     PreyExp[sp]=0;
     PredExp[sp]=0;
    }
    
     for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){	
         prey=v->RRR.rpar_PreyFrom[links];
         pred=v->RRR.rpar_PreyTo[links];
         
	       //v->RRR.rpar_DD[links] = 1. + (MHANDLE-1.)       * exp(HandleX[pred]);
         //v->RRR.rpar_VV[links] = 1. + (MSCRAMBLE-1.)     * exp(PredX[pred]+PreyX[prey]);
	       v->RRR.rpar_DD[links] = 1. +  exp(PredHandle[pred]+ PreyHandle[prey]);
         v->RRR.rpar_VV[links] = 1. +  exp(PredX[pred]     + PreyX[prey]);
         //v->RRR.rpar_HandleSwitch[links] = 1. + PredExp[pred] + PreyExp[prey];
        v->RRR.rpar_HandleSwitch[links] = exp(0.05 * (PredExp[pred] + PreyExp[prey]));
   }	
     

}
// -----------------------------------------------------------------------------

double SSQ_run(struct SimRun *v, double x[])
{
double ssq;
int i;

  // First set default parameters for sim
     path_to_rates(v);

  // Then apply x vector to the default non-juvenile parameters.
  
  // This applies X[1..NumLiving] to PB
     apply_vector_to_rates(v, x);
  
  // Then initialize stanzas
     initialize_stanzas(v);
     
  // Then apply parts of x vector that override stanza calcs.   
     
  // Then run the model and caluclate SSQ
     Adams_Basforth(v,0, fit_Years);
     ssq = SSQ(v);

  //if (!DISCARD_YEAR){ssq = SSQ();}
  //else             {ssq = BIGSSQ;}

  //printf ("A run: %g, %d\n",ssq, v->DISCARD_YEAR);
  return ssq;
 
}
//------------------------------------------------------------------------------

void *deriv_thread(void *threadarg)
{
   struct FitRun *F;
   double old, newSSQ;
   
   // FitRun structure is a SimRun structure with a couple of extra variables
   // (current SSQ, derivative step)
   
   // This line just copies the whole FitRun passed via a thread variable
   // into F (the fit run result structure)
      F = (struct FitRun *) threadarg;

   // If the control variable (from the control file) is 0, we're not fitting
   // to it, so we return a derivative of 0, otherwise we calculate the derivative
      if (F->v->fit_vec_control[F->dterm] > 0.0001) {

      // F->dterm is the index parameter number we're calculating the derivative for.
      // F->x[F->dterm] is the current value of that parameter.
      
      // Save the old paramter value
         old = F->x[F->dterm];
      // Temporarily increase F->x by the derivative step size.
         F->x[F->dterm] += DERIV_STEP;

      // With F->x bigger by DERIV_STEP, run the SSQ function
         //GIANT = 1; 
				 sprintf(F->v->RunID,"_%d_diff",F->dterm);
         newSSQ=SSQ_run(F->v, F->x);
         //GIANT=0;
      // Diagnostics 
         //int c;
			   //char str[10000]; 
         //sprintf(str,"Ugh,%d,%s,%.10g,%.10g,%.10g,%.10g,%.10g",
         // F->dterm , path_species[F->dterm%NUM_LIVING] 
				 //,newSSQ , F->baseSSQ , old , F->x[F->dterm] ,  
				 //(newSSQ - F->baseSSQ) / (F->x[F->dterm] - old));
				 //for (c=1; c<=fitN; c++){sprintf(str,"%s,%.10g",str, F->v->fit_SSQfit[c]);}
         // for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
			 	 // if(FORCE_BIO_FLAG[c]){
    	 	 //			 sprintf(str,"%s,%.10g",str,F->v->fit_SSQ[3][c]);                                }
         //  }
				 //sprintf(str,"%s,\n",str);
				 //cout << str;

      // Then derivative is new SSQ - base, divided by step size
      // F->df = (newSSQ - F->baseSSQ) / (F->x[F->dterm] - old);
			// KYA 8/5/08: setting the derivative in the above line was summing 
			// up the pieces with poor precision.  Now, the df is not calculated 
			// here, but the pieces of the derivative are saved in the run stats.  
			// This is returned to the main function to assemble the derivative 
			// vector.  For now, we just set df to 0.0 (while keeping the parts 
			// after the run).
         F->df = 0.0;
      // reset the saved old parameter value
         F->x[F->dterm] = old;
      }
      else {
         F->df = 0.0;
      }   
   pthread_exit(NULL);

}

//------------------------------------------------------------------------------

// Machine-dependent float truncation function, assumes fp locally declared
//#define TRUNC(ff) fp=(unsigned int *)&(ff);  fp[0]=(fp[0]&0xFFFFFC00);

void SSQ_deriv(struct SimRun *v, double x[],double df[])
{
   int r, rc, i, status, cyc, c, sp,t;
   int fun;
   //double *basefit;
   double *basedf;
   double ssq;
   //float baseforce, curforce;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];
   struct FitRun Runs[PRUNS];
   struct FitRun BaseRun;
   //unsigned int *fp;
   //unsigned int *fp1,*fp2;
   //char str[10000];
   char outFile[80]; 
   //FILE *qptr;
   double baseval, delval;
   
	 cout << "Starting Derivative..." << endl;
   
      //sprintf(outFile,"%s_FitBits%d.csv",OUTFOLDER,BigIt);
      //qptr=fopen(outFile,"w");
   
        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);   

        // basefit will hold the base (pre-step) SSQ parts for making the
        // derivative calculations.
           //basefit   =dvector(1,fitN+NUM_LIVING+NUM_DEAD);
           basedf    =dvector(1,NDIM);
           
					 BaseRun.v = new_SimRun();
           BaseRun.x = dvector(1,NDIM);
                  for (i=1; i<=NDIM; i++){
                      BaseRun.x[i]=x[i];
                      BaseRun.v->fit_vec_control[i] = v->fit_vec_control[i];
                  }
                  for (i=0; i<=4; i++){
                      for (sp=1; sp<=NUM_LIVING; sp++){
                          BaseRun.v->fit_ssq_control[i][sp] = v->fit_ssq_control[i][sp];
                      }
                  }   
           //GIANT = 1; sprintf(BaseRun.v->RunID,"_b_base");
           BaseRun.baseSSQ = SSQ_run(BaseRun.v, BaseRun.x);
           //GIANT = 0;
           
           //sprintf(outFile,"Rr_term0000");
           //output_run(BaseRun.v,outFile);

        // Allocate one ecosystem for each thread, and assign the starting vector
           for (r=0; r<PRUNS; r++){
               // Initialize the FitRun variables, one per thread with copy
               // of ecosystem and fitting vector.
                  Runs[r].v       = new_SimRun();
                  Runs[r].x       = dvector(1,NDIM);
               // Copy various controlling information into each threads' FitRun. 
                  for (i=1; i<=NDIM; i++){
                      Runs[r].x[i]=x[i];
                      Runs[r].v->fit_vec_control[i] = v->fit_vec_control[i];
                  }
                  for (i=0; i<=4; i++){
                      for (sp=1; sp<=NUM_LIVING; sp++){
                          Runs[r].v->fit_ssq_control[i][sp] = v->fit_ssq_control[i][sp];
                      }
                  }                
               // Do a base run for the first ecosystem to set the base SSQ and
							 // SSQ parts from which the derivative is calculated.  
							 // FitRuns other than 0, do a run to make sure everything is 
							 // initialized, zeroed, etc., though values aren't saved for this.
							    //GIANT = 1; sprintf(Runs[r].v->RunID,"_%d_base",r);
							    Runs[r].baseSSQ = SSQ_run(Runs[r].v, Runs[r].x);
                  //GIANT = 0;
                      // Diagnostics
                      // sprintf(str,"Ugh,%d,Base,%.10g,%.10g,%.10g,%.10g,%.10g",
                      // r,Runs[r].baseSSQ,Runs[r].baseSSQ ,0,0,0);
					            // for (c=1; c<=fitN; c++){sprintf(str,"%s,%.10g",str, Runs[r].v->fit_SSQfit[c]);}
                      // for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
			 		            // if(FORCE_BIO_FLAG[c]){
    	 				        // sprintf(str,"%s,%.10g",str,Runs[r].v->fit_SSQ[3][c]);
                      //          }
                      //  }
					            //  sprintf(str,"%s,\n",str);
					            //cout << str;
							 
							 // For the first run (thread 0) save values in baserun, for
							 // other runs nothing needs to be done (I think).							 
                   if (r==0){ 
                      for (i=1; i<=fitN; i++){
											     //basefit[i]=Runs[r].v->fit_SSQfit[i];
											}
                      for (i=1; i<=NUM_LIVING+NUM_DEAD; i++){
                           //basefit[i+fitN] = Runs[r].v->fit_SSQ[3][i];
                      }
                   }
                   else{ /* do nothing I think */ }

           }  // End of Initializing ecosystems for each thread loop
        //GIANT = 1;
        // Loop through all variables in chunks of PRUNS threads
           for (i=0; i<NDIM; i+=PRUNS){

          // Spawn the threads into deriv_thread routine
             for (r=0; r<PRUNS; r++){
             //  Mapping a thread to a variable, if statement for the last few if not divisible
                 fun = PRUNS-r-1;
                 Runs[r].dterm = i+fun+1;
                 if (Runs[r].dterm <=NDIM){
                    rc = pthread_create(&threads[r],&attr, deriv_thread, (void *) &Runs[r]);
                 }
             }
      
          // Wait for the threads to return
             for (r=0; r<PRUNS; r++){
                 if (Runs[r].dterm <=NDIM){   
                     rc=pthread_join(threads[r], (void **)&status);
                     if (rc){cout << r << "error" << endl;}
                 }
             }              
        
          // Now assemble the derivative pieces into a vector 
             for (r=0; r<PRUNS; r++){
                 if (Runs[r].dterm <=NDIM){ 
                      // DIAGNOSTICS by setting TESTFLAG
                         //if(TESTFLAG){
                  	     //   sprintf(outFile,"RunFlag%d",Runs[r].dterm );
                         //   output_run(Runs[r].v,outFile);
                         // }

                      // Save processor number for diagnostics
                         rr[Runs[r].dterm]=r;

                      // Now Load the derivative pieces into a single
                      // float vector DF, indexed by dterm
                         basedf[Runs[r].dterm] = 0.0;
                         
                      // recalculating SSQ here, does this help ???
                         ssq=SSQ(Runs[r].v);
                      
                      // Add derivative parts if control vector (fitting ON) 
											// for term is positive, otherwise derivative part is 0
                         if (Runs[r].v->fit_vec_control[Runs[r].dterm] > 0.0001){
                            // Fitting parts
                               for (c=1; c<=fitN; c++){
                                   baseval = BaseRun.v->fit_SSQfit[c]; //TRUNC(baseval);
                                   delval  = Runs[r].v->fit_SSQfit[c]; //TRUNC(delval);
                                   
                                   basedf[Runs[r].dterm] += (delval-baseval);
                                   dSdX[Runs[r].dterm][c] = (delval-baseval)/DERIV_STEP;
                                   //fp1 = (unsigned int *)&delval;
                                   //fp2 = (unsigned int *)&baseval;
                                   //fprintf(qptr,"%d,fit,%d,%d,%X,%X,%.20g,\n",r,Runs[r].dterm,c,
																	 //            fp1[0],
																	//						 fp2[0],
																	//						 delval-baseval);
																	 //cout << r << ",fit," << Runs[r].dterm << "," << c 
																	 //     << "," << Runs[r].v->fit_SSQfit[c] << "," << basefit[c] << endl;
				    								   }
				    								// Forcing parts 
                               for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
                                   baseval = BaseRun.v->fit_SSQ[3][c]; //TRUNC(baseval);
                                   delval  = Runs[r].v->fit_SSQ[3][c]; //TRUNC(delval);                               
                                   basedf[Runs[r].dterm] += (delval-baseval);
                                   dSdX[Runs[r].dterm][c+fitN] = (delval-baseval)/DERIV_STEP;
                                   //fp1 = (unsigned int *)&delval;
                                   //fp2 = (unsigned int *)&baseval;
                                   //fprintf(qptr,"%d,force,%d,%d,%X,%X,%.20g,\n",r,Runs[r].dterm,c,
																	 //            fp1[0],
																	//						 fp2[0],
																	//						 delval-baseval);                                 
																	 //cout << r << ",force," << Runs[r].dterm << "," << c 
																	 //     << "," << Runs[r].v->fit_SSQ[3][c] << "," << basefit[c+fitN] << endl;
                               }
                            // Divide whole derivative by Step Size     
                               basedf[Runs[r].dterm] /= DERIV_STEP;
                         }
                         
                         //cout << Runs[r].dterm << "," << basedf[Runs[r].dterm] << endl;
												  
                         df[Runs[r].dterm] =(double)basedf[Runs[r].dterm];
                         //TRUNC(df[Runs[r].dterm]);
                     // Temporary diagnostics
                         //if ((BigIt==0) && (Runs[r].dterm<=24)){
                         //  BaseRun.baseSSQ = SSQ_run(BaseRun.v, BaseRun.x);
                         //  sprintf(outFile,"Rr_term0000%d",Runs[r].dterm);
                         //  output_run(BaseRun.v,outFile);                        
                  	     //  sprintf(outFile,"Rr%d_term%d",r,Runs[r].dterm);
                         //  output_run(Runs[r].v,outFile);
                         //  if(Runs[r].dterm==24){exit(0);}
                         //}
                     // Save the derivative parts in a global dSdX matrix, this is
										 // entirely for printing diagnostics and should not be used for
										 // further calculations 
											   //for (c=1; c<=fitN; c++){
                         //     dSdX[Runs[r].dterm][c] = (Runs[r].v->fit_SSQfit[c]-basefit[c])/DERIV_STEP; //NOTE HARDCODED DERIVATIVE STEP
                         //}                    
                         //for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
                         //    dSdX[Runs[r].dterm][c+fitN] =
												//		    (Runs[r].v->fit_SSQ[3][c]-basefit[c+fitN])/DERIV_STEP;														    
												//				//if(FORCE_BIO_FLAG[c]){
												//				//  cout << Runs[r].dterm << "," << Runs[r].df << "," << c << "," << path_species[c]<< "," << Runs[r].v->fit_SSQ[3][c] << "," <<  basefit[c+fitN] << endl;
                        //        //}
                        // }
                      
											//YET MORE DIAGNOSTICS
											    //for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
                          //if(FORCE_BIO_FLAG[c]){
                 					//fprintf(tptr,"%d,%d,%d,%s,%s,%g,%g,%g,",BigIt,r,Runs[r].dterm,
                 					//                              path_species[Runs[r].dterm%NUM_LIVING],
													//	                             path_species[c],
													//	                             Runs[r].v->fit_SSQ[3][c],
							            //                               Runs[r].v->fit_ssq_control[3][c],
												  //																 Runs[r].v->fit_FORCEraw[c]); 
												  //	for (t=0; t<=fit_lastyear; t++){
                          //     fprintf(tptr,"%g,",v->fit_bioanom[c][t]);
                          // }
                          // fprintf(tptr,"\n");                     
                          //}
                          //}
                     
                 }  // END OF if dterm<NDIM loop

             // EVEN YET MORE DIAGNOSTICS
             //for (r=0; r<PRUNS; r++){
               /* if (Runs[r].dterm <=NDIM){
                 cout <<"PartSSQ," << i << "," << r << "," << Runs[r].dterm << "," << df[Runs[r].dterm] << ",";
                      cout << curforce << ",";
                      for (c=1; c<=fitN; c++){
                           cout << Runs[r].dSdx[c] << ",";
                      }
                 cout << endl;
                 cout <<"BasePartSSQ," << i << "," << r << "," << Runs[r].dterm << "," << df[Runs[r].dterm] << ",";
                      cout << baseforce << ",";
                      for (c=1; c<=fitN; c++){
                           cout << basefit[c] << ",";
                      }
                 cout << endl;
                 cout <<"DerivPart," << i << "," << r << "," << Runs[r].dterm << "," << df[Runs[r].dterm] << ",";
                      cout << (curforce-baseforce)/0.1 << ",";
                      for (c=1; c<=fitN; c++){
                           cout << dSdX[Runs[r].dterm][c] << ",";
                      }
                 cout << endl;                 
                 }*/
                 
             // Every so often print an update  
                if ((i+r+1)%NUM_LIVING==0){cout << "var: " << (i+r+1)/NUM_LIVING << " of " << NUM_VARS <<endl;}
             }  // END OF R LOOP           
             
          
          }  // END OF FULL DERIVATIVE LOOP

          // Print out vector
             //     for (i=0; i<=NDIM; i++){cout << i << " " << df[i] << endl;}
        //fclose(qptr);
        // Thread cleanup
           //free_dvector(basefit,1,fitN+NUM_LIVING+NUM_DEAD);
           free_dvector(basedf,1,NDIM);
           for (r=0; r<PRUNS; r++){
               //free_vector(Runs[r].dSdx,1,fitN);
               free_dvector(Runs[r].x,1,NDIM);
               //free(Runs[r].v);
               free_SimRun(Runs[r].v);
           } 
           free_SimRun(BaseRun.v);
           free_dvector(BaseRun.x,1,NDIM);
           
           pthread_attr_destroy(&attr);    
					            

}
#undef TRUNC
// -----------------------------------------------------------------------------
//#define RAND_RUNS 200
//#define RAND_RANGE 0.5
void *scatter_thread(void *threadarg)
{
   struct ScatterRun *F;
   int t, i, sp;
   int accept;
   //float newSSQ;
   double *randX;
   double *currentX;
   double *tmp;
   double SSQtest;
   double a_prob, p;
   //float baseSSQ;
   
   F       = (struct ScatterRun *) threadarg;
   currentX = dvector(1,NDIM);
   randX    = dvector(1,NDIM);
   accept   = 0;
   
   for (i=1; i<NDIM; i++){
       currentX[i] = F->x[i];
   }

   for (t=1; t<=SYSTEMS_TO_RUN; t++){
     for (i=1; i<=NDIM; i++){
       if (F->v->fit_vec_control[i] > 0.0001) {
          randX[i] = currentX[i] + (2*NOISE_RANGE*uniform(&F->v->rng)-NOISE_RANGE);
       }
       else{
          randX[i] = currentX[i];
       }
     }
     
     SSQtest      = SSQ_run(F->v, randX);

     a_prob = 1.0/(SSQtest/F->SSQlist[t-1]);
     p      = uniform(&F->v->rng);
     if ((p<a_prob) && (!F->v->DISCARD_YEAR)){
        accept++;
        tmp      = currentX;
        currentX = randX;
        randX    = tmp;
        //cout << F->v->thread_ID << " " << t <<" "<< accept <<" "<< (float)a_prob/(float)t 
        //     <<" "<< a_prob <<" "<< p <<" "<< SSQtest <<" "<< F->SSQlist[t-1] << endl;
     }  
     else {
        SSQtest = F->SSQlist[t-1];
     }      
     
     F->SSQlist[t]= SSQtest;
     for (i=1; i<=NDIM; i++){
         F->Xlist[i][t]=currentX[i];
     }         
     
     if (t%10 == 0){ 
          cout << F->v->thread_ID << " " << t <<" "<< accept <<" "<< 
               (double)accept/(double)t << endl;
     }
     
     for (sp=1; sp<=NUM_LIVING; sp++){
         //fprintf(fptr,"%d,%s,",sp,path_species[sp]);
         //for (i=0; i<NUM_VARS; i++){fprintf(fptr,"%g,",p[i*NUM_LIVING+sp]);}
         for (i=0; i<=4; i++){F->Slist[i*NUM_LIVING+sp][t] = F->v->fit_SSQ[i][sp];}
         //fprintf (fptr,"\n");
     }     
     
    //if (F->SSQlist[t] < BIGSSQ){
    //    cout << F->v->thread_ID << " " << t << " " << F->SSQlist[t] << endl;
    //}
    
   // end of time loop
      }
   
   free_dvector(randX,1,NDIM);
   free_dvector(currentX,1,NDIM);
   
   pthread_exit(NULL);

}
//#undef RAND_RANGE
//------------------------------------------------------------------------------
void SSQ_scatter(struct SimRun *v, double x[])
{
   int r, rc, i, status, cyc, c, sp,t,j;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];
   struct ScatterRun Runs[PRUNS];
   char outFile[80]; 
   FILE *dptr[NUM_LIVING+1], *sptr;
   const char* w = "w";
   double *Scount;
   cout << "Starting a Scramble..." << endl;
   
        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);   

        // Allocate one ecosystem for each thread, and assign the starting vector
           Scount =  dvector(1,NDIM);
           for (r=0; r<PRUNS; r++){
                //Runs[r].v       = (struct SimRun *)calloc(MAX_SYSTEMS,sizeof(struct SimRun));
                Runs[r].v       = new_SimRun();
                Runs[r].x       = dvector(1,NDIM);
                Runs[r].SSQlist = dvector(0,SYSTEMS_TO_RUN);
                Runs[r].Xlist   = dmatrix(1,NDIM,0,SYSTEMS_TO_RUN);
                Runs[r].Slist   = dmatrix(1,NDIM,0,SYSTEMS_TO_RUN);
                
                init_by_array(rseed+r, &Runs[r].v->rng); 
                Runs[r].v->thread_ID=rseed+r;
                
                for (i=1; i<=NDIM; i++){
                   Runs[r].x[i]=x[i];
                   Runs[r].v->fit_vec_control[i] = v->fit_vec_control[i];
                }
                for (i=0; i<=4; i++){
                  for (sp=1; sp<=NUM_LIVING; sp++){
                       Runs[r].v->fit_ssq_control[i][sp] = v->fit_ssq_control[i][sp];
                  }
                }                
                
                if (r==0){
                   Runs[r].baseSSQ    = SSQ_run(Runs[r].v, Runs[r].x);
                   Runs[r].SSQlist[0] = Runs[r].baseSSQ;
                   cout << Runs[r].baseSSQ << endl;
                   for (sp=1; sp<=NUM_LIVING; sp++){
                       for (i=0; i<=4; i++){
                           Runs[r].Slist[i*NUM_LIVING+sp][0] = 
                           Runs[r].v->fit_SSQ[i][sp];
                       }
                   }                    
                }
                else{
                   Runs[r].baseSSQ    = Runs[0].baseSSQ;
                   Runs[r].SSQlist[0] = Runs[r].baseSSQ;
                   cout << Runs[r].baseSSQ << endl;
                   for (sp=1; sp<=NUM_LIVING; sp++){
                       for (i=0; i<=4; i++){
                           Runs[r].Slist[i*NUM_LIVING+sp][0] = 
                           Runs[0].Slist[i*NUM_LIVING+sp][0];
                       }
                   } 
                }
           }
        
        // Loop through all variables in chunks of PRUNS threads
           //for (i=0; i<NDIM; i+=PRUNS){

          // Spawn the threads
             for (r=0; r<PRUNS; r++){
             //  Mapping a thread to a variable, if statement for the last few if not divisible
                 //Runs[r].dterm = i+r+1;
                 //if (Runs[r].dterm <=NDIM){
                    rc = pthread_create(&threads[r],&attr, scatter_thread, (void *) &Runs[r]);
                 //}
             }

          // Wait for the threads to return
             for (r=0; r<PRUNS; r++){
                 //if (Runs[r].dterm <=NDIM){   
                     rc=pthread_join(threads[r], (void **)&status);
                     if (rc){cout << r << "error" << endl;}
                // }
             //  Every so often print an update  
             
                 //if ((i+r+1)%NUM_LIVING==0){cout << "var: " << (i+r+1)/NUM_LIVING << " of " << NUM_VARS <<endl;}
             }                             
          //}
        // Thread cleanup
          for (r=0; r<PRUNS; r++){
               for (i=1; i<=SYSTEMS_TO_RUN; i++){
                   cout << r << " " << i << " " << Runs[r].SSQlist[i] << endl;
               }     
          }
          
        sprintf(outFile,"%s_sum_scramble.csv",OUTFOLDER,sp);
        sptr=fopen(outFile,w);
        for (sp=1; sp<=NUM_LIVING; sp++){
             for (j = 0; j<4; j++){
                 i = NUM_LIVING*j+sp;
                 if (Runs[0].Slist[i][0] > EPS){
                      Scount[i]=1;
                      fprintf(sptr,"%s_%s,",path_species[sp],FitNames[j]);
                 }
                 else {
                      Scount[i]=0;
                 }
            }
        }
        fprintf (sptr,"\n");
        
        for (sp=1; sp<=NUM_LIVING; sp++){            
            sprintf(outFile,"%s_sp_%d_scramble.csv",OUTFOLDER,sp);
		    		dptr[sp]=fopen(outFile,w);
            fprintf(dptr[sp],"r,t,SSQ,");
            for (j=0; j<NUM_VARS; j++){
                fprintf(dptr[sp],"%s_%s",path_species[sp],VarNames[j]);
            }
            fprintf (dptr[sp],"\n");
        }
        
        for (r=0; r<PRUNS; r++){
            for (t=0; t<=SYSTEMS_TO_RUN; t++){         
                 for (sp=1; sp<=NUM_LIVING; sp++){
                      fprintf(dptr[sp],"%d,%d,%g,%d,",r,t,Runs[r].SSQlist[t],sp);
                      for (j = 0; j<NUM_VARS; j++){
                           i = NUM_LIVING*j+sp;
                           fprintf(dptr[sp],"%g,",Runs[r].Xlist[i][t]);
                           if ((j<4) && (Scount[i]==1)){
                              fprintf(sptr,"%g,",Runs[r].Slist[i][t]);
                           }
                      }
                      fprintf(dptr[sp],"\n");
                }
                fprintf(sptr,"\n");
              }
         }
         for (sp=1; sp<=NUM_LIVING; sp++){ 
             fclose(dptr[sp]);
        }               
        fclose(sptr);
         //oldvals[i] = rpar_ActiveRespFrac[i];

                 
           free_dvector(Scount,1,NDIM);         
           for (r=0; r<PRUNS; r++){
               free_dvector(Runs[r].x,1,NDIM);
               free_dvector(Runs[r].SSQlist,0,SYSTEMS_TO_RUN);
               free_dmatrix(Runs[r].Xlist,1,NDIM,0,SYSTEMS_TO_RUN);
               free_dmatrix(Runs[r].Slist,1,NDIM,0,SYSTEMS_TO_RUN);                
               free_SimRun(Runs[r].v);
           } 
           pthread_attr_destroy(&attr);               

}
//#undef RAND_RUNS
// -----------------------------------------------------------------------------

int solve_dfp(struct SimRun *v, double *p, int its)
{
  int iter, i, sp, c;
	//float *p;
  double fret,*dp;
  char outFile[80]; 
  FILE *fptr;
  time_t start1, end1;
  
  //start1=time(NULL);

  // Allocating vectors:  
     //p=vector(1,NDIM);           // p is vector of parameters to fit (called x in later subroutines)
   	 dp=dvector(1,NDIM);          // dp contains partial derivatives dSSQ/dp for each p	
  	 dSdX=dmatrix(1,NDIM,1,fitN+NUM_LIVING+NUM_DEAD); // full partial derivative matrix for sensitivity analyses
	   
	// Initialize P vector to 0 (can use other starting points if desired) 
     //for (i=1; i<=NDIM; i++){p[i]=0.0;}
     
  // Set starting biomass to first biomass in time series
     //for (sp=1; sp<=NUM_LIVING; sp++){
     //     i=NUM_LIVING*1 + sp;
     //     p[i] = fit_bioStart[sp];
     //     cout << i << " " << path_species[sp] << " " << p[i] << endl;
     //}
  
  // Do a first run after starting vector is applied
     start1=time(NULL);
        fret = SSQ_run(v,p);
        printf ("Start  SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
        sprintf(outFile,"_iter0");
        output_run(v, outFile);
        outstep(0,p,v); /*bin_outstep(0,p,v);*/
     end1=time(NULL);
     cout << "One run in: " << (float)(end1-start1)/60. << " Min." << endl;
     cout << "Est. step time: " << ((float)(end1-start1)/60.) * NDIM / PRUNS << " Min." << endl;
     
     start1=time(NULL);
  // Now call the general minimization routine from numerical recepies
  // SSQ_run and SSQ_deriv are the functions that calculate SSQ and dSSQ/dP
     ITMAX = its;
     //dfpmin(v,p,NDIM,GTOL,&iter,&fret, SSQ_run, SSQ_deriv);
     frprmn(v,p,NDIM,GTOL,&iter,&fret, SSQ_run, SSQ_deriv);
  
	   printf("Iterations: %3d\n",iter);
	   printf("Final  SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
  // Now call the SSQ and derivative routine a final time to get diagnostics
  // (e.g. steepness of SSQ surface at final solution)
     
		 //fret = SSQ_run(v,p);
     //printf("Intermediate test SSQ %14.6g, \n",fret);
     //sprintf(outFile,"_iter999");
     //output_run(v, outFile);
     //outstep(999,p,v); bin_outstep(999,p,v);
     //SSQ_deriv(v,p,dp);
     
     fret = SSQ_run(v,p);
     printf("Last   SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
     sprintf(outFile,"_iter1000");
     output_run(v, outFile);
     outstep(1000,p,v); /*bin_outstep(1000,p,v);*/

  // Cleanup
     free_dmatrix(dSdX,1,NDIM,1,fitN+NUM_LIVING+NUM_DEAD);
     //free_vector(p,1,NDIM);
     free_dvector(dp,1,NDIM);
     end1=time(NULL);
     cout << "solved in: " << (float)(end1-start1)/60. << " Min." << endl; 
  
	return 0;
}

//#undef NDIM 
#undef GTOL

// -----------------------------------------------------------------------------

double SSQ(struct SimRun *v){
    int s, t, measure_month, pred;
    double ssq, diff, obs, est, sd, invvar; // ssqpart;
    double fitpart;
    double forcepart;
    //double pathpart, pathtot;
    double dietpart;
    double dietNsamples;
    //double diet_other_OBS[NUM_GROUPS], diet_other_EST[NUM_GROUPS];
		//unsigned int *fp;
    //int Jlist[NUM_LIVING+1];
        
    ssq=0.0;
    
    for (s=1; s<=NUM_LIVING+NUM_DEAD; s++){
            v->fit_SSQ[0][s]=0.0;
            v->fit_SSQ[1][s]=0.0;
            v->fit_SSQ[2][s]=0.0;
            v->fit_SSQ[3][s]=0.0;
    }
    
    //memset(Jlist, 0,(NUM_LIVING+1)*sizeof(int));  
    dietNsamples = 100;
    v->fit_diettot=0.0;
    
    for (s=1; s<=dietN; s++){
        dietpart=0.0;
        v->fit_diet_SSQ[s] = 0.0;
        for (t=0; t<=fit_Years; t++){
				    est = v->fit_diet_EST[s][t];
				    obs =    fit_diet_OBS[s][t];
				    //sd  =    1.0;
				    if (((est>0.0000001) && (est<0.9999999)) && ((obs>0.0000001) && (obs<0.9999999))){
						   diff        = -dietNsamples * obs * log(est);
						   dietpart    = diff;
						   //cout << s << " " << est << " " << obs << " " << -dietNsamples * obs * log(est) << endl; 
               v->fit_diet_SSQ[s] += dietpart;
               v->fit_diettot += dietpart;
						}
        
        }
         //cout << "TOT " << v->fit_diet_SSQ[s] << endl;
    }
    //exit(0);
    
    v->fit_fittot=0;
    for (s=1; s<=fitN; s++){
        fitpart=0.0;
		    //v->fit_SSQ[s]=0.0;
		    for (t=0; t<=fit_Years; t++){
				    est = v->fit_EST[s][t];
				    obs = fit_OBS[s][t];
            sd  = fit_SD[s][t];
				    if ((obs>0) && (sd>0)){
               sd =  sqrt(log(1.0+sd*sd/(obs*obs)));
            }
				    //forced = v->fit_bioanom[fit_ApplyTo[s]][t];
				    if ((est>0) && (obs>0)){
						   diff=(log(obs)-log(est))/sd;
						   //fittot     += diff*diff;
						   fitpart    += log(sd) + ONEHALFLOGTWOPI + 0.5 * diff*diff;
						   //v->fit_SSQ[s] += diff*diff;
						}
		
				}
			
				v->fit_SSQraw[s]=fitpart;
				//TRUNC(v->fit_SSQraw[s]);
        // DO WE TAKE PATH FITS OUT OF JUVENILES, ONES WITH TIMESERIES??? KYA 7/26/07 
        switch (fit_type[s]){
               case -100:
               case 0:
                      v->fit_SSQwt[s] = v->fit_ssq_control[0][fit_ApplyTo[s]];
                      fitpart *= v->fit_SSQwt[s];
                      v->fit_SSQ[0][fit_ApplyTo[s]] += fitpart; 
                      //Jlist[fit_ApplyTo[s]]=1;
               break;
               case 50:
                      v->fit_SSQwt[s] = v->fit_ssq_control[2][fit_ApplyTo[s]];
                      fitpart *= v->fit_SSQwt[s];
                      v->fit_SSQ[2][fit_ApplyTo[s]] += fitpart; 
               break;               
               case 6:
                      v->fit_SSQwt[s] = v->fit_ssq_control[1][fit_ApplyTo[s]];
                      fitpart *= v->fit_SSQwt[s];
                      v->fit_SSQ[1][fit_ApplyTo[s]] += fitpart;
               break;
               default:
               break;
        }
        
        v->fit_SSQfit[s]=fitpart;
        //TRUNC(v->fit_SSQfit[s]);
				//fit type 0 is fit to biomass, first ssq component in output
//				if (fit_type[s]==0){fitpart *= v->fit_ssq_control[0][fit_ApplyTo[s]];
//                            v->fit_SSQ[0][fit_ApplyTo[s]] += fitpart; 
//                            Jlist[fit_ApplyTo[s]]=1;
//                            }
//        //fit type 6 is fit to total catch, second ssq component in output
//				if (fit_type[s]==6){fitpart *= v->fit_ssq_control[1][fit_ApplyTo[s]];
//                            v->fit_SSQ[1][fit_ApplyTo[s]] = fitpart;
//                           }

				v->fit_fittot += fitpart;
		}
    
    v->fit_forcetot=0.0;
    for (s=1; s<=NUM_LIVING+NUM_DEAD; s++){
         forcepart=0.0;
         v->fit_FORCEraw[s]=0.0;
      	 for (t=0; t<=fit_Years; t++){
      	     v->fit_FORCEraw[s] += v->fit_bioanom[s][t] * v->fit_bioanom[s][t];
             //forcepart += v->fit_bioanom[s][t] * v->fit_bioanom[s][t];
             //forcetot  += v->fit_bioanom[s][t] * v->fit_bioanom[s][t];
         }
         //if (v->fit_bioanom[s][PATH_YEAR_START-fit_StartYear] >0.0)
         //   {Jlist[s]=2;}
         
         //SSQ component for the forced anomaly comes out fourth in output
         forcepart = v->fit_ssq_control[3][s] * v->fit_FORCEraw[s];
         v->fit_SSQ[3][s] = forcepart;
         //TRUNC(v->fit_SSQ[3][s]);
         v->fit_forcetot += forcepart;
    }
    
       
//    for (s=1; s<=juv_N; s++){
//        Jlist[juv_JuvNum[s]]=s;     
//    }

//    pathtot=0.0;

/*  
  // This path fitting (removed on 7/30/07) gave a high weight to path fits
  // (multiplied by number of years of fit, read into invvar).  This is now
  // moved into fit weights.  Reason for number of years is, arbitrary upweighting
  // (i.e., I dunno).
    invvar= (double)fit_Years;
    t=PATH_YEAR_START-fit_StartYear;
    for (s=1; s<=NUM_LIVING; s++){
         pathpart=0.0;
         //cout << path_species[s] << "," << t << "," << v->out_BB[s][t] << "," << path_BB[s];
         if (!Jlist[s]){
         //for (t=PATH_YEAR_START-fit_StartYear; t<=PATH_YEAR_END-fit_StartYear; t++){
              est = v->out_BB[s][t];
              obs =    path_BB[s];
              //cout << path_species[s] << "," << log(est) << "," << log(obs) << endl;
				  if ((est>0) && (obs>0)){
						   diff=(log(obs)-log(est));
						   //pathtot     += diff*diff*invvar;
						   pathpart    += diff*diff*invvar;
						   //cout << est << "," << obs << "," << invvar << "," << diff << "," << diff*diff << "," << pathtot << "," << pathpart;
					}
					}
					//SSQ component for fit to the path base comes out third in output
					pathpart *= v->fit_ssq_control[2][s];
					//v->fit_SSQ[2][s] = pathpart;
					pathtot += pathpart;
					//cout << "," << v->fit_SSQ[2][s] << endl;
    }
*/
    // fifth SSQ output is currently not used, 18th is forthcoming from Kerim
        
    //ssq = 1.0*forcetot + 1.0*pathtot + 1.0*fittot;
    //cout << "PATH:  "  << pathtot << endl;
    //cout << "FIT:   "  << fittot << endl;
    //cout << "FORCE: "  << forcetot << endl;
    // with PathPart moved
      ssq = 1.0*v->fit_forcetot + 1.0*v->fit_fittot + 1.0*v->fit_diettot;
      //TRUNC(ssq); 
    return ssq;
}
// -----------------------------------------------------------------------------

// The FOLLOWING MINIMIZATION ROUTINES are taken from Numerical Recepies in C
// without comment.  They are cut and paste except for a couple of print statements

#define FREEALL fclose(fptr); fclose(tptr); free_dvector(xi,1,n);free_dvector(pnew,1,n); \
free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n); \
free_dvector(dg,1,n);

void dfpmin(struct SimRun *v, double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(struct SimRun*, double []), void (*dfunc)(struct SimRun*, double [], double []))
{
	void lnsrch(struct SimRun *v, int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(struct SimRun*, double []));
	int check,i,its,j,r;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  // KYA Added
	double sss;  char outFile[80]; char ilabel[80];
  FILE *fptr;
  
  if (ITMAX == 0){return;}

     BigIt=0;

      sprintf(outFile,"%s_ForcePartials.csv",OUTFOLDER);
      tptr=fopen(outFile,"w");
      
      sprintf(outFile,"%s_fitDerivs.csv",OUTFOLDER);
      fptr=fopen(outFile,"w");
      fprintf(fptr,"Step,Species,var,x,Control,dSdx,");
      for (i=1; i<=fitN; i++){
          fprintf(fptr,"Fit_%s_%d,",path_species[fit_ApplyTo[i]],fit_type[i]);
      }
      for (i=1; i<=NUM_LIVING+NUM_DEAD; i++){
          if(FORCE_BIO_FLAG[i]){
             fprintf(fptr,"Force%s,",path_species[i]);
          }
      }      
      fprintf(fptr,"EndTot,\n");
      
	dg=dvector(1,n);
	g=dvector(1,n);
	hdg=dvector(1,n);
	hessin=dmatrix(1,n,1,n);
	pnew=dvector(1,n);
	xi=dvector(1,n);
	printf ("first\n");
	fp=(*func)(v,p);
	printf ("Second\n");
	(*dfunc)(v,p,g);
  	printf ("Third\n");
	  sprintf(ilabel,"%d",0);
	  outderiv(v,fptr, p, g, ilabel);
	//outstep(10000+its,p,v,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		lnsrch(v,n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		// KYA Added
    	 sss = SSQ_run(v,p);
    	 cout << "Iter: " << its << "SSQ: " << sss << endl;
	     sprintf(outFile,"RunFit_iter%d",its);
       output_run(v,outFile);
       outstep(its,p,v); /*bin_outstep(its,p,v);*/
       BigIt=its;
		   //if (its ==2){TESTFLAG=1;}else{TESTFLAG=0;}
		//exit(0);

		if (test < TOLX) {
    	 sprintf(ilabel,"TOLX exit");
	     outderiv(v,fptr, p, g, ilabel);
			FREEALL
			return;
		}

		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(v,p,g);
		   sprintf(ilabel,"%d",its);
	     outderiv(v,fptr, p, g, ilabel);
		
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
    	 sprintf(ilabel,"gtol exit");
	     outderiv(v,fptr, p, g, ilabel);
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=i;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin[j][i]=hessin[i][j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("too many iterations in dfpmin");
    	 sprintf(ilabel,"too many exit");
	     outderiv(v,fptr, p, g, ilabel);
	
	FREEALL
}
//#undef ITMAX
//#undef EPS
//#undef TOLX
#undef STPMX
#undef FREEALL
#undef NRANSI


//#define TOLX 1.0e-7

// -----------------------------------------------------------------------------   
void lnsrch(struct SimRun *v, int n, double xold[], double fold, double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(struct SimRun*, double []))
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;
   
	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;

for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(v,x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);

if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam > 0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		alam=FMAX(tmplam,0.1*alam);
	}
}
#undef ALF
#undef TOLX

//---------------------------------------------------------------

   
int ncom;
double *pcom,*xicom,(*nrfunc)(struct SimRun*, double []);

double f1dim(struct SimRun *v, double x)
{
	int j;
	double f,*xt;
   
	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(v,xt);
	free_dvector(xt,1,ncom);
	return f;
}
#undef NRANSI
//---------------------------------------------------------------

#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
      
double brent(struct SimRun *vv, double ax, double bx, double cx, 
             double tol,double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;
   
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f1dim(vv,x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {

r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));

fu=f1dim(vv,u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI

//---------------------------------------------------------------
#define NRANSI
#define GOLD 1.618034
#define GLIMIT 6.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
   
void mnbrak(struct SimRun *v, double *ax, double *bx, double *cx, double *fa, 
            double *fb, double *fc)
{
	double ulim,u,r,q,fu,dum;
   
	*fa=f1dim(v,*ax);
	*fb=f1dim(v,*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=f1dim(v,*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));

ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=f1dim(v,u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=f1dim(v,u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=f1dim(v,u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,f1dim(v,u))

}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=f1dim(v,u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=f1dim(v,u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI

//---------------------------------------------------------------
/* note #undef's at end of file */
#define NRANSI
#define TOL 2.0e-4
   
void linmin(struct SimRun *v, double p[], double xi[], int n, 
            double *fret, double (*func)(struct SimRun*, double []))
{

	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
   
	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(v,&ax,&xx,&bx,&fa,&fx,&fb);
	*fret=brent(v,ax,xx,bx,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];

}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}
#undef TOL
#undef NRANSI

//---------------------------------------------------------------------
/* note #undef's at end of file */
#define NRANSI
      
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);
   
void frprmn(struct SimRun *v, double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(struct SimRun*, double []), void (*dfunc)(struct SimRun*, double [], double []))
{

	int i,j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

  // KYA Added
	double sss;  char outFile[80]; char ilabel[80];
  FILE *fptr;
  
  if (ITMAX == 0){return;}

     BigIt=0;

      sprintf(outFile,"%s_ForcePartials.csv",OUTFOLDER);
      tptr=fopen(outFile,"w");
      
      sprintf(outFile,"%s_fitDerivs.csv",OUTFOLDER);
      fptr=fopen(outFile,"w");
      fprintf(fptr,"Step,Species,var,x,Control,dSdx,");
      for (i=1; i<=fitN; i++){
          fprintf(fptr,"Fit_%s_%d,",path_species[fit_ApplyTo[i]],fit_type[i]);
      }
      for (i=1; i<=NUM_LIVING+NUM_DEAD; i++){
          if(FORCE_BIO_FLAG[i]){
             fprintf(fptr,"Force%s,",path_species[i]);
          }
      }      
      fprintf(fptr,"EndTot,\n");
      
  its = 0;
    	 sss = SSQ_run(v,p);
    	 //cout << "Iter: " << its << "SSQ: " << sss << endl;
    	 printf("Iter%2d SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         its,sss,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
    	 
	     sprintf(outFile,"RunFit_iter%d",its);
       output_run(v,outFile);
       outstep(its,p,v); /*bin_outstep(its,p,v);*/
       BigIt=its;
			 	   
	g=dvector(1,n);
	h=dvector(1,n);
	xi=dvector(1,n);
	fp=(*func)(v,p);
	(*dfunc)(v,p,xi);
		  sprintf(ilabel,"%d",0);
	    outderiv(v,fptr, p, g, ilabel);
	
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(v,p,xi,n,fret,func);

		// KYA Added
    	 sss = SSQ_run(v,p);
    	 //cout << "Iter: " << its << "SSQ: " << sss << endl;
    	 printf("Iter%2d SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         its,sss,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
	     sprintf(outFile,"RunFit_iter%d",its);
       output_run(v,outFile);
       outstep(its,p,v); /*bin_outstep(its,p,v);*/
       BigIt=its;
		
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
		   sprintf(ilabel,"NORMAL Exit");
	     outderiv(v,fptr, p, g, ilabel);
FREEALL
			return;
		}
		fp= *fret;
		(*dfunc)(v,p,xi);
		   sprintf(ilabel,"%d",its);
	     outderiv(v,fptr, p, g, ilabel);

		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
		   sprintf(ilabel,"Impossible!");
	     outderiv(v,fptr, p, g, ilabel);		
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
    	 sprintf(ilabel,"too many exit");
	     outderiv(v,fptr, p, g, ilabel);
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI

