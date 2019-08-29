
// Load in all variables, basic packages, etc.
#include "ecosim_threads.h"   // Main headers (standard C packages) AND global variables

// Local includes that require global variables contained in ecosim.h
#include "simio.cpp"
#include "fitfunctions.cpp"
#include "fftest.cpp"
#include "simlabs.cpp" 
#include "sense.cpp"
//#include "charts.cpp"

#ifdef MPI_COMPILE
#include "MPI_sim.cpp"
#else
int MPI_test(int argc, char *argv[], char submode){
  cout << "Attempting to run mpi routines on non-mpi compile.  Exiting."<< endl;
  exit(1);
}
#endif

// MAIN CODE BODY

// ------------------------------------------------------------------------

int main(int argc, char **argv)
   {
   int error;

   // Check memory of structures for ballparking size
      cout << sizeof(struct SimRun) << endl;
      cout << sizeof(struct RatePar) << endl;
   
   // Read command line options
      error=commParse(argc, argv);
      if(error){exit(1);}

   // Variables needed to set for basic run
      //fit_Years=MAX_YEARS; // (Years to run if there's no fitting file)

   // Read in forcing functions (climate, fishing)
   // Note:  read_juveniles has to be first so the other files
   //        can apply climate and fitting to juvenile groups.
	 // Note:  These are all global variables and should be treated 
	 //        as const if at all possible, make local copies if needed
	 //        to change.      
      read_juveniles();
      read_climate();
      new_read_fitting();
      read_diets();

      read_guilds();      

   // Some Basic calculations that may be useful
      TL_From_FlowMat();


   // Now we have path data and forcing files loaded.  What next??
   // Probably one of the following run modes:
   
      cout << "Run Mode " << RUNMODE ;
      switch (RUNMODE){
      // MPI
         case 'M':
         MPI_test(argc,argv,SUBMODE);
         break;
      // Run fitting
         case 'f':
           cout << " is fitting." << endl;
           simple_fitted();
           break;
      // Run once with full/complete outputs
         case 'o':
           cout << " does one run." << endl;
           one_run();
           break;  
      // Generate random noise
         case 'r':
           cout << " is Noise Generation." << endl;     
           Noisy_Run();
           break;
      // Scramble from fitted
         case 's':
           cout << " is Scramble from fitted." << endl;     
           scramble_fitted();
           break;    
      // Generate a derivative field
         case 'd':
           cout << " is Derivative field from fitted." << endl;     
           point_derivative();
           break;
      // Run a long cycle (more than 200 years)  
         case 'l':
           cout << " is a long run." << endl;                
           long_run();
           break;              
      // This generates random ecosystems and tests/rejects them
         case 'e':
            cout << " is a random sense set." << endl;
            // KYA June 3 2008, random-ecosystem-series is depreciated.
            // You should call fromfitted and set up an input vector file
            // (otherwise all vectors will default to 0, e.g. handling time on)
            random_fromfitted_series();
            //random_ecosystem_series();
            break;
      // This loads already generated ecosystems and outputs them
         case 'i':
            cout << " loads a random sense set." << endl;
            loaded_ecosystem_series();
            break;
            
         case 'a':
           cout << " assessment series stats." << endl;
           loaded_ecosystem_series_assess();
           break;
           
      // This loads already generated ecosystems, perturbs and outputs results
         case 'p':
            cout << " perturbs a random sense set." << endl;
            perturb_loaded_series();
            break;
      // Kerim's Labs
          case 'K':
             cout << " does what Kerim wants it to." << endl;
             sim_labs();
             break;
         
      case '\0':
         cout << "No Run Mode Specified. " << endl;
         break;
      default:
         cout << "Run Mode " << RUNMODE << " unrecognized." << endl;
         break;
   }  //end RUNMODE switch
   
   // Generic Thread Cleanup command (in case of lingering threads)
      pthread_exit(NULL);
   return 0;
}



// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
int one_run(void)
{
struct SimRun *v;
double ssq;

int i, link, sp;
//float *fit_vector;
char outFile[80];					   
FILE *dmpptr;
const char* w = "w";
int rlookup[NUM_GROUPS+1][NUM_GROUPS+1];
          // Allocate memory for the SimRun
             //fit_vector=vector(0,NDIM);
             //v = (struct SimRun *)calloc(1,sizeof(struct SimRun)); 
                        	 //cout << "jo1" << endl;
             v=new_SimRun();  
  
             //for (i=0; i<=NDIM; i++){v->fit_vector[i]=0.0;}
					// Sets up baseline run
					           	 //cout << "jo2" << endl;
             base_system(v);  
					           	 //cout << "jo3" << endl;

      		// Run Adams-Basforth
      		   FLOWOUT = 0; // Setting flowout to 1 dumps flow files for each year during the base run
      		   //printf("hey\n");
					   Adams_Basforth(v, 0, fit_Years);
      		   //printf("hey\n");
             FLOWOUT = 0;
					           	 //cout << "jo4" << endl;
          // If it didn't crash in first run, calculate SSQ and save the system
           	 ssq = 0.0;
           	 //cout << "jo5" << endl;
             if (!v->DISCARD_YEAR){
                 ssq = SSQ(v);
                 sprintf(outFile,"_base");
                 output_run(v, outFile);
                 outstep(99,v->fit_vector,v);
             }
	           printf("PreLoadSSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         ssq,v->fit_fittot, v->fit_forcetot, v->fit_diettot);

         //  Load starting fit vector (or set to 0 if none existing)    
             load_instep(v);  
         //  This solves for best fit and returns best fitting rates in v
             //v->fit_vector[NUM_LIVING*1 + 2] = -5.0;  //S perm whales
             //v->fit_vector[NUM_LIVING*1 + 6] = -0.6;  // Jumpback 
             //v->fit_vector[NUM_LIVING*1 + 7] = -2.8;  // Fin
             //v->fit_vector[NUM_LIVING*1 + 8] = -4.6;  // Sei
             //v->fit_vector[NUM_LIVING*1 + 9] = -5.0;  // Right
             //solve_dfp(v,v->fit_vector);
                          
        //   Now do a final run with new rates
             // initialize_stanzas(v);
             // Adams_Basforth(v, 0, fit_Years);
             //ssq = 0.0;
             FLOWOUT=0; 
             //FLOWOUT=1;
             AGEOUT = 1;            
             ssq = SSQ_run(v, v->fit_vector);
             AGEOUT = 0;
             FLOWOUT=0;
             //deriv_field(v, v->fit_vector);
	           printf("Final  SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         ssq,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
             if (!v->DISCARD_YEAR){
                 //ssq = SSQ(v);
                 sprintf(outFile,"_single");
                 output_run(v, outFile);
                 //chart_run(v);

                 outstep(100,v->fit_vector,v);
             }
            


           link_dump(v);  
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
                
          //free_vector(fit_vector,1,NDIM);
          free_SimRun(v);
						 
return 0;
}
// ------------------------------------------------------------------------------

int long_run(void)
{
struct SimRun *v;
float ssq;
int i, link, sp, t;
//float *fit_vector;
char outFile[80];				   
FILE *dmpptr;
FILE *fptr;  
const char* w = "w";
int rlookup[NUM_GROUPS+1][NUM_GROUPS+1];

          // Allocate memory for a Sim Run
             v=new_SimRun();  
          // Load starting fit vector (or set to 0 if none existing)    
             load_instep(v);  
          // First set default parameters for sim
             path_to_rates(v);
          // Then apply x vector to the default non-juvenile parameters.  
             apply_vector_to_rates(v, v->fit_vector);
          // Then initialize stanzas
             initialize_stanzas(v);         
          // initialize the random seed
             init_by_array(rseed, &v->rng);
              
          // Then run the model and caluclate SSQ
             sprintf(outFile,"%s_long_bio.csv",OUTFOLDER);
             fptr=fopen(outFile,w);
             for (REPEAT = 0; REPEAT<=LONG_CYCLES; REPEAT++){
             
                 for (sp=0; sp<=NUM_GROUPS; sp++){
		                 for (t=0; t<=fit_Years*STEPS_PER_YEAR; t++){
                     				  // force_byprey[sp][t]=1.0;
				                      v->force_bymort[sp][t]=1.0 + 0.1 * gaussian(&v->rng);
				                      // force_byrecs[sp][t]=1.0;
				             }
		             } 
                 //output_run(v, outFile);
                 cout << "Long Step " << REPEAT << endl;
                 Adams_Basforth(v,0, fit_Years);
                    
            		// Save biomass to a file
            		   if (REPEAT == 0){
                     fprintf(fptr,"Year,");
                     for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                        if(v->discarded[sp]){fprintf(fptr,"dd%s,",path_species[sp]);}
                        else             {fprintf(fptr,"%s,",path_species[sp]);}
						         }   fprintf(fptr,"\nGuild,");
                     for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                          fprintf(fptr,"%s,",guildlist[path_species[sp]].c_str());
						         }   fprintf(fptr,"\n");						         
						         
                   }
                   for (t=0; t<fit_Years; t++){
                      fprintf(fptr,"%d,",t + fit_StartYear + REPEAT*fit_Years);
                         for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                             fprintf(fptr,"%g,",v->out_BB[sp][t]);}   
                         fprintf(fptr,"\n");
						       }                 
             }            
                 
						 fclose(fptr);

          // De-allocate the memory
             free_SimRun(v);
						 
return 0;
}
// ------------------------------------------------------------------------------


void link_dump(struct SimRun *v)
{
int links,pred,prey;
char outFile[80];					   
FILE *dmpptr;
const char* w = "w";

sprintf(outFile,"%s_links.csv",OUTFOLDER);
dmpptr=fopen(outFile,w);

    fprintf(dmpptr,"links,prey,pred,QQ,HandleSwitch,DD,VV\n");
    for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
		    prey = v->RRR.rpar_PreyFrom[links];
		    pred = v->RRR.rpar_PreyTo[links];

        fprintf(dmpptr,"%d,%d,%d,%g,%g,%g,%g,\n",
           links,
           prey,
           pred,
           v->RRR.rpar_QQ[links], 
           v->RRR.rpar_HandleSwitch[links],
           v->RRR.rpar_DD[links],
           v->RRR.rpar_VV[links]);
//        Q =   v->RRR.rpar_QQ[links] * pred_YY[pred] * pow(prey_YY[prey],v->RRR.rpar_HandleSwitch[links]) *
//				    ( v->RRR.rpar_DD[links] / ( v->RRR.rpar_DD[links] - 1.0 + 
//						                     pow(v->RRR.rpar_HandleSelf[pred] * prey_YY[prey]   + 
//																 (1. - v->RRR.rpar_HandleSelf[pred]) * HandleSuite[pred],
//                                                     v->RRR.rpar_HandleSwitch[links])) )*
//            ( v->RRR.rpar_VV[links] / ( v->RRR.rpar_VV[links] - 1.0 + 
//                                 v->RRR.rpar_ScrambleSelf[pred] * pred_YY[pred] + 
//						                     (1. - v->RRR.rpar_ScrambleSelf[pred]) * PredSuite[prey]) );
    }
    fclose(dmpptr);
}
//  -----------------------------------------------------------------------------
int point_derivative(void)
{
struct SimRun *v;
float ssq;
int i;
//float *fit_vector;
char outFile[80];	

          // Allocate memory for the SimRun
             v=new_SimRun();  
					// Sets up baseline run
             base_system(v);  
      		// Run Adams-Basforth
             load_instep(v);  
             FLOWOUT=1;
             ssq = SSQ_run(v, v->fit_vector);
             FLOWOUT=0;
          // If it didn't crash in first run, calculate SSQ and save the system
             cout << "loaded run with SSQ: " << ssq << endl;

          // Do the derivative
             deriv_field(v, v->fit_vector);             

             free_SimRun(v);
						 
return 0;
}

// ------------------------------------------------------------------------------
int scramble_fitted(void)
{
struct SimRun *v;
double ssq;
int i;
//float *fit_vector;
char outFile[80];					   

          // Allocate memory for the SimRun
             v=new_SimRun();  
					// Sets up baseline run
             base_system(v);  
      		// Run Adams-Basforth
      		   //FLOWOUT = 0; // Setting flowout to 1 dumps flow files for each year during the base run
					   //Adams_Basforth(v, 0, fit_Years);
             //FLOWOUT = 0;
          // If it didn't crash in first run, calculate SSQ and save the system
           	 //ssq = 0.0;
             //cout << "base run with SSQ: " << ssq << endl;
         //  Load starting fit vector (or set to 0 if none existing)  
             load_instep(v);  
             ssq = SSQ_run(v, v->fit_vector);
          // If it didn't crash in first run, calculate SSQ and save the system
             cout << "loaded run with SSQ: " << ssq << endl;

// PUT THIS NEXT LINE BACK IN TO TURN ON FITTING -- JUST OUT FOR DERIVS            
           //  deriv_field(v, v->fit_vector);

             SSQ_scatter(v,v->fit_vector);
                          
          free_SimRun(v);
						 
return 0;
}

//-------------------------------------------------------------------------------
int simple_fitted(void)
{
struct SimRun *v;
double ssq;
int i;
//float *fit_vector;
char outFile[80];					   

          // Allocate memory for the SimRun
             //fit_vector=vector(0,NDIM);
             //v = (struct SimRun *)calloc(1,sizeof(struct SimRun)); 
             v=new_SimRun();  

             //for (i=0; i<=NDIM; i++){v->fit_vector[i]=0.0;}
					// Sets up baseline run
             base_system(v);  
             //deriv_field(v, v->fit_vector);
      		// Run Adams-Basforth
      		   FLOWOUT = 0; // Setting flowout to 1 dumps flow files for each year during the base run
					   Adams_Basforth(v, 0, fit_Years);
             FLOWOUT = 0;
          // If it didn't crash in first run, calculate SSQ and save the system
           	 //ssq = 0.0;
             //if (!v->DISCARD_YEAR){
             //    ssq = SSQ(v);
             //    sprintf(outFile,"_base");
             //    output_run(v, outFile);
             //}
             //cout << "base run with SSQ: " << ssq << endl;
         //  Load starting fit vector (or set to 0 if none existing) 
				     //bin_instep(v); 
             load_instep(v);  
             ssq = SSQ_run(v, v->fit_vector);
          // If it didn't crash in first run, calculate SSQ and save the system
           	 //ssq = 0.0;
             //if (!v->DISCARD_YEAR){
             //    ssq = SSQ(v);
             //    sprintf(outFile,"_prep");
             //    output_run(v, outFile);
             //}
             printf ("Load   SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		                 ssq,v->fit_fittot, v->fit_forcetot, v->fit_diettot);


// PUT THIS NEXT LINE BACK IN TO TURN ON FITTING -- JUST OUT FOR DERIVS            
           //  deriv_field(v, v->fit_vector);

             solve_dfp(v,v->fit_vector,200);
                          
        //   Now do a final run with new rates
             // initialize_stanzas(v);
             // Adams_Basforth(v, 0, fit_Years);
             //ssq = 0.0;             
             ssq = SSQ_run(v, v->fit_vector);
             //deriv_field(v, v->fit_vector);
             if (!v->DISCARD_YEAR){
                 //ssq = SSQ(v);
                 sprintf(outFile,"_fit");
                 output_run(v, outFile);
             }
             cout << "final run with SSQ: " << ssq << endl;

          //free_vector(fit_vector,1,NDIM);
          free_SimRun(v);
						 
return 0;
}


// -----------------------------------------------------------------------------

int deriv_field(struct SimRun *v, double x[]){
int LL,sp,i;
float ff,pp;
int jj,kk;

    const char* w = "w";
	  FILE *fptr;   
    char outFile[80];

     // First set default parameters for sim
        path_to_rates(v);

     // Then apply x vector to the default non-juvenile parameters.
  
     // This applies X[1..NumLiving] to PB
        apply_vector_to_rates(v, x);
  
      // Then initialize stanzas
         initialize_stanzas(v);

         //LL=(NUM_GROUPS+1)*sizeof(float);
         //memset(v->deriv_dyt, 0,LL);
         //memset(v->deriv_DerivT, 0,LL);

   // Set Start Biomass and other state variables
	    for (sp=0; sp<=NUM_GROUPS; sp++){
			     v->state_BB[sp]    = v->RRR.rpar_B_BaseRef[sp];
			     //v->fit_bioanom[sp][0] = 0.0;
			     v->state_Ftime[sp] = 1.0;
			     //v->discarded[sp]   = 0;
			     //v->month_noise[sp] = 0.0;
			}
	 // Set State Variables for Stanzas
		 // for (i=1; i<=juv_N; i++){
		 //    if (v->stanzaBasePred[juv_JuvNum[i]]>0){
		 //    pred_YY[juv_JuvNum[i]] = v->state_Ftime[juv_JuvNum[i]] * 
		 //		        v->stanzaPred[juv_JuvNum[i]]/v->stanzaBasePred[juv_JuvNum[i]];
		 //    pred_YY[juv_AduNum[i]] = v->state_Ftime[juv_AduNum[i]] * 
		 //		        v->stanzaPred[juv_AduNum[i]]/v->stanzaBasePred[juv_AduNum[i]];
		 //    }
		 // }
      sprintf(outFile,"%s_deriv.csv",OUTFOLDER);
			fptr=fopen(outFile,w);
			fprintf(fptr,"Species,Bmult,Gain,Loss,\n");
      for (sp=0; sp<=NUM_LIVING; sp++){   
               for (i=1; i<=juv_N; i++){
                   {v->stanzaPred[juv_JuvNum[i]]=v->stanzaBasePred[juv_JuvNum[i]];}
                   {v->stanzaPred[juv_AduNum[i]]=v->stanzaBasePred[juv_AduNum[i]];}
               }   
          cout << sp << endl;
        for  (jj=-5; jj<=1; jj++){
          pp = pow(10.0,(float)jj);
          for (kk=2; kk<=19; kk++){  
               ff = pp * (float)kk;
               v->state_BB[sp] = v->RRR.rpar_B_BaseRef[sp] * ff;
               for (i=1; i<=juv_N; i++){
               if (sp == juv_JuvNum[i]){v->stanzaPred[sp]=ff*v->stanzaBasePred[sp];}
               if (sp == juv_AduNum[i]){v->stanzaPred[sp]=ff*v->stanzaBasePred[sp];}
               }
              deriv_master(v,0,0,0);
              fprintf(fptr,"%s,%g,%g,%g,\n",path_species[sp],ff,
                      (v->deriv_TotGain[sp] - v->deriv_LossPropToQ[sp])/v->state_BB[sp],
                       v->deriv_LossPropToB[sp]/v->state_BB[sp]
                      );
          }
         }
         v->state_BB[sp] = v->RRR.rpar_B_BaseRef[sp];
      }
      fclose(fptr);

}

// -----------------------------------------------------------------------------

int deriv_at_time(struct SimRun *v, int y, int m, int d)
{
int LL,sp,i;
float ff,pp;
int jj,kk;
float old_juv_pred[MAX_SPLIT+1];
float old_adu_pred[MAX_SPLIT+1];
float old_BB;
float old_dyt[NUM_GROUPS+1];
float old_dT[NUM_GROUPS+1];

    const char* w = "w";
	  FILE *fptr;   
    char outFile[80];

    for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
       old_dyt[i]=v->deriv_dyt[i];
       old_dT[i]=v->deriv_DerivT[i];
    }

    for (i=1; i<=juv_N; i++){
        old_juv_pred[i]=v->stanzaPred[juv_JuvNum[i]];
        old_adu_pred[i]=v->stanzaPred[juv_AduNum[i]];
    }

      sprintf(outFile,"%s_year%d_deriv.csv",OUTFOLDER,y+fit_StartYear);
			fptr=fopen(outFile,w);
			fprintf(fptr,"Species,Bmult,Gain,Loss,M2,M0,F,\n");

      for (sp=0; sp<=NUM_LIVING; sp++){  
          old_BB=v->state_BB[sp];
          for (i=1; i<=juv_N; i++){
                   {v->stanzaPred[juv_JuvNum[i]]=old_juv_pred[i];}
                   {v->stanzaPred[juv_AduNum[i]]=old_adu_pred[i];}
          }   
          cout << y+ fit_StartYear << " " << sp << endl;        
          
        for  (jj=-5; jj<=1; jj++){
          pp = pow(10.0,(float)jj);
          for (kk=2; kk<=19; kk++){  
               ff = pp * (float)kk;       
                 v->state_BB[sp] = v->RRR.rpar_B_BaseRef[sp] * ff;
              for (i=1; i<=juv_N; i++){
                   if (sp == juv_JuvNum[i]){v->stanzaPred[sp]=ff*v->stanzaBasePred[sp];}
                  if (sp == juv_AduNum[i]){v->stanzaPred[sp]=ff*v->stanzaBasePred[sp];}
               }

               deriv_master(v,y,m,d);

              fprintf(fptr,"%s,%g,%g,%g,%g,%g,%g,\n",path_species[sp],
                       ff,
                      (v->deriv_TotGain[sp] - v->deriv_LossPropToQ[sp])/v->state_BB[sp],
                       v->deriv_LossPropToB[sp]/v->state_BB[sp],
                       v->deriv_FoodLoss[sp]/v->state_BB[sp],
                       v->deriv_MzeroLoss[sp]/v->state_BB[sp], 
                       v->deriv_FishingLoss[sp]/v->state_BB[sp]
                      );
          }
         }
         v->state_BB[sp] = old_BB;
         //deriv_master(v,y,m,d);
      }
      
          for (i=1; i<=juv_N; i++){
                   {v->stanzaPred[juv_JuvNum[i]]=old_juv_pred[i];}
                   {v->stanzaPred[juv_AduNum[i]]=old_adu_pred[i];}
          }       
      //deriv_master(v,y,m,d);
      for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
          v->deriv_dyt[i]=old_dyt[i];
          v->deriv_DerivT[i]=old_dT[i];
      }

      fclose(fptr);

}

// -----------------------------------------------------------------------------
int Adams_Basforth (struct SimRun *v, int StartYear, int EndYear){
    
int y, m, d, c, j;
int sp, t, i, ageMo, s, link,prey,pred,links,gr;
int inbound;
double old_B, new_B, nn, ww, bb, pd;
double newbio; float ystep;
double bioanom, bioratio;
double out_CC[NUM_GROUPS+1];
double startyear_forcebio[NUM_GROUPS+1];
int Jlist[NUM_GROUPS+1];
int Alist[NUM_GROUPS+1];
unsigned int LL;
//const char* fname = "StanzaOut.csv";
char outFile[80];
const char* w = "w";
FILE *dptr[NUM_LIVING];
FILE *rptr[MAX_SPLIT+1];
    //fptr=fopen(fname,w);
   // Set derivative terms to 0

//unsigned int *fp;

memset(Jlist, 0,(NUM_LIVING+1)*sizeof(int));   
memset(Alist, 0,(NUM_LIVING+1)*sizeof(int)); 
      for (s=1; s<=juv_N; s++){
        Jlist[juv_JuvNum[s]]=s;     
        Alist[juv_AduNum[s]]=s; 
      } 

for (i=1; i<=juv_N; i++){
          for (t=0; t<=fit_Years; t++){
          v->out_RR[i][t]=0.0;
          v->out_EE[i][t]=0.0;
          v->out_SSB[i][t]=0.0;
          }
}

inbound=0;

if (AGEOUT){
             for (j=1; j<=juv_N; j++){
    		     // Save age structure to a file (only for one species numbered in Adams-Basforth)
     		         sprintf(outFile,"%s_juv%d_age.csv",OUTFOLDER,j);
                 rptr[j]=fopen(outFile,w);
                  fprintf(rptr[j],"Year,");
                 for (ageMo=v->firstMoJuv[j]; ageMo<v->lastMoAdu[j]; ageMo+=12){
                     fprintf(rptr[j],"NN%d,WW%d,",ageMo,ageMo);
                 }
                 fprintf(rptr[j],"NNplus,WWplus,\n");
             }
}
            
if (!REPEAT){
      LL=(NUM_GROUPS+1)*sizeof(double);
         memset(v->deriv_dyt, 0,LL);
         memset(v->deriv_DerivT, 0,LL);
                    
   // Set Start Biomass and other state variables
	    for (sp=0; sp<=NUM_GROUPS; sp++){
			     v->state_BB[sp]    = v->RRR.rpar_B_BaseRef[sp];
			     //v->fit_bioanom[sp][0] = 0.0;
           for (y=0; y<=fit_Years; y++){
             //bioMat[sp][y] = v->state_BB[sp];
             //v->out_BB[sp][y] = v->state_BB[sp];
             v->fit_bioanom[sp][y]=0.0;
           }
			     v->state_Ftime[sp] = 1.0;
			     v->discarded[sp]   = 0;
			     //v->month_noise[sp] = 0.0;
			}

        for (c=1; c<=fitN; c++){
             for (y=0; y<=fit_Years; y++){
                  v->fit_EST[c][y]=-1.0;
             }
        }

	 // Set State Variables for Stanzas
      SplitSetPred(v);
      
					// set to forced biomass
					// INITIAL SET NOT PENALIZED AS PART OF BIOANOM
  		           for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
  		               if ( (FORCED_BIO[sp][StartYear] >0.0))
  		                  {
                           //ystep    = ((float)(m)) / ((float)(STEPS_PER_YEAR));
                           newbio   = FORCED_BIO[sp][StartYear];// +  ystep * (FORCED_BIO[sp][y+1] - FORCED_BIO[sp][y]);
                           bioanom  = newbio - v->state_BB[sp];
                           bioratio = newbio/v->state_BB[sp];
                           //cout << sp << " " << y << " " << ystep << " " << FORCED_BIO[sp][y] << " " << FORCED_BIO[sp][y+1] << endl;
                           if  ( (v->RRR.rpar_NoIntegrate[sp] == sp) || (v->RRR.rpar_NoIntegrate[sp] == 0) ){
                                v->state_BB[sp] *= bioratio;
                           }
                           else {
                                 for (i = 1; i<= juv_N; i++){
                                      if (sp == juv_AduNum[i]){
                                          //for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++)
                                          for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++)
                                              {
                                              v->NageS[i][ageMo] *= bioratio;
                                              }
                                      }
                                 } 
                                
                           }
                                                                  
                        }
                }	
                SplitSetPred(v);

     deriv_master(v,0,0,0);
     deriv_master(v,0,0,0);

}
					           	 // cout << "jod1" << endl;
//GIANT_DUMP
 if (GIANT){giant_dump(v,v->RunID);}


// MAIN LOOP STARTS HERE     
      for (y=StartYear; y<EndYear; y++){

        LL=(NUM_GROUPS+1)*sizeof(double);
        memset(out_CC,0,LL); 
        
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
             //bioMat[sp][y] = v->state_BB[sp];
             //v->out_BB[sp][y] = v->state_BB[sp];
             v->fit_bioanom[sp][y]=0.0;
        }
        for (i=1; i<=juv_N; i++){
           v->out_SSB[i][y] = v->SpawnBio[i];
        }
        
        for (c=1; c<=fitN; c++){
             switch (fit_type[c]){
               case -100:
               case 50:
               case 0:
                  v->fit_EST[c][y] = v->state_BB[fit_ApplyTo[c]];
                  //TRUNC(v->fit_EST[c][y]);
               break;
               case 6:
               //   v->fit_EST[c][y] = out_CC[fit_ApplyTo[c]];
               break;
               default:
               break;
             }
        }


        float tmpB,tmpPred,ff;
        int ttt;
                
        if ((FLOWOUT)&&(y%1==0)){deriv_at_time(v, y, 0, 0);}
        
        for (m=0; m<STEPS_PER_YEAR; m++){
         
		      for (d=0; d<STEPS_PER_MONTH; d++){
		                 // printf("%d %d %d\n",y,m,d);
		          // Calculate Derivative for a given timestep
	               deriv_master(v,y,m,d);
	                           //printf("%d %d %d\n",y,m,d);
	                FILE *fptr;
	                char outfile[80];
	                
	                if ((FLOWOUT)&&(d==0)&&(m==MEASURE_MONTH)){
	                    sprintf(outfile,"%s_flow%d.csv",OUTFOLDER,y);
	                    fptr=fopen(outfile,w);
	                    fprintf (fptr,"Name,");
                      for (i=0; i<=NUM_GROUPS; i++){
                          fprintf(fptr,"%s,",path_species[i]);
                      }
                      fprintf(fptr,"\n");
                      for (sp=0; sp<=NUM_GROUPS; sp++){
                          fprintf(fptr,"%s,",path_species[sp]);
                          for (i=0; i<=NUM_GROUPS; i++){
                              fprintf(fptr,"%g,",v->deriv_ConsMat[sp][i]);
                          }
                          fprintf(fptr,"\n");
                      }
                      fclose(fptr);
                   }

	               									           	 // cout << "jod6 " << m << endl;                 

	               //cout << d << endl;
	               // Rate saving needs to be here to associate correct
	               // biomass with correct rates.
	               
	               for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
	                   v->monthly_BB[sp][y*STEPS_PER_YEAR+m] = v->state_BB[sp];
	                   v->monthly_PP[sp][y*STEPS_PER_YEAR+m] = v->deriv_TotGain[sp] - v->deriv_LossPropToQ[sp]; 
	               }
	               
                 if ((d==0)&&(m==MEASURE_MONTH)){
	                  for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
 	                       v->out_BB[sp][y] = v->state_BB[sp];
 	                          //TRUNC(v->out_BB[sp][y]);														 
					               v->out_MM[sp][y] =(v->deriv_MzeroLoss[sp]+v->deriv_FoodLoss[sp])/v->state_BB[sp];
					               v->out_M2[sp][y] =(v->deriv_FoodLoss[sp])/v->state_BB[sp];
					               v->out_QB[sp][y] = v->deriv_FoodGain[sp]/v->state_BB[sp];
					               v->out_PB[sp][y] = (v->deriv_FoodGain[sp]- v->deriv_UnAssimLoss[sp]  - v->deriv_ActiveRespLoss[sp])/v->state_BB[sp];
                         }

                    for (gr=1; gr<=NUM_GROUPS; gr++){
                        v->out_Effort[gr][y] = v->fish_Effort[gr];
                    }
                    
                    for (link=1; link<=NUM_LIVING; link++){
                         v->fit_diet_EST[link][y] = 1.0;
										}
                    for (link = NUM_LIVING+1; link<=dietN; link++){
                         pred = int(fit_diet_lookup[link]/1000); 
                         prey = fit_diet_lookup[link] - pred*1000;
                         if (prey<=NUM_LIVING+NUM_DEAD){
                         //cout << link << " " << pred << " " << prey << endl;     
                             v->fit_diet_EST[link][y]    = v->deriv_ConsMat[prey][pred]/v->deriv_FoodGain[pred];
                             v->fit_diet_EST[pred][y]   -= v->fit_diet_EST[link][y];
														  
														}                       
                         }	               
	               }
	            // Loop through species, applying Adams-Basforth or fast
	            // equilibrium depending on species.
	            
								 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                     
     					    // Adjust feeding time
	                   if (v->RRR.rpar_NoIntegrate[sp] < 0){pd = v->stanzaPred[sp];}
	                   else                         {pd = v->state_BB[sp];}
	                
									   //if (isnan(v->state_Ftime[sp])){
										 //   cout << "now" << sp << endl;
										 //}	                   
                     if((pd > 0) && (v->deriv_FoodGain[sp]>0)){
										      v->state_Ftime[sp] = 0.1 + 0.9 * v->state_Ftime[sp] * 
                                ((1.0 - v->RRR.rpar_FtimeAdj[sp]/(double)STEPS_PER_MONTH) 
																+ v->RRR.rpar_FtimeAdj[sp]/(double)STEPS_PER_MONTH * 
										            v->RRR.rpar_FtimeQBOpt[sp] / (v->deriv_FoodGain[sp] / pd));
                         }	

									   if (isnan(v->state_Ftime[sp])){
										    cout << "later" << sp <<" "<< v->RRR.rpar_FtimeQBOpt[sp] <<" "<< v->deriv_FoodGain[sp] <<" "<< pd << endl;
										 }	  
                     MAX_THRESHOLD(v->state_Ftime[sp], 2.0);
                     
                  // Biomass update for non-split groups
								     old_B=v->state_BB[sp];
								     new_B=v->state_BB[sp];
								     if (v->RRR.rpar_NoIntegrate[sp] == 0){
										    new_B        = (1.0 - SORWT) * v->deriv_biomeq[sp] 
												               + SORWT * v->state_BB[sp];
										    MIN_THRESHOLD(new_B,v->RRR.rpar_B_BaseRef[sp]*EPSILON);
										    MAX_THRESHOLD(new_B,v->RRR.rpar_B_BaseRef[sp]*BIGNUM);
										 }
										 else if (v->RRR.rpar_NoIntegrate[sp] == sp){
										  // First Timestep left out for now (is this a problem?)
											// if t=0 v->deriv_BB(i) = v->deriv_BB(i) + v->deriv_DerivT(i) * DELTA_T
                           new_B = v->state_BB[sp] + (DELTA_T / 2.0) * 
                                  (3.0 * (v->deriv_DerivT[sp]) - v->deriv_dyt[sp]); 
                           MIN_THRESHOLD(new_B,v->RRR.rpar_B_BaseRef[sp]*EPSILON);
                           MAX_THRESHOLD(new_B,v->RRR.rpar_B_BaseRef[sp]*BIGNUM);
										       }

                    if (isnan(new_B)){
                       cout << sp << " Y: " << y << " M: " << m << " D: " << d << " nan " << endl;
                       cout << v->deriv_FoodLoss[sp] << " FoodLoss" << endl;
                       cout << v->deriv_FoodGain[sp] << " FoodGain" << endl;
                       cout << v->deriv_UnAssimLoss[sp] << " UnAssimLoss" << endl;
                       cout << v->deriv_ActiveRespLoss[sp] << " ActiveRespLoss" << endl;   
                       cout << v->deriv_DetritalGain[sp] << " DetritalGain" << endl;
                       cout << v->deriv_FishingGain[sp] << " FishingGain" << endl;
                       cout << v->deriv_MzeroLoss[sp] << " MzeroLoss" << endl;
                       cout << v->deriv_FishingLoss[sp] << " FishingLoss" << endl;
                       cout << v->deriv_DetritalLoss[sp] << " DetritalLoss" << endl;
                       cout << v->deriv_FishingThru[sp] << " FishingThru" << endl;                    

                       cout << v->deriv_dyt[sp] << " dyt" << endl;  
											 cout << v->deriv_DerivT[sp] << " DerivT" << endl;  
											 cout << v->deriv_biomeq[sp] << " biomeq" << endl;                        
                       
                       for (prey=0; prey<=NUM_LIVING; prey++){
                            cout << prey << " ,"
                            << v->fit_vector[NUM_LIVING*0 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*1 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*2 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*3 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*4 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*5 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*6 + prey] << " ,"
                            << v->fit_vector[NUM_LIVING*6 + prey] << endl;
                       }
											 cout << v->RRR.rpar_ActiveRespFrac[sp]   << " ," 
							              << v->RRR.rpar_UnassimRespFrac[sp]  << endl;
							                            
							                            
                       for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
		                        prey = v->RRR.rpar_PreyFrom[links];
		                        pred = v->RRR.rpar_PreyTo[links];
		                        if (( prey==sp) || (pred==sp)){                  
		                        cout << prey << " " << pred << " " <<
		                        v->state_BB[prey] << " " << 
		                        v->state_BB[pred] << " " << 
		                        v->deriv_ConsMat[prey][pred] << " " << 
                            v->RRR.rpar_QQ[links] << " " << 
										        v->RRR.rpar_HandleSwitch[links] << " " << 
				                    v->RRR.rpar_DD[links] << " " <<  
						                v->RRR.rpar_HandleSelf[pred] << " " <<  
                            v->RRR.rpar_VV[links] << " " << 
                            v->RRR.rpar_ScrambleSelf[pred] << " " << 
                            endl;}
          
		                        }

                       exit(0);
										}
										    
										v->state_BB[sp] = new_B;
										//MIN_THRESHOLD(v->state_BB[sp],v->RRR.rpar_B_BaseRef[sp]*EPSILON);                                 
						     // sum up catch at every time step
								    if (fabs(v->state_BB[sp]/old_B-1.0) > EPSILON){
								       out_CC[sp] += ((v->deriv_FishingLoss[sp] * DELTA_T)/old_B) *
									                       (v->state_BB[sp] - old_B) /
										 	    							 log(v->state_BB[sp]/old_B);
									  }
									  else {out_CC[sp] += (v->deriv_FishingLoss[sp] * DELTA_T);
								    }		  

								 /*if (isnan(v->state_BB[sp]))
								    {cout << sp << " between " << y  << " " << m << " " << d << endl; 
										     cout << y   << "," <<  m << "," <<  d << "," << v->state_BB[sp] << ",";
		                     cout << sp << "," << v->deriv_dyt[sp]  << "," << v->deriv_DerivT[sp] << "," 
		                          << v->deriv_TotGain[sp] << "," << v->deriv_TotLoss[sp]<< "," <<endl;
		           
		              cout << sp << "," << v->deriv_FoodGain[sp]  << "," << v->deriv_DetritalGain[sp] << "," 
		                     << v->deriv_FishingGain[sp] << ","<<endl;		
		           
		                 cout << sp << "," << v->deriv_UnAssimLoss[sp]  << "," << v->deriv_ActiveRespLoss[sp] << "," 
		                 << v->deriv_FoodLoss[sp] << "," << v->deriv_MzeroLoss[sp] << "," << 
							          v->deriv_FishingLoss[sp] << "," << v->deriv_DetritalLoss[sp] << "," << endl;			        
		                    exit(0);
										}*/
									         
								 } // End of species loop
								

								    
 			        // Decide whether to stop running due to too high or low
 			           if (y<BURN_TIME){  // only do this in the first BURN_TIME years
 			                v->DISCARD_YEAR=0;
 			                // Not checking detritus because Discards and Offal go to 0
 			                // when starting with 0 catch.
 			                // PROTECT JUVENILES FROM DEATH
                      /*if (!Jlist[sp])*/{
                      for (sp=1; sp<=NUM_LIVING; sp++){
                         if (v->state_BB[sp] < v->RRR.rpar_B_BaseRef[sp] * LO_DISCARD){v->DISCARD_YEAR=y; 
							                                                          v->discarded[sp]=1; }
                         if (v->state_BB[sp] > v->RRR.rpar_B_BaseRef[sp] * HI_DISCARD){v->DISCARD_YEAR=y;
							                                                             v->discarded[sp]=1; v->discarded[sp]=1; }
	                    }
	                    }
				              // Check for NaN and infinity errors
                         for (sp=1; sp<=NUM_LIVING; sp++){
                             if (isinf(v->state_BB[sp]) || isnan(v->state_BB[sp])){
										             v->DISCARD_YEAR=y+1; v->discarded[sp]=1;
										             //cout << "boo!" <<endl;
										         }
								          }		
							  
				              if (v->DISCARD_YEAR){y=EndYear; m=STEPS_PER_YEAR; d=STEPS_PER_MONTH;}	
				         }
				         

				       
					}  // End of days loop

          // Stanza calculations
					   update_stanzas(v,y,m+1);
             SplitSetPred(v);
          								           	 // cout << "jod5 " << m << endl;  
          for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
              if (isinf(v->state_BB[sp]) || isnan(v->state_BB[sp])){
                 inbound++;
                 // Just show the first 10 crashed, often one crash drives the
                 // whole system.
                 if (inbound<=10){
                    cout << y << " " << m << " " << sp << " out of bounds" << endl;
                    }
                 //if (y<BURN_TIME){
                 //   v->DISCARD_YEAR=y+1; v->discarded[sp]=1;
                 //   y=EndYear; m=STEPS_PER_YEAR; d=STEPS_PER_MONTH;
                 //}
                       
              }
          }    
					// set to forced biomass
  		           for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
  		               //if ((  FORCED_BIO[sp][y+1] > 0.0) && (!(FORCED_BIO[sp][y] >0.0)))
                     //   {FORCED_BIO[sp][y] = v->state_BB[sp];
                     //   }
  		               if ((FORCED_BIO[sp][y+1]>0.0))

  		                  {
  		                     if (m==0){startyear_forcebio[sp] = v->state_BB[sp];}
                           ystep    = ((double)(m)) / ((double)(STEPS_PER_YEAR));
                           newbio   = startyear_forcebio[sp] +  ystep * (FORCED_BIO[sp][y+1] - startyear_forcebio[sp]);

                           bioanom  = newbio - v->state_BB[sp];
                           v->fit_bioanom[sp][y] += bioanom;
                           //TRUNC(v->fit_bioanom[sp][y]); // not ideal, rather do it once
                           bioratio = newbio/v->state_BB[sp];
                           //cout << sp << " " << y << " " << ystep << " " << FORCED_BIO[sp][y] << " " << FORCED_BIO[sp][y+1] << endl;
                           if  ( (v->RRR.rpar_NoIntegrate[sp] == sp) || (v->RRR.rpar_NoIntegrate[sp] == 0) ){
                                v->state_BB[sp] *= bioratio;
                                //cout << "Sboo!" << path_species[sp] << " " << v->RRR.rpar_NoIntegrate[sp] << endl;
                           }
                           else {
                                 //for (i = 1; i<= juv_N; i++){
                                 //cout << "Aboo!" << sp << endl;
                                 i = Alist[sp];
                                 if (i>0){
                                         //// cout << "jboo!" << sp << endl;
                                      /*if (sp == juv_AduNum[i])*/


                                          for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++)
                                              {
                                              v->NageS[i][ageMo] *= bioratio;
                                              }
                                              v->NageS[i][v->firstMoJuv[i]] *= bioratio;
                                              v->SpawnBio[i] *= bioratio;
                                      }
                                 i = Jlist[sp];
                                 if (i>0){ 
                                          //cout << "boo!" << sp << endl;
                                          for (ageMo=v->firstMoJuv[i]+1; ageMo<=v->lastMoJuv[i]; ageMo++)
                                              {
                                              v->NageS[i][ageMo] *= bioratio;
                                              }
                                              //v->NageS[i][v->firstMoJuv[i]] *= bioratio;
                                              //v->SpawnBio[i] *= bioratio;
                                 }																     
                           }                                    




                        }
                }	
								           	 // cout << "jod4 " << m << endl;  
          // Update stanza state vars
					   SplitSetPred(v);
					   for (i = 1; i<= juv_N; i++){
                v->out_EE[i][y] += v->NageS[i][v->firstMoJuv[i]];
                v->out_RR[i][y] += v->NageS[i][v->firstMoAdu[i]] * v->WageS[i][v->firstMoAdu[i]];
             }
 			       if (y<BURN_TIME){  // only do this in the first BURN_TIME years
					    for (i=1; i<=juv_N; i++){
					        //TURN OFF JUVENILE DEATH
					        sp=juv_JuvNum[i];

                       if (v->state_BB[sp] < v->RRR.rpar_B_BaseRef[sp] * LO_DISCARD){v->DISCARD_YEAR=y; 
							                                                        v->discarded[sp]=1;}
                       if (v->state_BB[sp] > v->RRR.rpar_B_BaseRef[sp] * HI_DISCARD){v->DISCARD_YEAR=y;
							                                                           v->discarded[sp]=1;}

                 sp=juv_AduNum[i];
                         if (v->state_BB[sp] < v->RRR.rpar_B_BaseRef[sp] * LO_DISCARD){v->DISCARD_YEAR=y; 
							                                                          v->discarded[sp]=1;}
                         if (v->state_BB[sp] > v->RRR.rpar_B_BaseRef[sp] * HI_DISCARD){v->DISCARD_YEAR=y;
							                                                             v->discarded[sp]=1;}							                                                             
					      //v->out_SpawnBio[i][y*STEPS_PER_YEAR+m+1]    = SpawnBio[i];
                //v->out_SpawnEnergy[i][y*STEPS_PER_YEAR+m+1] = SpawnEnergy[i];
                //v->out_eggs[i][y*STEPS_PER_YEAR+m+1]        = NageS[i][0];
                //v->out_Nrec[i][y*STEPS_PER_YEAR+m+1]        = NageS[i][v->firstMoAdu[i]];
                //v->out_WavgJuv[i][y*STEPS_PER_YEAR+m+1]     = v->state_BB[juv_JuvNum[i]]/v->state_NN[juv_JuvNum[i]];
                //v->out_WavgAdu[i][y*STEPS_PER_YEAR+m+1]     = v->state_BB[juv_AduNum[i]]/v->state_NN[juv_AduNum[i]];
                //v->out_GGJuv[i][y*STEPS_PER_YEAR+m+1]       = v->stanzaGGJuv[i];
                //v->out_GGAdu[i][y*STEPS_PER_YEAR+m+1]       = v->stanzaGGAdu[i];
                //v->out_stanzaPredJuv[i][y*STEPS_PER_YEAR+m+1] = v->stanzaPred[juv_JuvNum[i]];
                //v->out_stanzaPredAdu[i][y*STEPS_PER_YEAR+m+1] = v->stanzaPred[juv_AduNum[i]];
              }
             }
             if (v->DISCARD_YEAR){y=EndYear; m=STEPS_PER_YEAR; d=STEPS_PER_MONTH;}
       
       //cout << "gauss" << endl;
       // WHY IS THERE NOISE HERE?????
       //      for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
       //         v->month_noise[sp] = NOISE_SCALE * gaussian(&v->rng);
			 //	     }
				//cout << "month"  << endl;
								           	 // cout << "jod3 " << m << endl;  
        // Output juvenile details
           if (AGEOUT){
                for (j=1; j<=juv_N; j++){
    	  	     // Save age structure to a file (only for one species numbered in Adams-Basforth)
                    fprintf(rptr[j],"%d,",y);
                    for (ageMo=v->firstMoJuv[j]; ageMo<v->lastMoAdu[j]; ageMo+=12){
                        fprintf(rptr[j],"%g,%g,",v->NageS[j][ageMo], v->WageS[j][ageMo]);
                    }
                    fprintf(rptr[j],"%g,%g,\n",
                            v->NageS[j][v->lastMoAdu[j]], 
                            v->WageS[j][v->lastMoAdu[j]]);
                }
          }
      // cout << "jod3a " << m << endl;  
		    }  // End of months loop
  
		    for (c=1; c<=fitN; c++){
             switch (fit_type[c]){
               case 0:
               //   v->fit_EST[c][y] = v->state_BB[fit_ApplyTo[c]];
               break;
               case 6:
                  v->fit_EST[c][y] = out_CC[fit_ApplyTo[c]];
                  //TRUNC(v->fit_EST[c][y]);
               break;
               default:
               break;
             }
        }
        
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
             v->out_CC[sp][y] = out_CC[sp];
             v->out_TotF[sp][y] = v->out_CC[sp][y]/v->out_BB[sp][y];
        }
        
        
		    if (y==fit_lastyear){y==EndYear;}
		  
								           	 // cout << "jod2 " << y << endl;  
      }  // End of years loop
      
      // Close juvenile output files
      if (AGEOUT){
             for (j=1; j<=juv_N; j++){
                 fclose(rptr[j]);
             }
      }
      

}

//------------------------------------------------------------------------------
int base_system(struct SimRun *v){

int i;
     //cout << "umpapa" << endl;
     path_to_rates(v);
		 //cout << "humpapa" << endl;
		 // For changing vuls, for example
     //   for (i=1; i<=v->RRR.rpar_NumPredPreyLinks; i++){
		 //       v->RRR.rpar_XX[i] = 1. + exp(5.0*(uniform(&v->rng)-0.5));
		 //  }
		 		  
     initialize_stanzas(v);
    //cout << "dumpapa" << endl;
		 // Or for changing Ftime
		 //   for (i=1; i<=NUM_LIVING; i++){
		 //       v->RRR.rpar_FtimeAdj=uniform(&v->rng);
		 //    }
		 return 0;
}

// -----------------------------------------------------------------------------

int path_to_rates(struct SimRun *v){

    int sp, links, LL;
    float AA, VV;
    //struct SimRun v;
    
		// Energetics for Living and Dead Groups
    for (sp=1; sp<=NUM_LIVING + NUM_DEAD; sp++){
		    // Reference biomass for calculating YY
					 v->RRR.rpar_B_BaseRef[sp] = path_BB[sp]; 
				// Mzero proportional to (1-EE)
		       v->RRR.rpar_MzeroMort[sp] = path_PB[sp] * (1.0 - path_EE[sp]);
		    // Unassimilated is the proportion of CONSUMPTION that goes to
		    // detritus.  Note in Viz, this is written as proportion of
		    // respiration, so slightly different bookkeeping
           v->RRR.rpar_UnassimRespFrac[sp] = path_GS[sp];
        // Active respiration is proportion of CONSUMPTION that goes
        // to "heat"
        // TODO:  Passive respiration is left out here, as it should be
        // merged with bioenergetics of Von Bertalanffy
           if (path_QB[sp] > EPSILON ){
	  			    v->RRR.rpar_ActiveRespFrac[sp]    = 1.0 - (path_PB[sp]/path_QB[sp])
							                                 - v->RRR.rpar_UnassimRespFrac[sp];
				      v->RRR.rpar_FtimeAdj[sp]  = 0.0; 
              v->RRR.rpar_FtimeQBOpt[sp] = path_QB[sp];
							v->RRR.rpar_PBopt[sp] = path_PB[sp];						                                 
		    		}
				    else {
						  v->RRR.rpar_ActiveRespFrac[sp]  = 0.0;
				      v->RRR.rpar_FtimeAdj[sp]  = 0.0; 
              v->RRR.rpar_FtimeQBOpt[sp] = 1.0;
							v->RRR.rpar_PBopt[sp] = path_PB[sp];				    
				    }
				// NoIntegrate setting (done for QB, calculated this way so detritus
				// and PrimProd are considered PB/1.0)
				// TODO:  IS THIS LIMIT CORRECT??
				// CHANGING STEP METHOD:  OLD Changes NoIntegrate Based On Current
				// OPTION 1
           if (path_PB[sp]/(1.0 - v->RRR.rpar_ActiveRespFrac[sp] - v->RRR.rpar_UnassimRespFrac[sp]) >
					    2 * STEPS_PER_YEAR * STEPS_PER_MONTH)   
					         {v->RRR.rpar_NoIntegrate[sp]= 0;}
					 else    {v->RRR.rpar_NoIntegrate[sp]=sp;} 
				// OPTION 2:  Change based on Original Only: Note:  random systems keep it, change that??
        //   if (path_QB[sp] >
				//	    2 * STEPS_PER_YEAR * STEPS_PER_MONTH)   
				//	         {v->RRR.rpar_NoIntegrate[sp]= 0;}
				//	 else    {v->RRR.rpar_NoIntegrate[sp]=sp;}         	 
				// OPTION 3:  Change based on fixed Steps Per Month of 1
        //   if (path_PB[sp]/(1.0 - v->RRR.rpar_ActiveRespFrac[sp] - v->RRR.rpar_UnassimRespFrac[sp]) >
				//	    2 * STEPS_PER_YEAR * 1)   
				//	         {v->RRR.rpar_NoIntegrate[sp]= 0;}
				//	 else    {v->RRR.rpar_NoIntegrate[sp]=sp;} 				
				// OPTION 4: 2 and 3 combined
        //   if (path_QB[sp] >
				//	    2 * STEPS_PER_YEAR * 1)   
				//	         {v->RRR.rpar_NoIntegrate[sp]= 0;}
				//	 else    {v->RRR.rpar_NoIntegrate[sp]=sp;} 				
				
					 
					 
				// Feeding Time parameters

		}

    // Defaults for group 0 ('outside')
       v->RRR.rpar_B_BaseRef[0]=1.0;  
	   	 v->RRR.rpar_MzeroMort[0]=0.0;
       v->RRR.rpar_UnassimRespFrac[0]=0.0; 
		   v->RRR.rpar_ActiveRespFrac[0]=0.0;
       v->RRR.rpar_FtimeQBOpt[0] = 1.0;
			 v->RRR.rpar_PBopt[0] = 1.0;	
       v->RRR.rpar_HandleSelf[0]   = HandleSelfWt;
			 v->RRR.rpar_ScrambleSelf[0] = ScrambleSelfWt;
			 //v->RRR.rpar_HandleSwitch[0] = PREYSWITCH;

    // Fishing Effort defaults to 0 for non-gear, 1 for gear
       for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++)           {v->fish_Effort[sp]=0.0;}
       for (sp=NUM_LIVING+NUM_DEAD+1; sp<=NUM_GROUPS; sp++){v->fish_Effort[sp]=1.0;}

		// Now we initialize predator/prey pairs
		   v->RRR.rpar_NumPredPreyLinks = 0;
		   
		// Start by putting in primary production as a consumption of group 0
		// so we use vuls to set PP curve rather than maxPB term
       for (sp = 1; sp<=NUM_LIVING; sp++){
            // This test for Prim Prod assumes that there are no "mixed"
            // groups with both PP and consumption
            if ((path_PB[sp]>EPSILON) && (path_QB[sp]<=EPSILON)){
						    v->RRR.rpar_NumPredPreyLinks++;
						    v->RRR.rpar_PreyFrom[v->RRR.rpar_NumPredPreyLinks] = 0;
						    v->RRR.rpar_PreyTo[v->RRR.rpar_NumPredPreyLinks] = sp;
						    //v->RRR.rpar_XX[v->RRR.rpar_NumPredPreyLinks] = MDULL;
						    //v->RRR.rpar_HH[v->RRR.rpar_NumPredPreyLinks] = MSINGLEHANDLE;
                v->RRR.rpar_DD[v->RRR.rpar_NumPredPreyLinks] = MHANDLE;
								v->RRR.rpar_VV[v->RRR.rpar_NumPredPreyLinks] = MSCRAMBLE;						    
						    v->RRR.rpar_QQ[v->RRR.rpar_NumPredPreyLinks] = path_BB[sp]*path_PB[sp];						    

						}   

       }
    
    // Now we look through all predator prey pairs
        int pred, prey;
        for (prey=1; prey<=NUM_GROUPS; prey++){
            v->RRR.rpar_HandleSelf[prey]   = HandleSelfWt;
						v->RRR.rpar_ScrambleSelf[prey] = ScrambleSelfWt;
					  
				    for (pred=1; pred<=NUM_GROUPS; pred++){
						    if (path_DC[prey][pred]>EPSILON){
								   v->RRR.rpar_NumPredPreyLinks++;
						       v->RRR.rpar_PreyFrom[v->RRR.rpar_NumPredPreyLinks] = prey;
						       v->RRR.rpar_PreyTo[v->RRR.rpar_NumPredPreyLinks] = pred;
						       //v->RRR.rpar_XX[v->RRR.rpar_NumPredPreyLinks] = MDULL;
						       //v->RRR.rpar_HH[v->RRR.rpar_NumPredPreyLinks] = MSINGLEHANDLE;
						       v->RRR.rpar_HandleSwitch[v->RRR.rpar_NumPredPreyLinks] = PREYSWITCH;
						       v->RRR.rpar_DD[v->RRR.rpar_NumPredPreyLinks] = MHANDLE;
								   v->RRR.rpar_VV[v->RRR.rpar_NumPredPreyLinks] = MSCRAMBLE;	
						       v->RRR.rpar_QQ[v->RRR.rpar_NumPredPreyLinks] = path_QB[pred] *
									                                     path_BB[pred] *
									                                     path_DC[prey][pred];								   
								}
						}
				}
				
      
 // For Handling time and Hotspots
 
    LL=(NUM_GROUPS+1)*sizeof(float);
      memset(v->RRR.rpar_PredTotWeight, 0, LL);
			memset(v->RRR.rpar_PreyTotWeight, 0, LL);  
 
    for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
		    prey = v->RRR.rpar_PreyFrom[links];
		    pred = v->RRR.rpar_PreyTo[links];
                   VV = v->RRR.rpar_VV[links] * 
                        v->RRR.rpar_QQ[links] / 
                        path_BB[prey]; 
                   AA = (2 * v->RRR.rpar_QQ[links] * VV) / 
                        (VV * path_BB[pred] * path_BB[prey] - 
                        v->RRR.rpar_QQ[links] * path_BB[pred]); 
                  v->RRR.rpar_PredPredWeight[links] = AA * path_BB[pred]; 
                  v->RRR.rpar_PreyPreyWeight[links] = AA * path_BB[prey]; 
                  v->RRR.rpar_PredTotWeight[prey] += v->RRR.rpar_PredPredWeight[links]; 
                  v->RRR.rpar_PreyTotWeight[pred] += v->RRR.rpar_PreyPreyWeight[links];
     }
				int li;
				for (li=1; li <= v->RRR.rpar_NumPredPreyLinks; li++){

        //if (path_QB[RRR.rpar_PreyTo[li]]>EPSILON){
            //cout << li 
						//        << "," << RRR.rpar_PreyFrom[li] 
					//					<< "," << RRR.rpar_PreyTo[li]  
						//        << "," << RRR.rpar_PredPredWeight[li] 
						//				<< "," << RRR.rpar_PreyPreyWeight[li]
						//        << "," << RRR.rpar_PredTotWeight[RRR.rpar_PreyFrom[li]]
						//        << "," << RRR.rpar_PreyTotWeight[RRR.rpar_PreyTo[li]]
						//        << endl;
            v->RRR.rpar_PredPredWeight[li] = v->RRR.rpar_PredPredWeight[li]/ v->RRR.rpar_PredTotWeight[v->RRR.rpar_PreyFrom[li]]; 
            v->RRR.rpar_PreyPreyWeight[li] = v->RRR.rpar_PreyPreyWeight[li]/ v->RRR.rpar_PreyTotWeight[v->RRR.rpar_PreyTo[li]]; 
            }       
        //}   
        
     // GEAR FISHING
        int grn, gr, dest,det;
        v->RRR.rpar_NumFishingLinks = 0;
        for (grn=1; grn<=NUM_GEARS; grn++){
             //RRR.rpar_B_BaseRef[sp] = 1.0;
             gr=grn+NUM_LIVING+NUM_DEAD;
     		     for(prey=1; prey<=NUM_GROUPS; prey++){
				        if (path_Catch[grn][prey] > EPSILON * v->RRR.rpar_B_BaseRef[prey]){
				            dest=0;
				            v->RRR.rpar_NumFishingLinks++;
				            v->RRR.rpar_FishingFrom[v->RRR.rpar_NumFishingLinks]=prey;
				            v->RRR.rpar_FishingThrough[v->RRR.rpar_NumFishingLinks]=gr;
				            v->RRR.rpar_FishingTo[v->RRR.rpar_NumFishingLinks]=dest;
				            v->RRR.rpar_FishingQ[v->RRR.rpar_NumFishingLinks]= 
										     path_Catch[grn][prey]/v->RRR.rpar_B_BaseRef[prey];
								}
				        if (path_Discards[grn][prey] > EPSILON * v->RRR.rpar_B_BaseRef[prey]){
				        for (det=1; det<=NUM_DEAD; det++){
				            if (path_DetFate[det][gr] > EPSILON){
				                dest = det + NUM_LIVING;
				                v->RRR.rpar_NumFishingLinks++;
				                v->RRR.rpar_FishingFrom[v->RRR.rpar_NumFishingLinks]=prey;
				                v->RRR.rpar_FishingThrough[v->RRR.rpar_NumFishingLinks]=gr;
				                v->RRR.rpar_FishingTo[v->RRR.rpar_NumFishingLinks]=dest;
				                v->RRR.rpar_FishingQ[v->RRR.rpar_NumFishingLinks]= 
										         path_Discards[grn][prey]*path_DetFate[det][gr]/
														 v->RRR.rpar_B_BaseRef[prey];				                
										}
				        } 
				        
                }
            }
        }

   int i;
   for (i = 1; i<= juv_N; i++){
        v->RRR.rpar_drawn_K[i]    =  juv_VonBK[i];
	      v->RRR.rpar_drawn_AduZ[i] =  juv_aduEqAgeZ[i];
		    v->RRR.rpar_drawn_JuvZ[i] =  juv_juvEqAgeZ[i];
		    v->RRR.rpar_SpawnX[i]            = 10000.0;
		    //v->RRR.rpar_SpawnX[i]            = 2.0;
        v->RRR.rpar_SpawnEnergy[i]       = 1.0; 
        v->RRR.rpar_SpawnAllocR[i]       = 1.0;
        v->RRR.rpar_SpawnAllocG[i]       = 1.0;
   }
   
    //cout << v->RRR.rpar_NumPredPreyLinks << " " << NUM_PREDPREYLINKS << endl;
    //cout << v->RRR.rpar_NumFishingLinks << " " << NUM_CATCHLINKS << endl;
   // Set Start Biomass and other state variables
	    for (sp=0; sp<=NUM_GROUPS; sp++){
			     v->state_BB[sp]    = v->RRR.rpar_B_BaseRef[sp];
			     //cout << sp << "," << RRR.rpar_B_BaseRef[sp] << endl;
			     v->state_Ftime[sp] = 1.0;
			     v->discarded[sp]   = 0;
			}
     // Detrital links
     //   int RRR.rpar_NumDetLinks; 
     //   int RRR.rpar_DetFrom[NUM_DEAD*NUM_GROUPS+1]; 
     //   int RRR.rpar_DetTo[NUM_DEAD*NUM_GROUPS+1]; 
		 //cout << "finished path to rates." <<endl; 
}

// ----------------------------------------------------------------------------
// Initialize juvenile adult or "stanza" age structure in sim

//KYA Got rid of SplitNo and SplitWage as they were duplicates
//    for NageS and WageS

int initialize_stanzas(struct SimRun *v){
//  monthly age, species counter
    int ageMo, i, links, sp;
    
//  survival rate, previous survival rate
    double survRate_juv_month;
    double survRate_adu_month;
    //float survByMonth;
    //float prevSurv;   
   
//  sum biomass, sum juv biomass, sum to calculate consumption, consumption
    double BioPerEgg;  // Was SumBio
    double juvBio; 
    double sumK;
    double Qmult[NUM_GROUPS+1];
    double Qcheck[NUM_GROUPS+1];
    double juvCons;
    double juvQB;
    double aduQB;
    //float drawn_K, drawn_AduZ, drawn_JuvZ;  
		    
 // Set Q multiplier to 1 (used to rescale consumption rates)
    for (i = 0; i<= NUM_GROUPS; i++){Qmult[i]=1.0;}
    
 // loop over split species groups to initialize all parameters  
    for (i = 1; i<= juv_N; i++){
		// RANDOMIZATION
		   //v->RRR.rpar_drawn_K    =   (1. + RandScale *0.5   * ( uniform(&v->rng)- 0.5)) * juv_VonBK[i];
			 //v->RRR.rpar_drawn_AduZ =   (1. + RandScale *0.5   * ( uniform(&v->rng)- 0.5)) * juv_aduEqAgeZ[i];
			 //v->RRR.rpar_drawn_JuvZ =   (1. + RandScale *1.0   * ( uniform(&v->rng)- 0.5)) * juv_juvEqAgeZ[i];
			         
    // first calculate the number of months in each age group.  
	  // For end stanza (plus group), capped in EwE @400, this is 90% of 
		// relative max wt, uses generalized exponent D
       v->firstMoJuv[i] = 0;
       v->lastMoJuv[i]  = juv_RecAge[i]* STEPS_PER_YEAR - 1;
       v->firstMoAdu[i] = juv_RecAge[i]* STEPS_PER_YEAR;
       v->lastMoAdu[i]  = (int)(  log (1.0 - pow(double(0.9) , double(1.- juv_VonBD[i]))) / 
				                       ( -3. * v->RRR.rpar_drawn_K[i] * (1.- juv_VonBD[i]) / STEPS_PER_YEAR) 
														);
       if (v->lastMoAdu[i]> MAX_MONTHS_STANZA) {v->lastMoAdu[i]= MAX_MONTHS_STANZA;} 
      // cout << "alast Mo:" << v->lastMoJuv[i] << ":" << i << ":" << endl;        
     // Energy required to grow a unit in weight (scaled to Winf = 1)
          v->vBM[i] = 1.0 - 3.0*v->RRR.rpar_drawn_K[i]/STEPS_PER_YEAR;

				   
    // fill monthly vectors for each species, rel weight and consumption at age
		// Uses generalized vonB (exponent is d).    
       for (ageMo = 0; ageMo <= v->lastMoAdu[i]; ageMo++){
            v->WageS[i][ageMo] = pow(double(1.0 -
						                       exp(-3. * v->RRR.rpar_drawn_K[i] * (1. - juv_VonBD[i]) * ageMo/STEPS_PER_YEAR)) ,
																	 double(1./(1.- juv_VonBD[i])));
            v->WWa[i][ageMo] = pow(double(v->WageS[i][ageMo]), double(juv_VonBD[i])); 
						//if (i==8){cout << ageMo << " " << WageS[i][ageMo] << " " << WWa[i][ageMo] << endl;}
        }
       //cout << endl << "LastMo " << lastMoAdu[i] << " vbM " << vBM[i] << endl;       
    // Survival rates for juveniles and adults
    // KYA Switched survival rate out of BAB == assume all on adults
       survRate_juv_month = exp(-(v->RRR.rpar_drawn_JuvZ[i])/ STEPS_PER_YEAR);
       //survRate_juv_month = exp(-(v->RRR.rpar_PBopt[juv_JuvNum[i]])/ STEPS_PER_YEAR);
       survRate_adu_month = exp(-(v->RRR.rpar_drawn_AduZ[i])/ STEPS_PER_YEAR);
           
       //survRate_adu_month = exp(-(v->RRR.rpar_PBopt[juv_AduNum[i]])/ STEPS_PER_YEAR);    
    // For numbers, first we figure out numbers per 1 egg (biomass per egg)
			 BioPerEgg = 0;  // This is ADULT biomass per egg
			 // First month is set to 1, so this is N per egg
    			v->NageS[i][0] = 1.0;
			 // Months 1...LastAdu-1
    			for (ageMo = 1; ageMo < v->lastMoAdu[i]; ageMo++){
           if (ageMo <= v->firstMoAdu[i]) {    
							 v->NageS[i][ageMo] = v->NageS[i][ageMo-1] * survRate_juv_month;
					 }
					 else {
								v->NageS[i][ageMo] = v->NageS[i][ageMo-1] * survRate_adu_month;
								//BioPerEgg += NageS[i][ageMo] * WageS[i][ageMo];
					 }
					 		   				 
			 }
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "a ";} cout << endl;			 
			 // Make last age group into a plus group
			 if (juv_RecMonth[i]>0){
    			v->NageS[i][v->lastMoAdu[i]] = v->NageS[i][v->lastMoAdu[i]-1] * pow(double(survRate_adu_month),double(12.)) 
							                         / (1. - pow(double(survRate_adu_month),double(12.)));
		   }
			 else {
    			v->NageS[i][v->lastMoAdu[i]] = v->NageS[i][v->lastMoAdu[i]-1] * survRate_adu_month 
							                         / (1. - survRate_adu_month);			 					                         
			 }
			 BioPerEgg += v->NageS[i][v->lastMoAdu[i]] * v->WageS[i][ageMo];
					
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << " ";} cout << endl;					
			 // FIX FIX for SARAH THIS NEXT LINE ADDED TO FIX MISSING AGE IN ABOVE LOOP
    		for (ageMo = 0; ageMo < v->lastMoAdu[i]; ageMo++){
					  if (juv_RecMonth[i]>0){ 
					        if ((ageMo + juv_RecMonth[i]) % 12 != 0 ){
									   v->NageS[i][ageMo] = 0.;
									}
						}
	  				if (ageMo >= v->firstMoAdu[i]) {
							 BioPerEgg += v->NageS[i][ageMo] * v->WageS[i][ageMo];
					  }		 		   				 
			 }			 
			 					
    // Actual Eggs is Biomass of path/Biomass per recruit
       v->recruits[i] = v->RRR.rpar_B_BaseRef[juv_AduNum[i]] / BioPerEgg;

    // REVISIT THIS!!!! KYA Also took BAB out of here (don't know if this is a juv. or adu term)
       v->RzeroS[i] = v->recruits[i];
       //RzeroS[i] = recruits[i] * exp(juv_aduBAB[i]/ STEPS_PER_YEAR);;      

    // Now Scale up all numbers to match actual (Ecopath or random) biomass
       for (ageMo = 0; ageMo <= v->lastMoAdu[i]; ageMo++){
            v->NageS[i][ageMo] *= v->recruits[i];         
       }
         
    // calc biomass of juvenile group using adult group input
       juvBio = 0;
       for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){
           juvBio += v->NageS[i][ageMo] * v->WageS[i][ageMo];
       }    
          
    // First calculate adult assimilated consumption
        sumK = 0;
        for (ageMo= v->firstMoAdu[i]; ageMo<= v->lastMoAdu[i]; ageMo++){
            sumK += v->NageS[i][ageMo] * v->WWa[i][ageMo];
        } 
        
		// SumK is adjusted to QB, ratio between the two is A term in respiration
		// Actually, no, since Winf could scale it up, it could be above or below
		// Q/B
			 sumK  = v->RRR.rpar_FtimeQBOpt[juv_AduNum[i]] * v->RRR.rpar_B_BaseRef[juv_AduNum[i]] / sumK;
       aduQB =  v->RRR.rpar_FtimeQBOpt[juv_AduNum[i]];

    // Calculate juvenile assimilated consumption        
        juvCons = 0;
        for (ageMo=v->firstMoJuv[i]; ageMo<= v->lastMoJuv[i]; ageMo++){
            juvCons += v->NageS[i][ageMo] * v->WWa[i][ageMo];
        }    
    //  Scale Juvenile consumption by same A as adults
        //Qmult[juv_JuvNum[i]] = sumK;
        juvCons *= sumK; //Qmult[juv_JuvNum[i]];
        juvQB = juvCons / juvBio; 
        
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "b ";} cout << endl;
    //  Calculate spawning biomass as the amount of biomass over Wmat.
        v->SpawnBio[i] = 0;
    //    cout << "alast Mo:" << v->lastMoJuv[i] << ":" << i << ":" << endl;
        for (ageMo = 0; ageMo <= v->lastMoAdu[i]; ageMo++){ 
            // REC_CHANGE
            //if (WageS[i][ageMo] > WmatWinf[i]){
            //    SpawnBio[i]    += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]);
            //    SpawnEnergy[i] += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]) 
						//		                                  * (WWa[i][ageMo] - (WageS[i][ageMo] - WageS[i][ageMo-1]));
						// }
            if ( (v->WageS[i][ageMo] > Wmat001[i]) && (ageMo > Amat001[i])) {
                v->SpawnBio[i]  += v->WageS[i][ageMo] * v->NageS[i][ageMo] /
                                (1. + exp(- ((v->WageS[i][ageMo]-Wmat50[i])/WmatSpread[i]) 
																		      - ((ageMo          - Amat50[i])/AmatSpread[i])
																));                             
            }
            // END REC_CHANGE   
        }
     //   cout << "balast Mo:" << v->lastMoJuv[i] << ":" << i << ":" << endl;
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c ";} cout << endl;
	  //cout << "Kansas" << endl;
    // BaseSpawnCost = RzeroS[i] * 1/(1. - sumK)     
    // KYA Spawn X is Beverton-Holt.  To turn off Beverton-Holt, set
    // SpawnX to 10000 or something similarly hight.  2.0 is half-saturation.
    // 1.000001 or so is minimum.
       //v->RRR.rpar_SpawnX[i]            = 10000.0;
       //v->RRR.rpar_SpawnEnergy[i]       = 1.0; 
       //v->RRR.rpar_SpawnAllocR[i]       = 1.0;
       //v->RRR.rpar_SpawnAllocG[i]       = 1.0;
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c1 ";} cout << endl;
       //v->baseSpawnEnergy[i]   = v->SpawnEnergy[i];
       v->baseSpawnBio[i]      = v->SpawnBio[i];
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c2 ";} cout << endl;
			 // REC_CHANGE  
          v->EggsStanza[i]        = v->SpawnBio[i];
          //EggsStanza[i]        = SpawnEnergy[i];
       // END REC_CHANGE
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c3 ";} cout << endl;
       v->baseEggsStanza[i]    = v->EggsStanza[i];
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c4 ";} cout << endl;
       v->Rbase[i] = v->NageS[i][v->lastMoJuv[i]+1] * 
                                 v->WageS[i][v->lastMoJuv[i]+1];     
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << " " << i << " c5 ";} cout << endl;   
    //   cout << "aclast Mo:" << v->lastMoJuv[i] << ":" << i << ":" << endl; 
    // Set NoIntegrate flag (negative in both cases, different than split pools)
       v->RRR.rpar_NoIntegrate[juv_JuvNum[i]] = -juv_JuvNum[i];
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "c6 ";} cout << endl;
       v->RRR.rpar_NoIntegrate[juv_AduNum[i]] = -juv_AduNum[i];
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "d ";} cout << endl;       
    // leaving out code to rescale recruitment if recruitment is seasonal
    // calculates RscaleSplit[i] to give same annual avg as equal monthly rec
       v->RscaleSplit[i] = 1;
    //cout << "Kansas1" << endl;
    // Now Reset juvenile B and QB in line with age stucture

						//cout << path_species[juv_JuvNum[i]] << " NewJuvBio " << juvBio << " " << RRR.rpar_B_BaseRef[juv_JuvNum[i]] 
						//     << " " << juvBio/RRR.rpar_B_BaseRef[juv_JuvNum[i]] << endl;
            //cout << path_species[juv_JuvNum[i]] << " NewJuvQB " << juvQB << " " << path_QB[juv_JuvNum[i]] 
						//     << " " << juvQB/path_QB[juv_JuvNum[i]]  << endl;
						//cout << endl;
						
						Qmult[juv_JuvNum[i]] = juvCons / 
						     (v->RRR.rpar_FtimeQBOpt[juv_JuvNum[i]]*v->RRR.rpar_B_BaseRef[juv_JuvNum[i]]);
				    v->RRR.rpar_FtimeAdj[juv_JuvNum[i]]  = 0.0;
				    v->RRR.rpar_FtimeAdj[juv_AduNum[i]]  = 0.0;
            v->RRR.rpar_FtimeQBOpt[juv_JuvNum[i]] = juvQB;
            v->RRR.rpar_FtimeQBOpt[juv_AduNum[i]] = aduQB; 
            v->RRR.rpar_B_BaseRef[juv_JuvNum[i]]  = juvBio;
		//							cout << "dalast Mo:" << v->lastMoJuv[i] << ":" << i << ":" << endl;					
		//for (sp = 1; sp<= juv_N; sp++){cout << v->lastMoJuv[sp] << "e ";} cout << endl;
    }
    //    for (i = 1; i<= juv_N; i++){
     //cout << "Oklhaa12 " << i << " " 
    //      << v->firstMoJuv[i] << " " << v->lastMoJuv[i] << endl;    
    //}
    
    
//cout << "Kansas2" << endl;
//  Reset QQ in line with new comsumtion
    for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
        //cout << RRR.rpar_QQ[links] << " " <<Qmult[RRR.rpar_PreyTo[links]] << endl;
        v->RRR.rpar_QQ[links] *= Qmult[v->RRR.rpar_PreyTo[links]];
    }
    // check for consistency
    for (links=0; links<=NUM_GROUPS; links++){Qcheck[links]=0.0;}
    for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
        Qcheck[v->RRR.rpar_PreyTo[links]] += v->RRR.rpar_QQ[links];
    }    
    
    //for (i = 1; i<= juv_N; i++){
    //        cout << i << " " << Qcheck[juv_JuvNum[i]] << " " << path_QB[juv_JuvNum[i]]* path_BB[juv_JuvNum[i]]
	  //					     << " " << Qcheck[juv_JuvNum[i]]/(path_QB[juv_JuvNum[i]]* path_BB[juv_JuvNum[i]])  << endl;    
    //  }
    // Having looped through all species, set state levels before calling derivative
    //cout << "Kansas3" << endl;
    //for (i = 1; i<= juv_N; i++){
    // cout << "Oklhaa11 " << i << " " 
    //      << v->firstMoJuv[i] << " " << v->lastMoJuv[i] << endl;    
    //}
    
    SplitSetPred(v);
    //cout << "Kansas4" << endl;
    for (i = 1; i<= juv_N; i++){    
        v->stanzaBasePred[juv_JuvNum[i]] = v->stanzaPred[juv_JuvNum[i]];
        v->stanzaBasePred[juv_AduNum[i]] = v->stanzaPred[juv_AduNum[i]];        
    }
    
//  CALCULATE FIRST DERIVATIVE to get deriv_FoodGain
    //deriv_master(v,0,0,0);        

        //what is correct param for StartEatenBy??? deriv_FoodGain???
        //so do I need to call deriv_master before this can be set up?
// KYA if this is growth per month, shouldn't deriv_foodgain be per month??
    //cout << "Kansas5" << endl;
    for (i = 1; i<= juv_N; i++){
        for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){
            v->SplitAlpha[i][ageMo] = (v->WageS[i][ageMo+1] - 
                                   v->vBM[i]*v->WageS[i][ageMo])*
                                   v->stanzaBasePred[juv_JuvNum[i]] /        //BB should be stanzaPred
                                   (v->RRR.rpar_FtimeQBOpt[juv_JuvNum[i]] * 
																	  v->RRR.rpar_B_BaseRef[juv_JuvNum[i]]);         //deriv_FoodGain[juv_JuvNum[i]]
        }
        for (ageMo=v->firstMoAdu[i]; ageMo<v->lastMoAdu[i]; ageMo++){        
            v->SplitAlpha[i][ageMo] = (v->WageS[i][ageMo+1] - 
                                   v->vBM[i]*v->WageS[i][ageMo])*
                                   v->stanzaBasePred[juv_AduNum[i]] /       //BB should be stanzaPred
                                   (v->RRR.rpar_FtimeQBOpt[juv_AduNum[i]] *
																	  v->RRR.rpar_B_BaseRef[juv_AduNum[i]]);       //deriv_FoodGain[juv_AduNum[i]]
        }
        v->SplitAlpha[i][v->lastMoAdu[i]] = v->SplitAlpha[i][v->lastMoAdu[i]-1];
        

    }
   //   cout << "Kentucky" << endl;
   //cout << "finished initialize stanzas" << endl;
}

// ----------------------------------------------------------------------------
// SplitSetPred function called in sim stanza initialize and update

int SplitSetPred(struct SimRun *v){
//  monthly age, species counter
    int ageMo, i;
//  holding vars for eggs, survival, growth for stanzas   
    double Bt, pt, Nt;
//cout << "Oklha" << endl;
    //for (i = 1; i<= juv_N; i++){
    // cout << "Oklhaa10 " << i << " " 
    //      << v->firstMoJuv[i] << " " << v->lastMoJuv[i] << endl;    
    //}
    for (i = 1; i<= juv_N; i++){
    //cout << "Oklhaa " << i << endl;
        Bt = 1e-30;
        pt = 1e-30;
        Nt = 1e-30;
        //loop over juv timesteps
     //cout << "Oklhaa9 " << Bt << pt << Nt << i << " " 
     //     << v->firstMoJuv[i] << " " << v->lastMoJuv[i] << endl;
        for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){
            //cout << ageMo << endl;
            Bt = Bt + v->NageS[i][ageMo] * v->WageS[i][ageMo];
            pt = pt + v->NageS[i][ageMo] * v->WWa[i][ageMo];
            Nt = Nt + v->NageS[i][ageMo];
            
        }
        //cout << "Oklha5 " << i << endl;
        v->state_BB[juv_JuvNum[i]]   = Bt;
        v->stanzaPred[juv_JuvNum[i]] = pt;
        v->state_NN[juv_JuvNum[i]]   = Nt;
//cout << "Oklha1" << endl;        
        //loop over adult timesteps
        Bt = 1e-30;
        pt = 1e-30;
        Nt = 1e-30;       
        for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++){
            Bt = Bt + v->NageS[i][ageMo] * v->WageS[i][ageMo];
            pt = pt + v->NageS[i][ageMo] * v->WWa[i][ageMo];
            Nt = Nt + v->NageS[i][ageMo];
        }
//cout << "Oklha2" << endl; 
        v->state_BB[juv_AduNum[i]] = Bt;
        v->stanzaPred[juv_AduNum[i]] = pt;
        v->state_NN[juv_AduNum[i]] = Nt;    
        v->stanzaGGJuv[i] =  v->deriv_FoodGain[juv_JuvNum[i]]/v->stanzaPred[juv_JuvNum[i]]; 
        v->stanzaGGAdu[i] =  v->deriv_FoodGain[juv_AduNum[i]]/v->stanzaPred[juv_AduNum[i]]; 
//cout << "Oklha3" << i << endl; 
    }
//cout << "Oklha3" << endl; 

   
}


// ----------------------------------------------------------------------------
// Update juvenile adult or "stanza" age structure during sim run 
// on monthly timesteps (not hardwiring months, but recommended)

int update_stanzas(struct SimRun *v, int yr, int mon){

//  monthly age, species counter
    int ageMo, i;
//  holding vars for eggs, survival, growth, numbers for stanzas   
    double Su, Gf, Nt;

    double propRepro;
    
    //loop over split species groups to update n, wt, biomass in sim  
    for (i = 1; i<= juv_N; i++){
        v->SpawnBio[i] = 0;
        //v->SpawnEnergy[i] = 0;
     // loop over juv timesteps with juv survival and growth
        Su = exp(-v->deriv_LossPropToB[juv_JuvNum[i]]/STEPS_PER_YEAR/v->state_BB[juv_JuvNum[i]]);    
        Gf = v->deriv_FoodGain[juv_JuvNum[i]]/v->stanzaPred[juv_JuvNum[i]];  // BB should be stanzaPred
             
				v->stanzaGGJuv[i] = Gf;

         //v->RRR.rpar_SpawnAllocR[i]       = 1.0;
         //v->RRR.rpar_SpawnAllocG[i]       = 1.0;
						         
        for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){   
            v->NageS[i][ageMo] = v->NageS[i][ageMo] * Su;
            v->WageS[i][ageMo] = v->vBM[i] * v->WageS[i][ageMo] + Gf * v->SplitAlpha[i][ageMo];
            // REC_CHANGE
            //if (WageS[i][ageMo] > WmatWinf[i]) {
            //    SpawnBio[i]    += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]);
            //    SpawnEnergy[i] += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]) 
						//		                                  * (WWa[i][ageMo] - (WageS[i][ageMo] - WageS[i][ageMo-1]));
            //}
            if ( (v->WageS[i][ageMo] > Wmat001[i]) && (ageMo > Wmat001[i])) {
                v->SpawnBio[i]  += v->WageS[i][ageMo] * v->NageS[i][ageMo] /
                                (1. + exp(- ((v->WageS[i][ageMo]-Wmat50[i])/WmatSpread[i]) 
																		      - ((ageMo          -Amat50[i])/AmatSpread[i])
																));                             
               //v->SpawnEnergy[i] = 1.0;
            }            
            // END REC_CHANGE
        }
				      
     // loop over adult timesteps with adult survival and growth
        Su = exp(-v->deriv_LossPropToB[juv_AduNum[i]]/STEPS_PER_YEAR/v->state_BB[juv_AduNum[i]]);
        Gf = v->deriv_FoodGain[juv_AduNum[i]]/v->stanzaPred[juv_AduNum[i]];   //BB should be stanzaPred
        
        v->stanzaGGAdu[i] = Gf;

        for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++){
            v->NageS[i][ageMo] = v->NageS[i][ageMo] * Su;
            v->WageS[i][ageMo] = v->vBM[i] * v->WageS[i][ageMo] + Gf * v->SplitAlpha[i][ageMo];
            // REC_CHANGE
            //if (WageS[i][ageMo] > WmatWinf[i]) {
            //   SpawnBio[i]    += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]);
            //    SpawnEnergy[i] += NageS[i][ageMo] * (WageS[i][ageMo] - WmatWinf[i]) 
						//		                                  * (WWa[i][ageMo] - (WageS[i][ageMo] - WageS[i][ageMo-1]));
            //}
						if ( (v->WageS[i][ageMo] > Wmat001[i]) && (ageMo > Amat001[i])) {
                v->SpawnBio[i]  += v->WageS[i][ageMo] * v->NageS[i][ageMo] /
                                (1. + exp(- ((v->WageS[i][ageMo]-Wmat50[i])/WmatSpread[i]) 
																		      - ((ageMo          -Amat50[i])/AmatSpread[i])
																));                             
               //v->SpawnEnergy[i] = 1.0;
            }     
						       
            // END REC_CHANGE
        }
        
        // KYA This is a Beverton-Holt curve between Spawning Biomass and
        // Number of Eggs Produced
        // Old way:  EggsStanza[i] = SpawnBio[i];
        // REC_CHANGE : back to the old way
        
        v->EggsStanza[i] = v->SpawnBio[i] * v->RRR.rpar_SpawnEnergy[i] * v->RRR.rpar_SpawnX[i] /
				                             (v->RRR.rpar_SpawnX[i] - 1.0 + 
			  															 (v->SpawnBio[i] / v->baseSpawnBio[i]));
			  															 
			  v->EggsStanza[i] *= force_byrecs[juv_JuvNum[i]][yr*STEPS_PER_YEAR+mon];		
        if (FORCED_REC[i][yr]>0.0){													 
           //v->EggsStanza[i] = v->baseEggsStanza[i] * FORCED_REC[i][yr]; //original to force to rec time series
           v->EggsStanza[i] *= FORCED_REC[i][yr]; //Sarah's alternate for relative rec forcing
           //cout << yr << " " << i << " " << FORCED_REC[i][yr]<< endl;
        }	
			  if (juv_RecMonth[i]>0){
				   if (mon == juv_RecMonth[i]) {v->EggsStanza[i] *=1. ;}
				 else                          {v->EggsStanza[i] *=0. ;}
				}
        //EggsStanza[i] = SpawnEnergy[i] * SpawnX[i] /
				//                             (SpawnX[i] - 1.0 + 
			  //														 (SpawnEnergy[i] / baseSpawnEnergy[i]));        
        // END REC_CHANGE
        
        //now update n and wt looping backward over age
        Nt = v->NageS[i][v->lastMoAdu[i]] + v->NageS[i][v->lastMoAdu[i]-1];
        if (Nt == 0){Nt = 1e-30;}
        v->WageS[i][v->lastMoAdu[i]] = (v->WageS[i][v->lastMoAdu[i]]*v->NageS[i][v->lastMoAdu[i]] + 
                           v->WageS[i][v->lastMoAdu[i]-1]*v->NageS[i][v->lastMoAdu[i]-1])/ Nt;
        v->NageS[i][v->lastMoAdu[i]] = Nt;
        
        for (ageMo=v->lastMoAdu[i]-1; ageMo>v->firstMoJuv[i]; ageMo--) {
            v->NageS[i][ageMo] = v->NageS[i][ageMo-1];
            v->WageS[i][ageMo] = v->WageS[i][ageMo-1];
  
        }

        // finally apply number of eggs        
        if (v->baseEggsStanza[i] > 0){
            v->NageS[i][v->firstMoJuv[i]] = v->RscaleSplit[i] * v->RzeroS[i] * 
                                      pow(double(v->EggsStanza[i] / v->baseEggsStanza[i]), 
                                      double(juv_RecPower[i]));
        }
        v->WageS[i][v->firstMoJuv[i]] = 0;

    // Uses generalized vonB (exponent is d).    
    // ADDED FOR STABILITY 4/13/07 (Unlucky Friday)
       for (ageMo = 0; ageMo <= v->lastMoAdu[i]; ageMo++){
            v->WWa[i][ageMo] = pow(double(v->WageS[i][ageMo]), double(juv_VonBD[i])); 
        }        
        

    }

    //last update biomass and pred index info for all spp
    //SplitSetPred(v);

}

// ----------------------------------------------------------------------------
int deriv_master(struct SimRun *v, int y, int m, int d){
    int sp, links, prey, pred, i;
    double Master_Density, Q;
    double prey_YY[NUM_GROUPS+1];
    double pred_YY[NUM_GROUPS+1];
    double PredSuite[NUM_GROUPS+1];
	  double HandleSuite[NUM_GROUPS+1];  
	  
    // NOT THREAD SAFE!!!! static float TerminalF[NUM_GROUPS+1];

    unsigned int LL;
    // Some derivative parts need to be set to zero
        LL=(NUM_GROUPS+1)*sizeof(double);
        memset(v->deriv_FoodLoss,		  0,LL);
        memset(v->deriv_FoodGain,		  0,LL);
        memset(v->deriv_UnAssimLoss,	  0,LL);
        memset(v->deriv_ActiveRespLoss,0,LL);   
        memset(v->deriv_DetritalGain,  0,LL);
        memset(v->deriv_FishingGain,   0,LL);
        memset(v->deriv_MzeroLoss,     0,LL);
        memset(v->deriv_FishingLoss,   0,LL);
        memset(v->deriv_DetritalLoss,  0,LL);
        memset(v->deriv_FishingThru,   0,LL);
	      memset(PredSuite,           0,LL);
	      memset(HandleSuite,         0,LL);
				    
  			LL=(NUM_GROUPS+1)*(NUM_GROUPS+1)*sizeof(double);
			  //memset(flowMat,0,LL);
			  memset(v->deriv_ConsMat,0,LL);
    //  Set YY = B/B(ref) for functional response 
    for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
		     prey_YY[sp]=v->state_Ftime[sp] * v->state_BB[sp]/v->RRR.rpar_B_BaseRef[sp];
		     pred_YY[sp]=v->state_Ftime[sp] * v->state_BB[sp]/v->RRR.rpar_B_BaseRef[sp];
		     //if (isnan(prey_YY[sp]))
		     //  {cout << sp << " " << v->state_Ftime[sp] << " " << v->state_BB[sp] << " " << v->RRR.rpar_B_BaseRef[sp] << endl;}
		}    
		prey_YY[0]=1.0;  pred_YY[0]=1.0;
		
		for (i=1; i<=juv_N; i++){
		  if (v->stanzaBasePred[juv_JuvNum[i]]>0){
		    pred_YY[juv_JuvNum[i]] = v->state_Ftime[juv_JuvNum[i]] * 
				        v->stanzaPred[juv_JuvNum[i]]/v->stanzaBasePred[juv_JuvNum[i]];
		    pred_YY[juv_AduNum[i]] = v->state_Ftime[juv_AduNum[i]] * 
				        v->stanzaPred[juv_AduNum[i]]/v->stanzaBasePred[juv_AduNum[i]];
		    }
		}
    // add "mediation by search rate" KYA 7/8/08
       for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
           pred_YY[sp] *= force_bysearch[sp][y*STEPS_PER_YEAR+m]; 
       }
       
	   //if (MDull<0.9999){
	      for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
		          prey = v->RRR.rpar_PreyFrom[links];
		          pred = v->RRR.rpar_PreyTo[links];
              PredSuite[prey]   += pred_YY[pred] * v->RRR.rpar_PredPredWeight[links];
              HandleSuite[pred] += prey_YY[prey] * v->RRR.rpar_PreyPreyWeight[links];
		  }
	   //}
	//if (y<1){
  //for (pred=1; pred<=NUM_GROUPS; pred++){
	//     cout << pred << "," << PredSuite[pred] << "," <<HandleSuite[pred] <<endl;
  //}
  //}
  //else {exit(0);}
    // CONSUMPTION (includes Primary Production)
    //  Calculate functional responses for all predator/prey links
    //LL=1*sizeof(float);
    //memset(&v->out_PP,		  0,LL);
//    for (sp=0; sp<=NUM_LIVING+NUM_DEAD; sp++){
        //if (PREYSWITCH == 1){
           //cout << sp << " "<< v->RRR.rpar_HandleSwitch[sp] << endl;
           //prey_YY[sp]     = pow( prey_YY[sp],   v->RRR.rpar_HandleSwitch[sp]);
           //HandleSuite[sp] = pow(HandleSuite[sp],v->RRR.rpar_HandleSwitch[sp]);
        //}
        //if (PREDSWITCH == 1){
        //   pred_YY[sp]     *= pred_YY[sp];
        //   PredSuite[sp] *= PredSuite[sp];
        //}
        
//    }    
//                v->RRR.rpar_HandleSelf[sp]   = HandleSelfWt;
//						v->RRR.rpar_ScrambleSelf[sp] = ScrambleSelfWt;
//					  v->RRR.rpar_HandleSwitch[sp] = PREYSWITCH;
    for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
		    prey = v->RRR.rpar_PreyFrom[links];
		    pred = v->RRR.rpar_PreyTo[links];

        // This is EwE Classic (c)(r)(tm)(4.0beta)(all rights reserved)
        // Master_Density =v->state_Ftime[pred] * pred_YY[pred];		    
				// KYA: Version with Ftime in YY
				//   Master_Density = pred_YY[pred];	 

        // This is the Combined
           //Master_Density =  MDULL     * pred_YY[pred] + 
           //                  MSCRAMBLE * PredSuite[prey] + 
           //                  MHANDLE   * HandleSuite[pred];

        // Here is the final equation
    		//   Q = v->RRR.rpar_QQ[links] * v->RRR.rpar_XX[links] * pred_YY[pred] * prey_YY[prey] / 
        //       (Master_Density + ((v->RRR.rpar_XX[links] - 1.0) * 
				//	                   (1.0 + v->state_Ftime[prey]) / 2.0) );

        // IWC approach to combined Handling and Predator density.  
				// Does not include Ftime!  

        // Multiplicative version which is too weird for me.
        // Q =   v->RRR.rpar_QQ[links] * 
        //     ( v->RRR.rpar_XX[links] * pred_YY[pred] / ((v->RRR.rpar_XX[links] - 1.0) + pred_YY[pred])     )*
        //     ( v->RRR.rpar_HH[links] * prey_YY[prey] / ((v->RRR.rpar_HH[links] - 1.0) + prey_YY[prey])     )*
				// 		( v->RRR.rpar_DD[links]                 / ((v->RRR.rpar_DD[links] - 1.0) + HandleSuite[pred]) )*
				// 		( v->RRR.rpar_VV[links]                 / ((v->RRR.rpar_VV[links] - 1.0) + PredSuite[prey])   );
        
				// Additive version    
        Q =   v->RRR.rpar_QQ[links] * pred_YY[pred] * pow(prey_YY[prey],v->RRR.rpar_HandleSwitch[links]) *
				    ( v->RRR.rpar_DD[links] / ( v->RRR.rpar_DD[links] - 1.0 + 
						                     pow(v->RRR.rpar_HandleSelf[pred] * prey_YY[prey]   + 
																 (1. - v->RRR.rpar_HandleSelf[pred]) * HandleSuite[pred],
                                                     v->RRR.rpar_HandleSwitch[links])) )*
            ( v->RRR.rpar_VV[links] / ( v->RRR.rpar_VV[links] - 1.0 + 
                                 v->RRR.rpar_ScrambleSelf[pred] * pred_YY[pred] + 
						                     (1. - v->RRR.rpar_ScrambleSelf[pred]) * PredSuite[prey]) );
            
				Q *= force_byprey[prey][y*STEPS_PER_YEAR+m]; 
				//if (isnan(Q)){
				//   cout << pred << " " << prey << " " << Q << " " << prey_YY[prey] <<
				//   " " << pred_YY[pred] << " " << v->RRR.rpar_QQ[links] << " " << v->RRR.rpar_XX[links] << endl;
				//}
        v->deriv_FoodLoss[prey]           += Q;
        v->deriv_FoodGain[pred]           += Q;
        v->deriv_UnAssimLoss[pred]        += Q * v->RRR.rpar_UnassimRespFrac[pred]; 
        v->deriv_ActiveRespLoss[pred]     += Q * v->RRR.rpar_ActiveRespFrac[pred];  
        
		    if (m==MEASURE_MONTH){ v->out_LinksM[links][y] = Q/v->state_BB[prey];}
		    v->deriv_ConsMat[prey][pred] += Q;
		    //if (prey==0){v->out_PP += Q;}
		}
		
    // Species-Specific Rates for living species only
    for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
        // deriv_PassiveRespLoss[sp] = 
        v->deriv_MzeroLoss[sp] = v->RRR.rpar_MzeroMort[sp] * v->state_BB[sp];
        v->deriv_ConsMat[sp][0] += v->deriv_MzeroLoss[sp];
    }	   
		   
		int gr,dest;
		double caught, totQ;		

		// This sets EFFORT by time series of gear-target combinations
		// if -1 is an input value, uses TERMINAL F (last non-negative F)
		
		   for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
		       if (y+m+d == 0){v->fish_Effort[gr]=1.0;}
		       else           {v->fish_Effort[gr]=1.0;} // NOTE DEFAULT!  THIS CAN BE CHANGED TO 1.0
       // Added 7/8/08 for forced effort
           if (FORCED_EFFORT[gr][y] > -0.001) 
					    {v->fish_Effort[gr]=FORCED_EFFORT[gr][y];}

			     if ((FORCED_TARGET[gr]>0) && (FORCED_CATCH[gr][y]>-EPSILON)){
			        totQ = 0.0;
			        sp   = FORCED_TARGET[gr];
			        for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){
    					    if ((v->RRR.rpar_FishingThrough[links] == gr) && 
		    					    (v->RRR.rpar_FishingFrom[links]) == sp){
				    					totQ += v->RRR.rpar_FishingQ[links];
									}
						  }
						  v->fish_Effort[gr] = FORCED_CATCH[gr][y]/ 
							                  (totQ * v->state_BB[sp]);
							if (FORCED_CATCH[gr][y] >= v->state_BB[sp])
							   {v->fish_Effort[gr] = (1.0-EPSILON)*(v->state_BB[sp])/ 
				    			                  (totQ * v->state_BB[sp]);}
					 }
					 // By putting F after catch, Frates override absolute catch
			     if ((FORCED_FTARGET[gr]>0) && (FORCED_FRATE[gr][y]>-EPSILON)){
			        totQ = 0.0;
			        sp   = FORCED_FTARGET[gr];
			        for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){
    					    if ((v->RRR.rpar_FishingThrough[links] == gr) && 
		    					    (v->RRR.rpar_FishingFrom[links]) == sp){
				    					totQ += v->RRR.rpar_FishingQ[links];
									}
						  }
						  v->fish_Effort[gr] = FORCED_FRATE[gr][y]/totQ;
							//if (FORCED_CATCH[gr][y] >= v->state_BB[sp])
							//   {v->fish_Effort[gr] = (1.0-EPSILON)*(v->state_BB[sp])/ 
				    //			                  (totQ * v->state_BB[sp]);}
					 }					 
					 
				   //if ((y==0) && (m==0) && (d==0)){
					 //    cout << path_species[gr] << " " << FORCED_TARGET[gr] << " " << path_species[sp] << " " 
					 //	      << v->state_BB[sp] << " " << FORCED_CATCH[gr][y] << " "
					 //		      << v->fish_Effort[gr] << endl;
					 //} 
					 
			 }
			 					 
    // Now apply effort to determine catch
    for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){
				 prey = v->RRR.rpar_FishingFrom[links];
				 gr   = v->RRR.rpar_FishingThrough[links];
				 dest = v->RRR.rpar_FishingTo[links];
				 caught = v->RRR.rpar_FishingQ[links] * v->fish_Effort[gr] * v->state_BB[prey]; 
				 if (m==MEASURE_MONTH){ v->out_LinksF[links][y] = caught; }
				 // if (caught>=v->state_BB[prey]){caught=(1.0-EPSILON)*(v->state_BB[prey]);}
         v->deriv_FishingLoss[prey] += caught;
         v->deriv_FishingThru[gr]   += caught;
         v->deriv_FishingGain[dest] += caught;
         v->deriv_ConsMat[prey][gr] += caught;
		}		

    // Special "CLEAN" fisheries assuming q=1
       for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
            //desired_catch = FORCED_CATCH[sp][t];
            //caught = v->fish_Effort[sp] * v->state_BB[sp];
            caught = FORCED_CATCH[sp][y] + FORCED_FRATE[sp][y] * v->state_BB[sp];
            if (caught <= -EPSILON)
						   {caught = v->TerminalF[sp] * v->state_BB[sp];}
            if (caught>=v->state_BB[sp]){caught=(1.0-EPSILON)*(v->state_BB[sp]);}
            v->deriv_FishingLoss[sp] += caught;
            v->deriv_FishingThru[0]  += caught;
            v->deriv_FishingGain[0]  += caught;
            v->TerminalF[sp] = caught/v->state_BB[sp];
       }
       
		// DETRITUS (Danger!  Seems to work, but I don't know why!)
		int det;
		float TotDetOut;
    for (sp=1; sp<=NUM_LIVING; sp++){
        TotDetOut = v->deriv_MzeroLoss[sp] + v->deriv_UnAssimLoss[sp];
        for (det=1; det<=NUM_DEAD; det++){
            v->deriv_DetritalGain[det+NUM_LIVING] += TotDetOut * path_DetFate[det][sp];
        }
    }
    for (sp=1; sp<=NUM_DEAD; sp++){
        TotDetOut = v->deriv_MzeroLoss[sp+NUM_LIVING]; // v->state_BB[sp+NUM_LIVING];
        v->deriv_DetritalLoss[sp+NUM_LIVING] += TotDetOut; // ????????????????
        v->deriv_MzeroLoss[sp+NUM_LIVING] = 0;
        for (det=1; det<=NUM_DEAD; det++){
            //v->deriv_DetritalLoss[sp] -= TotDetOut * path_DetFate[det][sp];  // ??????
            v->deriv_DetritalGain[det+NUM_LIVING] += TotDetOut * path_DetFate[det][sp+NUM_LIVING];
        }
    }    

    // Noise Addition
    for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
       v->deriv_FoodLoss[i]  *= v->force_bymort[i][y*STEPS_PER_YEAR+m];
       v->deriv_MzeroLoss[i] *= v->force_bymort[i][y*STEPS_PER_YEAR+m];
    }
    
    // FINAL ADDITION OF PARTS
    double noise;
    
    for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
        //if (m==6){v->out_TotF[i][y] = v->deriv_FishingLoss[i]/v->state_BB[i];}
        v->deriv_dyt[i]=v->deriv_DerivT[i];
        v->deriv_TotGain[i] = v->deriv_FoodGain[i]     + /*v->deriv_PrimProdGain[i] + */ 
                           v->deriv_DetritalGain[i] + v->deriv_FishingGain[i];
        
        v->deriv_LossPropToQ[i] = v->deriv_UnAssimLoss[i]  + v->deriv_ActiveRespLoss[i];
        v->deriv_LossPropToB[i] = v->deriv_FoodLoss[i]     + v->deriv_MzeroLoss[i] +
                               v->deriv_FishingLoss[i]  + v->deriv_DetritalLoss[i]; //+ 
                               //v->deriv_PassiveRespLoss[i]
        
				noise=0;
				//noise=(v->deriv_LossPropToB[i]) * v->month_noise[i];            
        v->deriv_TotLoss[i] = v->deriv_LossPropToQ[i] + v->deriv_LossPropToB[i];    
        v->deriv_DerivT[i]  = v->deriv_TotGain[i] - v->deriv_TotLoss[i] + noise; // + Noise
        // v->deriv_DerivT[i] = v->deriv_DerivT[i] * (1! - FREQ_LOCKED[i])
        
		// Set biomeq for "fast equilibrium" of fast variables
				if (v->state_BB[i] > 0) {
          v->deriv_biomeq[i] = v->deriv_TotGain[i] / 
					                 (v->deriv_TotLoss[i] / v->state_BB[i]);
        }
    //if ((y==0) && (m==0) && (d==0)){
    //cout << y   << "," <<  m << "," <<  d << "," << v->state_BB[i] << "," << v->RRR.rpar_B_BaseRef[i] << ",";
		//cout << i << "," << deriv_dyt[i]  << "," << deriv_DerivT[i] << "," 
		//           << deriv_TotGain[i] << "," << deriv_TotLoss[i]<< ",";
		//           
		//cout << i << "," << deriv_FoodGain[i]  << "," << deriv_DetritalGain[i] << "," 
		//           << deriv_FishingGain[i] << ",";		
		//           
		//cout << i << "," << deriv_UnAssimLoss[i]  << "," << deriv_ActiveRespLoss[i] << "," 
		//                 << deriv_FoodLoss[i] << "," << deriv_MzeroLoss[i] << "," << 
		//					          deriv_FishingLoss[i] << "," << deriv_DetritalLoss[i] << "," << endl;
		//}			        
		           
    }
    return 0;
}













