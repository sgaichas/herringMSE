
#include <gsl/gsl_histogram.h>
//int random_from_fitted
// -----------------------------------------------------------------------------


//------------------------------------------------------------------------------
int random_fromfitted_series(void)
{ 
   int r, rc, status, i;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];

   struct SimRun *Runs[PRUNS];

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      
    // CREATING a SERIES
    //printf ("test1\n");
        for (r=0; r<PRUNS; r++){
               Runs[r]=new_SimRun();
               //cout << "and 1" << endl;
               //for (i=0; i<=NDIM; i++){Runs[r]->fit_vector[i]=0.0;}
               if (r==0){load_instep(Runs[r]);}
               else {for (i=0; i<=NDIM; i++){
                         Runs[r]->fit_vector[i]=Runs[0]->fit_vector[i];}
               }
               //cout << "and 2" << endl;
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
               Runs[r]->thread_ID=rseed+r;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, save_a_series, (void *) &Runs[r]);
               //cout << "and 3" << endl;
                   printf ("Process %d ready, sir.\n",r);
               rc = pthread_create(&threads[r],&attr, save_a_series, (void *) Runs[r]);
               //cout << "and 4" << endl;
        }
      
        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }
      
        // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
                free_SimRun(Runs[r]);
           }         
   return 0;

}

// -----------------------------------------------------------------------------

// KYA Depreciated in favor of random_fromfitted Jume 3 2008
/*
int random_ecosystem_series(void)
{
   int r, rc, status;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];

   struct SimRun *Runs[PRUNS];

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      
    // CREATING a SERIES
        for (r=0; r<PRUNS; r++){
               Runs[r]=new_SimRun();
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
               Runs[r]->thread_ID=rseed+r;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, save_a_series, (void *) &Runs[r]);
               rc = pthread_create(&threads[r],&attr, save_a_series, (void *) Runs[r]);
        }
      
        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }
      
        // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
                free_SimRun(Runs[r]);
           }         
   return 0;

}
*/
// -----------------------------------------------------------------------------

int loaded_ecosystem_series(void)
{
   int r, rc, status;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];

   struct SimRun *Runs[PRUNS];

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      
    // CREATING a SERIES
        for (r=0; r<PRUNS; r++){
               Runs[r]=new_SimRun();
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
               Runs[r]->thread_ID=rseed+r;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, load_a_series, (void *) &Runs[r]);
               rc = pthread_create(&threads[r],&attr, load_a_series, (void *) Runs[r]);
        }
      
        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }
      
        // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
                free_SimRun(Runs[r]);
           }           
   return 0;

}

// -----------------------------------------------------------------------------
// added June 4 2008 by SKG to call load_and_perturb_series
int perturb_loaded_series(void)
{
   int r, rc, status;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];

   struct SimRun *Runs[PRUNS];

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      
    // CREATING a SERIES
        for (r=0; r<PRUNS; r++){
               Runs[r]=new_SimRun();
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
               Runs[r]->thread_ID=rseed+r;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, load_a_series, (void *) &Runs[r]);
               rc = pthread_create(&threads[r],&attr, load_and_perturb_series, (void *) Runs[r]);
        }
      
        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }
      
        // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
                free_SimRun(Runs[r]);
           }           
   return 0;

} 

// -----------------------------------------------------------------------------// -----------------------------------------------------------------------------

void *save_a_series(void *threadarg)
//void run_a_series(struct SimRun *v)
{
struct SimRun *v;
int i, k, sp, kept, big_cycle, moreflag;
time_t start, end;
char outFile[80];
float ssq;
const char* w = "w";
const char* ss = "%";
char c;
char logFile[80]; FILE *logptr;
char endFile[80]; FILE *endptr;
char binFile[80]; FILE *binptr;

 float sq[MAX_SYSTEMS];
   int sn[MAX_SYSTEMS];
float *sv[MAX_SYSTEMS];
double *newvec;

struct RatePar *Systems;
struct RatePar sys;

      v = (struct SimRun *) threadarg;
      Systems = (struct RatePar *)calloc(MAX_SYSTEMS,sizeof(struct RatePar)); 
          if (Systems==NULL){
             cout << "allocation failure" << endl;
             pthread_exit(NULL);
          }  

      start = time(NULL);
      sprintf(logFile,"%s_Seed%d_log.csv",OUTFOLDER,v->thread_ID);
      logptr=fopen(logFile,w);

      sprintf(endFile,"%s_Seed%d_end.csv",OUTFOLDER,v->thread_ID);
      endptr=fopen(endFile,w);

      fprintf(endptr,"Seed,System,Year,ssq,");
      for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){     
           fprintf(endptr,"B_%s,M_%s,",path_species[sp],path_species[sp]);
      }
      fprintf(endptr,"\n");
           
      kept=0;
      big_cycle=0;

      newvec  = dvector(0,NDIM);
      for (k=0; k<MAX_SYSTEMS; k++){
          sv[k]  = vector(0,NDIM);
      }

      sprintf(binFile,"%s_Seed%d.bin",OUTFOLDER,v->thread_ID);
      binptr=fopen(binFile,"wb");
      
      for (i=0; i<=SYSTEMS_TO_RUN; i++){
          // zeroeth run is baseline, every other one is
					// from a generated a random system
					   if (i<1){path_to_rates(v);  
                      apply_vector_to_rates(v,v->fit_vector);
             //         cout << "PreyLinks=" << v->RRR.rpar_NumPredPreyLinks << endl;
                     }   // Sets up baseline run
                else {random_system(v,newvec);}   // Sets up random system
             initialize_stanzas(v);
      		// Run Adams-Basforth
					   //Adams_Basforth(v, 0, fit_Years);
             Adams_Basforth(v, 0, BURN_TIME);
          // If it didn't crash, calculate SSQ and save the system
           	 ssq = 0.0;
           	            	 
             if (!v->DISCARD_YEAR){
                 ssq = SSQ(v);
                 //if (kept>MAX_SYSTEMS-1){kept=MAX_SYSTEMS-1;}
                 sq[kept] = ssq;
                 sn[kept] = i;
                 memcpy(sv[kept], newvec,sizeof(float)*(NDIM+1));
                 save_system(&Systems[kept], &v->RRR);
                 kept++;		    
             }
             
             if ((kept==MAX_SYSTEMS) || (i==SYSTEMS_TO_RUN)){
								cout << v->thread_ID << " got max or is done!" << endl;
                fwrite(&kept,sizeof(int),1,binptr);       // Number of systems following
                for (k=0; k<kept; k++){
                    fwrite(ss,1,sizeof(char),binptr);     // Separator
                    fwrite(&sn[k],1,sizeof(int),binptr);
                    fwrite(&sq[k],1,sizeof(float),binptr);
                    fwrite( sv[k],1,sizeof(float)*(NDIM+1),binptr);
                    fwrite(&Systems[k],1,sizeof(struct RatePar),binptr);
                }
                if (i==SYSTEMS_TO_RUN){moreflag=0;} else {moreflag=1;}
                fwrite(&moreflag,1,sizeof(int),binptr);
								big_cycle++;  
								kept=0; 
						 }
						 
          // Logging the run
             fprintf (logptr, "%d,",i);
				     if (v->DISCARD_YEAR){ 
                fprintf (logptr, "discarded,0,%d,",v->DISCARD_YEAR); // cout << i << " : discarded ";
	   				      for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
				  		        if (v->discarded[sp]){fprintf (logptr,"%s,",path_species[sp]); //cout << path_species[sp] << ",";
                  }}
						    //cout << " in year " << v->DISCARD_YEAR << endl;
				     }
					   else {fprintf (logptr, "kept,%g,0,",ssq);// cout << i << " : kept SSQ = " << ssq << endl;
					   }
					   fprintf (logptr, "\n"); 

         //  Saving the end state
             if (!v->DISCARD_YEAR){
                //output_run(v,outFile);
//                 fprintf(endptr,"%d,%d,%d,%g,",v->thread_ID,i,fit_Years+fit_StartYear-1,ssq);
                 fprintf(endptr,"%d,%d,%d,%g,",v->thread_ID,i,BURN_TIME+fit_StartYear-1,ssq);
                 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){     
//                     fprintf(endptr,"%g,%g,",v->out_BB[sp][fit_Years-1],v->out_MM[sp][fit_Years-1]);
                      fprintf(endptr,"%g,%g,",v->out_BB[sp][BURN_TIME-2],v->out_MM[sp][BURN_TIME-2]);
                 }
                fprintf(endptr,"\n");             
             }
             else{
                 fprintf(endptr,"%d,%d,%d,%g,",v->thread_ID,i,v->DISCARD_YEAR+fit_StartYear,ssq);
                 for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){     
                     fprintf(endptr,"%g,%g,",v->out_BB[sp][v->DISCARD_YEAR],v->out_MM[sp][v->DISCARD_YEAR]);
                 }
                fprintf(endptr,"\n");             
             }

					// Every few hundred iterations, print a status report
             float remaining,passed;
            
             if ((i%200 == 10)){
                 end=time(NULL);
                 passed = (float)(end-start);
                 remaining = (float)(SYSTEMS_TO_RUN-i);
                 cout <<"seed: "<< v->thread_ID <<" run: "<< i << " kept: "<< kept + MAX_SYSTEMS*big_cycle
                      <<" in " << passed/60.0 << " min. "<< (remaining * passed/(float)i)/60.0
                      <<" min. left est. saved: " << (int)((float)SYSTEMS_TO_RUN*(float)(kept + MAX_SYSTEMS*big_cycle)/(float)i) 
                      << " " << endl;
             }
          
     }
     fclose(logptr);
     fclose(endptr);
     
     // Now we save the good ecosystems in a binary file
        //cout << "FINAL seed: " << v->thread_ID << " run: " << i << " kept: " << kept + MAX_SYSTEMS*big_cycle << " " << endl;

        fclose(binptr);

// Finally, clean up memory and exit thread
   for (k=0; k<MAX_SYSTEMS; k++){
        free_vector(sv[k],0,NDIM);
   }
   free_dvector(newvec,0,NDIM);
   free(Systems);
   pthread_exit(NULL);      

}

// ---------------------------------------------------------------------
int multiseries_stats()
{
int r,thread_ID;
char in_vec[80];
int cols, rows, iout,sp,t,k,i;
int series_index[MAX_COLS+1], series_type[MAX_COLS+1];
string sLine, ins;
double flout;

        for (r=0; r<PRUNS; r++){
                thread_ID=rseed+r;
             // Open this file and declare parser that breaks lines into fields
                sprintf(in_vec,"%s_fitVectorsSeed%d.csv",OUTFOLDER,thread_ID);
                ifstream infile(in_vec);
                if (!infile){cout<<"No file loaded: "<<in_vec<<endl;return 1;}
		            else {cout << "loaded from " << in_vec << endl;}

                CSVParser parser;

              // read the first line of the file in and pass it to the parser
              // then count the number of columns in the top line
                 getline(infile,sLine);  parser << sLine;
                 cols=0;  ins="in";  
								 while (ins != ""){parser >> ins; if (ins == "") continue; cols++;}
		             
              // read the rest of the lines and put in the forcing    
                 t=0;
                 while (!infile.eof()) {
                       t++;
                       getline(infile,sLine); if (sLine == "") continue; 
		                   parser << sLine;    
		                   
                       parser >> k;      // Ecosystem number
                       parser >> flout;  // SSQ
                       parser >> sp;     // Species Number
                       parser >> ins;    // Species name
                       for (i=0; i<NUM_VARS; i++){
		                       parser >> flout;

		                   }
                 }
                 // Close the file
                    infile.close();
             }         
}


//------------------------------------------------------------------------------
void *load_a_series(void *threadarg)
{
struct SimRun *v;
int i, sp, kept,link, moreflag, k,links,yy;
time_t start, end;
char outFile[80]; FILE *fptr;
float ssq;
const char* w = "w";
const char* ss = "%";
char c;
char logFile[80]; FILE *logptr;
char endFile[80]; FILE *endptr;
char binFile[80]; FILE *binptr;
char catFile[80]; FILE *catptr;
//char consFile[80]; FILE *consptr;  //added Oct 2011 by SKG for consumption outputs

FILE *dmpptr[NUM_LIVING+1];
 float sq[MAX_SYSTEMS];
   int sn[MAX_SYSTEMS];
float *sv[MAX_SYSTEMS];
float *newvec;
struct RatePar *Systems;
struct RatePar sys;

//struct gsl_histogram *params[NUM_LIVING+1][NUM_VARS];

      // Translate pointer for threads
         v = (struct SimRun *) threadarg;

     // Allocate space to store ecosystems when loaded
        Systems = (struct RatePar *)calloc(MAX_SYSTEMS,sizeof(struct RatePar)); 
        if (Systems==NULL){cout<<"allocation failure"<<endl; pthread_exit(NULL);}  

       newvec  = vector(0,NDIM);
       for (k=0; k<MAX_SYSTEMS; k++){
          sv[k]  = vector(0,NDIM);
       }
    // Initialize Histograms
       //for (sp=1; sp<=NUM_LIVING; sp++){
       //     for (i=0; i<NUM_VARS; i++){
			//			    params[sp][i] = gsl_histogram_alloc(20);
			//			    gsl_histogram_set_ranges_uniform(params[sp][i],-6.0,6.0);
			//			}
			 //}

    // Load all ecosystems
       cout << "Loading ecosystems..." << endl;
       // THIS BINARY FORMAT MUST MATCH INPUT
          sprintf(binFile,"%s_Seed%d.bin",INFOLDER,v->thread_ID);
          binptr=fopen(binFile,"rb");

       // Vector output file
          sprintf(outFile,"%s_fitVectorsSeed%d.csv",OUTFOLDER,v->thread_ID);
          fptr=fopen(outFile,"w");
          fprintf(fptr,"System,SSQ,Num,Species,");
          for (i=0; i<NUM_VARS; i++){fprintf (fptr,"p%d,",i);}
          fprintf (fptr,"\n");

       // State output file
          sprintf(endFile,"%s_Seed%d_end.csv",OUTFOLDER,v->thread_ID);
          endptr=fopen(endFile,w);
          fprintf(endptr,"Seed,System,Species,Year,");
          fprintf(endptr,"outBB,outCC,outMM,outM2,outPB,outQB,outTotF,");
          fprintf(endptr,"\n");
    
      //  Detailed Catch file
          sprintf(catFile,"%s_Seed%d_catch.csv",OUTFOLDER,v->thread_ID);
          catptr=fopen(catFile,w);
          fprintf(catptr,"Seed,System,Species,Gear,Destination,");
          //fprintf(catptr,"Catch,");
          for (yy=0; yy<=fit_Years; yy++){fprintf(catptr,"Y%d,",yy+fit_StartYear);}
					fprintf(catptr,"\n"); 
/*
      //  Detailed Consumption file added October 2011 by SKG
          sprintf(consFile,"%s_Seed%d_cons.csv",OUTFOLDER,v->thread_ID);
          consptr=fopen(consFile,w);
          fprintf(consptr,"Seed,System,Eaten,EatenBy,Year,Cons");
          //fprintf(consptr,"Seed,System,Eaten,EatenBy,");
          //for (yy=0; yy<=fit_Years; yy++){fprintf(consptr,"Y%d,",yy+fit_StartYear);}
		  fprintf(consptr,"\n"); 
*/					          					          
          moreflag=1;
          start = time(NULL);
          while (moreflag){
                // Get number kept in a batch
   								 fread(&kept,sizeof(int),1,binptr);
								// Read in that many systems      
                   for (k=0; k<kept; k++){
                     fread(&c,1,sizeof(char),binptr); if (c!='%'){cout << "alignment error" << endl;}
                     fread(&sn[k],1,sizeof(int),binptr);
                     fread(&sq[k],1,sizeof(float),binptr);
                     fread( sv[k],1,sizeof(float)*(NDIM+1),binptr);
                     fread(&Systems[k],1,sizeof(struct RatePar),binptr);
                   }
                   cout << kept << " ecosystems loaded." << endl; 

                // See if there's another batch to run
                   fread(&moreflag,1,sizeof(int),binptr);
                
                // OUTPUT VECTOR INFO
                   for (k=0; k<kept; k++){                                        
                     for (sp=1; sp<=NUM_LIVING; sp++){
                         fprintf(fptr,"%d,%g,%d,%s,",sn[k],sq[k],sp,path_species[sp]);
                         for (i=0; i<NUM_VARS; i++){
												     fprintf(fptr,"%g,",sv[k][i*NUM_LIVING+sp]);
                             //gsl_histogram_increment(params[sp][i],(double)(sv[k][i*NUM_LIVING+sp]));
                         }
                         fprintf (fptr,"\n");
                     }     
								}                 

                   for (k=0; k<kept; k++){    
                    // Load system rate parameters into the SimRun variable          
                       load_system(v, &Systems[k]);
                    // Run the ecosystem and calculate the SSQ
                       Adams_Basforth(v, 0, fit_Years);
                       //ssq = SSQ(v);
                    // Outputs

           
                       for (yy=0; yy<=fit_Years; yy++){
                            for (sp=1; sp<=NUM_LIVING; sp++){
                                fprintf(endptr,"%d,%d,%d,%d,",v->thread_ID,sn[k],sp,
                                                          yy+fit_StartYear);
                                fprintf(endptr,"%g,",v->out_BB[sp][yy]);
                                fprintf(endptr,"%g,",v->out_CC[sp][yy]);
                                fprintf(endptr,"%g,",v->out_MM[sp][yy]);
                                fprintf(endptr,"%g,",v->out_M2[sp][yy]);
                                fprintf(endptr,"%g,",v->out_PB[sp][yy]);
                                fprintf(endptr,"%g,",v->out_QB[sp][yy]);
                                fprintf(endptr,"%g,",v->out_TotF[sp][yy]);                                
                                fprintf(endptr,"\n");
                             }
                         }
                         
                         
                           for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){

                                fprintf(catptr,"%d,%d,%d,%d,%d,",v->thread_ID,sn[k],
                                                              v->RRR.rpar_FishingFrom[links],
                                                              v->RRR.rpar_FishingThrough[links],
                                                              v->RRR.rpar_FishingTo[links]);
                             for (yy=0; yy<=fit_Years; yy++){fprintf(catptr,"%g,",v->out_LinksF[links][yy]);}				                       
													   fprintf(catptr,"\n");
													 }
/* 
						//added by SKG October 2011 for consumption outputs
                          for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){

								for (yy=0; yy<=fit_Years; yy++){
									fprintf(consptr,"%d,%d,%d,%d,%d,",v->thread_ID,sn[k],
                                                              v->RRR.rpar_PreyFrom[links],
                                                              v->RRR.rpar_PreyTo[links],
															  yy+fit_StartYear);	
								    fprintf(consptr,"%g,",v->out_LinksM[links][yy] * v->out_BB[v->RRR.rpar_PreyFrom[links]][yy]);				                       
								    fprintf(consptr,"\n");
								    }
								}

													 //     prey = v->RRR.rpar_FishingFrom[links];
				                   //     gr   = v->RRR.rpar_FishingThrough[links];
				                   //     dest = v->RRR.rpar_FishingTo[links];
				                   //caught = v->RRR.rpar_FishingQ[links] * v->fish_Effort[gr] * v->state_BB[prey]; 
				                   //}
				                   
				                   //}if (m==MEASURE_MONTH){v->out_LinksF[links][y] = caught;}
*/				 
				   
                       float remaining,passed;
                       if ((k%20 == 10)){
													   end=time(NULL);
                             passed = (float)(end-start);
                             remaining = (float)(kept-k);
													
													   cout <<"seed: "<< v->thread_ID <<" " << k << 
													   " of "<< kept << "in batch, " << passed/60.0 << " min. "
													   << (remaining * passed/(float)i)/60.0 <<" left in batch" << endl;
											 }
                   }

          }  // end of moreflag
          fclose(endptr);
          fclose(fptr); 
          fclose(binptr);  
		  fclose(catptr);  
		  //fclose(consptr);  

     //  for (sp=1; sp<=NUM_LIVING; sp++){
     //       for (i=0; i<NUM_VARS; i++){
		 //			     gsl_histogram_free(params[sp][i]);
		 //			}
		 //	 }

      for (k=0; k<MAX_SYSTEMS; k++){
           free_vector(sv[k],0,NDIM);
      }
      free_vector(newvec,0,NDIM);
      free(Systems);
      pthread_exit(NULL);
      
/* 
    // Open the output files
//       path_to_rates(v);  // needed to set number of links to make titles (silly)
        // not anymore, and it resets M0 so you dont get the perturbed parameter out
       for (sp=1; sp<=NUM_LIVING; sp++){
           sprintf(endFile,"%s_Seed%d_species%d.csv",OUTFOLDER,v->thread_ID,sp);
           dmpptr[sp]=fopen(endFile,w);
           fprintf(dmpptr[sp],"Seed,System,Species,Year,ssqOrig,ssqCheck,");
           fprintf(dmpptr[sp],"PB,QB,F,M2,Mzero,Mtot,B,");
//           for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
//               if (v->RRR.rpar_PreyFrom[link] == sp){
//                   fprintf(dmpptr[sp],"Mby_%s,",path_species[v->RRR.rpar_PreyTo[link]]);
//               }
//           }
//           for (link=1; link<=NUM_LIVING+NUM_DEAD; link++){     
//               fprintf(dmpptr[sp],"B_%s,",path_species[link]);
//           }
//           fprintf(dmpptr[sp],"SystemPP,");
//                 for (link=1 ; link<=juv_N; link++){
//                      fprintf(dmpptr[sp],"RS_%s,",path_species[juv_AduNum[link]]);
//                 }
           fprintf(dmpptr[sp],"\n");
        }
     // Cycle through all ecosystems and run them, producing output
        start = time(NULL);
        for (i=0; i<kept; i++){    
          // Load system rate parameters into the SimRun variable          
             load_system(v, &Systems[i]);
          // Run the ecosystem and calculate the SSQ
             Adams_Basforth(v, 0, BURN_TIME);
             ssq = SSQ(v);
          // Outputs
             int yy;
             for (yy=0; yy<=BURN_TIME; yy++){
             for (sp=1; sp<=NUM_LIVING; sp++){
                 //fprintf(dmpptr[sp],"Seed,System,Species,Year,ssqOrig,ssqCheck,");
                   fprintf(dmpptr[sp],"%d,%d,%s,%d,%g,%g,",v->thread_ID,sn[i],path_species[sp],
                                                         yy+fit_StartYear-1,sq[i],ssq);
                 //fprintf(dmpptr[sp],"F,Mtot,Mzero,QB,PB");
                   fprintf(dmpptr[sp],"%g,%g,%g,%g,%g,%g,%g,",
                                                  v->out_PB[sp][yy],
                                                  v->out_QB[sp][yy],
                                                  v->out_TotF[sp][yy],
                                                  v->out_M2[sp][yy],
                                                  v->RRR.rpar_MzeroMort[sp],
                                                  v->out_MM[sp][yy],
                                                  v->out_BB[sp][yy]);

//                   fprintf(dmpptr[sp],"%g,%g,%g,%g,%g,",v->out_TotF[sp],
//                                                  v->out_MM[sp][yy],
//                                                  v->RRR.rpar_MzeroMort[sp],
//                                                  v->out_QB[sp],
//                                                  v->out_PB[sp]);
//                 for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
//                     if (v->RRR.rpar_PreyFrom[link] == sp){
//                         fprintf(dmpptr[sp],"%g,",v->out_LinksM[link]);
//                     }
//                 }
//                 for (link=1; link<=NUM_LIVING+NUM_DEAD; link++){     
//                     fprintf(dmpptr[sp],"%g,",v->out_BB[link][yy]);
//                 }
//                //fprintf(dmpptr[sp],"SystemPP,");
//                 fprintf(dmpptr[sp],"%g,",v->out_PP);
//                 for (link=1 ; link<=juv_N; link++){
//                      fprintf(dmpptr[sp],"%g,",v->RRR.rpar_SpawnEnergy[link]);
//                 }
                 fprintf(dmpptr[sp],"\n");
              }
              }
             //fprintf(endptr,"%d,%d,%d,%g,%d,%g,",v->thread_ID,i,sn[i],sq[i],fit_Years+fit_StartYear-1,ssq);
             //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){     
             //    fprintf(endptr,"%g,%g,",v->out_BB[sp][fit_Years-1],v->out_MM[sp][fit_Years-1]);
             //}
             //fprintf(endptr,"\n");  

					// Every few hundred iterations, print a status report
             float remaining,passed;
             if ((i%200 == 10)){
                 end=time(NULL);
                 passed = (float)(end-start);
                 remaining = (float)(kept-i);
                 cout <<"seed: "<< v->thread_ID <<" loaded: "<< i << " of "<< kept
                      <<" in " << passed/60.0 << " min. "<< (remaining * passed/(float)i)/60.0
                      <<" min. left " << endl;
             }                   
        }
  
   // Finally, clean up memory and exit thread
      for (sp=1; sp<=NUM_LIVING; sp++){
          fclose(dmpptr[sp]);
      }
*/    			
			            
}           

// -----------------------------------------------------------------------------
// added June 4 2008 by SKG to recreate Sense outputs with files by actor
// last-decadal averages of B (or other outputs) in base and perturbed systems
// and difference between base and perturbed for plotting
// note, not keeping full outputs to reduce file handling and runtime

void *load_and_perturb_series(void *threadarg)
{
struct SimRun *v;
int i, sp, spn, kept,link, moreflag, k, t;
time_t start, end;
char outFile[80]; FILE *fptr;
float ssq;
const char* w = "w";
const char* ss = "%";
char c;
char logFile[80]; FILE *logptr;
char endFile[80]; FILE *dmpptr[NUM_LIVING+1]; FILE *endptr;
char binFile[80]; FILE *binptr;


 float sq[MAX_SYSTEMS];
   int sn[MAX_SYSTEMS];
float *sv[MAX_SYSTEMS];
struct RatePar *Systems;
struct RatePar sys;

float baseavg[NUM_LIVING][NUM_LIVING+NUM_DEAD];
float experavg[NUM_LIVING][NUM_LIVING+NUM_DEAD];
//float basemort[NUM_LIVING+NUM_DEAD];

      // Translate pointer for threads
         v = (struct SimRun *) threadarg;
         //perturbed=new_SimRun();
         //perturbed->thread_ID = v->thread_ID;
     // Allocate space to store ecosystems when loaded
        Systems = (struct RatePar *)calloc(MAX_SYSTEMS,sizeof(struct RatePar)); 
        if (Systems==NULL){cout<<"allocation failure"<<endl; pthread_exit(NULL);}  

       for (k=0; k<MAX_SYSTEMS; k++){
          sv[k]  = vector(0,NDIM);
       }

    // Load all ecosystems
       cout << "Loading ecosystems..." << endl;
       // THIS BINARY FORMAT MUST MATCH INPUT
          sprintf(binFile,"%s_Seed%d.bin",INFOLDER,v->thread_ID);
          binptr=fopen(binFile,"rb");	          		//cout << binFile << " binfile open" << endl;

    moreflag=1;
    start = time(NULL);
    while (moreflag){
                // Get number kept in a batch
   				fread(&kept,sizeof(int),1,binptr);		//cout << kept << " kept" << endl;
					// Read in that many systems      
                   for (k=0; k<kept; k++){				//cout << k << " k" << endl;
                     fread(&c,1,sizeof(char),binptr); if (c!='%'){cout << "alignment error" << endl;}
                     fread(&sn[k],1,sizeof(int),binptr);
                     fread(&sq[k],1,sizeof(float),binptr);
                     fread( sv[k],1,sizeof(float)*(NDIM+1),binptr);
                     fread(&Systems[k],1,sizeof(struct RatePar),binptr);
                   }
                   cout << kept << " ecosystems loaded." << endl; 

                // See if there's another batch to run
                   fread(&moreflag,1,sizeof(int),binptr);

     // Open the output files by ACTOR (perturbed species)
       for (sp=1; sp<=NUM_LIVING; sp++){
           sprintf(endFile,"%s_Seed%d_actor%d.csv",OUTFOLDER,v->thread_ID,sp);
           dmpptr[sp]=fopen(endFile,w);
           fprintf(dmpptr[sp],"Seed,System,Experiment,Actor,LastDecade,");
           for (spn=1; spn<=NUM_LIVING+NUM_DEAD; spn++){     
              //comment out all but one of below for alternate output
//              fprintf(dmpptr[sp],"avB_%s,",path_species[spn]);
              fprintf(dmpptr[sp],"avP_%s,",path_species[spn]);
//              fprintf(dmpptr[sp],"avM_%s,",path_species[spn]);
           }
           fprintf(dmpptr[sp],"\n");
        }

     // Cycle through all ecosystems and run them, producing output
        //start = time(NULL);
        for (i=0; i<kept; i++){    
          // Load system rate parameters into the SimRun variable          
             load_system(v, &Systems[i]);
             
          // Run the base ecosystem 
             Adams_Basforth(v, 0, BURN_TIME);
             //ssq = SSQ(v);

             for (sp=1; sp<=NUM_LIVING; sp++){              
                 //fprintf(dmpptr[sp],"Seed,System,Experiment,Actor,LastDecade,");
                   fprintf(dmpptr[sp],"%d,%d,Base,%s,%d,",v->thread_ID,sn[i], path_species[sp],
                                                         (BURN_TIME/10));
                // Outputs               
                   for (spn=1; spn<=NUM_LIVING+NUM_DEAD; spn++){  
                       baseavg[sp][spn]=0.0;   
                       int yy;
                       for (yy=0; yy<BURN_TIME; yy++){
                            //calculate average P for last 10 years--for B replace with commented line  
                            if(yy>=(BURN_TIME-10)){
                                baseavg[sp][spn] += (0.1 * v->out_BB[spn][yy] * v->out_PB[spn][yy]); //baseavg[sp][spn] += (0.1 * v->out_BB[spn][yy]);
                            //cout << "baseavg " << spn << " " << baseavg[spn] << " outBB " << v->out_BB[spn][yy] << endl;
                            }                         
                       }
                       //fprintf(dmpptr[sp],"avB_%s,",path_species[spn]);
                       fprintf(dmpptr[sp],"%g,", baseavg[sp][spn]);
                  }
                  fprintf(dmpptr[sp],"\n");
             }     

                // Run the perturbed ecosystem for each species
             for (sp=1; sp<=NUM_LIVING; sp++){
						     load_system(v, &Systems[i]);              
                 //fprintf(dmpptr[sp],"Seed,System,Experiment,Actor,LastDecade,");
                   fprintf(dmpptr[sp],"%d,%d,M_down10p,%s,%d,",v->thread_ID,sn[i], path_species[sp],
                                                         (BURN_TIME/10));
                //perturbation: increase M of sp by 10%
                   //basemort[sp] = 0.0;
                   //basemort[sp] = v->RRR.rpar_MzeroMort[sp];
                   //v->RRR.rpar_MzeroMort[sp] += (0.1 * v->RRR.rpar_MzeroMort[sp]);          
                
                //SKG changed perturbation to 10% decrease in M = 10% increase in P 
                
                   for (t=0; t<=BURN_TIME*STEPS_PER_YEAR; t++){
                       v->force_bymort[sp][t] = 0.9;
                       if (sp>1){v->force_bymort[sp-1][t] = 1.0;}
                   }
                   Adams_Basforth(v, 0, BURN_TIME);
                   //ssq = SSQ(v);
                // Outputs
                   for (spn=1; spn<=NUM_LIVING+NUM_DEAD; spn++){     
                       experavg[sp][spn]=0.0;
                       int yy;
                       for (yy=0; yy<BURN_TIME; yy++){
                            //calculate average P for last 10 years--for B replace with commented line   
                            if(yy>=BURN_TIME-10){
                                experavg[sp][spn] += (0.1 * v->out_BB[spn][yy] * v->out_PB[spn][yy]); //experavg[sp][spn] += (0.1 * v->out_BB[spn][yy]);
                            }                         
                       }
                       //fprintf(dmpptr[sp],"avB_%s,",path_species[spn]);
                       fprintf(dmpptr[sp],"%g,", experavg[sp][spn]);
                  }
                  //v->RRR.rpar_MzeroMort[sp] = basemort[sp];
                  fprintf(dmpptr[sp],"\n");
            }

                // calculate difference between perturbed and base
             for (sp=1; sp<=NUM_LIVING; sp++){              
                 //fprintf(dmpptr[sp],"Seed,System,Experiment,Actor,LastDecade,");
                   fprintf(dmpptr[sp],"%d,%d,Diff,%s,%d,",v->thread_ID,sn[i], path_species[sp],
                                                         (BURN_TIME/10));
                // Outputs for average difference calculation relative to base                 
                   for (spn=1; spn<=NUM_LIVING+NUM_DEAD; spn++){     
                       //fprintf(dmpptr[sp],"avB_%s,",path_species[spn]);
                       fprintf(dmpptr[sp],"%g,", ((experavg[sp][spn]-baseavg[sp][spn])/baseavg[sp][spn]));
                  }
                  fprintf(dmpptr[sp],"\n");

              } // end species loop

					// Every few hundred iterations, print a status report
             float remaining,passed;
             if ((i%20 == 10)){
                 end=time(NULL);
                 passed = (float)(end-start);
                 remaining = (float)(kept-i);
                 cout <<"seed: "<< v->thread_ID <<" loaded: "<< i << " of "<< kept
                      <<" in " << passed/60.0 << " min. "<< (remaining * passed/(float)i)/60.0
                      <<" min. left " << endl;
             }     
                           
        }// end systems loop
  
     }  // end of moreflag
          
     fclose(binptr);  
	 for (sp=1; sp<=NUM_LIVING; sp++){
		 fclose(dmpptr[sp]);
		 }            

// The below loop causes a sigfault; although sv[k] needs to be specified, I don't use it so it can't be freed?
//      for (k=0; k<MAX_SYSTEMS; k++){
//           free_vector(sv[k],0,NDIM);
//      }
      free(Systems);
      pthread_exit(NULL);

}           


// -----------------------------------------------------------------------------
void load_system(struct SimRun *v, struct RatePar *r){

int sp;
     path_to_rates(v);

		 memcpy(&v->RRR, r, sizeof(struct RatePar) );



     // Set Start Biomass and other state variables
	      for (sp=0; sp<=NUM_GROUPS; sp++){

           //if (path_PB[sp]/(1.0 - v->RRR.rpar_ActiveRespFrac[sp] - v->RRR.rpar_UnassimRespFrac[sp]) >
					 //   2 * STEPS_PER_YEAR * STEPS_PER_MONTH)   
					 //        {v->RRR.rpar_NoIntegrate[sp]= 0;}
					 //else    {v->RRR.rpar_NoIntegrate[sp]=sp;} 

			     v->state_BB[sp]    = v->RRR.rpar_B_BaseRef[sp];
			     v->state_Ftime[sp] = 1.0;
			     v->discarded[sp]   = 0;
			  }		 
		 
     initialize_stanzas(v);

}

// -----------------------------------------------------------------------------

void save_system(struct RatePar *s, struct RatePar *r)
{
     memcpy(s, r, sizeof(struct RatePar) );
}

//  -----------------------------------------------------------------------------
int random_system(struct SimRun *v, double *newvec){
   int sp,gr,links,LL,goodP,i;
   float catchVary, AA, VV;
   float x;
   double GAMMA_ALPHA;
   float QB[NUM_GROUPS+1];
   float PB[NUM_GROUPS+1];
   double DCtot[NUM_GROUPS+1];
   double DC[NUM_GROUPS*NUM_GROUPS+1];
   float XX[NUM_GROUPS*NUM_GROUPS+1]; // Single
    // Implements the Ecosense routines to make a random ecosystem
    // from pedigree inputs.  Does not filter ecosystem at all.
    // NOTE:  THIS IS DIFFERENCE FROM SENSE, Some of the things I did
    //        IN SENSE MAY BE DISTRIBUTIONALLY INCORRECT.  REVIEW NEEDED.
    
    // First make sure we're starting with the Path system
       path_to_rates(v);
    
    // Load any starting vectors
    // DON'T NEED THIS AS IT'S APPLIED THROUGHOUT DOWN BELOW!
       //apply_vector_to_rates(v,v->fit_vector);
       
    // Loop through groups
       for (sp=1; sp<=NUM_LIVING; sp++){
			      v->RRR.rpar_B_BaseRef[sp] *= ((pedigree[1][sp]*(2.0*uniform(&v->rng)-1.0)*0.3)+1.0); // B pedigree
			        x = v->fit_vector[NUM_LIVING*1 + sp];
  			      v->RRR.rpar_MzeroMort[sp] *= (( 2.0 * exp(x) /(1+exp(x))));
              v->RRR.rpar_MzeroMort[sp] *= ((pedigree[2][sp]*(2.0*uniform(&v->rng)-1.0)*0.3)+1.0); // PB pedigree
//			        if (sp==33){cout << pedigree[2][sp] << " " << v->RRR.rpar_MzeroMort[sp] << endl;}
            goodP=1;
			      while (goodP){
			       // First pick PB
                x = v->fit_vector[NUM_LIVING*0 + sp];
						    PB[sp] = path_PB[sp] * (( 2.0 * exp(x) /(1+exp(x))));
                PB[sp] *= ((pedigree[2][sp]*(2.0*uniform(&v->rng)-1.0)*0.3)+1.0); // PB pedigree
//                if (sp==33){cout << pedigree[2][sp] << " " << PB[sp] << endl;}                
             // Then pick QB if QB > 0
						    if (path_QB[sp] > EPSILON ){
    			         QB[sp] = path_QB[sp];
                   QB[sp] *= ((pedigree[3][sp]*(2.0*uniform(&v->rng)-1.0)*0.3)+1.0); // QB pedigree
    			         if (PB[sp]/QB[sp] > (0.99 - path_GS[sp]))
    			              {goodP++; if (goodP>100){
    			                            //cout << goodP << " " << PB[sp] << " " << QB[sp] << endl;
    			                            goodP=0;
												              PB[sp]=path_PB[sp];
												              QB[sp]=path_QB[sp];    			              
												              v->RRR.rpar_UnassimRespFrac[sp] = path_GS[sp];
	  			                            v->RRR.rpar_ActiveRespFrac[sp]  = 1.0 - (PB[sp]/QB[sp])
							                                - v->RRR.rpar_UnassimRespFrac[sp];												              
												           }
												
												}
    			         else {goodP=0; 
									       v->RRR.rpar_UnassimRespFrac[sp] = path_GS[sp];
	  			                v->RRR.rpar_ActiveRespFrac[sp]  = 1.0 - (PB[sp]/QB[sp])
							                                - v->RRR.rpar_UnassimRespFrac[sp];	
										}
    			      }
    			   // Primary Production case (QB=0)   
    			      else {goodP=0; 
    			            QB[sp] = PB[sp];
								      v->RRR.rpar_UnassimRespFrac[sp] = path_GS[sp];
								      v->RRR.rpar_ActiveRespFrac[sp]  = 0;
								}
    			  }
    			  v->RRR.rpar_FtimeQBOpt[sp] = QB[sp];
    			  v->RRR.rpar_PBopt[sp] = PB[sp];
    			  //cout << goodP << " " << PB[sp] << " " << QB[sp] << endl;
			      //while (respVary)
			 }

// Gsmma distributed diets
// Keeping the number of predator/prey links the same
// and setting diets that are drawn too low to a minimum value.

// GAMMA_ALPHA is the variance setting for diets.
// Alternative is to set separately for each predator based on pedigree
// 0.1 is high variance, and as much as numerics can handle without
// making too many 0's.  100 gives a fairly low variance with a normal
// distribution.

//#define GAMMA_ALPHA 50.0
    // SET DCtot FOR EACH GROUP TO 0
       LL=(NUM_GROUPS+1)*sizeof(double);
       memset(DCtot,0,LL);
    // Go through each predator/prey link and choose a gamma variable
    // with alpha =DC * GAMMA_ALPHA, beta =1
	  	 for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
           //DC[links]= path_DC[v->RRR.rpar_PreyFrom[links]][v->RRR.rpar_PreyTo[links]];
			     GAMMA_ALPHA = 50. * (1. - pedigree[4][v->RRR.rpar_PreyTo[links]]);           	                                  
			     DC[links] = d_gamma(GAMMA_ALPHA *
					                     (double)path_DC[v->RRR.rpar_PreyFrom[links]][v->RRR.rpar_PreyTo[links]],
					                     1.0, &v->rng);
			     if (DC[links]<EPSILON){DC[links]=2. * EPSILON;}
			     DCtot[v->RRR.rpar_PreyTo[links]] += DC[links];
			   }
		// Normalize and get QQ by multiplying by 
			 for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
			     DC[links] /=  DCtot[v->RRR.rpar_PreyTo[links]];
			     v->RRR.rpar_QQ[links] = (float)DC[links] * QB[v->RRR.rpar_PreyTo[links]] *
			                             v->RRR.rpar_B_BaseRef[v->RRR.rpar_PreyTo[links]] ;
			 }			   
//#undef GAMMA_ALPHA			 

 // Rule for Varying XX
    float PreyX[NUM_LIVING+NUM_DEAD+1],      PredX[NUM_LIVING+NUM_DEAD+1];
    float PredHandle[NUM_LIVING+NUM_DEAD+1], PreyHandle[NUM_LIVING+NUM_DEAD+1];
    float PredExp[NUM_LIVING+NUM_DEAD+1],    PreyExp[NUM_LIVING+NUM_DEAD+1];

    float vulvar;
    vulvar=VULVAR;

     memcpy(newvec, v->fit_vector, sizeof(float)*(NDIM+1));
     for (sp=1; sp<=NUM_LIVING; sp++){
        PreyHandle[sp] = 2.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*2 + sp];  
        PredHandle[sp] = 2.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*3 + sp]; 
        PreyX[sp]      = 2.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*4 + sp];
        PredX[sp]      = 2.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*5 + sp];   
        PreyExp[sp]    = 1.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*6 + sp];
        PredExp[sp]    = 1.0*vulvar*(uniform(&v->rng)-0.5) + v->fit_vector[NUM_LIVING*7 + sp];      
        newvec[NUM_LIVING*2 + sp] = PreyHandle[sp];
        newvec[NUM_LIVING*3 + sp] = PredHandle[sp];
        newvec[NUM_LIVING*4 + sp] = PreyX[sp];
        newvec[NUM_LIVING*5 + sp] = PredX[sp];
        newvec[NUM_LIVING*6 + sp] = PreyExp[sp];
        newvec[NUM_LIVING*7 + sp] = PredExp[sp];
     }
       
       
     PreyX[0] = 0;     PredX[0] = 0;
     PreyHandle[0] =0; PredHandle[0] =0;
     PredExp[0]=0;     PreyExp[0]=0;
     for (sp=NUM_LIVING+1; sp<=NUM_LIVING+NUM_DEAD; sp++){
         PreyX[sp] = 0;     PredX[sp] = 0;
         PreyHandle[sp] =0; PredHandle[sp] =0;
         PreyExp[sp]=0;     PredExp[sp]=0;
     }

    int pred, prey;
     for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){	
         prey=v->RRR.rpar_PreyFrom[links];
         pred=v->RRR.rpar_PreyTo[links];         
	       v->RRR.rpar_DD[links] = 1000. +  exp(PredHandle[pred]+ PreyHandle[prey]);
         v->RRR.rpar_VV[links] = 1. +  exp(PredX[pred]     + PreyX[prey]);
         v->RRR.rpar_HandleSwitch[links] = exp(0.05 * (PredExp[pred] + PreyExp[prey]));
     }	
/*
	  	 for (links=1; links<=v->RRR.rpar_NumPredPreyLinks; links++){
            x = v->fit_vector[NUM_LIVING*3 + v->RRR.rpar_PreyFrom[links]] +
                v->fit_vector[NUM_LIVING*4 + v->RRR.rpar_PreyTo[links]];
						    v->RRR.rpar_VV[links] = 1. + (MSCRAMBLE-1.)     * exp(x + 4. * (2.0*uniform(&v->rng)-1.0));
            x = v->fit_vector[NUM_LIVING*2 + v->RRR.rpar_PreyTo[links]];
						    v->RRR.rpar_DD[links] = 1. + (MHANDLE-1.)       * exp(x + 4. * (2.0*uniform(&v->rng)-1.0));
			 }
*/

 // For Handling time and Hotspots
 
     LL=(NUM_GROUPS+1)*sizeof(float);
      memset(v->RRR.rpar_PredTotWeight, 0, LL);
			memset(v->RRR.rpar_PreyTotWeight, 0, LL); 
 
    //int pred, prey;
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
				for (li=1; li<=v->RRR.rpar_NumPredPreyLinks; li++){

        //if (path_QB[v->RRR.rpar_PreyTo[li]]>EPSILON){
            v->RRR.rpar_PredPredWeight[li] = v->RRR.rpar_PredPredWeight[li]/ v->RRR.rpar_PredTotWeight[v->RRR.rpar_PreyFrom[li]]; 
            v->RRR.rpar_PreyPreyWeight[li] = v->RRR.rpar_PreyPreyWeight[li]/ v->RRR.rpar_PreyTotWeight[v->RRR.rpar_PreyTo[li]]; 
        //  }       
        }     
     
	  // Rule for Varying Catch
	     for (links=1; links<=v->RRR.rpar_NumFishingLinks; links++){
	         gr = v->RRR.rpar_FishingThrough[links]-NUM_LIVING-NUM_DEAD;
	         sp = v->RRR.rpar_FishingFrom[links];
		       catchVary = v->RRR.rpar_FishingQ[links] * 
					             path_BB[sp];
           //catchVary *= ((pedigree[4+gr][sp]*(2*uniform(&v->rng)-1.0))+1.0);
					 v->RRR.rpar_FishingQ[links] =  catchVary /
					                         v->RRR.rpar_B_BaseRef[v->RRR.rpar_FishingFrom[links]]; 
			 }
          
       // Defaults for group 0 ('outside')
          v->RRR.rpar_B_BaseRef[0]=1.0;  
	   	    v->RRR.rpar_MzeroMort[0]=0.0;
          v->RRR.rpar_UnassimRespFrac[0]=0.0; 
		      v->RRR.rpar_ActiveRespFrac[0]=0.0;
          v->RRR.rpar_FtimeQBOpt[0] = 1.0;
			    v->RRR.rpar_PBopt[0] = 1.0;	
   
      // Set Start Biomass and other state variables
	    for (sp=0; sp<=NUM_GROUPS; sp++){
			     v->state_BB[sp]    = v->RRR.rpar_B_BaseRef[sp];
			     v->state_Ftime[sp] = 1.0;
			     v->discarded[sp]   = 0;
			}

   for (i = 1; i<= juv_N; i++){
   	   v->RRR.rpar_drawn_K[i]     =   (1. + 0.5   * ( uniform(&v->rng)- 0.5)) * juv_VonBK[i];
		   v->RRR.rpar_drawn_AduZ[i]  =   (1. + 0.5   * ( uniform(&v->rng)- 0.5)) * juv_aduEqAgeZ[i];
		   v->RRR.rpar_drawn_JuvZ[i]  =   (1. + 1.0   * ( uniform(&v->rng)- 0.5)) * juv_juvEqAgeZ[i];
		   v->RRR.rpar_SpawnX[i]            = 10000.0;
       v->RRR.rpar_SpawnEnergy[i]       = 1.0;   //   (1. + 1.0   * ( uniform(&v->rng)- 0.5)); 
       v->RRR.rpar_SpawnAllocR[i]       = 1.0;
       v->RRR.rpar_SpawnAllocG[i]       = 1.0;
   }
   
   //initialize_stanzas(v);
}

// ---------------------------------------------------------------------

void *load_a_series_assess(void *threadarg)
{
struct SimRun *v;
int i, sp, kept,link, moreflag, k;
time_t start, end;
char outFile[80]; FILE *fptr;
float ssq;
const char* w = "w";
const char* ss = "%";
char c;
char logFile[80]; FILE *logptr;
char endFile[80]; FILE *endptr;
char binFile[80]; FILE *binptr;
FILE *dmpptr[NUM_LIVING+1];
 float sq[MAX_SYSTEMS];
   int sn[MAX_SYSTEMS];
float *sv[MAX_SYSTEMS];
float *newvec;
struct RatePar *Systems;
struct RatePar sys;

//struct gsl_histogram *params[NUM_LIVING+1][NUM_VARS];

      // Translate pointer for threads
         v = (struct SimRun *) threadarg;

     // Allocate space to store ecosystems when loaded
        Systems = (struct RatePar *)calloc(MAX_SYSTEMS,sizeof(struct RatePar)); 
        if (Systems==NULL){cout<<"allocation failure"<<endl; pthread_exit(NULL);}  

       newvec  = vector(0,NDIM);
       for (k=0; k<MAX_SYSTEMS; k++){
          sv[k]  = vector(0,NDIM);
       }
       
       map<string,int>::iterator iter, iter2;
       
       int g;
       string guild,guild2;
			 char ccr;
       //int gi;
       
    // Initialize Histograms
       //for (sp=1; sp<=NUM_LIVING; sp++){
       //     for (i=0; i<NUM_VARS; i++){
			//			    params[sp][i] = gsl_histogram_alloc(20);
			//			    gsl_histogram_set_ranges_uniform(params[sp][i],-6.0,6.0);
			//			}
			 //}

    // Load all ecosystems
       cout << "Loading ecosystems..." << endl;
       // THIS BINARY FORMAT MUST MATCH INPUT
          sprintf(binFile,"%s_Seed%d.bin",OUTFOLDER,v->thread_ID);
          binptr=fopen(binFile,"rb");

       // Vector output file
          sprintf(outFile,"%s_fitVectorsSeed%d.csv",OUTFOLDER,v->thread_ID);
          fptr=fopen(outFile,"w");
          fprintf(fptr,"System,SSQ,Num,Species,");
          for (i=0; i<NUM_VARS; i++){fprintf (fptr,"p%d,",i);}
          fprintf (fptr,"\n");

       // State output file
          sprintf(endFile,"%s_Seed%d_end.csv",OUTFOLDER,v->thread_ID);
          endptr=fopen(endFile,w);
          fprintf(endptr,"Seed,System,guildNum,guildName,Year,");
          fprintf(endptr,"BB,PB,MM,FF,");
          fprintf(endptr,"\n");

          moreflag=1;
          start = time(NULL);
          while (moreflag){
                // Get number kept in a batch
   								 fread(&kept,sizeof(int),1,binptr);
								// Read in that many systems      
                   for (k=0; k<kept; k++){
                     fread(&c,1,sizeof(char),binptr); if (c!='%'){cout << "alignment error" << endl;}
                     fread(&sn[k],1,sizeof(int),binptr);
                     fread(&sq[k],1,sizeof(float),binptr);
                     fread( sv[k],1,sizeof(float)*(NDIM+1),binptr);
                     fread(&Systems[k],1,sizeof(struct RatePar),binptr);
                   }
                   cout << kept << " ecosystems loaded." << endl; 

                // See if there's another batch to run
                   fread(&moreflag,1,sizeof(int),binptr);
                
                // OUTPUT VECTOR INFO
                   for (k=0; k<kept; k++){                                        
                     for (sp=1; sp<=NUM_LIVING; sp++){
                         fprintf(fptr,"%d,%g,%d,%s,",sn[k],sq[k],sp,path_species[sp]);
                         for (i=0; i<NUM_VARS; i++){
												     fprintf(fptr,"%g,",sv[k][i*NUM_LIVING+sp]);
                             //gsl_histogram_increment(params[sp][i],(double)(sv[k][i*NUM_LIVING+sp]));
                         }
                         fprintf (fptr,"\n");
                     }     
								}                 

                   for (k=0; k<kept; k++){    
                    // Load system rate parameters into the SimRun variable          
                       load_system(v, &Systems[k]);
                    // Run the ecosystem and calculate the SSQ
                       Adams_Basforth(v, 0, fit_Years);
                       //ssq = SSQ(v);
                    // Outputs
                       int yy;

                      for (iter=guildnum.begin(); iter!=guildnum.end(); iter++){
				               guild = (*iter).first;   // Guild name
                       g     = (*iter).second;  // Guild number for filenames
                         double bio,pb,mm,ff;
                         
                         for (yy=0; yy<=fit_Years; yy++){
                            bio = 0.0; 
														pb  = 0.0;
														mm  = 0.0;
														ff  = 0.0; 
                            for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                                 if  (guildlist[PathNames[sp]] == guild){
                                      bio += double(v->out_BB[sp][yy]);
                                      pb  += v->out_BB[sp][yy] * v->out_PB[sp][yy];
                                      mm  += v->out_BB[sp][yy] * v->out_MM[sp][yy] ;
                                      ff  += v->out_BB[sp][yy] * v->out_TotF[sp][yy];
													       }
								            }
								            if (bio>0.0){
					                     pb /= bio;
					                     mm /= bio;
					                     ff /= bio;
					                  } else {pb=0.0; mm=0.0; ff=0.0;}
					                  
                          
                          fprintf(endptr,"%d,%d,%d,%s,%d,",v->thread_ID,sn[k],g,guild.c_str(),
                                                          yy+fit_StartYear);
                                fprintf(endptr,"%g,",bio);
                                fprintf(endptr,"%g,",pb);
                                fprintf(endptr,"%g,",mm);
                                fprintf(endptr,"%g,",ff);                              
                                fprintf(endptr,"\n");
                          }
                      }
                      
                       float remaining,passed;
                       if ((k%20 == 10)){
													   end=time(NULL);
                             passed = (float)(end-start);
                             remaining = (float)(kept-k);
													
													   cout <<"seed: "<< v->thread_ID <<" " << k << 
													   " of "<< kept << "in batch, " << passed/60.0 << " min. "
													   << (remaining * passed/(float)i)/60.0 <<" left in batch" << endl;
											 }
                   }

          }  // end of moreflag
          fclose(endptr);
          fclose(fptr); 
          fclose(binptr);    

     //  for (sp=1; sp<=NUM_LIVING; sp++){
     //       for (i=0; i<NUM_VARS; i++){
		 //			     gsl_histogram_free(params[sp][i]);
		 //			}
		 //	 }

      for (k=0; k<MAX_SYSTEMS; k++){
           free_vector(sv[k],0,NDIM);
      }
      free_vector(newvec,0,NDIM);
      free(Systems);
      pthread_exit(NULL);
      
/* 
    // Open the output files
//       path_to_rates(v);  // needed to set number of links to make titles (silly)
        // not anymore, and it resets M0 so you dont get the perturbed parameter out
       for (sp=1; sp<=NUM_LIVING; sp++){
           sprintf(endFile,"%s_Seed%d_species%d.csv",OUTFOLDER,v->thread_ID,sp);
           dmpptr[sp]=fopen(endFile,w);
           fprintf(dmpptr[sp],"Seed,System,Species,Year,ssqOrig,ssqCheck,");
           fprintf(dmpptr[sp],"PB,QB,F,M2,Mzero,Mtot,B,");
//           for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
//               if (v->RRR.rpar_PreyFrom[link] == sp){
//                   fprintf(dmpptr[sp],"Mby_%s,",path_species[v->RRR.rpar_PreyTo[link]]);
//               }
//           }
//           for (link=1; link<=NUM_LIVING+NUM_DEAD; link++){     
//               fprintf(dmpptr[sp],"B_%s,",path_species[link]);
//           }
//           fprintf(dmpptr[sp],"SystemPP,");
//                 for (link=1 ; link<=juv_N; link++){
//                      fprintf(dmpptr[sp],"RS_%s,",path_species[juv_AduNum[link]]);
//                 }
           fprintf(dmpptr[sp],"\n");
        }
     // Cycle through all ecosystems and run them, producing output
        start = time(NULL);
        for (i=0; i<kept; i++){    
          // Load system rate parameters into the SimRun variable          
             load_system(v, &Systems[i]);
          // Run the ecosystem and calculate the SSQ
             Adams_Basforth(v, 0, BURN_TIME);
             ssq = SSQ(v);
          // Outputs
             int yy;
             for (yy=0; yy<=BURN_TIME; yy++){
             for (sp=1; sp<=NUM_LIVING; sp++){
                 //fprintf(dmpptr[sp],"Seed,System,Species,Year,ssqOrig,ssqCheck,");
                   fprintf(dmpptr[sp],"%d,%d,%s,%d,%g,%g,",v->thread_ID,sn[i],path_species[sp],
                                                         yy+fit_StartYear-1,sq[i],ssq);
                 //fprintf(dmpptr[sp],"F,Mtot,Mzero,QB,PB");
                   fprintf(dmpptr[sp],"%g,%g,%g,%g,%g,%g,%g,",
                                                  v->out_PB[sp][yy],
                                                  v->out_QB[sp][yy],
                                                  v->out_TotF[sp][yy],
                                                  v->out_M2[sp][yy],
                                                  v->RRR.rpar_MzeroMort[sp],
                                                  v->out_MM[sp][yy],
                                                  v->out_BB[sp][yy]);

//                   fprintf(dmpptr[sp],"%g,%g,%g,%g,%g,",v->out_TotF[sp],
//                                                  v->out_MM[sp][yy],
//                                                  v->RRR.rpar_MzeroMort[sp],
//                                                  v->out_QB[sp],
//                                                  v->out_PB[sp]);
//                 for (link=1; link<=v->RRR.rpar_NumPredPreyLinks; link++){
//                     if (v->RRR.rpar_PreyFrom[link] == sp){
//                         fprintf(dmpptr[sp],"%g,",v->out_LinksM[link]);
//                     }
//                 }
//                 for (link=1; link<=NUM_LIVING+NUM_DEAD; link++){     
//                     fprintf(dmpptr[sp],"%g,",v->out_BB[link][yy]);
//                 }
//                //fprintf(dmpptr[sp],"SystemPP,");
//                 fprintf(dmpptr[sp],"%g,",v->out_PP);
//                 for (link=1 ; link<=juv_N; link++){
//                      fprintf(dmpptr[sp],"%g,",v->RRR.rpar_SpawnEnergy[link]);
//                 }
                 fprintf(dmpptr[sp],"\n");
              }
              }
             //fprintf(endptr,"%d,%d,%d,%g,%d,%g,",v->thread_ID,i,sn[i],sq[i],fit_Years+fit_StartYear-1,ssq);
             //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){     
             //    fprintf(endptr,"%g,%g,",v->out_BB[sp][fit_Years-1],v->out_MM[sp][fit_Years-1]);
             //}
             //fprintf(endptr,"\n");  

					// Every few hundred iterations, print a status report
             float remaining,passed;
             if ((i%200 == 10)){
                 end=time(NULL);
                 passed = (float)(end-start);
                 remaining = (float)(kept-i);
                 cout <<"seed: "<< v->thread_ID <<" loaded: "<< i << " of "<< kept
                      <<" in " << passed/60.0 << " min. "<< (remaining * passed/(float)i)/60.0
                      <<" min. left " << endl;
             }                   
        }
  
   // Finally, clean up memory and exit thread
      for (sp=1; sp<=NUM_LIVING; sp++){
          fclose(dmpptr[sp]);
      }
*/    			
			            
}    

//----------------------------------------------------------------------
int loaded_ecosystem_series_assess(void)
{
   int r, rc, status;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];

   struct SimRun *Runs[PRUNS];

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      
    // CREATING a SERIES
        for (r=0; r<PRUNS; r++){
               Runs[r]=new_SimRun();
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
               Runs[r]->thread_ID=rseed+r;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, load_a_series, (void *) &Runs[r]);
               rc = pthread_create(&threads[r],&attr, load_a_series_assess, (void *) Runs[r]);
        }
      
        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }
      
        // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
                free_SimRun(Runs[r]);
           }           
   return 0;

}
