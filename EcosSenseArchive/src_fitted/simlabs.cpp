

//-------------------------------------------------------------------
int TL_From_FlowMat(){}
/*
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

int TL_From_FlowMat(){//Public Function TL_From_FlowMat(EFF() As Single) As Single()
float **LHS;         // Dim LHS!()
float *SumDC;        //Dim TL!()
int up,pred,prey;    //Dim up&, pred&, prey&

up  = NUM_LIVING + NUM_DEAD;   //up = UBound(EFF, 1): 
LHS   = matrix(1,up,1,up);     //ReDim LHS(1 To up, 1 To up)
SumDC = vector(1,up);			     //ReDim SumDC(1 To up)

for (pred=1; pred<=up; pred++){ //For pred = 1 To up
     SumDC[pred] = 0.0;
     for (prey=1; prey<=up; prey++){ //For prey = 1 To up
         LHS[pred][prey]=0.0;
				 SumDC[pred] += path_DC[prey][pred];           //LHS(pred, prey) = 0!: 
          //SumDC(pred) = SumDC(pred) + EFF(pred, prey)
     }// Next prey
} //Next pred

for (pred=1; pred<=up; pred++){ //For pred = 1 To up
    for (prey=1; prey<=up; prey++){ //For prey = 1 To up
         if(SumDC[pred]>0.0){  //If SingleComp(SumDC(pred), 0) > 0 Then _
            LHS[pred][prey] =  -path_DC[prey][pred]/ SumDC[pred];  //LHS(pred, prey) = -EFF(pred, prey) / SumDC(pred)
         }
    } //Next prey
    if (SumDC[pred]>0.0){ //If SingleComp(SumDC(pred), 0) > 0 Then
        LHS[pred][pred] = 1.0-path_DC[pred][pred]/SumDC[pred]; //LHS(pred, pred) = 1 - EFF(pred, pred) / SumDC(pred)
    }    
    else{
        LHS[pred][pred]=1.0; //LHS(pred, pred) = 1
    }//End If
} //Next pred

gsl_matrix *matLHS;
gsl_vector *vecTL;
gsl_vector *b;
gsl_permutation *p;
int sign;

matLHS = gsl_matrix_calloc(up,up);
b      = gsl_vector_calloc(up);
vecTL  = gsl_vector_calloc(up);
p      = gsl_permutation_calloc(up);

for (pred=1; pred<=up; pred++){ 
    for (prey=1; prey<=up; prey++){ 
        gsl_matrix_set(matLHS,pred-1,prey-1,LHS[pred][prey]);
    }
}
gsl_vector_set_all (b, 1.0);

// Solve LHS * b = TL using LU decomposition
   gsl_linalg_LU_decomp(matLHS,p,&sign);
   gsl_linalg_LU_solve(matLHS,p,b,vecTL);

for (pred=0; pred<up; pred++){
    TL[pred+1] = gsl_vector_get(vecTL,pred);
    //cout << path_species[pred+1] << "," << TL[pred+1] << endl;
}
		   
gsl_matrix_free(matLHS);
gsl_vector_free(b);
gsl_vector_free(vecTL);
gsl_permutation_free(p);

free_matrix(LHS,1,up,1,up);
free_vector(SumDC,1,up);

cout << "TL Matrix Decomposed" << endl;

} 
*/
//--------------------------------------------------------------

int scramble_series(void)
{
struct SimRun *v, *scr;
int i, sp, t, species;
float magnitude;
int duration, mag, loops;
float **noise;
const char* w = "w";
FILE *fptr;   
char outFile[80];					   

          // Allocate memory for the SimRun
             v   = new_SimRun();
						 scr = new_SimRun(); 
						 
						 init_by_array(0, &v->rng); 
 	 
					species = 33;	 
					magnitude = 1.1;
					duration = 40;
					
					loops=2;
					
					noise = matrix(0,NUM_LIVING,0,loops*MAX_YEARS*STEPS_PER_YEAR);
          for (sp=1; sp<=NUM_LIVING; sp++){
              for (t=1; t<=loops*MAX_YEARS*STEPS_PER_YEAR; t++){
                   noise[sp][t]= 0.1 * gaussian(&v->rng);
              }
          }

					for (mag = 2; mag <=8; mag+=2){
            for (duration = 0; duration <=100; duration +=10){
             magnitude = ((float)mag) * 0.1;
					// Sets up baseline run
					   ScrambleSelfWt = 1.0;
                base_system(v);
             ScrambleSelfWt = 0.0;    
						    base_system(scr);
      		// Run Adams-Basforth

            sprintf(outFile,"%s_scramble_mag%ddur%d_upperbio.csv",OUTFOLDER,mag,duration);
            fptr=fopen(outFile,w);
						    
						for (REPEAT = 0; REPEAT<loops; REPEAT++){
						
                 for (sp=1; sp<=NUM_LIVING; sp++){
                     for (t=1; t<=MAX_YEARS*STEPS_PER_YEAR; t++){
                         v->force_bymort[sp][t]= 1.0 +
                            noise[sp][t+loops*MAX_YEARS*STEPS_PER_YEAR];
                         scr->force_bymort[sp][t]= 1.0 +
                            noise[sp][t+loops*MAX_YEARS*STEPS_PER_YEAR];
                     }
                 }
											  
						     if (REPEAT ==0){
						    //for (sp=1; sp<=NUM_LIVING; sp++){
						       for (t=1; t<=duration * STEPS_PER_YEAR; t++){
        						    v->force_bymort[species][t]   +=magnitude;
												scr->force_bymort[species][t] +=magnitude;
	                 }
							  //}
						}
						//else{
						//    //for (sp=1; sp<=NUM_LIVING; sp++){
						//       for (t=1; t<=duration * STEPS_PER_YEAR; t++){
        	  //					    v->force_bymort[species][t]=1.0;
        		//				    scr->force_bymort[species][t]=magnitude;
	          //       }
						//	  //}											    
						// }    
						 cout << "Long Step " << REPEAT << " " << mag << " " << duration << endl;
						 Adams_Basforth(v, 0, MAX_YEARS);
						 Adams_Basforth(scr, 0, MAX_YEARS);
             //cout << "Done Step " << REPEAT << endl;
            
								// Save biomass to a file
            		   if (REPEAT == 0){
                     fprintf(fptr,"Year,Species,Base,Scram,\n");
                     //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                     //   if(v->discarded[sp]){fprintf(fptr,"ddbase_%s,",path_species[sp]);}
                     //   else             {fprintf(fptr,"base_%s,",path_species[sp]);}
                     //   if(scr->discarded[sp]){fprintf(fptr,"scram_dd%s,",path_species[sp]);}
                     //   else             {fprintf(fptr,"scram_%s,",path_species[sp]);}     
						         //}   fprintf(fptr,"\n");
                   }
                   for (t=0; t<MAX_YEARS; t++){
                       for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
                          if (TL[sp]>3.0){
                             fprintf(fptr,"%d,%s,",t + REPEAT*MAX_YEARS,path_species[sp]);              
                             fprintf(fptr,"%g,%g,",v->out_BB[sp][t],scr->out_BB[sp][t]);   
                             fprintf(fptr,"\n");
                          }
                         }
						       }                 
             }            
                 
						 fclose(fptr);
          }         
					}
					
					free_matrix(noise,0,NUM_LIVING,0,loops*MAX_YEARS*STEPS_PER_YEAR);      
          free_SimRun(v);
          free_SimRun(scr);
						 
return 0;
}

//----------------------------------------------------------------

int read_biomass(){

int cols, rows, iout,sp,t,c,iin, gr,j,i;
string sLine, ins;
double flout; 

   string in_series_names[NUM_LIVING+MAX_COLS+1];
   int    in_series_index[NUM_LIVING+MAX_COLS+1];
   int    in_series_type[NUM_LIVING+MAX_COLS+1];
  // FORMAT OF CSV INPUT FILE MUST BE 
	// Col 1 should be timestep (numbering arbitrary)
	// Row 1:  Series Name-
	// Row 2:  Group Number (0 for prim prod)
  // Row 3:  Series Type 

	// ALL LINES after that:  FITTING VALUE BY YEARS
		
		 dietN=0;
     for (c=0; c<=NUM_PREDPREYLINKS; c++){
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
           dietN=0;
           //if (flout>0.){fit_lastyear=t;}
           for (c=1; c<cols; c++){
		            parser >> flout;
                switch (in_series_type[c]){
                       case 40:
                          dietN++;
                          fit_diet_names[dietN]     = in_series_names[c];
                          fit_diet_type[dietN]      = in_series_type[c];
													fit_diet_lookup[dietN]  = in_series_index[c];
													fit_diet_OBS[dietN][t]  = flout;
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

//-------------------------------------------------------------------

struct NoiseSeries
{
   int spec_len, start_time, big_t, noise_len, thread_ID;
   double cv, ranscale;
   double *amp, *var, *offset;
   struct randvec rng;
   double *noise;
};

// -----------------------------------------------------------------------------
void *noisy_series(void *threadarg)
{
struct NoiseSeries *s;

       s = (struct NoiseSeries *) threadarg;
       
       NoiseSeries(s->spec_len, s->amp, s->var, s->offset, s->noise, 
	                 s->cv, s->ranscale, s->start_time, s->big_t, 
                   s->noise_len, &s->rng);

       pthread_exit(NULL); 
}
// -----------------------------------------------------------------------------

int Noisy_Run(void)
{
//#define SPEC_LEN     16000
//#define BIG_T        192000
#define SPEC_LEN     1000
#define BIG_T        12000
#define NUM_SERIES   NUM_LIVING 
#define NOISE_LEN    12000 * 5 
#define START_TIME   0
const char* w = "w";
FILE *fptr;   
char outFile[80]; 
time_t start, end;

   double **allnoise;

   int r, rc, status, rseed;
   pthread_attr_t attr;
   pthread_t threads[PRUNS];
   int i,t,j;

   struct NoiseSeries *Runs[PRUNS];
   
        // Main Output Matrix
           allnoise = dmatrix(0,NUM_SERIES,0,NOISE_LEN);

        // Required thread maintenence
           pthread_attr_init(&attr);
           pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // CREATING a SERIES
        rseed=0;
        for (r=0; r<PRUNS; r++){
               Runs[r] = (struct NoiseSeries *)calloc(1,sizeof(struct NoiseSeries));
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               init_by_array(rseed+r, &Runs[r]->rng); 
            // Initialize unchanging values
               Runs[r]->spec_len   = SPEC_LEN;
               Runs[r]->start_time = START_TIME;
               Runs[r]->big_t      = BIG_T;
               Runs[r]->noise_len  = NOISE_LEN;
               Runs[r]->ranscale = 0.9999;  
               Runs[r]->cv       = 1.0;

               Runs[r]->amp    = dvector(0,SPEC_LEN);
               Runs[r]->var    = dvector(0,SPEC_LEN);
               Runs[r]->offset = dvector(0,SPEC_LEN);
               Runs[r]->noise  = dvector(0,NOISE_LEN);

               for (i=0; i<=SPEC_LEN; i++){
			              Runs[r]->amp[i] = 1.0;
			              Runs[r]->var[i] = 1.0; 
	             }               
        }

  start = time(NULL);
  for (i=0; i<=NUM_SERIES/PRUNS; i++){        
           for (r=0; r<PRUNS; r++){
            // offset changed with each run
               for (j=0; j<=SPEC_LEN; j++){Runs[r]->offset[j] = genrand_res53(&Runs[r]->rng);}
            // Run ID number for output files
               Runs[r]->thread_ID=i*PRUNS+r;
            // Split the run onto a new thread and run it!
               cout << "Thread " << i*PRUNS+r << endl;
               rc = pthread_create(&threads[r],&attr, noisy_series, (void *) Runs[r]);
           }

        // Wait for the threads to come together
           for (r=0; r<PRUNS; r++){   
               rc=pthread_join(threads[r], (void **)&status);
               if (rc){cout << r << "error" << endl;}
           }

        // Assemble the numbers
           for (r=0; r<PRUNS; r++){
              if (Runs[r]->thread_ID<=NUM_SERIES){
        				for (t=0; t<=NOISE_LEN; t++){
                    allnoise[Runs[r]->thread_ID][t] = (Runs[r]->noise[t]);
                }
              }
           }
  end=time(NULL);
  cout << "Running for " << (end-start)/60.0 << " min." << endl;
  cout << "End in " <<  (1.0 - ((float)(i+1)/((float)NUM_SERIES/(float)PRUNS))) * 
                        (((float)(end-start))/60.0)/
                               ((float)(i+1)/((float)NUM_SERIES/(float)PRUNS)) << " min." << endl;
  }
  
  // Outputs
     sprintf(outFile,"%s_noise.csv",OUTFOLDER);
     fptr=fopen(outFile,w);
     fprintf(fptr,"Month,");  
     for (i=0; i<=NUM_SERIES; i++){fprintf(fptr,"Species%d,",i);}
     fprintf(fptr,"\n");
     for (t=0; t<=NOISE_LEN; t++){
         fprintf(fptr,"%d,",t);
         for (i=0; i<=NUM_SERIES; i++){
         //noiseval = (exp(-cv*cv/2 + cv*allnoise[i][t]));
           fprintf(fptr,"%g,",allnoise[i][t]);
		     }
		     fprintf(fptr,"\n");
     }   
     fclose(fptr);

   // Cleanup
      free_dmatrix(allnoise,0,NUM_SERIES,0,NOISE_LEN);
  
   // Thread cleanup
           pthread_attr_destroy(&attr);
           // Note:  Exit of main thread is in function main
           for (r=0; r<PRUNS; r++){
               free_dvector(Runs[r]->amp,0,SPEC_LEN);
               free_dvector(Runs[r]->var,0,SPEC_LEN);
               free_dvector(Runs[r]->offset,0,SPEC_LEN);
               free_dvector(Runs[r]->noise,0,NOISE_LEN);
               free(Runs[r]);
           }           

return 0;
  
#undef START_TIME
#undef NOISE_LEN
#undef NUM_SERIES
#undef BIG_T
#undef SPEC_LEN

}
//---------------------------------------------------------------
/*
int Noisy_Run(struct SimRun *v)
{
#define SPEC_LEN 500
#define BIG_T 12000   
#define START_TIME 0
const char* w = "w";
FILE *fptr;   
char outFile[80]; 

double cv, ranscale;
double amp[SPEC_LEN+1], var[SPEC_LEN+1], offset[SPEC_LEN+1];
double noise[MAX_YEARS*STEPS_PER_YEAR+1];
double allnoise[10+1][MAX_YEARS*STEPS_PER_YEAR+1];
int i,t,j;


// O = DLLNoiseSeries(SPEC_LEN, amp(0), Var(0), offset(0), Final(0), Freq_CV(j), 0, 0, BIG_T, UntilMonth, RS.SeedVec(0))
   for (i=0; i<=SPEC_LEN; i++){
			 amp[i] = 1.0;
			 var[i] = 1.0; 
	 }

   //for (i=0; i<=10; i++){
        cout << "ranseries " << i << endl;
        //ranscale = (double)i * 0.1;
        ranscale = 0.99;

        //NoiseSeries(SPEC_LEN, amp, var, offset, &allnoise[i][0], 
	      //            cv, ranscale, START_TIME, BIG_T, MAX_YEARS*STEPS_PER_YEAR);
				for (j=0; j<=SPEC_LEN; j++){offset[j] = genrand_res53(&v->rng);}
        NoiseSeries(SPEC_LEN, amp, var, offset, noise, 
	                  cv, ranscale, START_TIME, BIG_T, MAX_YEARS*STEPS_PER_YEAR);				     
 	 //}
   for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++)
       {force_byprey[0][t]= (float)(exp(-cv*cv/2 + cv*noise[t]));}

  cv       = 0.2;        
  for (i=1; i<=NUM_LIVING; i++){
        for (j=0; j<=SPEC_LEN; j++){offset[j] = genrand_res53(&v->rng);}
        cout << "ranseries " << i << endl;
        NoiseSeries(v, SPEC_LEN, amp, var, offset, noise, 
	                  1.0, ranscale, START_TIME, BIG_T, MAX_YEARS*STEPS_PER_YEAR);	
				for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++)
               {force_bymort[i][t]= (float)(exp(-cv*cv/2 + cv*noise[t]));}						
	}

sprintf(outFile,"%s/noise.csv",OUTFOLDER);
fptr=fopen(outFile,w);
for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++){
     fprintf(fptr,"%d,",t);
     for (i=0; i<=NUM_LIVING; i++){
         fprintf(fptr,"%g,",force_bymort[i][t]);
		 }
		 fprintf(fptr,"\n");
}
fclose(fptr);

#undef START_TIME
#undef BIG_T
#undef SPEC_LEN

}
*/


//---------------------------------------------------------------
void niche_mat(){

const char* w = "w";
FILE *fptr;   
char outFile[80];
double flow;
int group,cat;
int sp,step,prey,sp1,sp2,g;
int group_ID[NUM_GROUPS+1];
double group_MAT[NUM_GROUPS+1][NUM_GROUPS+1];
double group_MOR[NUM_GROUPS+1][NUM_GROUPS+1];

double group_TOT[NUM_GROUPS+1];
//int group_SIM[NUM_GROUPS+1][NUM_GROUPS+1];
int SIM_GROUPS;


for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
    group_ID[sp]=sp;
		}


//SIM_GROUPS=NUM_LIVING+NUM_DEAD;

for (step=1; step<NUM_LIVING+NUM_DEAD; step++){


   //fclose(fptr);
   
   SIM_GROUPS=0;
   for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
       if (group_ID[sp]>SIM_GROUPS){SIM_GROUPS=group_ID[sp];}
   }
   
   
    for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
         group_TOT[sp]=0.0;
         for (prey=0; prey<=NUM_LIVING+NUM_DEAD; prey++){
              group_MAT[sp][prey]=0.0;
              group_MOR[sp][prey]=0.0;
          }
    }
    
    for (sp=1; sp<=NUM_LIVING; sp++){
         for (prey=1; prey<=NUM_LIVING+NUM_DEAD; prey++){
             group_MAT[group_ID[sp]][group_ID[prey]] +=
                path_BB[sp]*path_PB[sp]*path_DC[prey][sp];
             group_TOT[group_ID[sp]] += path_BB[sp]*path_PB[sp]*path_DC[prey][sp];
             //cout << group_MAT[group_ID[sp]][group_ID[prey]] << endl;
				 }
				 if ((path_PB[sp]>EPSILON) && (path_QB[sp]<=EPSILON)){
				    group_MAT[group_ID[sp]][0] += path_BB[sp]*path_PB[sp];
				    group_TOT[group_ID[sp]]    += path_BB[sp]*path_PB[sp];
				 }     
    }

    for (sp=1; sp<=NUM_LIVING; sp++){
         for (prey=1; prey<=NUM_LIVING+NUM_DEAD; prey++){
             //group_MOR[group_ID[sp]][group_ID[prey]] +=
             //   path_QB[prey]*path_BB[prey]*path_DC[sp][prey];
             //group_TOT[group_ID[sp]] += path_QB[prey]*path_BB[prey]*path_DC[sp][prey];
             //cout << group_MOR[group_ID[sp]][group_ID[prey]] << endl;
				 }	     
    }

     for (sp=1; sp<=SIM_GROUPS; sp++){
         if (group_TOT[sp]>0.0){
				    for (prey=0; prey<=SIM_GROUPS; prey++){
				      group_MAT[sp][prey] /= group_TOT[sp];
				      //group_MOR[sp][prey] /= group_TOT[sp];
				    }
				 }
		 }

    double simNom, simDenom, simVal, simMax;
    int max1, max2;
    simMax = -1.0;
    for (sp1=1; sp1<=SIM_GROUPS-1; sp1++){
       for (sp2=sp1+1; sp2<=SIM_GROUPS; sp2++){
           simNom=0.0; simDenom=0.0;
					 for (prey=0; prey<=SIM_GROUPS; prey++){
					     //cout <<group_MAT[sp1][prey] <<" "<< group_MAT[sp2][prey] <<" "<< max2<<" "<<endl;
					     simNom += group_MAT[sp1][prey] * group_MAT[sp2][prey];
					     simDenom += group_MAT[sp1][prey] * group_MAT[sp1][prey] + 
					                 group_MAT[sp2][prey] * group_MAT[sp2][prey];

					     //simNom += group_MOR[sp1][prey] * group_MOR[sp2][prey];
					     //simDenom += group_MOR[sp1][prey] * group_MOR[sp1][prey] + 
					     //            group_MOR[sp2][prey] * group_MOR[sp2][prey];					                 
					                 
					 }
					 if (simDenom>0.0){
					    simVal =  simNom/(simDenom/2.0);
					 }
					 else {
					    simVal = 0.0;
					 }
					 if (simVal>=simMax){
					     simMax=simVal;
					     max1=sp1;
					     max2=sp2;
					 } 
		   }
    }
    
    
   
   //cout <<simMax <<" "<< max1 <<" "<< max2<<" "<<endl;
   //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
   //    if (group_ID[sp]==max2){cout << path_species[sp] << ",";}
   //}
   //cout << "joins " << endl;
   //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
   //    if (group_ID[sp]==max1){cout << path_species[sp] << ",";}
   //}
	 //cout << endl;

    int LostID, JoinedID; double LostVal, lostMaxTL;
    LostVal=simMax;
    JoinedID=max1;    
    lostMaxTL=-1.0;
    for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
        if ((group_ID[sp]==max2)&& (TL[sp]>lostMaxTL)){
				   lostMaxTL = TL[sp];
           LostID    = sp;
        }
        if (group_ID[sp]==max2){group_ID[sp]=max1;}
        if (group_ID[sp]>max2) {group_ID[sp]--;}
    }
    //SIM_GROUPS--;
   SIM_GROUPS=0;
   for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
       if (group_ID[sp]>SIM_GROUPS){SIM_GROUPS=group_ID[sp];}
   }
   
   int group_sig[NUM_GROUPS+1], group_count[NUM_GROUPS+1];
   double group_TL[NUM_GROUPS+1],group_inflow[NUM_GROUPS+1],maxTL;
   
   for (sp=1; sp<=SIM_GROUPS; sp++){
       maxTL=-1.0;
       group_TL[sp]=0.0; group_inflow[sp]=0.0; group_count[sp]=0;
       for (prey=1; prey<=NUM_LIVING+NUM_DEAD; prey++){
           if (group_ID[prey]==sp){
              if (TL[prey]>maxTL){maxTL=TL[prey]; group_sig[sp]=prey;}
              group_TL[sp]     += TL[prey]*path_BB[sp]*path_PB[sp];
              group_inflow[sp] += path_BB[sp]*path_PB[sp];
              group_count[sp]  ++;
					 }
       }
       if (group_inflow[sp]>0.0){group_TL[sp]/=group_inflow[sp];}
   }
   
   double flow, flowBetween[NUM_GROUPS+1],flowWithin[NUM_GROUPS+1];
   for (sp=1; sp<=SIM_GROUPS; sp++){
       flowBetween[sp]=0.0; flowWithin[sp]=0.0;
   }
   for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
       for (prey=1; prey<=NUM_LIVING+NUM_DEAD; prey++){
           flow = path_BB[sp]*path_PB[sp]*path_DC[prey][sp] +
                  path_QB[prey]*path_BB[prey]*path_DC[sp][prey];
           if (group_ID[sp]==group_ID[prey]){
              flowWithin[group_ID[sp]]+=flow;
           }
           else{
              flowBetween[group_ID[sp]]+=flow;
           }           
        }    
		}				 
		
	if (step==1){cout << "Groups,Sim,Combined,Gnum,HighTL,TL,Count,Inflow,Within,Between" << endl;}
   for (sp=1; sp<=SIM_GROUPS; sp++){
       cout << SIM_GROUPS << "," << LostVal;
       if (sp==JoinedID){cout << "," << path_species[LostID];}
                   else {cout << "," << "NA";}       
			 cout << "," << sp << "," << path_species[group_sig[sp]] 
			      << "," <<  group_TL[sp] << "," << group_count[sp] << "," << group_inflow[sp]  
			      << "," << flowWithin[sp] << "," << flowBetween[sp] << endl;
	 }
   
/*
char aWord[100];
char* xlabels[NUM_GROUPS+1];
int n;

n=0;
for (g=1; g<=SIM_GROUPS; g++){
     for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
          if (group_ID[sp]==g){
             strcpy(aWord, path_species[sp]);
						 int len = strlen(aWord)+1;
						 char* newSpace = new char[len];
						 strcpy(newSpace, aWord);
						 xlabels[sp]=newSpace;   
					}
     }
}

*/
   sprintf(outFile,"%s_%.3d_njump.csv",OUTFOLDER,SIM_GROUPS);
   fptr=fopen(outFile,w);
   fprintf(fptr,"Type,Step,Similarity,Group,GroupNum,PreyNum,Prey,PreyPer\n");
   
    for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
				 if ((path_PB[sp]>EPSILON) && (path_QB[sp]<=EPSILON)){
            fprintf(fptr,"Diet,%d,%f,%s,%d,%d,%s,%g\n",step,simMax,path_species[sp],group_ID[sp],
						                                0,"Zero",path_BB[sp]*path_PB[sp]);
				 }
				 //else{
         //   fprintf(fptr,"%d,%f,%s,%d,%d,%s,%g\n",step,simMax,path_species[sp],group_ID[sp],
				 //	                                0,"Zero",0.0);				 
				 //}
				 
        for (prey=1; prey<=NUM_LIVING+NUM_DEAD; prey++){
           //cout << "bleah " << sp << " " << prey << endl;
           if (path_DC[prey][sp]>1e-8){
            fprintf(fptr,"Diet,%d,%f,%s,%d,%d,%s,%g\n",step,simMax,path_species[sp],group_ID[sp],
						                                group_ID[prey],path_species[prey],
																						path_BB[sp]*path_PB[sp]*path_DC[prey][sp]);
					}
          if (path_DC[sp][prey]>1e-8){
            fprintf(fptr,"Mort,%d,%f,%s,%d,%d,%s,%g\n",step,simMax,path_species[sp],group_ID[sp],
						                                group_ID[prey],path_species[prey],
																						path_QB[prey]*path_BB[prey]*path_DC[sp][prey]);
					}					
        }
    }
    //cout << "bah!" << endl;
    //for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
     //   fprintf(fptr,"%d,",group_ID[sp]);
    //}
    //fprintf(fptr,"\n");    
    
fclose(fptr);    
}		  


exit(0);
		    
		  sprintf(outFile,"%s_nichemat.csv",OUTFOLDER);
      fptr=fopen(outFile,w);
      fprintf(fptr,"Groups,");
      fprintf(fptr,"PrimProd,");
      for (group=1; group<=NUM_LIVING+NUM_DEAD; group++){
          fprintf(fptr,"Diet%s,",path_species[group]);
      }
      for (group=1; group<=NUM_LIVING; group++){
          fprintf(fptr,"Mort%s,",path_species[group]);
      }      
      fprintf(fptr,"Mzero,\n");
      
      for (group=1; group<=NUM_LIVING+NUM_DEAD; group++){
          fprintf(fptr,"%s,",path_species[group]);
          if ((path_PB[group]>EPSILON) && (path_QB[group]<=EPSILON)){
             flow=path_BB[group]*path_PB[group];
          }
          else {flow=0.0;
					}
          fprintf(fptr,"%g,",flow);
          
          for (cat=1; cat<=NUM_LIVING+NUM_DEAD; cat++){
              flow = path_BB[group]*path_PB[group]*path_DC[cat][group];
              fprintf(fptr,"%g,",flow);
          }
          for (cat=1; cat<=NUM_LIVING; cat++){
              flow = path_BB[cat]*path_QB[cat]*path_DC[group][cat];
							fprintf(fptr,"%g,",flow);
				  }
				  flow = path_BB[group]*path_PB[group]*(1.0-path_EE[group]);
				  fprintf(fptr,"%g,",flow);
			fprintf(fptr,"\n");	 
		  }
			 
			fclose(fptr);		 
}

//---------------------------------------------------------------


int sim_labs(){
    // Noisy_Run(); //exit(0);
    //scramble_series();
      //exit(0);
      niche_mat();
}
//------------------------------------------------------------
