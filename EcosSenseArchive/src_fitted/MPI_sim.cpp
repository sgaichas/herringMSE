
#include <fftw3.h>


void getAmpPhase(fftw_complex *vals, int NN){

double A,B;
double R,phi;
int i;

for (i=0; i<NN; i++){
   A = vals[i][0];  B = vals[i][1];

   if (B==0.0){
      if (A<0.0) {phi=-3.14159265358979/2.0;}
      else       {phi=+3.14159265358979/2.0;}     
   }
   else{
      phi=atan(A/B);
   }

   R=sqrt(A*A + B*B);

   vals[i][0]=R;
   vals[i][1]=phi;

}
}
// -----------------------------------------------------------------------------
int noise_normalize(double *final, int TLEN, double cv){

	int t;
	double X, X2, MM, SD, TL;


	X=0.0; X2=0.0;
	for (t=0; t<TLEN; t++){
	     X  += final[t];
		 X2 += final[t] * final[t];
	}
	TL=(double)TLEN;
	MM   = X/TL;
	SD  = TL * X2 - X * X;
	SD  = sqrt(SD / (TL * (TL - 1.0) ) );
	if  (SD <= 0.0) {SD=1;}
	for (t=0; t<TLEN; t++){
		 final[t] = cv * (final[t]-MM)/SD;
	}

}

// -----------------------------------------------------------------------------
int MPI_noise(int myid,int numprocs){

#define TWOPI        6.283185307179586476925286 
//#define SPEC_LEN     1000
#define BIG_T        12000
#define NUM_SERIES   NUM_LIVING

#define CYCLES       11 

#define NOISE_LEN    24000+2400 
#define START_TIME   2400
#define SAMPLE_LEN   24000
#define SPEC_LEN     2000

   struct SimRun *v;

   MPI_Status status;
   int tag;
   
   int i, iters,p,t,f,S,sp,key_species;
   double cv, ranscale;
   double *amp, *offset;
   struct randvec rng;
   double *noise, *freq;
   double *allnoise;
   double *back, *allback;
   double *allfreqs;
   double *real, *allreal;
   double *imag, *allimag;
   double fq, phi;

   double **longB_out, **longB_amp, **longB_phi;
   double **longP_out, **longP_amp, **longP_phi;
   double *outamp, *outphi;
   int TIME_LEN;
   

   const char* w = "w";
   FILE *fptr, *cptr;   
   char outFile[80];	
   
   double *in;
   fftw_complex *out;
   fftw_plan plan, planback;

    // Allocate memory for a Sim Run
             v=new_SimRun();
   
   // Allocate storage for calculations and results
      //iters= (NUM_SERIES/numprocs) + 1;
      iters=2;
      
      allnoise=dvector(0,iters*numprocs*NOISE_LEN);
      allback=dvector(0,iters*numprocs*NOISE_LEN);
      allfreqs=dvector(0,iters*numprocs*SAMPLE_LEN);
      allreal=dvector(0,iters*numprocs*SAMPLE_LEN);
      allimag=dvector(0,iters*numprocs*SAMPLE_LEN);

      TIME_LEN  = (CYCLES-1)*fit_Years*STEPS_PER_YEAR;
      longB_out = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));
      longB_amp = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));
      longB_phi = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));  
      
      longP_out = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));
      longP_amp = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));
      longP_phi = dmatrix(0,(NUM_LIVING+1),0,(TIME_LEN+1));        
      
          
      noise     = dvector(0,TIME_LEN);
      outamp    = dvector(0,TIME_LEN);
      outphi    = dvector(0,TIME_LEN);
      
      //longB_res=dvector(fit_Years*STEPS_PER_YEAR*(CYCLES-1));
      
      amp    = dvector(0,SPEC_LEN); 
      offset = dvector(0,SPEC_LEN);               

      back   = dvector(0,NOISE_LEN*iters);
      freq   = dvector(0,SAMPLE_LEN*iters);
      real   = dvector(0,SAMPLE_LEN*iters);
      imag   = dvector(0,SAMPLE_LEN*iters);
      init_by_array(rseed+myid, &rng); 

  // Pointers needed for calculating spectra via FFT
     //in   = noise;
     in =   (double*) fftw_malloc(sizeof(double) * TIME_LEN);
     out  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TIME_LEN);
     plan = fftw_plan_dft_r2c_1d(TIME_LEN, in, out, FFTW_PRESERVE_INPUT);     

     //in = back;
     planback = fftw_plan_dft_c2r_1d(SAMPLE_LEN, out, in, FFTW_PRESERVE_INPUT);

  // MAIN SPECTRUM GENERATION LOOP
     for (i=0; i<iters; i++){

       // SET FUNCTIONS OF THE NOISE 
          ranscale = 0.9999;  
          cv       = 1.0;
          for (t=0; t<=SPEC_LEN; t++){
			        amp[t] = 1.0;
              offset[t] = (genrand_res53(&rng)-0.5)*(double)BIG_T; 
	        }

      // GENERATE THE NOISE

		     for (t=0; t<TIME_LEN; t++){
		      
          // FOURIER NOISE (NoiseSeries)
          //   noise[i*NOISE_LEN+t] = 0.0;
        	//   //for (S=1; S<=SPEC_LEN; S++){                 
		      //   //    fq=double(S)/(double(BIG_T));   
          //     for (S=4; S<=200*4; S++){
          //       fq=1.0/double(S);
			    //       noise[i*NOISE_LEN+t] += amp[S] * (  
				  //           (1.0-ranscale*genrand_res53(&rng)) * sin(TWOPI * (fq * (double)(t+START_TIME) + offset[S])) +
				  //	         (1.0-ranscale*genrand_res53(&rng)) * cos(TWOPI * (fq * (double)(t+START_TIME) + offset[S]))						                  
					//	     );
		      //   }
	        
	        S = myid*iters+i;
	        
	        fq  =  1.0/(5.0*(12.0*(double)(myid*iters+i)+1.0));
	        //phi = (double)i * 0.2 - 0.2; 
	        phi = 0; 
          
          noise[t] = sin(TWOPI * (fq * (double)t + phi));
          
          } // end of t noise generation loop

      // NORMALIZE NOISE
         noise_normalize(&noise[TIME_LEN],TIME_LEN,cv);
        
          // Load starting fit vector (or set to 0 if none existing)    
             //load_instep(v);  
          // First set default parameters for sim
             path_to_rates(v);
          // Then apply x vector to the default non-juvenile parameters.  
             //apply_vector_to_rates(v, v->fit_vector);
          // Then initialize stanzas
             initialize_stanzas(v);  

             for (REPEAT = 0; REPEAT<CYCLES; REPEAT++){   
                 //for (sp=0; sp<=NUM_GROUPS; sp++){
		                 for (t=0; t<=fit_Years*STEPS_PER_YEAR; t++){
		                      sp=121; v->force_bymort[sp][t]=1.+0.99*noise[REPEAT*fit_Years*STEPS_PER_YEAR+t];
		                      sp=122; v->force_bymort[sp][t]=1.+0.99*noise[REPEAT*fit_Years*STEPS_PER_YEAR+t];
		                      sp=123; v->force_bymort[sp][t]=1.+0.99*noise[REPEAT*fit_Years*STEPS_PER_YEAR+t];
		                      sp=124; v->force_bymort[sp][t]=1.+0.99*noise[REPEAT*fit_Years*STEPS_PER_YEAR+t];                          		                      
				             }
		             //} 
                 //output_run(v, outFile);
                 cout << "Long Step " << fit_Years << " " << REPEAT << endl;
                 Adams_Basforth(v,0, fit_Years);
                 cout << "Made it " << REPEAT << endl;
                 // Save biomass to a file
                 if (REPEAT){
                    for (sp=0; sp<=NUM_LIVING; sp++){
		                    for (t=0; t<fit_Years*STEPS_PER_YEAR; t++){
		                    //cout << sp <<","<< t << endl;
                             longB_out[sp][(REPEAT-1)*fit_Years*STEPS_PER_YEAR + t] =  v->monthly_BB[sp][t];
                             longP_out[sp][(REPEAT-1)*fit_Years*STEPS_PER_YEAR + t] =  v->monthly_PP[sp][t];
                        }   
						        }
                  }                
             } 


memcpy(in,noise,TIME_LEN*sizeof(double));
fftw_execute(plan); // RUN FFT
//getAmpPhase(out,TIME_LEN);

    for (f=0; f<TIME_LEN; f++){
        longB_amp[0][f]=out[f][0]; 
        longB_phi[0][f]=out[f][1];
        longP_amp[0][f]=out[f][0]; 
        longP_phi[0][f]=out[f][1];
    }
        
for (sp=1; sp<=NUM_LIVING; sp++){

    memcpy(in,longB_out[sp],TIME_LEN*sizeof(double));
    fftw_execute(plan); // RUN FFT
    //getAmpPhase(out,TIME_LEN);
    for (f=0; f<TIME_LEN; f++){
        longB_amp[sp][f]=out[f][0]; 
        longB_phi[sp][f]=out[f][1];
    }

    memcpy(in,longP_out[sp],TIME_LEN*sizeof(double));
    fftw_execute(plan); // RUN FFT
    //getAmpPhase(out,TIME_LEN);
    for (f=0; f<TIME_LEN; f++){
        longP_amp[sp][f]=out[f][0]; 
        longP_phi[sp][f]=out[f][1];
    }    
        
}



        sprintf(outFile,"%s_freq_series_%d_%d.csv",OUTFOLDER,i,myid);
        fptr=fopen(outFile,w);        
           fprintf(fptr,"t,noiseR,noiseI");
           for (sp=1; sp<=NUM_LIVING; sp++){
                fprintf(fptr,"R%g_%s,I%g_%s,",1.0/fq,path_species[sp],1.0/fq,path_species[sp]);
                }
           fprintf(fptr,"\n");  		    
                    
            for (t=0; t<TIME_LEN; t++){
                    fprintf(fptr,"%d,",t);
                    for (sp=0; sp<=NUM_LIVING; sp++){
                            fprintf(fptr,"%g,", longB_amp[sp][t],longB_phi[sp][t]);
                        }
                        fprintf(fptr,"\n");   
						        }        
        
        fclose(fptr);


/*        
        sprintf(outFile,"%s_time_series_%d_%d.csv",OUTFOLDER,i,myid);
        fptr=fopen(outFile,w);        
           fprintf(fptr,"t,noise,");
           for (sp=0; sp<=NUM_LIVING; sp++){
                fprintf(fptr,"%g_%s,",1.0/fq,path_species[sp]);
                }
           fprintf(fptr,"\n");  		    
                    
            for (t=0; t<TIME_LEN; t++){
                    fprintf(fptr,"%d,%g,",t,noise[t]);
                    for (sp=0; sp<=NUM_LIVING; sp++){
                            fprintf(fptr,"%g,", longB_out[sp][t]);
                        }
                        fprintf(fptr,"\n");   
						        }        
        
        fclose(fptr);
*/       
  cout << "Made it through" << endl;     
/*        
      // Successivly calculate frequencies         
         //in = &noise[(i+1)*NOISE_LEN-SAMPLE_LEN];

         //plan = fftw_plan_dft_r2c_1d(SAMPLE_LEN, in, out, FFTW_PRESERVE_INPUT);  
         memcpy(in,&noise[(i+1)*NOISE_LEN-SAMPLE_LEN],SAMPLE_LEN*sizeof(double));
         fftw_execute(plan); // RUN FFT
         
      // Store the power spectrum          
          for (f=0; f<SAMPLE_LEN; f++){
              freq[i*SAMPLE_LEN+f] =  out[f][0]*out[f][0] + out[f][1]*out[f][1];
              real[i*SAMPLE_LEN+f] =  out[f][0];
              imag[i*SAMPLE_LEN+f] =  out[f][1];   
          }
          
          
      // Try the backwards step
         //in =  &back[(i+1)*NOISE_LEN-SAMPLE_LEN];
         //planback = fftw_plan_dft_c2r_1d(SAMPLE_LEN, out, in, FFTW_PRESERVE_INPUT);
         fftw_execute(planback);   
         memset(&back[i*NOISE_LEN],0,NOISE_LEN*sizeof(double));
         memcpy(&back[(i+1)*NOISE_LEN-SAMPLE_LEN],in,SAMPLE_LEN*sizeof(double));
*/
     }
     exit(0);
     
     // Clean up FFT details
        fftw_destroy_plan(plan);
        fftw_destroy_plan(planback); 
        fftw_free(in);
        fftw_free(out);
     
/*     
     // Interprocess communication : spreading all results among processors
        // First collect all results into a single matrix on process 0
           tag=1;            
           if (myid>0){
              MPI_Send(noise, iters*NOISE_LEN, MPI_DOUBLE,0,tag+0,MPI_COMM_WORLD); 
              MPI_Send(freq,  iters*SAMPLE_LEN, MPI_DOUBLE,0,tag+1,MPI_COMM_WORLD); 
              MPI_Send(real,  iters*SAMPLE_LEN, MPI_DOUBLE,0,tag+2,MPI_COMM_WORLD); 
              MPI_Send(imag,  iters*SAMPLE_LEN, MPI_DOUBLE,0,tag+3,MPI_COMM_WORLD);               
              MPI_Send(back,  iters*NOISE_LEN, MPI_DOUBLE,0,tag+4,MPI_COMM_WORLD); 
           }              
           else{
               memcpy(&allnoise[0],noise,iters*NOISE_LEN*sizeof(double));
               memcpy(&allfreqs[0],freq,iters*SAMPLE_LEN*sizeof(double));
               memcpy(&allreal[0],real,iters*SAMPLE_LEN*sizeof(double));
               memcpy(&allimag[0],imag,iters*SAMPLE_LEN*sizeof(double));               
               memcpy(&allback[0],back,iters*NOISE_LEN*sizeof(double));
               for (p=1; p<numprocs; p++){
                  MPI_Recv(&allnoise[p*iters*NOISE_LEN],iters*NOISE_LEN,MPI_DOUBLE,p,tag,MPI_COMM_WORLD,&status);
                  MPI_Recv(&allfreqs[p*iters*SAMPLE_LEN],iters*SAMPLE_LEN,MPI_DOUBLE,p,tag+1,MPI_COMM_WORLD,&status);
                  MPI_Recv(&allreal[p*iters*SAMPLE_LEN],iters*SAMPLE_LEN,MPI_DOUBLE,p,tag+2,MPI_COMM_WORLD,&status);
                  MPI_Recv(&allimag[p*iters*SAMPLE_LEN],iters*SAMPLE_LEN,MPI_DOUBLE,p,tag+3,MPI_COMM_WORLD,&status);                  
                  MPI_Recv(&allback[p*iters*NOISE_LEN],iters*NOISE_LEN,MPI_DOUBLE,p,tag+4,MPI_COMM_WORLD,&status);
               }         
            }
        
         // Now broadcast matrices in process 0 to all other processes 
            MPI_Bcast(&allnoise[0],iters*numprocs*NOISE_LEN,MPI_DOUBLE,0,MPI_COMM_WORLD);    
            MPI_Bcast(&allfreqs[0],iters*numprocs*SAMPLE_LEN,MPI_DOUBLE,0,MPI_COMM_WORLD); 
            MPI_Bcast(&allreal[0],iters*numprocs*SAMPLE_LEN,MPI_DOUBLE,0,MPI_COMM_WORLD); 
            MPI_Bcast(&allimag[0],iters*numprocs*SAMPLE_LEN,MPI_DOUBLE,0,MPI_COMM_WORLD);           
            MPI_Bcast(&allback[0],iters*numprocs*NOISE_LEN,MPI_DOUBLE,0,MPI_COMM_WORLD); 
                      
  // Outputting results to file (process 0 only)        
     if (myid==0){
        sprintf(outFile,"%s_time_proc%d.csv",OUTFOLDER,myid);
        fptr=fopen(outFile,w);
        sprintf(outFile,"%s_freq_proc%d.csv",OUTFOLDER,myid);
        cptr=fopen(outFile,w);        

              for (p=0; p<5; p++){
                   for (i=0; i<iters; i++){
	                     fq  =  1.0/(10.0*((double)p+1.0));
	                     fprintf(fptr,"%g,%g,",fq,fq);
	                     fprintf(cptr,"P%g,R%g,I%g,",fq,fq,fq);
                   }
              }
              fprintf(fptr,"\n");
              fprintf(cptr,"\n");
              
              for (p=0; p<5; p++){
                   for (i=0; i<iters; i++){
	                     phi = (double)i * 0.2 - 0.2; 
	                    fprintf(fptr,"%g,%g,",phi,phi);
	                    fprintf(cptr,"P%g,R%g,I%g,",phi,phi,phi);
                   }
              }
              fprintf(fptr,"\n"); 
              fprintf(cptr,"\n");
              
         for (t=0; t<NOISE_LEN; t++){
              for (p=0; p<5; p++){
                   for (i=0; i<iters; i++){
                        fprintf(fptr,"%g,%g,",allnoise[(p*iters+i)*NOISE_LEN + t],
                                              allback[(p*iters+i)*NOISE_LEN + t]
                        );
                   }
              }
              fprintf(fptr,"\n");
         }

         for (t=0; t<SAMPLE_LEN/2+1; t++){
              fq=12.0*(double)(t)/(double)(SAMPLE_LEN);
              //fprintf(cptr,"%g,%g,",fq,1.0/fq);
              for (p=0; p<5; p++){
                   for (i=0; i<iters; i++){
                        fprintf(cptr,"%g,%g,%g,",allfreqs[(p*iters+i)*SAMPLE_LEN + t],
                                                 allreal[(p*iters+i)*SAMPLE_LEN + t],
                                                 allimag[(p*iters+i)*SAMPLE_LEN + t]
                        );
                   }
              }
              fprintf(cptr,"\n");
         }
         
         fclose(fptr);
         fclose(cptr);         
     }    
*/

  // General cleanup
      free_dvector(amp,0,SPEC_LEN);
      free_dvector(offset,0,SPEC_LEN);
      free_dvector(noise,0,NOISE_LEN*iters);
      free_dvector(back,0,NOISE_LEN*iters); 
      free_dvector(freq,0,SAMPLE_LEN*iters);                   
      free_dvector(real,0,SAMPLE_LEN*iters);
      free_dvector(imag,0,SAMPLE_LEN*iters);
      free_dvector(allnoise,0,iters*numprocs*NOISE_LEN);  
      free_dvector(allback,0,iters*numprocs*NOISE_LEN); 
      free_dvector(allfreqs,0,iters*numprocs*SAMPLE_LEN);           
      free_dvector(allreal,0,iters*numprocs*SAMPLE_LEN);
      free_dvector(allimag,0,iters*numprocs*SAMPLE_LEN);
      
  // De-allocate the memory
             free_SimRun(v);
             
#undef START_TIME
#undef NOISE_LEN
#undef NUM_SERIES
#undef BIG_T
#undef SPEC_LEN
#undef TWOPI
}

//-----------------------------------------------------------------------------------

int MPI_deriv(SimRun *v, double x[], double df[], int myid, int numprocs, int runswitch){

int i,sp, c, p;
double *basedf;
double old, newSSQ, ssq;
struct FitRun BaseRun;
struct FitRun DiffRun;
double baseval, delval;
double starttime, endtime;
double *df_list;
int    *df_proc, df_num;
int    tag,fun;
int    running;
int    *phonein;
  char outFile[80]; 
  FILE *fptr;
  
MPI_Status status;

    phonein   = ivector(0,numprocs);
    memset(phonein,0,(numprocs+1)*sizeof(int));

    basedf    = dvector(1,NDIM);

    df_list   = dvector(0,NDIM+1);
    df_proc   = ivector(0,NDIM+1);


    BaseRun.v = new_SimRun();
    BaseRun.x = dvector(1,NDIM);

    DiffRun.v  = new_SimRun();
    DiffRun.x  = dvector(1,NDIM);

running   = runswitch;
starttime = MPI_Wtime(); 

if (myid==0){cout << " Starting derivative " << runswitch << endl;}

//cout << myid << "_1_" << running << endl;

while(1){

    //cout << myid << "_2_" << running << endl;      
    MPI_Bcast(&x[1], NDIM,  MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    //cout << myid << "_3_" << running << endl; 
    MPI_Bcast(&running, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    //cout << myid << "_4_" << running << endl; 
    if (!running){break;}

    for (i=1; i<=NDIM; i++){
         BaseRun.x[i]=x[i];
         BaseRun.v->fit_vec_control[i] = v->fit_vec_control[i];
    }
    for (i=0; i<=4; i++){
        for (sp=1; sp<=NUM_LIVING; sp++){
            BaseRun.v->fit_ssq_control[i][sp] = v->fit_ssq_control[i][sp];
        }
    }   
    BaseRun.baseSSQ = SSQ_run(BaseRun.v, BaseRun.x);

    for (i=1; i<=NDIM; i++){
        DiffRun.x[i]=x[i];
        DiffRun.v->fit_vec_control[i] = v->fit_vec_control[i];
    }
    for (i=0; i<=4; i++){
        for (sp=1; sp<=NUM_LIVING; sp++){
            DiffRun.v->fit_ssq_control[i][sp] = v->fit_ssq_control[i][sp];
        }
    }    
   	DiffRun.baseSSQ = SSQ_run(DiffRun.v, DiffRun.x);

    //cout << myid << "_5_" << running << endl; 

    df_num    = 0;
    for (i=0; i<NDIM; i+=numprocs){
            //cout << myid << "_6_" << running << endl; 
        fun = numprocs - myid - 1;
        DiffRun.dterm = i+fun+1;
        if (DiffRun.dterm <=NDIM){
            // From deriv_thread
            if (DiffRun.v->fit_vec_control[DiffRun.dterm] > 0.0001) {
         
				      old = DiffRun.x[DiffRun.dterm];
              DiffRun.x[DiffRun.dterm] += DERIV_STEP;
				      newSSQ=SSQ_run(DiffRun.v, DiffRun.x);

			        // KYA 8/5/08: setting the derivative in the above line was summing 
			        // up the pieces with poor precision.  Now, the df is not calculated 
			        // here, but the pieces of the derivative are saved in the run stats.  
			        // This is returned to the main function to assemble the derivative 
			        // vector.  For now, we just set df to 0.0 (while keeping the parts 
			        // after the run).
              DiffRun.df = 0.0;
      
              DiffRun.x[DiffRun.dterm] = old;
            }
            else {
              DiffRun.df = 0.0;
            }   
        }
        
          // Now assemble the derivative pieces into a vector 
             /// CHANGED FOR MPI for (r=0; r<PRUNS; r++){
                 if (DiffRun.dterm <=NDIM){ 

                      // Save processor number for diagnostics
                         // rr[DiffRun.dterm]=myid;

                      // Now Load the derivative pieces into a single
                      // float vector DF, indexed by dterm
                         basedf[DiffRun.dterm] = 0.0;
                         
                      // recalculating SSQ here, does this help ??? MPI REMOVED
                         ssq=SSQ(DiffRun.v);
                      
                      // Add derivative parts if control vector (fitting ON) 
											// for term is positive, otherwise derivative part is 0
                         if (DiffRun.v->fit_vec_control[DiffRun.dterm] > 0.0001){
                            // Fitting parts
                               for (c=1; c<=fitN; c++){
                                   baseval = BaseRun.v->fit_SSQfit[c]; //TRUNC(baseval);
                                   delval  = DiffRun.v->fit_SSQfit[c]; //TRUNC(delval);
                                   
                                   basedf[DiffRun.dterm] += (delval-baseval);
                                   //dSdX[DiffRun.dterm][c] = (delval-baseval)/DERIV_STEP;

				    								   }
				    								// Forcing parts 
                               for (c=1; c<=NUM_LIVING+NUM_DEAD; c++){
                                   baseval = BaseRun.v->fit_SSQ[3][c]; //TRUNC(baseval);
                                   delval  = DiffRun.v->fit_SSQ[3][c]; //TRUNC(delval);                               
                                   basedf[DiffRun.dterm] += (delval-baseval);
                                   //dSdX[DiffRun.dterm][c+fitN] = (delval-baseval)/DERIV_STEP;

                               }
                            // Divide whole derivative by Step Size     
                               basedf[DiffRun.dterm] /= DERIV_STEP;
                         }
                         												  
                      // df[DiffRun.dterm] =(double)basedf[DiffRun.dterm];
                      // The derivative parts are saved in a list of part and
											// reference to step number; a practical way to gather
											// them up later                         
                         df_list[df_num] = basedf[DiffRun.dterm];
                         df_proc[df_num] = DiffRun.dterm;
                         df_num++;
                         
                                         
                 }  // END OF if dterm<NDIM loop
    }  // END OF i loop through variables in steps of numprocs 
    
    tag = 1;
        //cout << myid << "_7_" << running << endl; 
    if (myid>0){
       MPI_Send(&df_num,      1, MPI_INTEGER, 0, tag+0, MPI_COMM_WORLD);
       MPI_Send(df_list, df_num,  MPI_DOUBLE, 0, tag+1, MPI_COMM_WORLD); 
       MPI_Send(df_proc, df_num, MPI_INTEGER, 0, tag+2, MPI_COMM_WORLD);
       MPI_Send(&myid,        1, MPI_INTEGER, 0, tag+3, MPI_COMM_WORLD);
           //cout << myid << "_8_" << running << endl; 
    }              
    else{
        for (sp=0; sp<df_num; sp++){ df[df_proc[sp]] = df_list[sp]; }
        for (p=1; p<numprocs; p++){
             MPI_Recv(&df_num,      1, MPI_INTEGER, p, tag+0, MPI_COMM_WORLD, &status);
             MPI_Recv(df_list, df_num,  MPI_DOUBLE, p, tag+1, MPI_COMM_WORLD, &status); 
             MPI_Recv(df_proc, df_num, MPI_INTEGER, p, tag+2, MPI_COMM_WORLD, &status);
             for (sp=0; sp<df_num; sp++){df[df_proc[sp]] = df_list[sp];}
						 MPI_Recv(&phonein[p],  1, MPI_INTEGER, p, tag+3, MPI_COMM_WORLD, &status); 
						     //cout << myid << "_9_" << p << "_" << running << endl;         
				}
        //cout << myid << "_9a_" << running << endl; 
		    sprintf(outFile,"%s_PartialTest_%d_%d.csv",OUTFOLDER,myid,runswitch);
        //cout << myid << "_9b_" << running << endl;
        fptr=fopen(outFile,"w");
        //cout << myid << "_9c_" << running << endl;
        for (i=1; i<=NDIM; i++){fprintf(fptr,"%d,%g\n",i,df[i]);}
        //cout << myid << "_9d_" << running << endl;
        fclose(fptr);
        
        break;				
								          
    }

    //cout << myid << "_10_" << running << endl; 
    // Process 0 now broadcasts the new vector and the 'running' status to the
    // other processes.  Broadcasting a status of 0 ends the loop
    //MPI_Bcast(&x[1], NDIM, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    //MPI_Bcast(&running, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
		//cout << myid << " is on " << running << endl;    

}

   // Clean-up of memory
      free_dvector(df_list,0,NDIM+1);
      free_ivector(df_proc,0,NDIM+1);
      
      free_dvector(basedf,1,NDIM);

      free_dvector(DiffRun.x,1,NDIM);
      free_SimRun(DiffRun.v);
           
      free_SimRun(BaseRun.v);
      free_dvector(BaseRun.x,1,NDIM);
  		//cout << myid << " is retired " << endl;
  		//errflag=0;
			//if (myid==0){
			   //cout << "pieces reporting: "; 
			//   for(p=1; p<numprocs; p++) {cout << phonein[p] << " ";}
			   //cout << endl;
			//} 
      free_ivector(phonein,0,numprocs);

			endtime = MPI_Wtime(); 
			if (myid==0){cout << "deriv " << runswitch << " time " << endtime-starttime << endl;}
			if (myid==1){cout << "final deriv time " << endtime-starttime << endl;}
			//if (myid>0){cout << myid << " escape 1" << endl;} 
}

//-----------------------------------------------------------------------------------
// The FOLLOWING MINIMIZATION ROUTINES are taken from Numerical Recepies in C
// without comment.  They are cut and paste except for a couple of print statements
#define GTOL      0.01                   // Min %change of a parameter at final solution step
#define EPS       3.0e-8
#define TOLX      (4*EPS)
#define STPMX     1.0          //past runs at 0.1, recent at 0.5
#define ALF       1.0e-4
#define BIGSSQ    1000000.0
#define DERIV_STEP 0.1
// -----------------------------------------------------------------------------   
void MPI_lnsrch(struct SimRun *v, int n, double xold[], double fold, double g[], double p[], double x[], double *f, double stpmax, int *check)
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
		*f=SSQ_run(v,x);
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


//---------------------------------------------------------------


#define FREEALL MPI_deriv(v,p,g,myid,numprocs,0); free_dvector(xi,1,n);free_dvector(pnew,1,n); \
free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n); \
free_dvector(dg,1,n);

void MPI_dfpmin(struct SimRun *v, double p[], int n, double gtol, int *iter, double *fret, int myid, int numprocs)
{
//	void lnsrch(struct SimRun *v, int n, double xold[], double fold, double g[], double p[], double x[],
//		 double *f, double stpmax, int *check, double (*func)(struct SimRun*, double []));
	int check,i,its,j,r;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double *dg,*g,*hdg,**hessin,*pnew,*xi;
  int MPI_GO;
  // KYA Added
	double sss;  char outFile[80]; char ilabel[80];
  FILE *fptr;
  
  if (ITMAX == 0){return;}

     BigIt=0;

      //sprintf(outFile,"%s_ForcePartials.csv",OUTFOLDER);
      //tptr=fopen(outFile,"w");
      //
      //sprintf(outFile,"%s_fitDerivs.csv",OUTFOLDER);
      //fptr=fopen(outFile,"w");
      //fprintf(fptr,"Step,Species,var,x,Control,dSdx,");
      //for (i=1; i<=fitN; i++){
      //    fprintf(fptr,"Fit_%s_%d,",path_species[fit_ApplyTo[i]],fit_type[i]);
     // }
      //for (i=1; i<=NUM_LIVING+NUM_DEAD; i++){
      //    if(FORCE_BIO_FLAG[i]){
      //       fprintf(fptr,"Force%s,",path_species[i]);
      //    }
      //}      
      //fprintf(fptr,"EndTot,\n");
    

	//printf ("%d first\n",myid);

   g  = dvector(1,n);

   MPI_deriv(v,p,g,myid,numprocs,1);
   //if (myid>0){cout << myid << " escape 2" << endl;} 
   if (myid!=0){
     //if (myid>0){cout << myid << " escape 3" << endl;} 
     free_dvector(g,1,n); 
     //if (myid>0){cout << myid << " escape 4" << endl;} 
     return;
   } // No processes other than 0 go beyond here	
  
  //printf ("%d second\n",myid);
	dg     = dvector(1,n);
	hdg    = dvector(1,n);
	hessin = dmatrix(1,n,1,n);
	pnew   = dvector(1,n);
	xi     = dvector(1,n);

  fp=SSQ_run(v,p);
	
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
		
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		MPI_lnsrch(v,n,p,fp,g,xi,pnew,fret,stpmax,&check);
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
    	 cout << "Iter: " << its << " SSQ: " << sss << endl;
	     sprintf(outFile,"RunFit_iter%d",its);
       output_run(v,outFile);
       //outstep(its,p,v); /*bin_outstep(its,p,v);*/
       BigIt=its;
		   //if (its ==2){TESTFLAG=1;}else{TESTFLAG=0;}
		//exit(0);

		if (test < TOLX) {
    	 printf("TOLX exit\n");
	     //outderiv(v,fptr, p, g, ilabel);
			FREEALL
			return;
		}

		for (i=1;i<=n;i++) dg[i]=g[i];
       MPI_deriv(v,p,g,myid,numprocs,its+1);
		   //sprintf(ilabel,"%d",its);
	     //outderiv(v,fptr, p, g, ilabel);
		
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
    	 printf("gtol exit\n");
	     //outderiv(v,fptr, p, g, ilabel);
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
	//nrerror("too many iterations in dfpmin");
  printf("too many exit\n");
	     //outderiv(v,fptr, p, g, ilabel);
	
	FREEALL
}
//#undef ITMAX
//#undef EPS
//#undef TOLX


//#define TOLX 1.0e-7


//-----------------------------------------------------------------------------------
int MPI_solve_dfp(struct SimRun *v, double *p, int its, int myid, int numprocs)
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
  	 //dSdX=dmatrix(1,NDIM,1,fitN+NUM_LIVING+NUM_DEAD); // full partial derivative matrix for sensitivity analyses
	   
	// Initialize P vector to 0 (can use other starting points if desired) 
     //for (i=1; i<=NDIM; i++){p[i]=0.0;}
     
  // Set starting biomass to first biomass in time series
     //for (sp=1; sp<=NUM_LIVING; sp++){
     //     i=NUM_LIVING*1 + sp;
     //     p[i] = fit_bioStart[sp];
     //     cout << i << " " << path_species[sp] << " " << p[i] << endl;
     //}
  
  // Do a first run after starting vector is applied
     //start1=time(NULL);
     if (myid == 0){
        fret = SSQ_run(v,p);
		    printf ("Start  SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
        sprintf(outFile,"_iter0");
        output_run(v, outFile);
        outstep(0,p,v); /*bin_outstep(0,p,v);*/
      }
     //end1=time(NULL);
     //cout << "One run in: " << (float)(end1-start1)/60. << " Min." << endl;
     //cout << "Est. step time: " << ((float)(end1-start1)/60.) * NDIM / PRUNS << " Min." << endl;
     
     //start1=time(NULL);
  // Now call the general minimization routine from numerical recepies
  // SSQ_run and SSQ_deriv are the functions that calculate SSQ and dSSQ/dP
     ITMAX = its;
     MPI_dfpmin(v,p,NDIM,GTOL,&iter,&fret, myid, numprocs);
     //frprmn(v,p,NDIM,GTOL,&iter,&fret, SSQ_run, SSQ_deriv);
     //if (myid>0){cout << myid << " escape 5" << endl;} 
	   //printf("Iterations: %3d\n",iter);
	   if (myid==0){
	   printf("Final  SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
		  
  // Now call the SSQ and derivative routine a final time to get diagnostics
  // (e.g. steepness of SSQ surface at final solution)
     
		 //fret = SSQ_run(v,p);
     //printf("Intermediate test SSQ %14.6g, \n",fret);
     //sprintf(outFile,"_iter999");
     //output_run(v, outFile);
     //outstep(999,p,v); /*bin_outstep(999,p,v);*/
     //SSQ_deriv(v,p,dp);
     
     fret = SSQ_run(v,p);
     printf("Last   SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		         fret,v->fit_fittot, v->fit_forcetot, v->fit_diettot);
     sprintf(outFile,"_iter1000");
     output_run(v, outFile);
     outstep(1000,p,v); /*bin_outstep(1000,p,v);*/
     }
     
  // Cleanup
       //if (myid>0){cout << myid << " escape 6" << endl;} 
     //free_dmatrix(dSdX,1,NDIM,1,fitN+NUM_LIVING+NUM_DEAD);
     //free_vector(p,1,NDIM);
     free_dvector(dp,1,NDIM);
     //end1=time(NULL);
     //cout << "solved in: " << (float)(end1-start1)/60. << " Min." << endl; 
       //if (myid>0){cout << myid << " escape 7" << endl;} 
	return 0;
}

#undef GTOL                         // Min %change of a parameter at final solution step
#undef EPS       
#undef TOLX      
#undef STPMX               //past runs at 0.1, recent at 0.5
#undef ALF       
#undef BIGSSQ    
#undef DERIV_STEP 
//-----------------------------------------------------------------------------------
int MPI_simple_fitted(int myid, int numprocs)
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
      		   FLOWOUT = 0; // Setting flowout to 1 dumps flow files for each year during the base run
					   Adams_Basforth(v, 0, fit_Years);
             FLOWOUT = 0;
         //  Load starting fit vector (or set to 0 if none existing) 
             load_instep(v);  
             ssq = SSQ_run(v, v->fit_vector);
          // If it didn't crash in first run, calculate SSQ and save the system
             printf ("%d Load   SSQ:%9.2f, fit:%9.2f force:%9.2f diet:%9.2f\n",
		                 myid, ssq,v->fit_fittot, v->fit_forcetot, v->fit_diettot);

         // MAIN SOLVER ROUTINE
             MPI_solve_dfp(v,v->fit_vector,200,myid,numprocs);
             //if (myid>0){cout << myid << " escape 8" << endl;}              
        //   Now do a final run with new rates
				   if (myid==0){           
             ssq = SSQ_run(v, v->fit_vector);
             if (!v->DISCARD_YEAR){
                 sprintf(outFile,"_fit");
                 output_run(v, outFile);
                 rate_output(v);
             }
             cout << myid << " final output run with SSQ: " << ssq << endl;
           }
       // Cleanup
          free_SimRun(v);
						 
return 0;
}


//-----------------------------------------------------------------------------------

int MPI_test(int argc, char *argv[], char submode){

int myid, numprocs;
double starttime, endtime;

   struct SimRun *Runs;

    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 
    
    //fftw_test();
    //return(0);
    
    starttime = MPI_Wtime();
    
    switch(submode){ 
    
    case 'a':
    // CREATING a SERIES
       Runs=new_SimRun();
               //cout << "and 1" << endl;
               //for (i=0; i<=NDIM; i++){Runs[r]->fit_vector[i]=0.0;}
       load_instep(Runs);
               //cout << "and 2" << endl;
            // Put seed in random number generator
            // Seed should be a number between 0 and 511, higher numbers will repeat seeds
               //init_by_array(rseed+r, &Runs[r].rng); 
       init_by_array(rseed+myid, &Runs->rng); 
            // Run ID number for output files
               //Runs[r].thread_ID=rseed+r;
          
       printf("%d %f %f %f %f\n",myid,uniform(&Runs->rng),uniform(&Runs->rng),uniform(&Runs->rng),uniform(&Runs->rng));
       Runs->thread_ID=rseed+myid;
            // Split the run onto a new thread and run it!
               //rc = pthread_create(&threads[r],&attr, save_a_series, (void *) &Runs[r]);
               //cout << "and 3" << endl;
       save_a_series((void *) Runs);        

       free_SimRun(Runs);
       
       printf("%d %f",myid,endtime-starttime);

       break;
    
		case 'f':
		   MPI_simple_fitted(myid, numprocs);
		   break;
		   
    case 'b':
       MPI_noise(myid,numprocs);   
       
       break;
       
    default:
       break;
    } // end submode switch
  
    endtime = MPI_Wtime();
    MPI_Finalize(); 
    if (myid==0){cout << "parallel " << myid <<" done in " << endtime-starttime<< " seconds!" << endl;}
    return 0; 

}
