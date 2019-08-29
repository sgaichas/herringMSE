//
// ECOSIM_THREADS by KERIM AYDIN (Kerim.Aydin@noaa.gov)
//
// Project started out as "pure C" version of Ecosim, but then added C++ for 
// convenient IO libraries and string handling such as reading .csv files.
// Other than IO, "math" should be almost entirely C, with C data structures.
//
// ECOSYSTEM-SPECIFIC variables:  Some ecosystem-specific variables are
// input at compile time, including the Ecopath data and the number of
// processors for the machine this will be run on.  These, and other
// aspects of the numerics that the user may want to change can be found
// in ecosim_UserInputs.h
//
// Project requires a "complete" Ecopath model as input, formatted as a
// set of C input arrays in a header file.  Visual Basic routines for
// converting a "standard" Ecopath .eii file to such a header file are
// available from the author.  
//
// Note that a lot of the data structures in this software are particularly
// complex due to explicit parallel programming/threading for running on
// machines with multiple processors.  A single-threaded version is much
// simpler (but not under development right now).

// CODE STUCTURE       This main project consists of 6 files:
//
// ecosim_UserInputs.h Contains some ecosystem-specific variables that
//                     the user must set before compiling.
// ecosim_threads.h:   Contains the data structures, function definitions,
//                     and global variables.  Also has some basic full functions
//                     for allocating and freeing structures.
// ecosim_threads.cpp: Contains most of the basics of running ecosim,
//                     function main, and many controls.
// simio.cpp:          Contains input/output routines including command
//                     line switches, reading .csv files, and writing
//                     outputs.
// fitfunctions.cpp:   Contains the routines for fitting ecosim runs to
//                     data.                       
// All of the above are sewn together with #include statements in ecosim_threads.cpp.

// In addition, there are three stand-alone sets of functions which are compiled
// separately and then linked (all from open source code with copyrights indicated
// in the files):
// 1.  gauss_threads.h and gauss_threads.cpp are routines for generating particularly
//     strong pseudo-random numbers in a thread-safe manner (e.g separate threads
//     can have separate pseudo-random sequences without interferences on state).
// 2.  csvparser.h and csvparser.cpp implements a simple C++ class for reading
//     csv files.
// 3.  nrutil.h and nrutil.cpp contains a minimal set of functions from
//     Numerical Recepies in C (Press et al. 2nd edition), primarily the
//     c vector and matrix types.
       
// C header files (required standard libraries)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <pthread.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <sys/file.h>
#include <dirent.h>

#ifdef MPI_COMPILE
#include "mpi.h"
#endif

// C++ header files (required standard libraries)
#include <iostream>
#include <fstream>
#include <string>
#include <map>

// Local includes that are standalone packages
#include "csvparser.h"           // Package for reading in CSV files
#include "gauss_threads.h"       // Random Number generators
#define NRANSI                   // Numerical recipies needs
#include "nrutil.h"              // Numerical recipies support functions
//#include "charts.h"              // Chart Drawing software

using namespace std;

// DATA, FORCING AND FITTING INPUTS
#include "ecosim_UserInputs.h"

// Defaults for some parameters that can be changed by command line arguments
   int SYSTEMS_TO_RUN  = 0;
   int BURN_TIME       = 0;
   float NOISE_RANGE   = 0.0;
   float VULVAR       = 6.0;
   char RUNMODE        = '\0';
   char SUBMODE        = '\0';       
   unsigned long rseed = 0;

   float AREA_MULT     = 1.0;
   // Multiplicative way of specifying func. resp.
   //float MDULL          = 2.0;
   //float MSINGLEHANDLE  = 1000.0;

   // Currently used
   float MSCRAMBLE      = 2.0;
   float MHANDLE        = 1000.0;

   // Switching Power 
   //int PREDSWITCH        = 0;
   float PREYSWITCH        = 1.0;

   // 1.0 = no overlap, 0.0 = complete overlap
   float ScrambleSelfWt = 1.0;
   float HandleSelfWt   = 1.0;

   //float NOISE_SCALE   = 0.0;
   int REPEAT            = 0;
   int LONG_CYCLES       = 6;
	 // Some fixed definitions, shouldn't need to tweak
#define STEPS_PER_YEAR     12                   // Steps per year between stanzas ("months")
#define MIN_REC_FACTOR     6.906754779          // This is ln(1/0.001 - 1) used to set min. logistic matuirty to 0.001
#define MIN_THRESHOLD(x,y) if ((x)<(y))(x)=(y)  // Set x to min y if too small
#define MAX_THRESHOLD(x,y) if ((x)>(y))(x)=(y)  // Set x to max y if too big
#define DELTA_T            1.0/(STEPS_PER_YEAR*STEPS_PER_MONTH) // timestep in years
#define SORWT              0.5                  // Fast Equilibrium smoother
#define EPSILON            1E-8                 // Test threshold for "too close to zero"
#define BIGNUM             1E+8                 // Test threshold for "too big"

// The following variables must be FIXED for all ecosystems.  In other words,
// they should be loaded once and then not changed.  These should be similar
// to AD-model-builder input statics (as opposed to parameters)

// File names
     char OUTFOLDER[80];
     char INFOLDER[80];
     char CLIMATE_FILE[80];
     char JUV_FILE[80];
     char FITTING_FILE[80];
     char DIET_FILE[80];
     char FIT_VECTOR_FILE[80];
     char FIT_CONTROL_FILE[80];
     
// TIME SERIES FORCING FUNCTIONS
   // Force prey to be consumed more (including prey 0 for PrimProd))
      float force_byprey[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      float force_bymort[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      float force_byrecs[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      // add "mediation by search rate" KYA 7/8/08
         float force_bysearch[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
   // Variables prefaced by fish_ are Fishing related

// Guild variables
       string PathNames[NUM_GROUPS+1];
       map<string,string> guildlist;
       map<string,string> trophiclist;
       map<string,string> dattype;
       map<string,int>    spcolor;
       map<string,int>    guildcount;
       map<string,int>    guildnum;
       int NGuilds;
       string GuildOrder[NUM_GROUPS];
       
// Global storage NOT THREAD SAFE!!!
   //float flowMat[NUM_GROUPS+1][NUM_GROUPS+1];
   //float bioMat[NUM_GROUPS+1][MAX_YEARS+1];
   //float catMat[NUM_GROUPS+1][MAX_YEARS+1];
   int   FLOWOUT = 0;
   int   AGEOUT  = 0;
   int   GIANT   = 0;

// For SENSE ROUTINES
   int   fit_StartYear;
   int   fit_Years = MAX_YEARS;
   // Public FitYears&, FitMonths&
   int fitN;
   int fit_lastyear;
   string fit_names[NUM_LIVING+MAX_COLS+1];
   int    fit_type[NUM_LIVING+MAX_COLS+1];
   int    fit_ApplyTo[NUM_LIVING+MAX_COLS+1];
   int    fit_ApplyFrom[NUM_LIVING+MAX_COLS+1];
   int    fit_pedigree[NUM_LIVING+MAX_COLS+1];
   float  fit_Norm[NUM_LIVING+MAX_COLS+1];
   float  fit_SD[NUM_LIVING+MAX_COLS+1][MAX_YEARS+1];
   float  fit_OBS[NUM_LIVING+MAX_COLS+1][MAX_YEARS+1];
   float  fit_RAW[NUM_LIVING+MAX_COLS+1][MAX_YEARS+1];

   float  fit_bioStart[NUM_GROUPS+1];  
   
   // DIET fitting 12/4/07
      //fit_diet_tofrom[NUM_PREDPREYLINKS];
      int dietN;
      string fit_diet_names[NUM_LIVING+NUM_PREDPREYLINKS+1];
      int fit_diet_type[NUM_LIVING+NUM_PREDPREYLINKS+1];
      int fit_diet_lookup[NUM_LIVING+NUM_PREDPREYLINKS+1];
      double fit_diet_OBS[NUM_LIVING+NUM_PREDPREYLINKS+1][MAX_YEARS+1];
         
   //int bioN, bioApplyTo[MAX_COLS+1]; 
	 //float bioSSQpart[MAX_COLS+1];
   //float bioOBS[MAX_COLS+1][MAX_YEARS+1], bioSD[MAX_COLS+1][MAX_YEARS+1], 
   //      bioEST[MAX_COLS+1][MAX_YEARS+1]; 

   //int catN, catApplyTo[MAX_COLS+1]; 
	 //float catSSQpart[MAX_COLS+1];
   //float catOBS[MAX_COLS+1][MAX_YEARS+1], catSD[MAX_COLS+1][MAX_YEARS+1], 
   //      catEST[MAX_COLS+1][MAX_YEARS+1];
				          
   //float FORCED_EFFORT[NUM_GROUPS+1][MAX_YEARS+1];
   float FORCED_CATCH[NUM_GROUPS+1][MAX_YEARS+1];
   float FORCED_BIO[NUM_GROUPS+1][MAX_YEARS+1]; 
   float FORCED_BIO_RAW[NUM_GROUPS+1][MAX_YEARS+1]; 
   float FORCED_CATCH_RAW[NUM_GROUPS+1][MAX_YEARS+1];
   int   FORCE_BIO_FLAG[NUM_GROUPS+1];
	 float FORCED_FRATE[NUM_GROUPS+1][MAX_YEARS+1];
	 float FORCED_EFFORT[NUM_GROUPS+1][MAX_YEARS+1];  // Added 7/8/08 KYA
	 int FORCED_TARGET[NUM_GROUPS+1];
	 int FORCED_FTARGET[NUM_GROUPS+1];
	 float FORCED_REC[MAX_SPLIT+1][MAX_YEARS+1];
	 int    FORCE_PEDIGREE[NUM_GROUPS+1];
   // Public RecOBS!(), RecWTS!(), RecEST!(), RecApplyTo&(), recSSQparts!(), recSSQ!, recN&, recRevLookup!()
   // Public CatOBS!(), CatWTS!(), CatEST!(), CatApplyTo&(), catSSQparts!(), catSSQ!, catN&, catRevLookup!()
   // Public FspOBS!(), FspWTS!(), FspEST!(), FspApplyTo&(), FspSSQparts!(), FspSSQ!, FspN&, FspRevLookup!()
   // Public NumParams&
   // Public FORCED_BIO!(), FORCED_RECRUITS!(), FORCED_CATCH!(), NFORCEDCATCH&, FORCED_LINKS&(), START_CATCH!(), START_F!()
   float TL[NUM_GROUPS+1];
   
// Variables prefaced by juv_ are juvenile/adult input parameters
   int juv_N;
   string juv_Names[NUM_GROUPS+1];
   int juv_JuvNum[NUM_GROUPS+1];
	 int juv_AduNum[NUM_GROUPS+1];
	 int juv_RecAge[NUM_GROUPS+1];
	 float juv_VonBK[NUM_GROUPS+1];
	 float juv_RecPower[NUM_GROUPS+1];
	 float juv_aduEqAgeZ[NUM_GROUPS+1];
	 float juv_juvEqAgeZ[NUM_GROUPS+1];
	 float juv_Wmat[NUM_GROUPS+1];
	 float juv_Winf[NUM_GROUPS+1];
   float juv_VonBD[NUM_GROUPS+1];
   int juv_RecMonth[NUM_GROUPS+1];
    // REC_CHANGE
       //float WmatWinf[MAX_SPLIT+1];
       float Wmat50[MAX_SPLIT+1];
       float Wmat001[MAX_SPLIT+1];
       float WmatSpread[MAX_SPLIT+1];
       float Amat50[MAX_SPLIT+1];
       float Amat001[MAX_SPLIT+1];
       float AmatSpread[MAX_SPLIT+1];
    // END REC_CHANGE

// The RatePar structure contains all the parameters which may be fitted, 
// randomly generated, or otherwise changed to create new ecosystems.
// WARNING WARNING.  IF YOU CHANGE THE ORDER OF THESE OR ADD NEW ONES,
// ALL OLD SAVED ECOSYSTEMS WON'T LOAD ANYMORE
struct RatePar{
   // Variables prefaced by rpar_ are Ecosim rate paramters (for example, 
   // mortalities) calculated from Ecopath (for Sense, sometimes with noise).
  // or needed for sim running (vuls, feedtimeadjust)
     // Juvenile Numbers generated randomly
        float rpar_drawn_K[MAX_SPLIT+1];
        float rpar_drawn_AduZ[MAX_SPLIT+1];
        float rpar_drawn_JuvZ[MAX_SPLIT+1];
        float rpar_SpawnEnergy[MAX_SPLIT+1];
        float rpar_SpawnX[MAX_SPLIT+1];
        float rpar_SpawnAllocR[MAX_SPLIT+1];
        float rpar_SpawnAllocG[MAX_SPLIT+1];
     // Integration method
        int rpar_NoIntegrate[NUM_GROUPS+1];
     // Predator prey links - number of them, and from (prey) and to (predator)
        int rpar_NumPredPreyLinks; // Long
        int rpar_PreyFrom[NUM_PREDPREYLINKS+1]; // Long
        int rpar_PreyTo[NUM_PREDPREYLINKS+1]; // Long
        float rpar_QQ[NUM_PREDPREYLINKS+1];
        //float rpar_XX[NUM_PREDPREYLINKS+1]; // Single
        //float rpar_HH[NUM_PREDPREYLINKS+1];
        float rpar_VV[NUM_PREDPREYLINKS+1];
        float rpar_DD[NUM_PREDPREYLINKS+1];
        float rpar_HandleSelf[NUM_GROUPS+1];
        float rpar_ScrambleSelf[NUM_GROUPS+1];
        float rpar_HandleSwitch[NUM_PREDPREYLINKS+1];
        // rpar_PredAval(); // Single
        // rpar_PredVval(); // Single\
        // rpar_ZZ(); // Single
        // rpar_QX(); // Single
        // rpar_Xminus1(); // Single
        // rpar_PredPrey_ReverseLookup(); // Long
        // 'Primary production links
        // rpar_NumPrimProdLinks; // Long
        // rpar_PrimProdTo(); // Long
     // Detrital links
     // UNUSED???
        //int rpar_NumDetLinks; // Long
        //int rpar_DetFrom[NUM_DEAD*NUM_GROUPS+1]; // Long
        //int rpar_DetTo[NUM_DEAD*NUM_GROUPS+1]; // Long
     // Fishing links
        int rpar_NumFishingLinks; // Long
        int rpar_FishingFrom[NUM_CATCHLINKS*NUM_DEAD+1]; // Long
        int rpar_FishingThrough[NUM_CATCHLINKS*NUM_DEAD+1]; // Long
        int rpar_FishingTo[NUM_CATCHLINKS*NUM_DEAD+1]; // Long
        float rpar_FishingQ[NUM_CATCHLINKS*NUM_DEAD+1]; // Single
     // BASICS (Bioenergetics, foraging) (one for each living species)
        float rpar_B_BaseRef[NUM_GROUPS+1]; // Single
        float rpar_PBopt[NUM_GROUPS+1];
        // rpar_QB_BaseRef(); // Single
        float rpar_UnassimRespFrac[NUM_GROUPS+1]; // Single
        float rpar_ActiveRespFrac[NUM_GROUPS+1]; // Single
        // float rpar_PassiveRespMort[NUM_GROUPS+1]; // Single
        float rpar_MzeroMort[NUM_GROUPS+1]; // Single
        // float rpar_BioAcc[NUM_GROUPS+1]; // Single
        // float rpar_Bioen_A[NUM_GROUPS+1]; // Single
        // float rpar_Bioen_K[NUM_GROUPS+1]; // Single
        // float rpar_Bioen_Z[NUM_GROUPS+1]; // Single
        // float rpar_BioenBaseP[NUM_GROUPS+1]; // Single
        // float rpar_BioenPredX[NUM_GROUPS+1]; // Single
        // float rpar_BioenSatur[NUM_GROUPS+1]; // Single
        // float rpar_BioenUnassR[NUM_GROUPS+1]; // Single
        // float rpar_BioenActR[NUM_GROUPS+1]; // Single
        float rpar_FtimeAdj[NUM_GROUPS+1]; // Single
        // float rpar_Handling(); // Single
        float rpar_FtimeQBOpt[NUM_GROUPS+1]; // Single
     // 'Used for Scramble and Handling Time
        float rpar_PredTotWeight[NUM_GROUPS+1]; // Single
        float rpar_PreyTotWeight[NUM_GROUPS+1]; // Single
        float rpar_PredPredWeight[NUM_PREDPREYLINKS+1]; // Single
        float rpar_PreyPreyWeight[NUM_PREDPREYLINKS+1]; // Single
     // rpar_FishRateBase(); // Single
     // rpar_CatchDestFrac(); // Single
};

// Variables prefaced by state_ represent the current state of the system
// (e.g. biomass, numbers, current feeding time)
//  these are the stanza parameters and vectors from splitNo to survByMonth:
//  ratio of weight at maturity to 
//  weight at infinity, vonBertMonthly, vonBertannual, baserecruits,
//  baserec with BA, Carl's thing, used in sim
//  vector of relative weight at age by month, splitWage^2/3   
//  vector of monthly survival rate, monthly numbers
//  WageS and NageS are initialized equal to splitWage and splitNo, then 
//  change in sim, used to update eggs (recruitment)

// a SimRum contains a vector of rate parameters (split out for saving)
// and then the state of the system at any given moment, and the derivatives.
struct SimRun{

    struct randvec rng;   // The random number state holder for this run
    int    thread_ID;     // An ID number for convenience
    char   RunID[80];
    FILE  *dump;          // A File to dump results to
    struct RatePar RRR;   // The RATE PARAMETERS for the run
    
    //float fit_EST[MAX_COLS+1][MAX_YEARS+1];
    double *fit_vector;

    double **fit_EST;
    double **fit_diet_EST; //[NUM_PREDPREYLINKS+1][MAX_YEARS+1];
    //float fit_SSQ[MAX_COLS+1];
    double **fit_SSQ;
    double *fit_diet_SSQ;		
    double **fit_bioanom;

    double  *fit_vec_control;
    double **fit_ssq_control;
    //float fit_bioanom[NUM_GROUPS+1][MAX_YEARS+1];
    //float 
    //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
    double **out_BB;
    double **monthly_BB;
    double **monthly_PP;
    
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
    float **out_LinksF; //out_LinksM[NUM_PREDPREYLINKS+1];
    float **out_TotF;   //out_TotF[NUM_GROUPS+1];
    float **out_Effort;
    float out_PP;

    double **force_bymort; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];

    double  fit_SSQraw[NUM_GROUPS+MAX_COLS+1];
    double  fit_SSQfit[NUM_GROUPS+MAX_COLS+1];
    double  fit_SSQwt[NUM_GROUPS+MAX_COLS+1];
    double  fit_FORCEraw[NUM_GROUPS+1];    

    double fit_fittot, fit_diettot, fit_forcetot;
    
    double state_BB[NUM_GROUPS+1];
    double state_Ftime[NUM_GROUPS+1];
    //float NageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    //float WageS[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];  
    double **NageS, **WageS;
    double state_NN[NUM_GROUPS+1];
    double stanzaPred[NUM_GROUPS+1];
    double stanzaBasePred[NUM_GROUPS+1];
    
		double vBM[MAX_SPLIT+1];
    double recruits[MAX_SPLIT+1];
		double RzeroS[MAX_SPLIT+1];
    //float WWa[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    double **WWa;
    double stanzaGGJuv[MAX_SPLIT+1];
    double stanzaGGAdu[MAX_SPLIT+1];
    double baseEggsStanza[MAX_SPLIT+1];
    double EggsStanza[MAX_SPLIT+1];
    double SpawnBio[MAX_SPLIT+1];
    double baseSpawnBio[MAX_SPLIT+1];
    double Rbase[MAX_SPLIT+1];
    
    //float baseSpawnEnergy[MAX_SPLIT+1];
    
//  hardcoding only two stanzas for now; these dimension them by months
    int firstMoJuv[MAX_SPLIT+1], lastMoJuv[MAX_SPLIT+1]; 
    int firstMoAdu[MAX_SPLIT+1], lastMoAdu[MAX_SPLIT+1];
    //float SplitAlpha[MAX_SPLIT+1][MAX_MONTHS_STANZA+1];
    double **SplitAlpha;
    double RscaleSplit[MAX_SPLIT+1];  //not used right now

// OUTPUTS
   //float out_BB[NUM_GROUPS+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_BB[NUM_GROUPS+1][MAX_YEARS+1];
   //float out_nn[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_ww[MAX_MONTHS_STANZA+1][MAX_YEARS*STEPS_PER_YEAR+1];
   //float out_CC[NUM_GROUPS+1][MAX_YEARS+1];
   
   double fish_Effort[NUM_GROUPS+1];
   double TerminalF[NUM_GROUPS+1];

// FOR RECORDING SENSE RESULTS
   int DISCARD_YEAR;
	 int discarded[NUM_GROUPS+1];
	 //float month_noise[NUM_GROUPS + 1];

// Variables prefaced by deriv_ are derivative parts (e.g. current consumption)
   double deriv_TotGain[NUM_GROUPS + 1];
   double deriv_TotLoss[NUM_GROUPS + 1];
   double deriv_LossPropToB[NUM_GROUPS + 1];
   double deriv_LossPropToQ[NUM_GROUPS + 1];
   double deriv_ConsMat[NUM_GROUPS + 1][NUM_GROUPS + 1];
   double deriv_DerivT[NUM_GROUPS + 1];
   double deriv_dyt[NUM_GROUPS + 1];
   double deriv_biomeq[NUM_GROUPS + 1];

   double deriv_FoodGain[NUM_GROUPS + 1];
   double deriv_DetritalGain[NUM_GROUPS + 1];
   double deriv_FoodLoss[NUM_GROUPS + 1];
   double deriv_UnAssimLoss[NUM_GROUPS + 1];
   double deriv_ActiveRespLoss[NUM_GROUPS + 1];
   //float deriv_PassiveRespLoss[NUM_GROUPS + 1];
   double deriv_MzeroLoss[NUM_GROUPS + 1];
   double deriv_DetritalLoss[NUM_GROUPS + 1];
   double deriv_TotDetOut[NUM_GROUPS + 1];
   
   double deriv_FishingLoss[NUM_GROUPS + 1];
   double deriv_FishingGain[NUM_GROUPS + 1];
   double deriv_FishingThru[NUM_GROUPS + 1];
};  

struct FitRun{
   struct SimRun *v;
   double  baseSSQ;
   double  *x;
   //float  *dSdx;
   double  df;
   int    dterm;
};

struct ScatterRun{
   struct SimRun *v;
   double  baseSSQ;
   double *SSQlist;
   double  **Xlist;
   double  **Slist;
   double  *x;
};

// List of functions in ecosim.cpp (excluding main)

// SimIO routines
   int read_climate();
   int read_fitting();
   int read_juveniles();
   int read_diets();
   int output_run(struct SimRun *v, char* fbase);

// Run controls through thread management
   int simple_fitted(void);
   int scramble_fitted(void);
   int point_derivative(void);
   // KYA Depreciated in favor of random_fromfitted Jume 3 2008 
      // int random_ecosystem_series(void);
   int loaded_ecosystem_series(void);
   int perturb_loaded_series(void);
   int random_fromfitted_series(void);
   void *save_a_series(void *threadarg);
   void *load_a_series(void *threadarg);
   void *load_and_perturb_series(void *threadarg);
   //void *random_from_fit_a_series(void *threadarg);
   
// Initialization routines
   void load_system(struct SimRun *v, struct RatePar *r);
   void save_system(struct RatePar *s, struct RatePar *r);
   int path_to_rates(struct SimRun *v);
   int base_system(struct SimRun *v);
   int random_system(struct SimRun *v, double *newvec);
   int initialize_stanzas(struct SimRun *v); 
   int long_run(void);
   int one_run(void);   
   int deriv_field(struct SimRun *v, double x[]);
   void link_dump(struct SimRun *v); 
// Run routines
   int Adams_Basforth(struct SimRun *v, int StartYear, int EndYear);
   int SplitSetPred(struct SimRun *v);
   int update_stanzas(struct SimRun *v, int yr, int mon);
   int deriv_master(struct SimRun *v, int y, int m, int d);

// Fitting routines
   void load_instep(struct SimRun *v);
   void outstep(int st, float *p, struct SimRun *v);
      
   double SSQ(struct SimRun *v);
   void SSQ_scatter(struct SimRun *v, double x[]);
   int solve_dfp(struct SimRun *v, double *p, int its);
   
   void lnsrch(struct SimRun *v, int n, double xold[], double fold, double g[], double p[], double x[],
	             double *f, double stpmax, int *check, double (*func)(struct SimRun*, double []));
   void dfpmin(struct SimRun *v, double p[], int n, double gtol, int *iter, double *fret,
	             double(*func)(struct SimRun*, double []), void (*dfunc)(struct SimRun*, double [], double []));
   void frprmn(struct SimRun *v, double p[], int n, double ftol, int *iter, double *fret,
	             double (*func)(struct SimRun*, double []), void (*dfunc)(struct SimRun*, double [], double []));

// Frequency routines
   int Noisy_Run();

// -----------------------------------------------------------------------------

struct SimRun* new_SimRun(){
   struct SimRun *v;
   int sp,t, i;

   v = (struct SimRun *)calloc(1,sizeof(struct SimRun));

    v->out_BB      = dmatrix(0,NUM_GROUPS,0,fit_Years);
    v->monthly_BB  = dmatrix(0,NUM_GROUPS,0,fit_Years*STEPS_PER_YEAR);
    v->monthly_PP  = dmatrix(0,NUM_GROUPS,0,fit_Years*STEPS_PER_YEAR);

    v->out_CC      = matrix(0,NUM_GROUPS,0,fit_Years);
    v->out_RR      = matrix(0,juv_N,0,fit_Years);
    v->out_EE      = matrix(0,juv_N,0,fit_Years);
    v->out_SSB     = matrix(0,juv_N,0,fit_Years);
 
    v->out_MM      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
    v->out_M2      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
    v->out_PB      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
    v->out_QB      = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
    v->out_LinksM  = matrix(0,NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
    v->out_LinksF  = matrix(0,NUM_CATCHLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
    v->out_TotF    = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
    v->out_Effort  = matrix(0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
      
    v->fit_vector  = dvector(0,NDIM);

    v->fit_EST     = dmatrix(0,fitN,0,fit_Years);
    v->fit_diet_EST= dmatrix(0,NUM_LIVING+NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
    v->fit_diet_SSQ= dvector(0,NUM_LIVING+NUM_PREDPREYLINKS);
    v->fit_SSQ     = dmatrix(0,4,0,NUM_GROUPS);
    v->fit_bioanom = dmatrix(0,NUM_GROUPS,0,fit_Years);

    v->fit_vec_control = dvector(0,NDIM);
    v->fit_ssq_control = dmatrix(0,4,0,NUM_GROUPS);

    v->NageS       = dmatrix(0,juv_N,0,MAX_MONTHS_STANZA);
    v->WageS       = dmatrix(0,juv_N,0,MAX_MONTHS_STANZA);
    v->WWa         = dmatrix(0,juv_N,0,MAX_MONTHS_STANZA);
    v->SplitAlpha  = dmatrix(0,juv_N,0,MAX_MONTHS_STANZA);
    
    v->force_bymort = dmatrix(0,NUM_GROUPS,0,MAX_YEARS*STEPS_PER_YEAR);

     for (sp=0; sp<=NUM_GROUPS; sp++){
		     for (t=0; t<=MAX_YEARS*STEPS_PER_YEAR; t++){
				  //force_byprey[sp][t]=1.0;
				  v->force_bymort[sp][t]=force_bymort[sp][t];
				  //force_byrecs[sp][t]=1.0;
				 }
		 }

     for (i=0; i<=NDIM; i++){v->fit_vec_control[i]=1.0;}
     for (sp=0; sp<=NUM_GROUPS; sp++){
          for (i=0; i<=4; i++){
              v->fit_ssq_control[i][sp] = 1.0;
          }
     }

   return v;
    
}

void free_SimRun(struct SimRun *v){

     free_dmatrix(v->out_BB,0,NUM_GROUPS,0,fit_Years);
     free_dmatrix(v->monthly_BB,0,NUM_GROUPS,0,fit_Years*STEPS_PER_YEAR);
     free_dmatrix(v->monthly_PP,0,NUM_GROUPS,0,fit_Years*STEPS_PER_YEAR);

     free_matrix(v->out_CC,0,NUM_GROUPS,0,fit_Years);
     free_matrix(v->out_RR,0,juv_N,0,fit_Years);
     free_matrix(v->out_EE,0,juv_N,0,fit_Years);    
     free_matrix(v->out_SSB,0,juv_N,0,fit_Years); 
 
     free_matrix(v->out_MM,0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
     free_matrix(v->out_M2,0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
     free_matrix(v->out_PB,0,NUM_GROUPS,0,MAX_YEARS);  //out_PB[NUM_GROUPS+1];
     free_matrix(v->out_QB,0,NUM_GROUPS,0,MAX_YEARS);  //out_QB[NUM_GROUPS+1];
     free_matrix(v->out_LinksM,0,NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
     free_matrix(v->out_LinksF,0,NUM_CATCHLINKS,0,MAX_YEARS); 
     free_matrix(v->out_TotF,0,NUM_GROUPS,0,MAX_YEARS);   //out_TotF[NUM_GROUPS+1];
     free_matrix(v->out_Effort,0,NUM_GROUPS,0,MAX_YEARS); //[NUM_GROUPS+1][MAX_YEARS+1];
     
     free_dvector(v->fit_vector,0,NDIM);

     free_dmatrix(v->fit_EST,0,fitN,0,fit_Years);
     free_dmatrix(v->fit_diet_EST,0,NUM_LIVING+NUM_PREDPREYLINKS,0,MAX_YEARS); //out_LinksM[NUM_PREDPREYLINKS+1];
     free_dmatrix(v->fit_SSQ,0,4,0,NUM_GROUPS);
     free_dvector(v->fit_diet_SSQ,0,NUM_LIVING+NUM_PREDPREYLINKS);
     free_dmatrix(v->fit_bioanom,0,NUM_GROUPS,0,fit_Years);
 
     free_dvector(v->fit_vec_control,0,NDIM);
     free_dmatrix(v->fit_ssq_control,0,4,0,NUM_GROUPS);
        
     free_dmatrix(v->NageS,0,juv_N,0,MAX_MONTHS_STANZA);
     free_dmatrix(v->WageS,0,juv_N,0,MAX_MONTHS_STANZA);
     free_dmatrix(v->WWa,0,juv_N,0,MAX_MONTHS_STANZA);
     free_dmatrix(v->SplitAlpha,0,juv_N,0,MAX_MONTHS_STANZA);
     
     free_dmatrix(v->force_bymort,0,NUM_GROUPS,0,MAX_YEARS*STEPS_PER_YEAR);     

   free(v);
}
