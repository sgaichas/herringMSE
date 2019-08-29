
// Header file that contains Ecopath model generated from .eii file
#include "../GOM/GOM_mammaldietupdate.h"
//#include "../GOM/GOM.h"  //original GOM EMAX model
//#include "../GOAsense2012/GOA_PathData.h"
//#include "../EBSsense2012/EBS_PathData_nocrabs.h"
//#include "../AIcons/2011_AI_stanza_ped.h"

// "Year window" of Ecopath model input above, to synch with "year"
// column in input forcing/fitting variables
#define PATH_YEAR_START    2000
#define PATH_YEAR_END      2003

// Definitions for parallel processes
// Number of processes to generate (usually # of processors on computer)
#define PRUNS          8   
// Max number of ecosystems in memory at one time (may need to adjust
// due to memory constraints)
#define MAX_SYSTEMS    500	

// Some definitions which might be changed sometimes
#define NUM_VARS   8                    // Fitting parameter types
#define NUM_ROWS   NUM_LIVING           // Number of vals for each parameter type
#define NDIM       NUM_VARS*NUM_ROWS    // Number of parameters in maximization vector
#define MAX_YEARS          200       // Maximum number of years for Ecosim run
#define MAX_SPLIT          30        // Max number of juvenile/adult pairs (set to save memory)
#define MAX_COLS           255      // Max number of columns in a csv file.
#define STEPS_PER_MONTH    1        // Integration Steps per month (speed (5) vs. accuracy (10 to 30))
#define MEASURE_MONTH      6        // Month of the year for measuring fits
//#define MAX_MONTHS_STANZA  288      // Max number age pools within stanza calcs (months)
#define MAX_MONTHS_STANZA  400      // Max number age pools within stanza calcs (months)
#define LO_DISCARD         1e-3     // Discard Limit for sense
#define HI_DISCARD         1e+3     // Discard Limit for sense
