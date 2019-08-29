// Automatically generated parameters from Viz to C for ecosystem

// Numbers of each type of group
#define NUM_GROUPS 32
#define NUM_LIVING 29
#define NUM_DEAD 2
#define NUM_GEARS 1
//Numbers of Links
#define NUM_PREDPREYLINKS 350
#define NUM_CATCHLINKS 123

// Variable names path_<variable> are from Ecopath balance
// Group 0 is for outside the ecosystem (input or output)
// NUMGROUPS+1 is required to range array from 0 to NUM_GROUPS
// Note that P of detritus is amount entering detritus pool
// and if no B is set in 'path for detritus, P/B is set to 0.5 so B=2*inflow
char *path_species[NUM_GROUPS+1] = { "Outside", "Phytoplankton- Primary Producers", "Bacteria", "Microzooplankton", "Small copepods", "Large Copepods", "Gelatinous Zooplankton", "Micronekton", "Macrobenthos- polychaetes", "Macrobenthos- crustaceans", "Macrobenthos- molluscs", "Macrobenthos- other", "Megabenthos- filterers", "Megabenthos- other", "Shrimp et al.", "Larval-juv fish- all", "Small Pelagics- commercial", "Small Pelagics- other", "Small Pelagics- squid", "Small Pelagics- anadromous", "Medium Pelagics- (piscivores & other)", "Demersals- benthivores", "Demersals- omnivores", "Demersals- piscivores", "Sharks- pelagics", "HMS", "Pinnipeds", "Baleen Whales", "Odontocetes", "Sea Birds", "Discard", "Detritus-POC", "Fishery"};
float path_BB[NUM_GROUPS+1] = { 1, 22.126, 5.484, 4.885, 10.403, 11.955, 1.283, 4.874, 18.942, 4.04, 9.866, 24.936, 2.879, 3.505, 0.396, 0.207, 5.714, 1.24, 0.275, 0.153, 0.0229, 2.981, 0.4, 4.006, 0.00296, 0.00587, 0.063, 0.602, 0.0336, 0.0035, 0.442, 81.333, 0};
float path_PB[NUM_GROUPS+1] = { 1, 163.143, 91.25, 72, 30.918, 35, 35, 14.25, 2.55, 3.3, 2.24, 2.04, 0.864, 1.68, 2, 15, 0.52, 0.42, 1.4, 0.437, 0.649, 0.459, 0.54, 0.55, 0.15, 0.5, 0.067, 0.042, 0.04, 0.275, 55.84456, 29.90675, 0};
float path_QB[NUM_GROUPS+1] = { 1, 0, 380.208, 242.424, 127.75, 109.5, 146, 36.5, 17.5, 21, 13.72, 11.777, 10, 11.03, 5, 45, 1.882, 2, 2, 2, 1.428, 0.9, 0.9, 1.014, 0.623, 2.362, 4.85, 2.3, 8.5, 5.362, 0, 0, 0};
float path_EE[NUM_GROUPS+1] = { 1, 0.8799849, 0.8828096, 0.9400325, 0.9124651, 0.7455723, 0.9142609, 0.8789812, 0.9042419, 0.8891937, 0.857506, 0.8912097, 0.8695364, 0.8498231, 0.8038035, 0.8922101, 0.8667946, 0.8881668, 0.8681036, 0.9139099, 0.9032695, 0.918685, 0.8572816, 0.7878628, 0.8402626, 0.8127174, 0.2868933, 1.682192E-02, 0.5582092, 0.1244432, 0.3113834, 0.9953184, 0};
float path_BA[NUM_GROUPS+1] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float path_GS[NUM_GROUPS+1] = { 0, 0, 0.2, 0.1, 0.25, 0.25, 0.35, 0.25, 0.5, 0.5, 0.6, 0.5, 0.7, 0.3, 0.3, 0.15, 0.15, 0.35, 0.15, 0.15, 0.15, 0.3, 0.35, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.15, 0, 0, 0};
float path_DtImp[NUM_GROUPS+1] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Pedigree row 0 is placeholder
// Pedigree row 1 is biomass, row 2 PB, row 3 QB, row 4 diet, all others gears
float pedigree[4+NUM_GEARS+1][NUM_LIVING+NUM_DEAD+1] = { 
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0.1,0.8,0.8,0.3,0.3,0.5,0.5,0.3,0.3,0.3,0.3,0.3,0.5,0.1,0.5,0.1,0.3,0.3,0.3,0.3,0.1,0.1,0.1,0.5,0.5,0.3,0.3,0.3,0.3,0,0,
         0,0,0.6,0.6,0.5,0.5,0.6,0.6,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.1,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.4,0.2,0.5,0.5,0.5,0.5,0,0,
         0,0,0.6,0.6,0.5,0.5,0.6,0.6,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.1,0.3,0.3,0.3,0.3,0.1,0.1,0.1,0.4,0.4,0.5,0.5,0.5,0.5,0,0,
         0,0,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.5,0.1,0.3,0.3,0.3,0.1,0.1,0.1,0.1,0.3,0.3,0.1,0.5,0.5,0.3,0,0,
         0,0,0,0,0,0,0.7,0,0.7,0.7,0.7,0.7,0.3,0.3,0.1,0.7,0.1,0.5,0.1,0.3,0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.5,0.5,0
        };
// Order for DC is path_DC[PREY][PRED]
// Order for DetFate is path_DET[Detritus][NonDetritus]
// Order for Catch is path_Catch[Fishery][Target]
// Order for Discards is path_Discards[Fishery][Target]
float path_DC[NUM_GROUPS+1][NUM_GROUPS+1] = { 
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0.431,0.212,0.724,0.6311606,8.600464E-02,0.1821822,0.1240447,0.1350932,0.400104,0.201177,0.791,0,0.0751503,5.894106E-02,1.101101E-02,0.166093,0,0.012,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0.197,0,0,2.000108E-02,0,0.2951063,0.1370946,0.1950507,0.201177,0.119,0.1417364,0.4478958,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0.097,0.113,4.294289E-02,0.0500027,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,3.200001E-02,0.1098539,0.3300178,0.2212212,0,5.003452E-03,0,0,0,0,0,0.4305694,0.1191191,8.104539E-02,0,5.600001E-02,0,0,0,0,0,0,0,0.0310062,0,0,0,0,0,
         0,0,0,0,0,0.12783,0.3600194,0.3663664,0,2.701864E-02,0,0,0,0,0,0.2487512,0.4344345,0.6763788,0.1007984,0.9020001,0,0,0,0,3.306613E-02,0,0,0.5211042,0,3.398811E-02,0,0,0,
         0,0,0,0,0,2.396812E-02,4.900265E-02,0,0,0,0,0,0,0,0,0,1.901902E-02,2.201233E-02,0,0,0.0200056,9.018036E-03,8.308309E-02,4.304305E-02,0,0.094,0,0.005001,9.989182E-03,0,0,0,0,
         0,0,0,0,0,0,0,4.604604E-02,5.902125E-02,5.403729E-02,3.400885E-02,6.705901E-02,0,0,0.1513026,6.793207E-02,0.2172172,4.402466E-02,0.501996,0.02,0,6.112225E-02,3.303304E-02,2.202202E-02,0,0,8.003202E-02,0.325065,0.0329643,0.1409507,0,0,0,
         0,0,0,0,0,0,0,0,2.700973E-02,0.1350932,0,6.605813E-02,0,7.885332E-02,0,2.097902E-02,1.101101E-02,0,0,0.002,0,0.1472946,0.1661662,1.101101E-02,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,3.295617E-04,0,0,5.001801E-03,3.302278E-02,0,1.000881E-02,0,6.587746E-02,8.016032E-03,1.298701E-02,5.505506E-02,7.204035E-04,0.1506986,0.002,1.200336E-02,0.1843687,0.1661662,1.101101E-02,1.102204E-02,0,0,0.0460092,0,0,0,0,0,
         0,0,0,0,0,0,0,0,6.402305E-04,2.201519E-02,4.00104E-03,2.201938E-02,0,0.2375581,0,4.995005E-03,8.008009E-03,7.204035E-04,0,0,0,0.1232465,0.1661662,3.303304E-02,0,0,0,0.0140028,0,0,0,0,0,
         0,0,0,0,0,1.997343E-03,0,0,2.100756E-02,0.1350932,2.200572E-02,3.603171E-02,0,0.2375581,0.1052104,2.097902E-02,1.101101E-02,2.001121E-03,1.696607E-02,0,1.200336E-02,0.1843687,0.1111111,0.1091091,1.102204E-02,0,0,0.0260052,0,0,0,0,0,
         0,0,0,0,0,0,0,0,1.00036E-03,5.003452E-03,3.00078E-03,1.000881E-03,0,9.981433E-03,0,0,0,0,0,0,0,5.210421E-02,4.004004E-02,1.101101E-02,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,1.00036E-03,4.002762E-03,7.401924E-04,6.005285E-03,0,4.192202E-02,0,0,0,0,0,0,0.1940544,9.619239E-02,0.1111111,2.202202E-02,0,0,0,0.0040008,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.002004E-03,0,0,0,2.395209E-02,0,1.700476E-02,1.703407E-02,3.003003E-03,0.1211211,0,0,0,0,0,0.0219923,0,0,0,
         0,0,0,0,0,0,4.000216E-03,0,0,0,0,0,0,0,0,7.492507E-02,0.1081081,7.003923E-03,0.1676647,6.000001E-03,1.200336E-02,0,1.101101E-02,1.101101E-02,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,3.600194E-04,0,0,0,0,0,0,0,0,0,0,0,7.984033E-03,0,0.3641019,4.809619E-02,4.404405E-02,0.2702703,0.2344689,0.18,0.4831933,0.0070014,0.2727047,0.3118908,0,0,0,
         0,0,0,0,0,0,4.500243E-04,0,0,0,0,0,0,0,0,0,0,0,7.984033E-03,0,3.701036E-02,2.004008E-03,2.002002E-03,5.605606E-02,5.611223E-02,0.701,0.1610644,0.0030006,0.2467328,0.2679062,0,0,0,
         0,0,0,0,0,0,9.200497E-05,0,0,0,0,0,0,0,0,0,0,0,2.095808E-02,0,2.600728E-02,4.008017E-03,2.002002E-03,4.404405E-02,0.1783567,0.025,0,0.0030006,0.3306419,3.198881E-02,0,0,0,
         0,0,0,0,0,0,4.400238E-05,0,0,0,0,0,0,0,0,0,0,0,9.980041E-04,0,1.900532E-02,0,0,2.002002E-03,2.204409E-02,0,5.802321E-02,8.0016E-04,6.992428E-03,0.019993,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.202017E-04,0,0,0,3.106212E-02,0,6.002402E-04,0,8.291021E-05,3.498775E-04,0,0,0,
         0,0,0,0,0,0,0,0,0,2.001381E-03,0,2.001762E-03,0,2.99443E-03,0,0,2.002002E-03,0,0,0,0.1220342,3.707415E-02,1.201201E-02,1.101101E-02,5.611223E-02,0,8.003202E-02,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,3.10214E-04,0,1.201057E-04,0,8.584032E-04,0,0,2.002002E-03,0,0,0,4.201176E-02,6.012024E-03,9.009009E-03,4.004004E-03,8.917836E-02,0,5.702281E-02,0,1.598269E-02,1.299545E-02,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,1.000881E-03,0,1.996287E-03,0,0,2.002002E-03,0,0,0,0.1220342,4.008017E-03,6.006006E-03,0.2182182,7.815631E-02,0,8.003202E-02,0,0.0819113,2.898985E-02,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.306613E-02,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.102204E-02,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.306613E-02,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.102204E-02,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.204409E-02,0,0,0,1.997837E-03,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.306613E-02,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,6.002161E-03,7.004833E-03,6.00156E-03,9.007927E-03,0,3.892759E-02,0.0751503,0,0,0,0,0,0,1.202405E-02,2.302302E-02,0,0,0,0,0,0,0.1289549,0,0,0,
         0,0,0.569,0.494,0.131,6.191765E-02,0.1000054,0.1841842,0.4601657,0.2982058,0.3350871,0.3773321,0.09,0.1417364,0.1362725,5.894106E-02,0,0,0,0,0,1.202405E-02,1.101101E-02,0,5.611223E-02,0,0,0.0140028,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        };
float path_DetFate[NUM_DEAD+1][NUM_GROUPS+1] = { 
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0,
         0,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.001,0,0
        };
float path_Catch[NUM_GEARS+1][NUM_GROUPS+1] = { 
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0.0917,0.345,0.0622,0,0.874,0.0000151,0.0138,0.00561,0.0101,0.142,0.00718,0.301,0.00024,0.00182,0,0,0,0,0,0,0
        };
float path_Discards[NUM_GEARS+1][NUM_GROUPS+1] = { 
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,0,0,0,0,6.36E-07,0,0.000135,0.0000183,0.000189,0.000724,0.0302,0.104,0.0195,2.65E-09,0.131,0.000156,0.000562,0.0167,0.00303,0.043,0.00228,0.0904,0.0000721,0.000545,0.00115,0.000405,0.000139,0.0000588,0,0,0
        };
