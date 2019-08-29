

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


 
