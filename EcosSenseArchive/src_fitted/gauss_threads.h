
#define N 624
struct randvec
{
   unsigned long state[N]; /* the array for the state vector  */
   int left; //= 1;
   int initf; //= 0;
   unsigned long *next;
};
#undef N

/* initializes state[N] with a seed */
   void init_genrand(unsigned long s, struct randvec *v);
/* initialize by an array with init_key of length key_length*/
   void init_by_array(unsigned long key_offset, struct randvec *v);
/* generates a random number on [0,0xffffffff]-interval */   
   unsigned long genrand_int32(struct randvec *v);
/* generates a random number on [0,0x7fffffff]-interval */   
   long genrand_int31(struct randvec *v);
/* generates a random number on [0,1]-real-interval */   
   double genrand_real1(struct randvec *v);
/* generates a random number on [0,1)-real-interval */ 
   double uniform(struct randvec *v);
/* generates a random number on (0,1)-real-interval */
   double genrand_real3(struct randvec *v); 
/* generates a random number on [0,1) with 53-bit resolution*/
   double genrand_res53(struct randvec *v); 
// Gaussian distribution
   double gaussian(struct randvec *v);
// Gamma distribution (float and double versions)
   float gamma (float ia, float ib, struct randvec *v);
   double d_gamma (double ia, double ib,struct randvec *v);
// Noise series of summed random sine and cosine waves
   int NoiseSeries(int NS, double *amp, double *var, double *offset, double *final,
	   							 double cv, double randscale, int StartTime, int MAXT, int TLEN,
                    struct randvec *v);
								
