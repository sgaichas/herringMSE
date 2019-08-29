
// DIFFERENT DISTRIBUTIONS ARE AT THE END (Kerim Aydin)
/* 
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed) 
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

//#include <stdio.h>
#include <math.h>
#include <string.h>  // Added for memeset function
#include "gauss_threads.h"
#include "RandomSeeds.h"

#define TWOPI 6.283185307179586476925286

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

//static unsigned long state[N]; /* the array for the state vector  */
//static int left = 1;
//static int initf = 0;
//static unsigned long *next;

/* initializes state[N] with a seed */
void init_genrand(unsigned long s, struct randvec *v)
{
    int j;
    v->state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        v->state[j] = (1812433253UL * (v->state[j-1] ^ (v->state[j-1] >> 30)) + j); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array v->state[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        v->state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    v->left = 1; v->initf = 1;
}

/* initialize by an array with init_key of length key_length*/
/* init_key is the array for initializing keys */
/* key_length is its length */
#define key_length 8
void init_by_array(unsigned long key_offset, struct randvec *v)
//unsigned long init_key[], key_length;
{
    int i, j, k;
    int ko;
    init_genrand(19650218UL, v);
    i=1; j=0;
    ko = (key_offset%512)*8;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        v->state[i] = (v->state[i] ^ ((v->state[i-1] ^ (v->state[i-1] >> 30)) * 1664525UL))
          + MainRandSeed[ko+j] + j; /* non linear */
        v->state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { v->state[0] = v->state[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        v->state[i] = (v->state[i] ^ ((v->state[i-1] ^ (v->state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        v->state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { v->state[0] = v->state[N-1]; i=1; }
    }

    v->state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    v->left = 1; v->initf = 1;
}
#undef key_length

static void next_state(struct randvec *v)
{
    unsigned long *p=v->state;
    int j;

    /* if init_genrand() has not been called, */
    /* a default initial seed is used         */
    if (v->initf==0) init_genrand(5489UL, v);

    v->left = N;
    v->next = v->state;
    
    for (j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);

    for (j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], v->state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(struct randvec *v)
{
    unsigned long y;

    if (--v->left == 0) next_state(v);
    y = *v->next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(randvec *v)
{
    unsigned long y;

    if (--v->left == 0) next_state(v);
    y = *v->next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(struct randvec *v)
{
    unsigned long y;

    if (--v->left == 0) next_state(v);
    y = *v->next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double uniform(struct randvec *v)
{
    unsigned long y;

    if (--v->left == 0) next_state(v);
    y = *v->next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(struct randvec *v)
{
    unsigned long y;

    if (--v->left == 0) next_state(v);
    y = *v->next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(struct randvec *v) 
{ 
    unsigned long a=genrand_int32(v)>>5, b=genrand_int32(v)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

// ----------------------------------------------------------------------
// Gaussian function, rad uniform should be (0..1] so 1-genrand_real2
// (Value of 0 gives error in log)

double gaussian(struct randvec *v){
       double rad,theta,val;
					rad   =  sqrt (-2.0 * log(1.0-uniform(v)));
					theta =  TWOPI * uniform(v);
					return (rad * cos (theta));       
}

// -----------------------------------------------------------------------

#define EPSILON 0.00000001   
#define eee 2.718281828

float gamma (float ia, float ib, struct randvec *v)
{

int int_ia, j;
float x,del,xfrac;
float v1,v2;

float v0, nu;
    if ((ia<=EPSILON) || (ib<=EPSILON)) { return 0.0; }
    
    int_ia=int(ia);  
		del = ia - (float)int_ia;
    x=0.0;
    // Gamma (n, 1) ~ sum(j)[-ln Uk] where Uk ~ uniform(0..1]
    for (j=1; j<=int_ia; j++){
		    x -= log(1.0-uniform(v));     		    
		}
		
    xfrac = 0.0;
    v0 = (eee/(eee + del));
    if (del >= EPSILON){
       do {
          v1 = (1.0-uniform(v));
          v2 = (1.0-uniform(v));
          if ( v1<= v0){
				     xfrac = exp( (1.0/del) * log (v1/v0));
					 	 nu    = v2 * exp( (del-1.0) * log (xfrac));				 
				  }
				  else{
				     xfrac = 1.0 - log ((v1 - v0)/(1-v0));
				     nu    = v2 * exp (-xfrac);
				  }
				 }
				 while (nu > exp(-xfrac) * exp( (del-1.0) * log (xfrac)));
/*		   do {
			    do {
				     do {
				 	      v1=(1.0-uniform()); //ran1(idum);
					      v2=2.0*(1.0-uniform())-1.0;
				        } 
						 while (v1*v1+v2*v2 > 1.0);
				     y=v2/v1;
				     am=ia-1;
				     s=sqrt(2.0*am+1.0);
				     xfrac=s*y+am;
			    } 
					while (xfrac <= 0.0);
			   e=(1.0+y*y)*exp(am*log(xfrac/am)-s*y);
       } 
			 while ((1.0-uniform()) > e);
*/  
    }
    return (1.0/ib) * (xfrac + x);

}

#undef eee
#undef EPSILON
// -------------------------------------------------
#define EPSILON 1e-16   
#define eee 2.71828182845904523536

double d_gamma (double ia, double ib, struct randvec *v)
{

int int_ia, j;
double x,del,xfrac;
double v1,v2;
double v0, nu;

    if ((ia<=EPSILON) || (ib<=EPSILON)) { return 0.0; }
    
    int_ia=int(ia);  
		del = ia - (double)int_ia;
    x=0.0;
    // Gamma (n, 1) ~ sum(j)[-ln Uk] where Uk ~ uniform(0..1]
    for (j=1; j<=int_ia; j++){
		    x -= log(1.0-genrand_res53(v));     		    
		}
		
    xfrac = 0.0;
    v0 = (eee/(eee + del));
    if (del >= EPSILON){
       do {
          v1 = (1.0-genrand_res53(v));
          v2 = (1.0-genrand_res53(v));
          if ( v1<= v0){
				     xfrac = exp( (1.0/del) * log (v1/v0));
					 	 nu    = v2 * exp( (del-1.0) * log (xfrac));				 
				  }
				  else{
				     xfrac = 1.0 - log ((v1 - v0)/(1-v0));
				     nu    = v2 * exp (-xfrac);
				  }
				 }
				 while (nu > exp(-xfrac) * exp( (del-1.0) * log (xfrac)));

    }
    return (1.0/ib) * (xfrac + x);

}

#undef eee
#undef EPSILON

// -------------------------------------------------

#define TWOPI 6.283185307179586476925286 
int NoiseSeries(int NS,
								double *amp,
								double *var,
								double *offset,
								double *final,
								double cv,
								double randscale,
								int StartTime,
								int MAXT,
								int TLEN,
                struct randvec *v)
{
	int S,t;
	double freq;
	//float *outvec;
	double X, X2, MM, SD, TL;
  struct randvec rng;   // The random number state holder for this run
	memset(final, 0, (1+TLEN)*sizeof(double));


  
	for (S=1; S<=NS; S++){                 
		freq=double(S)/(double(MAXT));   
		for (t=1; t<=TLEN; t++){
			final[t] += amp[S] * var[S] * (  
				    (1.0-randscale*genrand_res53(v)) * sin(TWOPI * ( freq * (double)(t+StartTime) + offset[S])) +
						(1.0-randscale*genrand_res53(v)) * cos(TWOPI * ( freq * (double)(t+StartTime) + offset[S]))						                  
						);
		}
	}
 
	X=0.0; X2=0.0;
	for (t=1; t<=TLEN; t++){
	     X  += final[t];
		 X2 += final[t] * final[t];
	}
	TL=(double)TLEN;
	MM   = X/TL;
	SD  = TL * X2 - X * X;
	SD  = sqrt(SD / (TL * (TL - 1.0) ) );
	if  (SD <= 0.0) {SD=1;}
	for (t=1; t<=TLEN; t++){
		 final[t] = cv * (final[t]-MM)/SD;
	}

return 0;
}
#undef TWOPI

