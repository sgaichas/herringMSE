
//#include <fftw3.h>
int fftw_test()
     {
         double *in;
//         fftw_complex *out;
//        fftw_plan p;
         int i, N;
         
         N=1000;
         //in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//         in  = (double*) fftw_malloc(sizeof(double) * N);
//         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
         //p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//         p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_PRESERVE_INPUT);

         for (i=0; i<N; i++){in[i]=sin(double(i*2*3.14159265358979)/10.0)+sin(double(i*2*3.14159265358979)/75.0);}

//         fftw_execute(p); /* repeat as needed */
         
//         for (i=0; i<N; i++){
//             cout << i << "," << in[i] << "," << out[i][0] << "," << out[i][1] 
//             <<"," << out[i][0]*out[i][0] + out[i][1] * out[i][1] << endl; 
//         }
//         fftw_destroy_plan(p);
//         fftw_free(in); fftw_free(out);
    }
