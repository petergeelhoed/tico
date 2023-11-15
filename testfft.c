#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

#include "mylib.h"

int main (int argc, char *argv[])
{
    int NN = 20;
    int evalue = 1;
    int ipeak[NN];
    
    fftw_complex *peak = fftw_alloc_complex(NN);
    fftw_complex *peak2 = fftw_alloc_complex(NN);
    fftw_complex *corr = fftw_alloc_complex(NN);
    fftw_complex *tmp = fftw_alloc_complex(NN);
    fftw_complex *filter = fftw_alloc_complex(NN);
    fftw_complex *filterFFT = makeFilter(evalue, NN);

    fftw_plan reversefilter = fftw_plan_dft_1d(NN, filterFFT, filter, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan forwardpeak = fftw_plan_dft_1d(NN, peak, tmp, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reversecorr = fftw_plan_dft_1d(NN, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fprintf(stderr,"======filter FFT=====\n");
    for (int j = 0; j < NN; j++) {fprintf(stderr,"%6d %f %f\n",j,filterFFT[j][0],filterFFT[j][0]);}
    fftw_execute(reversefilter);

    for (int j = 0; j < NN; j++)
    { 
        peak2[(j+1)%NN][0] = ((j==4)+(j==7));
        peak2[j][1] = 0;
        peak[j][0] = ((j==4)+(j==7));
        peak[j][1] = 0;
        ipeak[j] = peak[j][0];
    }

    fftw_complex *cor = crosscor(NN,peak,peak2);
    fprintf(stderr,"========croscor ===\n");
    for (int j = 0; j < NN; j++) fprintf(stderr,"%6d %f\n",j,cor[j][0]);
exit(0);
    fprintf(stderr,"========filter ===\n");
    for (int j = 0; j < NN; j++) fprintf(stderr,"%6d %f\n",j,filter[j][0]);
    fprintf(stderr,"=======filter norm====\n");
//    normalise(NN,filter);
    for (int j = 0; j < NN; j++) {fprintf(stderr,"%6d %f\n",j,filter[j][0]); }
    fprintf(stderr,"======peak FFT=====\n");
    fftw_execute(forwardpeak);
    for (int j = 0; j < NN; j++) {fprintf(stderr,"%6d %f %f\n",j,peak[j][0],peak[j][0]);}
    fprintf(stderr,"====corr == peak =======\n");
    fftw_execute(reversecorr);
    for (int j = 0; j < NN; j++) {corr[j][0] /= NN;}
    for (int j = 0; j < NN; j++) {fprintf(stderr,"%6d %f == %f \n",j,corr[j][0],peak[j][0]);}


    fftw_complex *conv = convolute(NN,ipeak,filterFFT);
    fprintf(stderr,"====convolution== filter == ipeak =====\n");
    for (int j = 0; j < NN; j++) {fprintf(stderr,"%6d %f == %f == %d \n",j,conv[j][0],filterFFT[j][0], ipeak[j]);}

exit (0);
}
