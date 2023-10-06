#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

#include "mylib.h"

int main (int argc, char *argv[])
{
    int buffer_frames = 20;
    int evalue = 4;
    
    fftw_complex *filterFFT = makeFilter(evalue, buffer_frames);

    int total[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) total[j] = 0;



    int peak[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) peak[j] = (int)(100*sin((float)(j+3)));;
    int der[buffer_frames];
    for (int j = 0; j < buffer_frames; j++) der[j] = (int)(100*sin((float)j));
 
    for (int j = 0; j < buffer_frames; j++) fprintf(stderr,"%6d        %6d %6d %6d %f\n",j,total[j],der[j],peak[j],filterFFT[j][0]);

    int val = 0;


    int maxpos = fftfit(
            der,
            total,
            peak,
            &val,
            filterFFT,
            buffer_frames);

//    total should contain der * 20 * maxcor^2
//    maxpos is peak position (so 10 for no shift)
//    val has masval *16

    for (int j = 0; j < buffer_frames; j++) fprintf(stderr,"%6d %6d %6d %6d %6d %6d\n",j,-maxpos+buffer_frames/2,total[j],der[j],peak[j],val);

printf("mkfilte\n");

    fftw_complex *in = fftw_alloc_complex(buffer_frames);
    fftw_complex *out = fftw_alloc_complex(buffer_frames);
    for (int j = 0; j < buffer_frames; j++) 
    {
        in[j][0] = (j%5==1);
        in[j][1]=0;
    }

    fftw_plan forward = fftw_plan_dft_1d(buffer_frames, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d(buffer_frames, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(forward);
    for (int j = 0; j < buffer_frames; j++) fprintf(stderr, "%3d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f \n",j,in[j][0],in[j][1],
            out[j][0],out[j][1],
            filterFFT[j][0],filterFFT[j][1]);

    for (int j=0; j < buffer_frames ; j++)
    {
        out[j][0] = (out[j][0]*filterFFT[j][0] - out[j][1]*filterFFT[j][1])/buffer_frames;
        out[j][1] = (out[j][0]*filterFFT[j][1] + out[j][1]*filterFFT[j][0])/buffer_frames;
    }
    fftw_execute(reverse);
    double sum= 0.0;
    for (int j = 0; j < buffer_frames; j++) 
    {
        fprintf(stderr,"%3d %12.4f %12.4f\n",j,in[j][0],in[j][1]);
        sum+=in[j][0];
    }
        fprintf(stderr,"SUM %12.4f\n",sum);


    exit (0);
}
