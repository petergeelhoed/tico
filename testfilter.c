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

    for (int j = 0; j < buffer_frames; j++) 
    {
        float a = filterFFT[j][0];
        float b = a;
       fprintf(stderr,"%3d %12.4f %12.4f\n",j,a,b);
    }
 //   filterFFT[0][0] = 0;
    exit (0);
}
