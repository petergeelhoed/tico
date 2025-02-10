#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

int main()
{
    unsigned int NN = 20;
    unsigned int evalue = 2;

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    for (unsigned int j = 0; j < NN; j++)
    {
        //        float a = filterFFT[j][0];
        //       fprintf(stderr,"%3d %12.4f %12.4f\n",j,a,filterFFT[j][1]);
    }
    int blah[NN];
    int orig[NN];
    for (unsigned int j = 0; j < NN; j++)
    {
        blah[j] = (j == 10) * 10000; //(int)1000*sin((float)j/2.);
        orig[j] = blah[j];
    }

    struct myarr input = {blah, 0, NN};
    fftw_complex* out = convolute(input, filterFFT);

    for (unsigned int j = 0; j < NN; j++)
    {
        printf("%3d %12d %12d %12f %12f\n",
               j,
               blah[j],
               orig[j],
               out[j][0],
               out[j][1]);
    }

    exit(0);
}
