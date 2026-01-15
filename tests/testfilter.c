#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

int main()
{
    unsigned int ArrayLength = 20;
    unsigned int evalue = 2;

    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        //        float a = filterFFT[j][0];
        //       fprintf(stderr,"%3d %12.4f %12.4f\n",j,a,filterFFT[j][1]);
    }
    int blah[ArrayLength];
    int orig[ArrayLength];
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        blah[j] = (j == 10) * 10000; //(int)1000*sin((float)j/2.);
        orig[j] = blah[j];
    }

    struct myarr input = {blah, 0, ArrayLength};
    fftw_complex* out = convolute(input, filterFFT);

    for (unsigned int j = 0; j < ArrayLength; j++)
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
