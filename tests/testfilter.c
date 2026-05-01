#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

// Magic number constants
#define FILTER_ARRAY_LENGTH 20
#define FILTER_EVALUE 2
#define FILTER_BLAH_INDEX 10
#define FILTER_BLAH_VALUE 10000

int main(void)
{
    unsigned int evalue = FILTER_EVALUE;

    fftw_complex* filterFFT = makeFilter(evalue, FILTER_ARRAY_LENGTH);

    for (unsigned int j = 0; j < FILTER_ARRAY_LENGTH; j++)
    {
        //        float a = filterFFT[j][0];
        //       fprintf(stderr,"%3d %12.4f %12.4f\n",j,a,filterFFT[j][1]);
    }
    int blah[FILTER_ARRAY_LENGTH];
    int orig[FILTER_ARRAY_LENGTH];
    for (unsigned int j = 0; j < FILTER_ARRAY_LENGTH; j++)
    {
        blah[j] = (j == FILTER_BLAH_INDEX) *
                  FILTER_BLAH_VALUE; //(int)1000*sin((float)j/2.);
        orig[j] = blah[j];
    }

    struct myarr input = {blah, 0, FILTER_ARRAY_LENGTH};
    fftw_complex* out = convolute(input, filterFFT);

    for (unsigned int j = 0; j < FILTER_ARRAY_LENGTH; j++)
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
