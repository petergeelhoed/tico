#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"
#include "mysync.h"

#define ArrayLength 20
int main(void)
{
    unsigned int evalue = 1;
    int ipeak[ArrayLength];

    fftw_complex* peak = fftw_alloc_complex(ArrayLength);
    fftw_complex* peak2 = fftw_alloc_complex(ArrayLength);
    fftw_complex* corr = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmp = fftw_alloc_complex(ArrayLength);
    fftw_complex* filter = fftw_alloc_complex(ArrayLength);
    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);

    fftw_plan reversefilter = fftw_plan_dft_1d(
        (int)ArrayLength, filterFFT, filter, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan forwardpeak = fftw_plan_dft_1d(
        (int)ArrayLength, peak, tmp, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reversecorr = fftw_plan_dft_1d(
        (int)ArrayLength, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fprintf(stderr, "======filter FFT=====\n");
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        fprintf(stderr, "%6d %f %f\n", j, filterFFT[j][0], filterFFT[j][0]);
    }
    fftw_execute(reversefilter);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        peak2[(j + 1) % ArrayLength][0] = ((j == 4) + (j == 7));
        peak2[j][1] = 0;
        peak[j][0] = ((j == 4) + (j == 7));
        peak[j][1] = 0;
        ipeak[j] = (int)peak[j][0];
    }

    fftw_complex* cor = crosscor(ArrayLength, peak, peak2);
    fprintf(stderr, "========croscor ===\n");
    for (unsigned int j = 0; j < ArrayLength; j++)
        fprintf(stderr, "%6d %f\n", j, cor[j][0]);
    exit(0);
    fprintf(stderr, "========filter ===\n");
    for (unsigned int j = 0; j < ArrayLength; j++)
        fprintf(stderr, "%6d %f\n", j, filter[j][0]);
    fprintf(stderr, "=======filter norm====\n");
    //    normalise(ArrayLength,filter);
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        fprintf(stderr, "%6d %f\n", j, filter[j][0]);
    }
    fprintf(stderr, "======peak FFT=====\n");
    fftw_execute(forwardpeak);
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        fprintf(stderr, "%6d %f %f\n", j, peak[j][0], peak[j][0]);
    }
    fprintf(stderr, "====corr == peak =======\n");
    fftw_execute(reversecorr);
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        corr[j][0] /= ArrayLength;
    }
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        fprintf(stderr, "%6d %f == %f \n", j, corr[j][0], peak[j][0]);
    }

    struct myarr input = {ipeak, 0, ArrayLength};
    fftw_complex* conv = convolute(input, filterFFT);
    fprintf(stderr, "====convolution== filter == ipeak =====\n");
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        fprintf(stderr,
                "%6d %f == %f == %d \n",
                j,
                conv[j][0],
                filterFFT[j][0],
                ipeak[j]);
    }

    exit(0);
}
