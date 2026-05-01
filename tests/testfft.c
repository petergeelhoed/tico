// Magic number constants
#define FFT_ARRAY_LENGTH 20
#define FFT_EVALUE 1
#define FFT_INDEX_1 4
#define FFT_INDEX_2 7

#include "myfft.h"
#include "mylib.h"
#include "mysync.h"
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    unsigned int evalue = FFT_EVALUE;
    int ipeak[FFT_ARRAY_LENGTH];

    fftw_complex* peak = fftw_alloc_complex(FFT_ARRAY_LENGTH);
    fftw_complex* peak2 = fftw_alloc_complex(FFT_ARRAY_LENGTH);
    fftw_complex* corr = fftw_alloc_complex(FFT_ARRAY_LENGTH);
    fftw_complex* tmp = fftw_alloc_complex(FFT_ARRAY_LENGTH);
    fftw_complex* filter = fftw_alloc_complex(FFT_ARRAY_LENGTH);
    fftw_complex* filterFFT = makeFilter(evalue, FFT_ARRAY_LENGTH);

    // Removed unused variable 'reversefilter' to silence warnings
    fftw_plan forwardpeak = fftw_plan_dft_1d(FFT_ARRAY_LENGTH,
                                             peak,
                                             tmp,
                                             FFTW_FORWARD,
                                             FFTW_ESTIMATE);
    fftw_plan reversecorr = fftw_plan_dft_1d(FFT_ARRAY_LENGTH,
                                             tmp,
                                             corr,
                                             FFTW_BACKWARD,
                                             FFTW_ESTIMATE);

    (void)fprintf(stderr, "======filter FFT=====\n");
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        peak2[(j + 1) % FFT_ARRAY_LENGTH][0] =
            ((j == FFT_INDEX_1) + (j == FFT_INDEX_2));
        peak2[j][1] = 0;
        peak[j][0] = ((j == FFT_INDEX_1) + (j == FFT_INDEX_2));
        peak[j][1] = 0;
        ipeak[j] = (int)peak[j][0];
    }

    fftw_complex* cor = crosscor(FFT_ARRAY_LENGTH, peak, peak2);
    (void)fprintf(stderr, "========croscor ===\n");
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr, "%6d %f\n", j, cor[j][0]);
    }
    (void)fprintf(stderr, "========filter ===\n");
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr, "%6d %f\n", j, filter[j][0]);
    }
    (void)fprintf(stderr, "=======filter norm====\n");
    //    normalise(FFT_ARRAY_LENGTH,filter);
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr, "%6d %f\n", j, filter[j][0]);
    }
    (void)fprintf(stderr, "======peak FFT=====\n");
    fftw_execute(forwardpeak);
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr, "%6d %f %f\n", j, peak[j][0], peak[j][0]);
    }
    (void)fprintf(stderr, "====corr == peak =======\n");
    fftw_execute(reversecorr);
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        corr[j][0] /= FFT_ARRAY_LENGTH;
    }
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr, "%6d %f == %f \n", j, corr[j][0], peak[j][0]);
    }

    struct myarr input = {ipeak, 0, FFT_ARRAY_LENGTH};
    fftw_complex* conv = convolute(input, filterFFT);
    (void)fprintf(stderr, "====convolution== filter == ipeak =====\n");
    for (unsigned int j = 0; j < FFT_ARRAY_LENGTH; j++)
    {
        (void)fprintf(stderr,
                      "%6d %f == %f == %d \n",
                      j,
                      conv[j][0],
                      filterFFT[j][0],
                      ipeak[j]);
    }

    exit(0);
}
