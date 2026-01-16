#include "crosscorint.h"

#include "mydefs.h"
#include "myfft.h"

#include <math.h>
#include <stdlib.h>

// Perform cross-correlation using FFT
void crosscorint(unsigned int ArrayLength,
                 const int* array,
                 const int* ref,
                 int* cross)
{
    // Allocate memory for FFTW complex arrays
    fftw_complex* tmparr = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmpref = fftw_alloc_complex(ArrayLength);
    if (!tmparr || !tmpref)
    {
        (void)fprintf(stderr, "Memory allocation failed in crosscorint\n");
        if (tmparr)
        {
            fftw_free(tmparr);
        }
        if (tmpref)
        {
            fftw_free(tmpref);
        }
        fftw_cleanup();
        exit(ERROR_ALLOCATE_MEM);
    }

    // Initialize FFTW complex arrays with input data
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        tmparr[j][0] = array[j];
        tmparr[j][1] = 0.0;
        tmpref[j][0] = ref[j];
        tmpref[j][1] = 0.0;
    }

    // Perform cross-correlation
    fftw_complex* coor = crosscor(ArrayLength, tmparr, tmpref);
    if (!coor)
    {
        (void)fprintf(stderr, "Cross-correlation failed in crosscorint\n");
        fftw_free(tmpref);
        fftw_free(tmparr);
        fftw_cleanup();
        exit(ERROR_ALLOCATE_MEM);
    }

    // Store results in the output array
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        cross[j] = (int)round(coor[j][0] * ArrayLength);
    }

    // Free allocated memory
    fftw_free(tmpref);
    fftw_free(tmparr);
    fftw_free(coor);
    fftw_cleanup();
}
