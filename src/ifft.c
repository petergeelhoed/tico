#include <ctype.h>
#include <errno.h>
#include <fftw3.h>
#include <limits.h>
#include <locale.h> // Optional: setlocale for decimal point
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
#include "parseargs.h"

int main(void)
{
    if (setlocale(LC_NUMERIC, "C") == NULL)
    {
        (void)fprintf(stderr, "Cannot set LC_NUMERIC to C\n");
    };

    unsigned int bufferLength = INIT_N;
    unsigned int index = 0;

    double* tmpy = calloc(bufferLength, sizeof(double));
    double* tmpx = calloc(bufferLength, sizeof(double));
    if (!tmpx || !tmpy)
    {
        (void)fprintf(stderr, "Memory allocation failed (initial buffers)\n");
        free(tmpx);
        free(tmpy);
        return -2;
    }

    double triplet[3];
    for (;;)
    {
        int nrDoubles = getDoublesFromStdin(3, triplet);
        if (nrDoubles < 0)
        {
            break;
        }
        if (nrDoubles < 3)
        {
            continue;
        }

        tmpx[index] = triplet[1];
        tmpy[index] = triplet[2];
        index++;

        if (index == bufferLength)
        {
            unsigned int newN = (bufferLength * 3) / 2;

            double* newTmpy = realloc(tmpy, newN * sizeof(double));
            if (!newTmpy)
            {
                (void)fprintf(stderr, "Memory allocation failed (grow tmpy)\n");
                free(tmpx);
                free(tmpy);
                return -2;
            }
            tmpy = newTmpy;

            double* newTmpx = realloc(tmpx, newN * sizeof(double));
            if (!newTmpx)
            {
                (void)fprintf(stderr, "Memory allocation failed (grow tmpx)\n");
                free(tmpy);
                free(tmpx);
                return -2;
            }
            tmpx = newTmpx;

            bufferLength = newN;
        }
    }

    if (index == 0)
    {
        (void)fprintf(stderr, "No valid input parsed from stdin\n");
        free(tmpy);
        free(tmpx);
        return -1;
    }

    unsigned int arrayLength = index;

    fftw_complex* input = fftw_alloc_complex(arrayLength);
    if (!input)
    {
        (void)fprintf(stderr, "fftw_alloc_complex(input) failed\n");
        free(tmpy);
        free(tmpx);
        return -2;
    }
    for (index = 0; index < arrayLength; index++)
    {
        input[index][0] = tmpx[index];
        input[index][1] = tmpy[index];
    }
    free(tmpy);
    free(tmpx);

    fftw_complex* out = fftw_alloc_complex(arrayLength);
    if (!out)
    {
        (void)fprintf(stderr, "fftw_alloc_complex(out) failed\n");
        fftw_free(input);
        return -2;
    }

    fftw_plan plan = fftw_plan_dft_1d(
        (int)arrayLength, input, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (!plan)
    {
        (void)fprintf(stderr, "fftw_plan_dft_1d failed\n");
        fftw_free(input);
        fftw_free(out);
        return -2;
    }

    fftw_execute(plan);

    for (index = 0; index < arrayLength; index++)
    {
        printf("%u %g %g\n", index, out[index][0], out[index][1]);
    }

    fftw_destroy_plan(plan);
    fftw_free(input);
    fftw_free(out);
    fftw_cleanup();
    return 0;
}
