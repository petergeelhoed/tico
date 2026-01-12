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

int main(void)
{
    if (setlocale(LC_NUMERIC, "C") == NULL)
    {
        (void)fprintf(stderr, "Cannot set LC_NUMERIC to C\n");
    };

    unsigned int iN = INIT_N;
    unsigned int i = 0;

    double* tmpy = calloc(iN, sizeof(double));
    double* tmpx = calloc(iN, sizeof(double));
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
        int k = getDoublesFromStdin(3, triplet);
        if (k < 0)
        {
            break;
        }
        if (k < 3)
        {
            continue;
        }

        tmpx[i] = triplet[1];
        tmpy[i] = triplet[2];
        i++;

        if (i == iN)
        {
            unsigned int newN = (iN * 3) / 2;

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

            iN = newN;
        }
    }

    if (i == 0)
    {
        (void)fprintf(stderr, "No valid input parsed from stdin\n");
        free(tmpy);
        free(tmpx);
        return -1;
    }

    unsigned int N = i;

    fftw_complex* in = fftw_alloc_complex(N);
    if (!in)
    {
        (void)fprintf(stderr, "fftw_alloc_complex(in) failed\n");
        free(tmpy);
        free(tmpx);
        return -2;
    }
    for (i = 0; i < N; i++)
    {
        in[i][0] = tmpx[i];
        in[i][1] = tmpy[i];
    }
    free(tmpy);
    free(tmpx);

    fftw_complex* out = fftw_alloc_complex(N);
    if (!out)
    {
        (void)fprintf(stderr, "fftw_alloc_complex(out) failed\n");
        fftw_free(in);
        return -2;
    }

    fftw_plan p =
        fftw_plan_dft_1d((int)N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    if (!p)
    {
        (void)fprintf(stderr, "fftw_plan_dft_1d failed\n");
        fftw_free(in);
        fftw_free(out);
        return -2;
    }

    fftw_execute(p);

    for (i = 0; i < N; i++)
    {
        printf("%u %g %g\n", i, out[i][0], out[i][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    return 0;
}
