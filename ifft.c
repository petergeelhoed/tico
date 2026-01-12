#include <errno.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"

int main()
{
    unsigned int iN = INIT_N;
    unsigned int i = 0;

    double* tmpy = calloc(iN, sizeof(double));
    double* tmpx = calloc(iN, sizeof(double));
    double real;
    double img;
    int ret = 0;
    int index;
    char line[LINESIZE];
    char* endptr;
    char* next_start;
    while (fgets(line, sizeof(line), stdin))
    {
        errno = 0;
        index = strtol(line, &endptr, DECIMAL);
        (void)index;

        if (line == endptr)
        {
            continue; // No number found on this line
        }
        next_start = endptr;
        real = strtod(next_start, &endptr);
        if (next_start == endptr)
        {
            continue; // Second number not found
        }
        next_start = endptr;
        img = strtod(next_start, &endptr);
        if (next_start == endptr)
        {
            continue; // Third number not found
        }

        tmpy[i] = img;
        tmpx[i] = real;

        i++;
        if (i == iN)
        {
            iN *= 3;
            iN /= 2;
            double* tmp2 = realloc(tmpy, iN * sizeof(double));
            if (tmp2)
            {
                tmpy = tmp2;
            }
            else
            {
                (void)fprintf(stderr, "Memory allocation failed");
                return -2;
            }
            tmp2 = realloc(tmpx, iN * sizeof(double));
            if (tmp2)
            {
                tmpx = tmp2;
            }
            else
            {
                (void)fprintf(stderr, "Memory allocation failed");
                return -2;
            }
        }
    }
    if (ret == 0)
    {
        (void)fprintf(stderr, "Failed to parse double\n");
        return -1;
    }

    unsigned int N = i;

    fftw_complex* in = fftw_alloc_complex(N);
    for (i = 0; i < N; i++)
    {
        in[i][0] = tmpx[i];
        in[i][1] = tmpy[i];
    }
    free(tmpy);
    free(tmpx);

    fftw_complex* out = fftw_alloc_complex(N);

    /* forward Fourier transform, save the result in 'out' */
    fftw_plan p =
        fftw_plan_dft_1d((int)(N), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    for (i = 0; i < N; i++)
    {
        printf("%d %g %g\n", i, out[i][0], out[i][1]);
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    return 0;
}
