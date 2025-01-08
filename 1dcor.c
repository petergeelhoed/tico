#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

int main()
{

    unsigned int NN = 1000;
    double d;
    double d1;
    double* peak = calloc(NN, sizeof(double));
    double* peak2 = calloc(NN, sizeof(double));

    unsigned int n = 0;

    while (scanf("%lf", &d) != EOF && scanf("%lf", &d1) != EOF)
    {
        peak[n] = d;
        peak2[n] = d1;
        n++;
        if (n == NN)
        {
            NN = NN * 3 / 2;
            double* p = realloc(peak, NN * sizeof(double));
            double* q = realloc(peak2, NN * sizeof(double));

            if (p == NULL || q == NULL)
            {
                printf("out of memory");
            }
            else
            {
                peak = p;
                peak2 = q;
            }
        }
    }
    fftw_complex* in = fftw_alloc_complex(n);
    fftw_complex* in2 = fftw_alloc_complex(n);
    for (unsigned int i = 0; i < n; ++i)
    {
        in[i][0] = peak[i];
        in[i][1] = 0.0;
        in2[i][0] = peak2[i];
        in2[i][1] = 0.0;
    }

    fftw_complex* cross = crosscor(n, in, in2);
    for (unsigned int i = 0; i < n; i++)
    {
        printf("%d %lf %lf %lf\n", i, peak[i], peak2[i], cross[i][0]);
    }
    free(peak);
    free(peak2);
    fftw_free(in);
    fftw_free(in2);
    fftw_free(cross);

    fftw_cleanup();
    exit(0);
}
