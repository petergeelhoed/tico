#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

int main()
{
    unsigned int ArrayLength = DEFAULT_NN;
    double readValue0;
    double readValue1;
    double* peak = calloc(ArrayLength, sizeof(double));
    double* peak2 = calloc(ArrayLength, sizeof(double));

    unsigned int length = 0;
    char line[LINESIZE];
    char* endptr;
    char* next_start;

    // Read lines until EOF or error
    while (fgets(line, sizeof(line), stdin))
    {
        errno = 0;

        // 1. Convert first number
        readValue0 = strtod(line, &endptr);
        if (line == endptr)
        {
            continue; // No number found on this line
        }
        // 2. Convert second number starting where the first one ended
        next_start = endptr;
        readValue1 = strtod(next_start, &endptr);
        if (next_start == endptr)
        {
            continue; // Second number not found
        }
        peak[length] = readValue0;
        peak2[length] = readValue1;
        length++;
        if (length == ArrayLength)
        {
            ArrayLength = ArrayLength * 3 / 2;
            double* tmpArray0 = realloc(peak, ArrayLength * sizeof(double));
            double* tmpArray1 = realloc(peak2, ArrayLength * sizeof(double));

            if (tmpArray0 == NULL || tmpArray1 == NULL)
            {
                printf("out of memory");
            }
            else
            {
                peak = tmpArray0;
                peak2 = tmpArray1;
            }
        }
    }
    fftw_complex* in1 = fftw_alloc_complex(length);
    fftw_complex* in2 = fftw_alloc_complex(length);
    for (unsigned int i = 0; i < length; ++i)
    {
        in1[i][0] = peak[i];
        in1[i][1] = 0.0;
        in2[i][0] = peak2[i];
        in2[i][1] = 0.0;
    }

    fftw_complex* cross = crosscor(length, in1, in2);
    for (unsigned int i = 0; i < length; i++)
    {
        printf("%d %lf %lf %lf\n", i, peak[i], peak2[i], cross[i][0]);
    }
    free(peak);
    free(peak2);
    fftw_free(in1);
    fftw_free(in2);
    fftw_free(cross);

    fftw_cleanup();
    exit(0);
}
