
#include "mylib.h"
#include "mymath.h"
#include "parseargs.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const size_t initial_size = 1000000;

    // Prefer heap to avoid large stack frames
    double* xarr = (double*)malloc(initial_size * sizeof(double));
    double* yvec = (double*)malloc(initial_size * sizeof(double));
    double* wvec = (double*)malloc(initial_size * sizeof(double));
    if (!xarr || !yvec || !wvec)
    {
        (void)fprintf(stderr, "memory allocation failed\n");
        free(xarr);
        free(yvec);
        free(wvec);
        return 1;
    }

    size_t arrayLength = 0;
    double arr[3];

    // Read one line at a time, up to 3 doubles per line
    for (;;)
    {
        int nr_doubles = getDoublesFromStdin(3, arr);
        if (nr_doubles < 0)
        {
            // EOF or read error -> stop
            break;
        }
        if (nr_doubles == 0)
        {
            // Blank line or no numbers -> skip to next line
            continue;
        }
        if (nr_doubles < 3)
        {
            // Line did not have 3 doubles -> depending on your policy:
            // either skip, or fill defaults, or break. Here we skip.
            (void)fprintf(stderr,
                          "Warning: line had only %d double(s); skipping\n",
                          nr_doubles);
            continue;
        }

        if (arrayLength >= initial_size)
        {
            (void)fprintf(
                stderr,
                "Capacity reached (%zu). Remaining input will be ignored.\n",
                initial_size);
            break;
        }

        xarr[arrayLength] = arr[0];
        yvec[arrayLength] = arr[1];
        wvec[arrayLength] = arr[2];
        arrayLength++;
    }

    // Perform linear regression
    unsigned int predictors = 1;     // number of predictors
    unsigned int terms = predictors; // printed terms count control

    double coeffs[2] = {0.0, 0.0};
    matlinreg(coeffs, xarr, (unsigned int)arrayLength, predictors, yvec, wvec);

    for (unsigned int i = 0; i < terms + 1; i++)
    {
        (void)printf("%10.6g ", coeffs[i]);
    }
    (void)printf("\n");

    free(xarr);
    free(yvec);
    free(wvec);
    return 0;
}
