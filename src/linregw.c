
#include "mymath.h"
#include "parseargs.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const size_t initialSize = 1000000;

    // Prefer heap to avoid large stack frames
    double* xarr = (double*)malloc(initialSize * sizeof(double));
    double* yvec = (double*)malloc(initialSize * sizeof(double));
    double* wvec = (double*)malloc(initialSize * sizeof(double));
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
        int nrDoubles = getDoublesFromStdin(3, arr);
        if (nrDoubles < 0)
        {
            // EOF or read error -> stop
            break;
        }
        if (nrDoubles == 0)
        {
            // Blank line or no numbers -> skip to next line
            continue;
        }
        if (nrDoubles < 3)
        {
            // Line did not have 3 doubles -> depending on your policy:
            // either skip, or fill defaults, or break. Here we skip.
            (void)fprintf(stderr,
                          "Warning: line had only %d double(s); skipping\n",
                          nrDoubles);
            continue;
        }

        if (arrayLength >= initialSize)
        {
            (void)fprintf(
                stderr,
                "Capacity reached (%zu). Remaining input will be ignored.\n",
                initialSize);
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
