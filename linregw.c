
#include "mylib.h"
#include "mymath.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const size_t n = 1000000;

    // Prefer heap to avoid large stack frames
    double* xarr = (double*)malloc(n * sizeof(double));
    double* yvec = (double*)malloc(n * sizeof(double));
    double* wvec = (double*)malloc(n * sizeof(double));
    if (!xarr || !yvec || !wvec)
    {
        (void)fprintf(stderr, "Memory allocation failed\n");
        free(xarr);
        free(yvec);
        free(wvec);
        return 1;
    }

    size_t N = 0;
    double arr[3];

    // Read one line at a time, up to 3 doubles per line
    for (;;)
    {
        int k = getDoublesFromStdin(3, arr);
        if (k < 0)
        {
            // EOF or read error -> stop
            break;
        }
        if (k == 0)
        {
            // Blank line or no numbers -> skip to next line
            continue;
        }
        if (k < 3)
        {
            // Line did not have 3 doubles -> depending on your policy:
            // either skip, or fill defaults, or break. Here we skip.
            (void)fprintf(
                stderr, "Warning: line had only %d double(s); skipping\n", k);
            continue;
        }

        if (N >= n)
        {
            fprintf(
                stderr,
                "Capacity reached (%zu). Remaining input will be ignored.\n",
                n);
            break;
        }

        xarr[N] = arr[0];
        yvec[N] = arr[1];
        wvec[N] = arr[2];
        N++;
    }

    // Perform linear regression
    unsigned int M = 1; // number of predictors
    unsigned int T = M; // printed terms count control

    double coeffs[2] = {0.0, 0.0};
    matlinreg(coeffs, xarr, (unsigned int)N, M, yvec, wvec);

    for (unsigned int i = 0; i < T + 1; i++)
    {
        printf("%10.6g ", coeffs[i]);
    }
    printf("\n");

    free(xarr);
    free(yvec);
    free(wvec);
    return 0;
}
