#include "mymath.h"
#include <stdio.h>
#include <stdlib.h>

// Function to perform linear regression with weighted inputs
int main()
{
    unsigned int n = 1000000;
    double* xarr = calloc(n, sizeof(double));
    double* yvec = calloc(n, sizeof(double));
    double* wvec = calloc(n, sizeof(double));

    if (!(xarr && yvec && wvec))
    {
        fprintf(stderr, "Cannot allocate memory\n");
        exit(EXIT_FAILURE);
    }

    unsigned int N = 0;
    double* x = xarr;
    double* y = yvec;
    double* w = wvec;
    int keepRunning = 1;

    // Read input data for x, y, and weights
    while (keepRunning)
    {
        keepRunning = ((scanf("%lf", x) == 1) && (scanf("%lf", y) == 1) &&
                       (scanf("%lf", w) == 1));
        x++;
        y++;
        w++;
        if (keepRunning)
            N++;
    }

    unsigned int M = 1;
    unsigned int T = M;

    // Perform linear regression
    double* coeffs = matlinreg(xarr, N, M, yvec, wvec);

    for (unsigned int i = 0; i < T + 1; i++)
    {
        printf("%10.6g ", coeffs[i]);
    }
    printf("\n");

    free(xarr);
    free(yvec);
    free(wvec);
    free(coeffs);

    return 0;
}
