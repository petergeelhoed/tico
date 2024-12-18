#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main()
{
    unsigned int n = 1000000;
    double* xarr = calloc(n, sizeof(double));
    double* yvec = calloc(n, sizeof(double));
    double* wvec = calloc(n, sizeof(double));

    if (!(xarr && yvec && wvec))
    {

        printf("Cannot allocate memory\n");
        exit(-1);
    }
    unsigned int N = 0;
    double* x = xarr;
    double* y = yvec;
    double* w = wvec;
    int keepRunning = 1;
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

    double* tmp = matlinreg(xarr, N, M, yvec, wvec);
    for (unsigned int i = 0; i < T + 1; i++)
    {
        printf("%10.6g ", tmp[i]);
    }
    printf("\n");
    exit(0);
}
