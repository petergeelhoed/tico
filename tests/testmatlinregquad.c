#define MATQ_M 2
#define MATQ_N 4
#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main(void)
{
    // 1+2x-0.5x^2
    double xarr[MATQ_M * MATQ_N] = // NOLINTBEGIN[readability-magic-numbers]
        {1.0,
         1.0,
         2.0,
         4.0,
         3.0,
         9.0,
         4.0,
         16.0}; // NOLINTEND[readability-magic-numbers]

    unsigned int matT = 1;
    unsigned int matS = 4;
    double yvec[4] =          // NOLINTBEGIN[readability-magic-numbers]
        {2.5, 3.0, 2.5, 1.0}; // NOLINTEND[readability-magic-numbers]
    double wvec[4] =          // NOLINTBEGIN[readability-magic-numbers]
        {1.0, 1.0, 2.0, 1.0}; // NOLINTEND[readability-magic-numbers]

    for (unsigned int j = 0; j < MATQ_N; j++)
    {
        for (unsigned int i = 0; i < MATQ_M; i++)
        {
            printf("%8.1f ", xarr[i + j * MATQ_M]);
        }
        printf("\n");
    }
    printf("\n");

    for (unsigned int j = 0; j < matS; j++)
    {
        for (unsigned int i = 0; i < matT; i++)
        {
            printf("%8.1fw%3.1f ", yvec[i + j * matT], wvec[i + j * matT]);
        }
        printf("\n");
    }

    printf("\n");
    double tmp[MATQ_M + 1];
    matlinreg(tmp, xarr, MATQ_N, MATQ_M, yvec, wvec);
    printf("\n");

    for (unsigned int i = 0; i < MATQ_M + 1; i++)
    {
        printf("%8.3f ", tmp[i]);
    }
    printf("\n");
    exit(0);
}
