#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"
// Magic number constants
#define TRANSPONE_M 3
#define TRANSPONE_N 2
int main(void)
{
    double xarr[TRANSPONE_N * TRANSPONE_M];

    for (unsigned int j = 0; j < TRANSPONE_N; j++)
    {
        for (unsigned int i = 0; i < TRANSPONE_M; i++)
        {
            xarr[i + j * TRANSPONE_M] = i + j * TRANSPONE_M;
            printf("%12.4f ", xarr[i + j * TRANSPONE_M]);
        }
        printf("\n");
    }
    printf("\n");
    transpone(xarr, TRANSPONE_N, TRANSPONE_M);
    printf("\n");

    for (unsigned int j = 0; j < TRANSPONE_M; j++)
    {
        for (unsigned int i = 0; i < TRANSPONE_N; i++)
        {
            printf("%12.4f ", xarr[i + j * TRANSPONE_N]);
        }
        printf("\n");
    }

    exit(0);
}
