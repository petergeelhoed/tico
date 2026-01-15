#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main(void)
{
    unsigned int M = 3;
    unsigned int N = 3;
    double xarr[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    unsigned int T = 2;
    unsigned int S = 3;
    double yvec[6] = {1.0, 2.0, 1.0, 2.0, 1.0, 2.0};

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            printf("%8.1f ", xarr[i + j * M]);
        }
        printf("\n");
    }
    printf("\n");

    for (unsigned int j = 0; j < S; j++)
    {
        for (unsigned int i = 0; i < T; i++)
        {
            printf("%8.1f ", yvec[i + j * T]);
        }
        printf("\n");
    }

    printf("\n");
    double* tmp = mulmat(xarr, N, M, yvec, S, T);
    printf("\n");

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < T; i++)
        {
            printf("%8.3f ", tmp[i + j * T]);
        }
        printf("\n");
    }
    exit(0);
}
