#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main()
{
    unsigned int M = 2;
    unsigned int N = 4;
    // 1+2x-0.5x^2
    double xarr[8] = {1.0,1.0, 2.0,4.0, 3.0,9.0, 4.0,16.0};


    unsigned int T = 1;
    unsigned int S = 4;
    double yvec[4] = {2.5,3.0,2.5,1.0};
    double wvec[4] = {1.0,1.0,2.0,1.0};

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
            printf("%8.1fw%3.1f ", yvec[i + j * T] , wvec[i + j * T]);
        }
        printf("\n");
    }

    printf("\n");
    double* tmp = matlinreg(xarr,N,M, yvec,wvec);
    printf("\n");

        for (unsigned int i = 0; i < M+1; i++)
        {
            printf("%8.3f ", tmp[i]);
        }
        printf("\n");
    exit(0);
}
