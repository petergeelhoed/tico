#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main()
{
    unsigned int M = 1;
    unsigned int N = 4;
    double xarr[4] = {1.0, 2.0, 3.0, 4.0};


    unsigned int T = 1;
    unsigned int S = 4;
    double yvec[4] = {3.0,5.0,8.0,9.0};
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

        for (unsigned int i = 0; i < T+1; i++)
        {
            printf("%8.3f ", tmp[i]);
        }
        printf("\n");
    exit(0);
}
