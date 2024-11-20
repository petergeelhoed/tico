#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main()
{
    unsigned int M = 3;
    unsigned int N = 2;
    double xarr[N*M];

    for (unsigned int j = 0 ; j < N ; j++)
    {
        for (unsigned int i = 0 ; i < M ; i++)
        {
            xarr[i+j*M] = i+j*M;
            printf("%12.4f ", xarr[i+j*M]);
        }
        printf("\n");
    }
    printf("\n");
    transpone(xarr,N,M);
    printf("\n");

    for (unsigned int j = 0 ; j < M ; j++)
    {
        for (unsigned int i = 0 ; i < N ; i++)
        {
            printf("%12.4f ", xarr[i+j*N]);
        }
        printf("\n");
    }

    exit(0);
}
