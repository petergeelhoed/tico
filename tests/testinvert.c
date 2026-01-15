#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main()
{
    unsigned int M = 3;
    unsigned int N = 3;
    double xarr[9] = {1.0, 2.0, -1.0, 2.0, 1.0, 2.0, -1.0, 2.0, 1.0};

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            printf("%8.1f ", xarr[i + j * M]);
        }
        printf("\n");
    }
    printf("\n");
    invert(xarr, N, M);
    printf("\n");

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            printf("%8.3f ", xarr[i + j * M]);
        }
        printf("\n");
    }
    exit(0);
}
