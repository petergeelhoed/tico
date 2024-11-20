#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>



void transpone(double *arr, unsigned int N, unsigned int M)
{
    double *tmp = calloc(N*M, sizeof(double));
    unsigned int k = 0;
    for (unsigned int i=0 ; i < M ; ++i )
    {
        for (unsigned int j=0 ; j < N ; ++j )
        {
            tmp[k++] = arr[i+j*M];
        }
    }
    memcpy(arr, tmp, sizeof(double)*M*N);
    free(tmp);
}
