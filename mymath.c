#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void transpone(double* arr, unsigned int N, unsigned int M)
{
    double* tmp = calloc(N * M, sizeof(double));
    unsigned int k = 0;
    for (unsigned int i = 0; i < M; ++i)
    {
        for (unsigned int j = 0; j < N; ++j)
        {
            tmp[k++] = arr[i + j * M];
        }
    }
    memcpy(arr, tmp, sizeof(double) * M * N);
    free(tmp);
}

void invert(double* arr, unsigned int N, unsigned int M)
{
    // make matrix twice as wide
    double* tmp = calloc(N * M * 2, sizeof(double));
    unsigned M2 = M * 2;
    for (unsigned int j = 0; j < N; ++j)
    {
        for (unsigned int i = 0; i < M; ++i)
        {
            tmp[i + j * M2] = arr[i + j * N];
        }
    }
    // make second part diagonal 1
    for (unsigned int i = 0; i < M; ++i)
    {
        tmp[i + M + i * M2] = 1.;
    }

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int k = 0; k < N; k++)
        {
            if (j != k)
            {
                double f = tmp[j + k * M2] / tmp[j + j * M2];
                for (unsigned int i = 0; i < M2; i++)
                {
                    tmp[i + k * M2] = tmp[i + k * M2] - tmp[i + j * M2] * f;
                }
            }
        }
    }
    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            arr[i + j * M] = tmp[i + M + j * M2] / tmp[j + j * M2];
        }
    }
    free(tmp);
}

double* mulmat(double* arr,
               unsigned int N,
               unsigned int M,
               double* vec,
               unsigned int S,
               unsigned int T)
{
    if (M != S)
        exit(-1);
    double* tmp = calloc(N * T, sizeof(double));

    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int l = 0; l < T; l++)
        {
            for (unsigned int i = 0; i < M; i++)
            {
                tmp[l + j * T] += arr[i + M * j] * vec[i * T + l];
            }
        }
    }
    return tmp;
}
