#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mymath.h"

struct mymat
{
    double* mat;
    unsigned int M;
    unsigned int N;
};

void printmat(double* arr, unsigned int N, unsigned int M)
{
    printf("\n");
    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            printf("%8.2g", arr[i + M * j]);
        }
        printf("\n");
    }
    printf("\n");
}

void transpone(double* arr, unsigned int N, unsigned int M)
{
    double* tmp = (double*)calloc(N * M, sizeof(double));
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in transpone\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < M; ++i)
    {
        for (unsigned int j = 0; j < N; ++j)
        {
            tmp[j + i * N] = arr[i + j * M];
        }
    }
    memcpy(arr, tmp, sizeof(double) * M * N);
    free(tmp);
}

void invert(double* arr, unsigned int N, unsigned int M)
{
    // make matrix twice as wide
    double* tmp = (double*)calloc(N * M * 2, sizeof(double));
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in invert\n");
        exit(EXIT_FAILURE);
    }
    unsigned int M2 = M * 2;
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
        tmp[i + M + i * M2] = 1.0;
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
                    tmp[i + k * M2] -= tmp[i + j * M2] * f;
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
    {
        fprintf(stderr, "Matrix multiplication dimension mismatch\n");
        exit(EXIT_FAILURE);
    }
    double* tmp = (double*)calloc(N * T, sizeof(double));
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in mulmat\n");
        exit(EXIT_FAILURE);
    }
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

void matlinreg(double coeffs[2],
               double* arr,
               unsigned int N,
               unsigned int M,
               double* vec,
               double* weight)
{
    double* xarr = (double*)calloc((M + 1) * N, sizeof(double));
    double* xarrT = (double*)calloc((M + 1) * N, sizeof(double));
    if (xarrT == NULL || xarr == NULL)
    {
        fprintf(stderr, "Memory allocation failed in matlinreg\n");
        free(xarr);
        free(xarrT);
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = 0; j < N; j++)
    {
        xarr[j * (M + 1)] = 1.0;
        for (unsigned int i = 0; i < M; i++)
        {
            xarr[i + 1 + (M + 1) * j] = arr[i + j * M];
        }
    }

    memcpy(xarrT, xarr, (M + 1) * N * sizeof(double));
    transpone(xarrT, N, M + 1);

    // multiply by inverse of weigths
    for (unsigned int j = 0; j < M + 1; j++)
    {
        for (unsigned int i = 0; i < N; i++)
        {
            xarrT[i + N * j] *= weight[i] * weight[i];
        }
    }

    double* xtwx = mulmat(xarrT, M + 1, N, xarr, N, M + 1);
    if (xtwx != NULL)
    {
        invert(xtwx, M + 1, M + 1);

        memcpy(xarrT, xarr, (M + 1) * N * sizeof(double));
        transpone(xarrT, N, M + 1);
        double* pipe = mulmat(xtwx, M + 1, M + 1, xarrT, M + 1, N);
        if (pipe != NULL)
        {
            for (unsigned int j = 0; j < M + 1; j++)
            {
                for (unsigned int i = 0; i < N; i++)
                {
                    pipe[i + N * j] *= weight[i] * weight[i];
                }
            }

            double* cffs = mulmat(pipe, M + 1, N, vec, N, 1);
            coeffs[0] = cffs[0];
            coeffs[1] = cffs[1];
            free(cffs);
            free(pipe);
        }
    }
    free(xtwx);
    free(xarr);
    free(xarrT);
}

void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               const unsigned int npeaks)
{
    unsigned int fitwindow = (i > npeaks) ? npeaks : i;

    if (i >= fitwindow && maxvals->arrd != NULL && maxes->arr != NULL && subpos->arrd!= NULL)
    {
        unsigned int m = 0;

        double* x = (double*)calloc(fitwindow, sizeof(double));
        double* y = (double*)calloc(fitwindow, sizeof(double));
        double* w = (double*)calloc(fitwindow, sizeof(double));
        if (x == NULL || y == NULL || w == NULL)
        {
            fprintf(stderr, "Memory allocation failed in fitNpeaks\n");
            free(x);
            free(y);
            free(w);
            exit(EXIT_FAILURE);
        }
        for (unsigned int k = 0; k < fitwindow; k++)
        {
            y[m] = (double)maxes->arr[i - k] + subpos->arrd[i - k];
            x[m] = (double)k;
            w[m] = maxvals->arrd[i - k] * maxvals->arrd[i - k];
            m++;
        }
        if (m > 1)
        {
            double coeffs[2] = {0.0, 0.0};
            matlinreg(coeffs, x, m, 1, y, w);
            *a = coeffs[0];
            *b = coeffs[1];
        }
        free(x);
        free(y);
        free(w);
    }
}
