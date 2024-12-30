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
    double* tmp = calloc(N * M, sizeof(double));
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in transpone\n");
        return;
    }
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
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in invert\n");
        return;
    }
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
    {
        fprintf(stderr, "Matrix multiplication dimension mismatch\n");
        exit(-1);
    }
    double* tmp = calloc(N * T, sizeof(double));
    if (tmp == NULL)
    {
        fprintf(stderr, "Memory allocation failed in mulmat\n");
        return NULL;
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

// Perform linear regression
double* matlinreg(
    double* arr, unsigned int N, unsigned int M, double* vec, double* weight)
{
    double* coeffs = calloc(M + 1, sizeof(double));
    double* xarr = calloc((M + 1) * N, sizeof(double));
    double* xarrT = calloc((M + 1) * N, sizeof(double));
    if (xarrT == NULL || coeffs == NULL || xarr == NULL)
    {
        fprintf(stderr, "Memory allocation failed in matlinreg\n");
        free(coeffs);
        free(xarr);
        free(xarrT);
        return NULL;
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
    invert(xtwx, M + 1, M + 1);

    memcpy(xarrT, xarr, (M + 1) * N * sizeof(double));
    transpone(xarrT, N, M + 1);
    double* pipe = mulmat(xtwx, M + 1, M + 1, xarrT, M + 1, N);

    for (unsigned int j = 0; j < M + 1; j++)
    {
        for (unsigned int i = 0; i < N; i++)
        {
            pipe[i + N * j] *= weight[i] * weight[i];
        }
    }

    coeffs = mulmat(pipe, M + 1, N, vec, N, 1);

    free(pipe);
    free(xtwx);
    free(xarr);
    free(xarrT);

    return coeffs;
}

void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const unsigned int npeaks)
{
    unsigned int fitwindow = i > npeaks ? npeaks : i;

    if (i >= fitwindow)
    {
        unsigned int m = 0;

        double* x = calloc(fitwindow, sizeof(double));
        double* y = calloc(fitwindow, sizeof(double));
        double* w = calloc(fitwindow, sizeof(double));
        if (x == NULL || y == NULL || w == NULL)
        {
            fprintf(stderr, "Memory allocation failed in fitNpeaks\n");
            free(x);
            free(y);
            free(w);
            return;
        }
        for (unsigned int k = 0; k < fitwindow; k++)
        {
            y[m] = (double)maxes->arr[i - k];
            x[m] = (double)k;
            w[m] = maxvals->arrd[i - k] * maxvals->arrd[i - k];
            m++;
        }
        if (m > 1)
        {
            double* coeffs = matlinreg(x, m, 1, y, w);
            if (coeffs != NULL)
            {
                *a = coeffs[0];
                *b = coeffs[1];
                free(coeffs);
            }
        }
        free(x);
        free(y);
        free(w);
    }
}
