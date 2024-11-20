#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void printmat(double* arr, unsigned int N, unsigned int M)
{
        printf("\n");
    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int i = 0; i < M; i++)
        {
            printf("%8.2f",arr[i+M*j]);
        }
        printf("\n");
    }
        printf("\n");
}

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

double* matlinreg(double* arr,
               unsigned int N,
               unsigned int M,
               double* vec,
               double* weight)
{
    double* coeffs = calloc(M+1, sizeof(double));
    double* xarr = calloc((M+1)*N, sizeof(double));
    //printmat(arr, N, M);
    for (unsigned int j = 0; j < N; j++)
    {
        xarr[j*(M+1)] = 1.0;
        for (unsigned int i = 0; i < M; i++)
        {
            xarr[i+1+(M+1)*j] = arr[i+j*M];
        }
    }
    //printmat(xarr, N, M+1);


    double* xarrT = calloc((M+1)*N, sizeof(double));
    memcpy(xarrT, xarr,(M+1)*N*sizeof(double));
    transpone(xarrT,N,M+1);
    
    //printmat(xarrT,M+1,N);


    // multiply by inverse of weigths
    for (unsigned int j = 0; j < M+1; j++)
    {
        for (unsigned int i = 0; i < N; i++)
        {
            xarrT[i+N*j] *= weight[i];
        }
    }
    //printmat(xarrT,M+1,N);
    //printmat(xarr,N,M+1);
    double* xtwx =  mulmat(xarrT,M+1,N,xarr, N, M+1);
    //printmat(xtwx,M+1,M+1);
    invert (xtwx,M+1,M+1);
    //printmat(xtwx,M+1,M+1);
    
    // this was overwritten
    memcpy(xarrT, xarr,(M+1)*N*sizeof(double));
    transpone(xarrT,N,M+1);
    double *pipe = mulmat(xtwx,M+1,M+1, xarrT,M+1,N);
    //printmat(pipe,M+1,N);
    for (unsigned int j = 0; j < M+1; j++)
    {
        for (unsigned int i = 0; i < N; i++)
        {
            pipe[i+N*j] *= weight[i];
        }
    }
    //printmat(pipe,M+1,N);

    coeffs = mulmat(pipe,M+1,N, vec,N,1);

    free(pipe);
    free(xtwx);
    free(xarr);
    free(xarrT);

    //printmat(coeffs,M+1,1);

    return coeffs;
}
void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const int* maxvals,
               const int* maxes,
               const int cvalue,
               const unsigned int npeaks)
{
    unsigned int fitwindow = i > npeaks ? npeaks : i;

    if (i >= fitwindow)
    {
        unsigned int m = 0;

        double* x = calloc(fitwindow, sizeof(double));
        double* y = calloc(fitwindow, sizeof(double));
        double* w = calloc(fitwindow, sizeof(double));
        for (unsigned int k = 0; k < fitwindow; k++)
        {
            if (maxvals[i - k] > cvalue)
            {
                y[m] = (double)maxes[i-k];
                x[m] = (double)k;
                w[m] = 1.0;
                m++;
            }
        }
        if (m > 1)
        {
            double *coeffs = matlinreg(x,m,1,y,w);
            *a = coeffs[0];
            *b = coeffs[1];
        }
    }
}
