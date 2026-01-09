#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "myarr.h"
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
        (void)fprintf(stderr, "Memory allocation failed in transpone\n");
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
        (void)fprintf(stderr, "Memory allocation failed in invert\n");
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

double* mulmat(const double* matrix,
               unsigned int N,
               unsigned int M,
               const double* vector,
               unsigned int S,
               unsigned int T)
{
    if (M != S)
    {
        (void)fprintf(stderr, "Matrix multiplication dimension mismatch\n");
        exit(EXIT_FAILURE);
    }
    double* tmp = (double*)calloc(N * T, sizeof(double));
    if (tmp == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in mulmat\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = 0; j < N; j++)
    {
        for (unsigned int l = 0; l < T; l++)
        {
            for (unsigned int i = 0; i < M; i++)
            {
                tmp[l + j * T] += matrix[i + M * j] * vector[i * T + l];
            }
        }
    }
    return tmp;
}

void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int N,
               unsigned int M,
               double* vec,
               const double* weight)
{
    double* xarr = (double*)calloc((M + 1) * N, sizeof(double));
    double* xarrT = (double*)calloc((M + 1) * N, sizeof(double));
    if (xarrT == NULL || xarr == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in matlinreg\n");
        free(xarr);
        free(xarrT);
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = 0; j < N; j++)
    {
        xarr[j * (M + 1)] = 1.0;
        for (unsigned int i = 0; i < M; i++)
        {
            xarr[i + 1 + (M + 1) * j] = xmat[i + j * M];
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

void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int N,
                const double* vec,
                const double* weight)
{
    double Sw = 0.0;
    double Swx = 0.0;
    double Swy = 0.0;
    double Swxx = 0.0;
    double Swxy = 0.0;

    for (unsigned int i = 0; i < N; i++)
    {
        double x = xmat[i];
        double y = vec[i];
        double w = weight[i] * weight[i];
        Sw += w;
        Swx += w * x;
        Swy += w * y;
        Swxx += w * x * x;
        Swxy += w * x * y;
    }
    double denom = Sw * Swxx - Swx * Swx;
    coeffs[0] = 0.0;
    coeffs[1] = 0.0;
    if (denom != 0.0)
    {
        coeffs[1] = (Sw * Swxy - Swx * Swy) / denom;
        coeffs[0] = (Swy - coeffs[1] * Swx) / Sw;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }
}

// Accumulate sufficient stats for weighted y = a + b x
static int fastlinreg_sufficient_stats(
    long double* Sw,
    long double* Swx,
    long double* Swy,
    long double* Swxx,
    long double* Swxy,
    long double* Sw2,
    const unsigned int i,
    const unsigned int fitwindow,
    const struct myarr* maxvals,
    const struct myarr* maxes,
    const struct myarr* subpos,
    int skip_outliers,
    double a,
    double b,
    double stdev_threshold // absolute threshold; compare with residual
)
{
    *Sw = *Swx = *Swy = *Swxx = *Swxy = *Sw2 = 0.0L;

    for (unsigned int k = 0; k < fitwindow; ++k)
    {
        unsigned int idx = i - k;

        // x = k
        double x = (double)k;

        // y = maxes->arr[idx] + subpos->arrd[idx]
        double y = (double)maxes->arr[idx] + subpos->arrd[idx];

        // w = maxvals->arrd[idx]
        double w = maxvals->arrd[idx];
        if (!(w > 0.0))
        {
            continue;
        } // skip non-positive weights

        if (skip_outliers)
        {
            double r = y - (a + b * x);
            if (fabs(r) > stdev_threshold)
            {
                continue;
            }
        }

        long double wl = (long double)w;
        long double xl = (long double)x;
        long double yl = (long double)y;

        *Sw += wl;
        *Swx += wl * xl;
        *Swy += wl * yl;
        *Swxx += wl * xl * xl;
        *Swxy += wl * xl * yl;
        *Sw2 += wl * wl;
    }

    // Return success if we have at least two distinct weighted points
    return (*Sw > 0.0L) ? 1 : 0;
}

static int solve_weighted_line(double* a,
                               double* b,
                               long double Sw,
                               long double Swx,
                               long double Swy,
                               long double Swxx,
                               long double Swxy)
{
    long double denom = Sw * Swxx - Swx * Swx;
    const long double eps = 1e-18L;
    if (fabsl(denom) < eps || Sw <= 0.0L)
    {
        return 0; // degenerate
    }
    long double bL = (Sw * Swxy - Swx * Swy) / denom;
    long double aL = (Swy - bL * Swx) / Sw;
    *a = (double)aL;
    *b = (double)bL;
    return 1;
}

// Compute weighted residual SSE for given (a,b) over the window
static long double compute_weighted_SSE(double a,
                                        double b,
                                        const unsigned int i,
                                        const unsigned int fitwindow,
                                        const struct myarr* maxvals,
                                        const struct myarr* maxes,
                                        const struct myarr* subpos,
                                        long double* Sw_out,
                                        long double* Sw2_out)
{
    long double SSE = 0.0L;
    long double Sw = 0.0L;
    long double Sw2 = 0.0L;

    for (unsigned int k = 0; k < fitwindow; ++k)
    {
        unsigned int idx = i - k;
        double x = (double)k;
        double y = (double)maxes->arr[idx] + subpos->arrd[idx];
        double w = maxvals->arrd[idx];
        if (!(w > 0.0))
        {
            continue;
        }

        double r = y - (a + b * x);
        long double wl = (long double)w;
        SSE += wl * (long double)(r * r);
        Sw += wl;
        Sw2 += wl * wl;
    }
    if (Sw_out)
    {
        *Sw_out = Sw;
    }
    if (Sw2_out)
    {
        *Sw2_out = Sw2;
    }
    return SSE;
}

void fitNpeaks(double* a,
               double* b,
               const unsigned int i,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               const unsigned int npeaks,
               const double SDthreshold)
{
    unsigned int fitwindow = (i > npeaks) ? npeaks : i;

    if (fitwindow > 1 && i >= fitwindow && maxvals && maxvals->arrd && maxes &&
        maxes->arr && subpos && subpos->arrd)
    {
        // Pass 1: fit using all points (positive weights)
        long double Sw;
        long double Swx;
        long double Swy;
        long double Swxx;
        long double Swxy;
        long double Sw2;
        int ok1 = fastlinreg_sufficient_stats(&Sw,
                                              &Swx,
                                              &Swy,
                                              &Swxx,
                                              &Swxy,
                                              &Sw2,
                                              i,
                                              fitwindow,
                                              maxvals,
                                              maxes,
                                              subpos,
                                              /*skip_outliers=*/0,
                                              /*a=*/0.0,
                                              /*b=*/0.0,
                                              /*stdev_threshold=*/0.0);
        if (!ok1)
        {
            // No valid data
            *a = 0.0;
            *b = 0.0;
            return;
        }

        double a1 = 0.0;
        double b1 = 0.0;
        if (!solve_weighted_line(&a1, &b1, Sw, Swx, Swy, Swxx, Swxy))
        {
            (void)fprintf(stderr, "Degenerate data in initial fit\n");
            *a = 0.0;
            *b = 0.0;
            return;
        }

        // Compute weighted residual SSE and sigma
        long double Sw_res;
        long double Sw2_res;
        long double SSE = compute_weighted_SSE(
            a1, b1, i, fitwindow, maxvals, maxes, subpos, &Sw_res, &Sw2_res);

        const long double two = 2.0L;
        // degrees of freedom
        long double df = Sw_res - two;
        double sigma = 0.0;
        if (df > 0.0L && SSE >= 0.0L)
        {
            sigma = sqrt((double)(SSE / df));
        }

        double thresh = SDthreshold * sigma;

        // Pass 2: re-fit excluding outliers |residual| > thresh
        ok1 = fastlinreg_sufficient_stats(&Sw,
                                          &Swx,
                                          &Swy,
                                          &Swxx,
                                          &Swxy,
                                          &Sw2,
                                          i,
                                          fitwindow,
                                          maxvals,
                                          maxes,
                                          subpos,
                                          /*skip_outliers=*/(thresh > 0.0),
                                          a1,
                                          b1,
                                          thresh);

        double a2 = a1;
        double b2 = b1;
        if (ok1 && solve_weighted_line(&a2, &b2, Sw, Swx, Swy, Swxx, Swxy))
        {
            *a = a2;
            *b = b2;
        }
        else
        {
            // Fall back to initial fit if second pass degenerates
            *a = a1;
            *b = b1;
        }
    }
}
