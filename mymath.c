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
    unsigned int Ncols;
    unsigned int Nrows;
};

void printmat(double* arr, unsigned int Nrows, unsigned int Ncols)
{
    printf("\n");
    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int i = 0; i < Ncols; i++)
        {
            printf("%8.2g", arr[i + Ncols * j]);
        }
        printf("\n");
    }
    printf("\n");
}

void transpone(double* arr, unsigned int Nrows, unsigned int Ncols)
{
    double* tmp = (double*)calloc(Nrows * Ncols, sizeof(double));
    if (tmp == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in transpone\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < Ncols; ++i)
    {
        for (unsigned int j = 0; j < Nrows; ++j)
        {
            tmp[j + i * Nrows] = arr[i + j * Ncols];
        }
    }
    memcpy(arr, tmp, sizeof(double) * Ncols * Nrows);
    free(tmp);
}

void invert(double* arr, unsigned int Nrows, unsigned int Ncols)
{
    // make matrix twice as wide
    double* tmp = (double*)calloc(Nrows * Ncols * 2, sizeof(double));
    if (tmp == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in invert\n");
        exit(EXIT_FAILURE);
    }
    unsigned int doubleNcols = Ncols * 2;
    for (unsigned int j = 0; j < Nrows; ++j)
    {
        for (unsigned int i = 0; i < Ncols; ++i)
        {
            tmp[i + j * doubleNcols] = arr[i + j * Nrows];
        }
    }
    // make second part diagonal 1
    for (unsigned int i = 0; i < Ncols; ++i)
    {
        tmp[i + Ncols + i * doubleNcols] = 1.0;
    }

    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int k = 0; k < Nrows; k++)
        {
            if (j != k)
            {
                double f = tmp[j + k * doubleNcols] / tmp[j + j * doubleNcols];
                for (unsigned int i = 0; i < doubleNcols; i++)
                {
                    tmp[i + k * doubleNcols] -= tmp[i + j * doubleNcols] * f;
                }
            }
        }
    }
    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int i = 0; i < Ncols; i++)
        {
            arr[i + j * Ncols] =
                tmp[i + Ncols + j * doubleNcols] / tmp[j + j * doubleNcols];
        }
    }
    free(tmp);
}

double* mulmat(const double* matrix0,
               unsigned int Nrows,
               unsigned int Ncols,
               const double* matrix1,
               unsigned int Mrows,
               unsigned int Mcols)
{
    if (Ncols != Mrows)
    {
        (void)fprintf(stderr, "Matrix multiplication dimension mismatch\n");
        exit(EXIT_FAILURE);
    }
    double* tmp = (double*)calloc(Nrows * Mcols, sizeof(double));
    if (tmp == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in mulmat\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int l = 0; l < Mcols; l++)
        {
            for (unsigned int i = 0; i < Ncols; i++)
            {
                tmp[l + j * Mcols] +=
                    matrix0[i + Ncols * j] * matrix1[i * Mcols + l];
            }
        }
    }
    return tmp;
}

void matlinreg(double coeffs[2],
               const double* xmat,
               unsigned int Nrows,
               unsigned int Ncols,
               double* vec,
               const double* weight)
{
    double* xarr = (double*)calloc((Ncols + 1) * Nrows, sizeof(double));
    double* xarrT = (double*)calloc((Ncols + 1) * Nrows, sizeof(double));
    if (xarrT == NULL || xarr == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in matlinreg\n");
        free(xarr);
        free(xarrT);
        exit(EXIT_FAILURE);
    }
    for (unsigned int j = 0; j < Nrows; j++)
    {
        xarr[j * (Ncols + 1)] = 1.0;
        for (unsigned int i = 0; i < Ncols; i++)
        {
            xarr[i + 1 + (Ncols + 1) * j] = xmat[i + j * Ncols];
        }
    }

    memcpy(xarrT, xarr, (Ncols + 1) * Nrows * sizeof(double));
    transpone(xarrT, Nrows, Ncols + 1);

    // multiply by inverse of weigths
    for (unsigned int j = 0; j < Ncols + 1; j++)
    {
        for (unsigned int i = 0; i < Nrows; i++)
        {
            xarrT[i + Nrows * j] *= weight[i] * weight[i];
        }
    }

    double* xtwx = mulmat(xarrT, Ncols + 1, Nrows, xarr, Nrows, Ncols + 1);
    if (xtwx != NULL)
    {
        invert(xtwx, Ncols + 1, Ncols + 1);

        memcpy(xarrT, xarr, (Ncols + 1) * Nrows * sizeof(double));
        transpone(xarrT, Nrows, Ncols + 1);
        double* pipe =
            mulmat(xtwx, Ncols + 1, Ncols + 1, xarrT, Ncols + 1, Nrows);
        if (pipe != NULL)
        {
            for (unsigned int j = 0; j < Ncols + 1; j++)
            {
                for (unsigned int i = 0; i < Nrows; i++)
                {
                    pipe[i + Nrows * j] *= weight[i] * weight[i];
                }
            }

            double* cffs = mulmat(pipe, Ncols + 1, Nrows, vec, Nrows, 1);
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
                unsigned int Npoints,
                const double* vec,
                const double* weight_arr)
{
    double Sum_w = 0.0;
    double Sum_wx = 0.0;
    double Sum_wy = 0.0;
    double Sum_wxx = 0.0;
    double Sum_wxy = 0.0;

    for (unsigned int i = 0; i < Npoints; i++)
    {
        const double x_val = xmat[i];
        const double y_val = vec[i];
        const double weight = weight_arr[i] * weight_arr[i];
        Sum_w += weight;
        Sum_wx += weight * x_val;
        Sum_wy += weight * y_val;
        Sum_wxx += weight * x_val * x_val;
        Sum_wxy += weight * x_val * y_val;
    }
    double denom = Sum_w * Sum_wxx - Sum_wx * Sum_wx;
    coeffs[0] = 0.0;
    coeffs[1] = 0.0;
    if (denom != 0.0)
    {
        coeffs[1] = (Sum_w * Sum_wxy - Sum_wx * Sum_wy) / denom;
        coeffs[0] = (Sum_wy - coeffs[1] * Sum_wx) / Sum_w;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }
}

// Accumulate sufficient stats for weighted y = par_a + par_b x
static int fastlinreg_sufficient_stats(
    long double* Sum_w,
    long double* Sum_wx,
    long double* Sum_wy,
    long double* Sum_wxx,
    long double* Sum_wxy,
    long double* Sum_w2,
    const unsigned int cur_pos,
    const unsigned int fitwindow,
    const struct myarr* maxvals,
    const struct myarr* maxes,
    const struct myarr* subpos,
    int skip_outliers,
    double par_a,
    double par_b,
    double stdev_threshold // absolute threshold; compare with residual
)
{
    *Sum_w = *Sum_wx = *Sum_wy = *Sum_wxx = *Sum_wxy = *Sum_w2 = 0.0L;

    for (unsigned int k = 0; k < fitwindow; ++k)
    {
        unsigned int idx = cur_pos - k;

        const double x_val = (double)k;
        const double y_val = (double)maxes->arr[idx] + subpos->arrd[idx];
        const double weight = maxvals->arrd[idx];
        if (!(weight > 0.0))
        {
            continue;
        } // skip non-positive weights

        if (skip_outliers)
        {
            double deviation = y_val - (par_a + par_b * x_val);
            if (fabs(deviation) > stdev_threshold)
            {
                continue;
            }
        }

        long double w_l = (long double)weight;
        long double x_l = (long double)x_val;
        long double y_l = (long double)y_val;

        *Sum_w += w_l;
        *Sum_wx += w_l * x_l;
        *Sum_wy += w_l * y_l;
        *Sum_wxx += w_l * x_l * x_l;
        *Sum_wxy += w_l * x_l * y_l;
        *Sum_w2 += w_l * w_l;
    }

    // Return success if we have at least two distinct weighted points
    return (*Sum_w > 0.0L) ? 1 : 0;
}

static int solve_weighted_line(double* par_a,
                               double* par_b,
                               long double Sum_w,
                               long double Sum_wx,
                               long double Sum_wy,
                               long double Sum_wxx,
                               long double Sum_wxy)
{
    long double denom = Sum_w * Sum_wxx - Sum_wx * Sum_wx;
    const long double eps = 1e-18L;
    if (fabsl(denom) < eps || Sum_w <= 0.0L)
    {
        return 0; // degenerate
    }
    long double b_L = (Sum_w * Sum_wxy - Sum_wx * Sum_wy) / denom;
    long double a_L = (Sum_wy - b_L * Sum_wx) / Sum_w;
    *par_a = (double)a_L;
    *par_b = (double)b_L;
    return 1;
}

// Compute weighted residual SSE for given (par_a,par_b) over the window
static long double compute_weighted_SSE(double par_a,
                                        double par_b,
                                        const unsigned int cur_pos,
                                        const unsigned int fitwindow,
                                        const struct myarr* maxvals,
                                        const struct myarr* maxes,
                                        const struct myarr* subpos,
                                        long double* Sum_w_out,
                                        long double* Sum_w2_out)
{
    long double SSE = 0.0L;
    long double Sum_w = 0.0L;
    long double Sum_w2 = 0.0L;

    for (unsigned int k = 0; k < fitwindow; ++k)
    {
        unsigned int idx = cur_pos - k;
        double x_val = (double)k;
        double y_val = (double)maxes->arr[idx] + subpos->arrd[idx];
        double weight = maxvals->arrd[idx];
        if (!(weight > 0.0))
        {
            continue;
        }

        double deviation = y_val - (par_a + par_b * x_val);
        long double w_l = (long double)weight;
        SSE += w_l * (long double)(deviation * deviation);
        Sum_w += w_l;
        Sum_w2 += w_l * w_l;
    }
    if (Sum_w_out)
    {
        *Sum_w_out = Sum_w;
    }
    if (Sum_w2_out)
    {
        *Sum_w2_out = Sum_w2;
    }
    return SSE;
}

void fitNpeaks(double* par_a,
               double* par_b,
               const unsigned int cur_pos,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               const unsigned int npeaks,
               const double SDthreshold)
{
    unsigned int fitwindow = (cur_pos > npeaks) ? npeaks : cur_pos;

    if (fitwindow > 1 && cur_pos >= fitwindow && maxvals && maxvals->arrd &&
        maxes && maxes->arr && subpos && subpos->arrd)
    {
        // Pass 1: fit using all points (positive weights)
        long double Sum_w;
        long double Sum_wx;
        long double Sum_wy;
        long double Sum_wxx;
        long double Sum_wxy;
        long double Sum_w2;
        int ok1 = fastlinreg_sufficient_stats(&Sum_w,
                                              &Sum_wx,
                                              &Sum_wy,
                                              &Sum_wxx,
                                              &Sum_wxy,
                                              &Sum_w2,
                                              cur_pos,
                                              fitwindow,
                                              maxvals,
                                              maxes,
                                              subpos,
                                              /*skip_outliers=*/0,
                                              /*par_a=*/0.0,
                                              /*par_b=*/0.0,
                                              /*stdev_threshold=*/0.0);
        if (!ok1)
        {
            // No valid data
            *par_a = 0.0;
            *par_b = 0.0;
            return;
        }

        double a_1 = 0.0;
        double b_1 = 0.0;
        if (!solve_weighted_line(
                &a_1, &b_1, Sum_w, Sum_wx, Sum_wy, Sum_wxx, Sum_wxy))
        {
            (void)fprintf(stderr, "Degenerate data in initial fit\n");
            *par_a = 0.0;
            *par_b = 0.0;
            return;
        }

        // Compute weighted residual SSE and sigma
        long double Sum_w_res;
        long double Sum_w2_res;
        long double SSE = compute_weighted_SSE(a_1,
                                               b_1,
                                               cur_pos,
                                               fitwindow,
                                               maxvals,
                                               maxes,
                                               subpos,
                                               &Sum_w_res,
                                               &Sum_w2_res);

        const long double two = 2.0L;
        // degrees of freedom
        long double degFree = Sum_w_res - two;
        double sigma = 0.0;
        if (degFree > 0.0L && SSE >= 0.0L)
        {
            sigma = sqrt((double)(SSE / degFree));
        }

        double thresh = SDthreshold * sigma;

        // Pass 2: re-fit excluding outliers |residual| > thresh
        ok1 = fastlinreg_sufficient_stats(&Sum_w,
                                          &Sum_wx,
                                          &Sum_wy,
                                          &Sum_wxx,
                                          &Sum_wxy,
                                          &Sum_w2,
                                          cur_pos,
                                          fitwindow,
                                          maxvals,
                                          maxes,
                                          subpos,
                                          /*skip_outliers=*/(thresh > 0.0),
                                          a_1,
                                          b_1,
                                          thresh);

        double a_2 = a_1;
        double b_2 = b_1;
        if (ok1 && solve_weighted_line(
                       &a_2, &b_2, Sum_w, Sum_wx, Sum_wy, Sum_wxx, Sum_wxy))
        {
            *par_a = a_2;
            *par_b = b_2;
        }
        else
        {
            // Fall back to initial fit if second pass degenerates
            *par_a = a_1;
            *par_b = b_1;
        }
    }
}
