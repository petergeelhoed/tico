#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "myarr.h"
#include "mydefs.h"
#include "mymath.h"

struct mymat
{
    double* mat;
    unsigned int Ncols;
    unsigned int Nrows;
};

size_t getmaxpos(const int* array, size_t ArrayLength)
{
    int maxTick = -INT_MAX;
    size_t maxIndex = 0;
    for (size_t j = 0; j < ArrayLength; j++)
    {
        if (array[j] > maxTick)
        {
            maxTick = array[j];
            maxIndex = j;
        }
    }
    return maxIndex;
}

void linreg(const double* xarr,
            const double* yarr,
            size_t ArrayLength,
            double* intercept,
            double* slope,
            double* stdev)
{
    double sumX = 0;
    double sumY = 0;
    double sumXx = 0;
    double sumXy = 0;
    double sumYy = 0;
    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        sumY += yarr[i];
        sumXx += xarr[i] * xarr[i];
        sumX += xarr[i];
        sumXy += xarr[i] * yarr[i];
        sumYy += yarr[i] * yarr[i];
    }

    *intercept = (sumY * sumXx - sumX * sumXy) /
                 ((double)ArrayLength * sumXx - sumX * sumX);
    *slope = ((double)ArrayLength * sumXy - sumX * sumY) /
             ((double)ArrayLength * sumXx - sumX * sumX);
    *stdev = sqrt((sumYy - 2 * (*intercept) * sumY - 2 * (*slope) * sumXy +
                   2 * (*intercept) * (*slope) * sumX +
                   (*intercept) * (*intercept) * (double)ArrayLength +
                   (*slope) * (*slope) * sumXx) /
                  (double)ArrayLength);
}

void transpone(double* arr, unsigned int Nrows, unsigned int Ncols)
{
    double* transposeBuffer = (double*)calloc(Nrows * Ncols, sizeof(double));
    if (transposeBuffer == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in transpone\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < Ncols; ++i)
    {
        for (unsigned int j = 0; j < Nrows; ++j)
        {
            transposeBuffer[j + i * Nrows] = arr[i + j * Ncols];
        }
    }
    memcpy(arr, transposeBuffer, sizeof(double) * Ncols * Nrows);
    free(transposeBuffer);
}

void invert(double* arr, unsigned int Nrows, unsigned int Ncols)
{
    // make matrix twice as wide
    double* augmentedMatrix =
        (double*)calloc(Nrows * Ncols * 2, sizeof(double));
    if (augmentedMatrix == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in invert\n");
        exit(EXIT_FAILURE);
    }
    unsigned int doubleNcols = Ncols * 2;
    for (unsigned int j = 0; j < Nrows; ++j)
    {
        for (unsigned int i = 0; i < Ncols; ++i)
        {
            augmentedMatrix[i + j * doubleNcols] = arr[i + j * Nrows];
        }
    }
    // make second part diagonal 1
    for (unsigned int i = 0; i < Ncols; ++i)
    {
        augmentedMatrix[i + Ncols + i * doubleNcols] = 1.0;
    }

    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int k = 0; k < Nrows; k++)
        {
            if (j != k)
            {
                double factor = augmentedMatrix[j + k * doubleNcols] /
                                augmentedMatrix[j + j * doubleNcols];
                for (unsigned int i = 0; i < doubleNcols; i++)
                {
                    augmentedMatrix[i + k * doubleNcols] -=
                        augmentedMatrix[i + j * doubleNcols] * factor;
                }
            }
        }
    }
    for (unsigned int j = 0; j < Nrows; j++)
    {
        for (unsigned int i = 0; i < Ncols; i++)
        {
            arr[i + j * Ncols] = augmentedMatrix[i + Ncols + j * doubleNcols] /
                                 augmentedMatrix[j + j * doubleNcols];
        }
    }
    free(augmentedMatrix);
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
    double* resultMatrix = (double*)calloc(Nrows * Mcols, sizeof(double));
    if (resultMatrix == NULL)
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
                resultMatrix[l + j * Mcols] +=
                    matrix0[i + Ncols * j] * matrix1[i * Mcols + l];
            }
        }
    }
    return resultMatrix;
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
    // multiply by inverse of weights
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
        free(xtwx);
    }
    free(xarr);
    free(xarrT);
}

void fastlinreg(double coeffs[2],
                const double* xmat,
                unsigned int Npoints,
                const double* vec,
                const double* weightArr)
{
    double sumWeight = 0.0;
    double sumWeightX = 0.0;
    double sumWeightY = 0.0;
    double sumWeightXX = 0.0;
    double sumWeightXY = 0.0;

    for (unsigned int i = 0; i < Npoints; i++)
    {
        const double xVal = xmat[i];
        const double yVal = vec[i];
        const double weight = weightArr[i] * weightArr[i];
        sumWeight += weight;
        sumWeightX += weight * xVal;
        sumWeightY += weight * yVal;
        sumWeightXX += weight * xVal * xVal;
        sumWeightXY += weight * xVal * yVal;
    }
    double denom = sumWeight * sumWeightXX - sumWeightX * sumWeightX;
    coeffs[0] = 0.0;
    coeffs[1] = 0.0;
    if (fabs(denom) > DOUBLE_LIMIT)
    {
        coeffs[1] = (sumWeight * sumWeightXY - sumWeightX * sumWeightY) / denom;
        coeffs[0] = (sumWeightY - coeffs[1] * sumWeightX) / sumWeight;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }
}

// Accumulate sufficient stats for weighted y = intercept + slope x
static int fastlinreg_sufficient_stats(
    long double* Sum_w,
    long double* Sum_wx,
    long double* Sum_wy,
    long double* Sum_wxx,
    long double* Sum_wxy,
    long double* Sum_w2,
    const unsigned int curPos,
    const unsigned int fitwindow,
    const struct myarr* maxvals,
    const struct myarr* maxes,
    const struct myarr* subpos,
    int skipOutliers,
    double intercept,
    double slope,
    double stdevThreshold // absolute threshold; compare with residual
)
{
    *Sum_w = *Sum_wx = *Sum_wy = *Sum_wxx = *Sum_wxy = *Sum_w2 = 0.0L;

    for (unsigned int k = 0; k < fitwindow; ++k)
    {
        unsigned int idx = curPos - k;

        const double xVal = (double)k;
        const double yVal = (double)maxes->arr[idx] + subpos->arrd[idx];
        const double weight = maxvals->arrd[idx];
        if (!(weight > 0.0))
        {
            continue;
        } // skip non-positive weights

        if (skipOutliers)
        {
            double deviation = yVal - (intercept + slope * xVal);
            if (fabs(deviation) > stdevThreshold)
            {
                continue;
            }
        }

        long double weightLong = (long double)weight;
        long double xValueLong = (long double)xVal;
        long double yValueLong = (long double)yVal;

        *Sum_w += weightLong;
        *Sum_wx += weightLong * xValueLong;
        *Sum_wy += weightLong * yValueLong;
        *Sum_wxx += weightLong * xValueLong * xValueLong;
        *Sum_wxy += weightLong * xValueLong * yValueLong;
        *Sum_w2 += weightLong * weightLong;
    }

    // Return success if we have at least two distinct weighted points
    return (*Sum_w > 0.0L) ? 1 : 0;
}

static int solve_weighted_line(double* intercept,
                               double* slope,
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
    *intercept = (double)a_L;
    *slope = (double)b_L;
    return 1;
}

// Compute weighted residual SSE for given (intercept,slope) over the window
static long double compute_weighted_SSE(double intercept,
                                        double slope,
                                        const unsigned int curPos,
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
        unsigned int idx = curPos - k;
        double xVal = (double)k;
        double yVal = (double)maxes->arr[idx] + subpos->arrd[idx];
        double weight = maxvals->arrd[idx];
        if (!(weight > 0.0))
        {
            continue;
        }

        double deviation = yVal - (intercept + slope * xVal);
        long double weightLong = (long double)weight;
        SSE += weightLong * (long double)(deviation * deviation);
        Sum_w += weightLong;
        Sum_w2 += weightLong * weightLong;
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

void fitNpeaks(double* intercept,
               double* slope,
               const unsigned int curPos,
               const struct myarr* maxvals,
               const struct myarr* maxes,
               const struct myarr* subpos,
               const unsigned int npeaks,
               const double SDthreshold)
{
    unsigned int fitwindow = (curPos > npeaks) ? npeaks : curPos;

    if (fitwindow > 1 && curPos >= fitwindow && maxvals && maxvals->arrd &&
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
                                              curPos,
                                              fitwindow,
                                              maxvals,
                                              maxes,
                                              subpos,
                                              /*skipOutliers=*/0,
                                              /*intercept=*/0.0,
                                              /*slope=*/0.0,
                                              /*stdevThreshold=*/0.0);
        if (!ok1)
        {
            // No valid data
            *intercept = 0.0;
            *slope = 0.0;
            return;
        }

        double interceptPass1 = 0.0;
        double slopePass1 = 0.0;
        if (!solve_weighted_line(&interceptPass1,
                                 &slopePass1,
                                 Sum_w,
                                 Sum_wx,
                                 Sum_wy,
                                 Sum_wxx,
                                 Sum_wxy))
        {
            (void)fprintf(stderr, "Degenerate data in initial fit\n");
            *intercept = 0.0;
            *slope = 0.0;
            return;
        }

        // Compute weighted residual SSE and sigma
        long double Sum_w_res;
        long double Sum_w2_res;
        long double SSE = compute_weighted_SSE(interceptPass1,
                                               slopePass1,
                                               curPos,
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
                                          curPos,
                                          fitwindow,
                                          maxvals,
                                          maxes,
                                          subpos,
                                          /*skipOutliers=*/(thresh > 0.0),
                                          interceptPass1,
                                          slopePass1,
                                          thresh);

        double interceptPass2 = interceptPass1;
        double slopePass2 = slopePass1;
        if (ok1 && solve_weighted_line(&interceptPass2,
                                       &slopePass2,
                                       Sum_w,
                                       Sum_wx,
                                       Sum_wy,
                                       Sum_wxx,
                                       Sum_wxy))
        {
            *intercept = interceptPass2;
            *slope = slopePass2;
        }
        else
        {
            // Fall back to initial fit if second pass degenerates
            *intercept = interceptPass1;
            *slope = slopePass1;
        }
    }
}

int nearlyEqual(double number0, double number1)
{
    if (isnan(number0) || isnan(number1))
    {
        return 0; // NaNs are never equal
    }
    if (isinf(number0) || isinf(number1))
    {
        return 0;
    }

    const double absEps = DOUBLE_LIMIT;
    const double relEps = DOUBLE_LIMIT;

    const double diff = fabs(number0 - number1);
    if (diff <= absEps)
    {
        return 1;
    }

    const double maxab = fmax(fabs(number0), fabs(number1));
    return diff <= relEps * maxab;
}

int shiftHalf(size_t value, size_t ArrayLength)
{
    return ((int)value + (int)ArrayLength / 2) % (int)(ArrayLength) -
           (int)(ArrayLength / 2);
}

// mods an int with a signed int, but makes sure the result is positive
size_t modSigned(int value, size_t ArrayLength)
{
    return (size_t)(((value % (int)ArrayLength) + (int)ArrayLength) %
                    (int)ArrayLength);
}
