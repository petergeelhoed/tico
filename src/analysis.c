#include "analysis.h"

#include "config.h"
#include "crosscorint.h"
#include "mydefs.h"
#include "mymath.h"
#include "mysync.h"
#include "parseargs.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Calculates total statistics and prints unweighted raw rate.
 *
 * @param count Number of samples.
 * @param maxpos Array of maximum positions.
 * @param ArrayLength Length of the array.
 * @param threshold Threshold value for filtering.
 * @param rate Sampling rate.
 */
static void calculateTotal(size_t count,
                           double* maxpos,
                           size_t ArrayLength,
                           double threshold,
                           double rate)
{
    double slope = 0.0;
    double intercept = 0.0;
    double stdev = 0.0;
    double* xarr = calloc(count, sizeof(double));
    if (xarr == NULL)
    {
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < count; ++i)
    {
        xarr[i] = (double)i;
    }

    linreg(xarr, maxpos, count, &intercept, &slope, &stdev);

    (void)fprintf(stderr,
                  "unweighted raw rate: %f s/d, %zu samples σ=%.2gms\n",
                  -slope * SECS_DAY / (double)ArrayLength,
                  count,
                  stdev * THOUSAND / rate);
    unsigned int maxIndex = 0;

    double deviation;

    for (unsigned int j = 0; j < 3; ++j)
    {
        for (size_t i = 0; i < count; ++i)
        {
            deviation =
                fabs((maxpos[i] - (intercept + xarr[i] * slope)) / stdev);
            if (deviation < threshold)
            {
                maxpos[maxIndex] = maxpos[i];
                xarr[maxIndex] = xarr[i];
                maxIndex++;
            }
        }
        count = maxIndex;
        maxIndex = 0;

        linreg(xarr, maxpos, count, &intercept, &slope, &stdev);

        (void)fprintf(stderr,
                      "after %.1fc3 (%.2gms) removal: %.2f s/d, %zu samples\n",
                      threshold,
                      stdev * THOUSAND / rate,
                      -slope * SECS_DAY / (double)ArrayLength,
                      count);
        threshold /= 2;
    }

    free(xarr);
}

void calculateTotalFromFile(size_t count,
                            FILE* rawfile,
                            size_t ArrayLength,
                            double threshold,
                            double rate)
{
    errno = 0;
    if (fseek(rawfile, 0, SEEK_SET) == -1)
    {
        (void)fprintf(stderr, "fseek failed with %d\n", errno);
        return;
    }
    unsigned int index = 0;
    double* all = (double*)calloc(count, sizeof(double));
    if (all)
    {
        size_t bufsize = BUFFER_SIZE;
        char* buf = (char*)malloc(bufsize * sizeof(char));
        if (buf == NULL)
        {
            (void)fprintf(stderr, "Cannot allocate memory for totalFromFile\n");
            free(all);
            return;
        }
        while (getline(&buf, &bufsize, rawfile) > 0 && index < count)
        {
            if (buf[0] != '#')
            {
                all[index++] = getDouble(buf);
            }
        }
        free(buf);
        calculateTotal(count, all, ArrayLength, threshold, rate);
        free(all);
    }
    else
    {
        (void)fprintf(stderr, "Cannot allocate memory for totalFromFile\n");
    }
}

double getBeatError(const struct myarr* totaltick, double rate, int verbose)
{
    size_t ArrayLength = totaltick->ArrayLength;
    int* cross = malloc((ArrayLength / 2) * sizeof(int));
    if (cross == NULL)
    {
        (void)fprintf(stderr, "Cannot allocate memory for getBeatError\n");
        exit(EXIT_FAILURE);
    }
    crosscorint(ArrayLength / 2,
                totaltick->arr,
                totaltick->arr + ArrayLength / 2,
                cross);
    if (verbose)
    {
        syncwrite(cross, ArrayLength / 2, "beaterror");
        syncwrite(totaltick->arr, ArrayLength / 2, "t1");
        syncwrite(totaltick->arr + ArrayLength / 2, ArrayLength / 2, "t2");
    }
    size_t postick = getmaxpos(cross, ArrayLength / 2);
    free(cross);
    return shiftHalf(postick, ArrayLength / 2) * THOUSAND / rate;
}
