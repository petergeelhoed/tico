#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
#include "mymath.h"
#include "mysync.h"
#include "parseargs.h"
#include "resources.h"

/* Prints header on line or at the top */
void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds)
{
    if (everyline)
    {
        static char tmp[DOUBLE_BUF];

        (void)sprintf(tmp, "%4.2f", beatError);
        (void)sprintf(tmp + BEAT_WIDTH - 1, "ms%+5.1f", fittedRate);
        (void)sprintf(tmp + BEAT_WIDTH + RATE_WIDTH - 2, "s/d");
        tmp[EVERY_WIDTH] = '\0';
        (void)fprintf(stderr, "%s", tmp);
    }
    else
    {
        (void)fprintf(
            stderr,
            "\033[s\033[2;0H\033[0K%8.2fms   %9.1fs/d   %12.2fs\033[u",
            beatError,
            fittedRate,
            seconds);
    }
}

void printspaces(int maxpos,
                 double hexvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double avg_pos,
                 unsigned int correlationThreshold)
{
    while (maxpos < (int)mod)
    {
        maxpos += (int)mod;
    }
    while (avg_pos < (double)mod)
    {
        avg_pos += (double)mod;
    }
    const unsigned int default_columns = 80;
    columns = columns > MAX_COLUMNS ? default_columns : columns;
    size_t width = (size_t)modSigned(maxpos, mod) * columns / mod;
    size_t widtha =
        (size_t)modSigned((int)lround(avg_pos), mod) * columns / mod;

    char spaces[MAX_COLUMNS];
    memset(spaces, ' ', width);
    spaces[width] = '\0';
    if (widtha < width)
    {
        spaces[widtha] = '|';
    }

    (void)fprintf(stderr,
                  "%s%s%X\033[0m",
                  spaces,
                  (unsigned int)hexvalue < correlationThreshold ? "\033[31m"
                                                                : "\033[32m",
                  (int)hexvalue);

    memset(spaces, ' ', width);
    if (widtha > width)
    {
        spaces[widtha - width - 1] = '|';
        spaces[widtha - width] = '\0';
        (void)fprintf(stderr, "%s", spaces);
    }
    (void)fprintf(stderr, "\n");
}

void writefile(FILE* filePtr, int* array, unsigned int ArrayLength)
{
    if (filePtr)
    {
        for (unsigned int j = 0; j < ArrayLength; j++)
        {
            (void)fprintf(filePtr, "%d\n", array[j]);
        }
    }
}

void writefileDouble(FILE* filePtr, double* array, unsigned int ArrayLength)
{
    if (filePtr)
    {
        for (unsigned int j = 0; j < ArrayLength; j++)
        {
            (void)fprintf(filePtr, "%f\n", array[j]);
        }
    }
}

static void calculateTotal(unsigned int count,
                           double* maxpos,
                           unsigned int ArrayLength,
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

    for (unsigned int i = 0; i < count; ++i)
    {
        xarr[i] = (double)i;
    }

    linreg(xarr, maxpos, count, &intercept, &slope, &stdev);

    /*
       intercept /= ArrayLength*ArrayLength;
       slope /= ArrayLength;
       stdev /= rate;
     */

    (void)fprintf(stderr,
                  "unweighted raw rate: %f s/d, %d samples σ=%.2gms\n",
                  -slope * SECS_DAY / ArrayLength,
                  count,
                  stdev * THOUSAND / rate);
    unsigned int maxIndex = 0;

    double deviation;

    for (unsigned int j = 0; j < 3; ++j)
    {
        for (unsigned int i = 0; i < count; ++i)
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
                      "after %.1fσ (%.2gms) removal: %.2f s/d, %d samples\n",
                      threshold,
                      stdev * THOUSAND / rate,
                      -slope * SECS_DAY / ArrayLength,
                      count);
        threshold /= 2;
    }

    free(xarr);
}

static void calculateTotalFromFile(unsigned int count,
                                   FILE* rawfile,
                                   unsigned int ArrayLength,
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
        size_t bufsize = BUF_SIZE;
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

void print_finals(CapConfig* cfg,
                  AppResources* res,
                  unsigned int ArrayLength,
                  unsigned int totalTickTock,
                  int toothshift)
{
    if (cfg->fpposition)
    {
        calculateTotalFromFile(totalTickTock,
                               cfg->fpposition,
                               ArrayLength,
                               cfg->SDthreshold,
                               cfg->rate);
    }
    if (cfg->fpmaxcor)
    {
        printTOD(cfg->fpmaxcor);
    }

    if (cfg->fptotal)
    {
        for (unsigned int t = 0; t < cfg->teeth; ++t)
        {
            struct myarr* tmp = res->teethArray[t];
            if (tmp != NULL)
            {
                struct myarr cumulativeTick = *tmp;
                toothshift = getshift(*res->teethArray[0], cumulativeTick);
                for (unsigned int j = 0; j < ArrayLength; ++j)
                {
                    (void)fprintf(cfg->fptotal,
                                  "%d %d %u %d\n",
                                  (int)j + toothshift,
                                  cumulativeTick.arr[j],
                                  t,
                                  toothshift);
                }
            }
            (void)fprintf(cfg->fptotal, "\n\n");
        }
    }
}

double getBeatError(const struct myarr* totaltick, double rate, int verbose)
{
    unsigned int ArrayLength = totaltick->ArrayLength;
    int* cross = malloc(ArrayLength / 2 * sizeof(int));
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
    unsigned int postick = getmaxpos(cross, ArrayLength / 2);
    free(cross);
    return shiftHalf(postick, ArrayLength / 2) * THOUSAND / rate;
}

void fillReference(FILE* fpDefPeak, struct myarr* reference, unsigned int teeth)
{
    if (fpDefPeak != NULL)
    {
        int arr[4];
        for (unsigned int t = 0; t < teeth; t++)
        {
            for (unsigned int j = 0; j < reference->ArrayLength; j++)
            {
                int int_count = getIntsFromStdin(4, arr);
                if (int_count < 0)
                {
                    break;
                }
                if (int_count < 3)
                {
                    (void)fprintf(
                        stderr,
                        "not enough values in -D <default peak file>\n 4 "
                        "columns required, %u samples and %u teeth\n",
                        reference->ArrayLength,
                        teeth);
                    exit(EXIT_FAILURE);
                }
                int value = arr[1];

                reference->arr[((int)j + (int)reference->ArrayLength) %
                               (int)reference->ArrayLength] = value;
            }
        }
        (void)fclose(fpDefPeak);
        fpDefPeak = NULL;
    }
    else
    {
        const int peakheight[3] = {100000, 80000, 60000};
        const int peakpos[3] = {0, 400, 800};

        for (int i = 0; i < 3; i++)
        {
            reference
                ->arr[reference->ArrayLength / 4 - (unsigned int)peakpos[i]] =
                peakheight[i];
            reference->arr[3 * reference->ArrayLength / 4 -
                           (unsigned int)peakpos[i]] = peakheight[i];
        }
    }
}
