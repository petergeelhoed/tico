#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
#include "mysync.h"
#include "parseargs.h"

/* Prints header on line or at the top */
void printheader(double fittedRate,
                 unsigned int everyline,
                 double beatError,
                 double seconds)
{
    if (everyline)
    {
        char line[EVERY_WIDTH + 1];
        memset(line, ' ', EVERY_WIDTH);
        (void)snprintf(line, BEAT_WIDTH, "%4.2f", beatError);
        (void)snprintf(
            line + BEAT_WIDTH - 1, RATE_WIDTH + 1, "ms%+5.1f", fittedRate);
        (void)sprintf(line + BEAT_WIDTH + RATE_WIDTH - 2, "s/d");
        (void)fprintf(stderr, "%s", line);
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

void linreg(const double* xarr,
            const double* yarr,
            unsigned int ArrayLength,
            double* par_a,
            double* par_b,
            double* par_s)
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_xx = 0;
    double sum_xy = 0;
    double sum_yy = 0;
    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        sum_y += yarr[i];
        sum_xx += xarr[i] * xarr[i];
        sum_x += xarr[i];
        sum_xy += xarr[i] * yarr[i];
        sum_yy += yarr[i] * yarr[i];
    }

    *par_a = (sum_y * sum_xx - sum_x * sum_xy) /
             (ArrayLength * sum_xx - sum_x * sum_x);
    *par_b = (ArrayLength * sum_xy - sum_x * sum_y) /
             (ArrayLength * sum_xx - sum_x * sum_x);
    *par_s = sqrt((sum_yy - 2 * (*par_a) * sum_y - 2 * (*par_b) * sum_xy +
                   2 * (*par_a) * (*par_b) * sum_x +
                   (*par_a) * (*par_a) * ArrayLength +
                   (*par_b) * (*par_b) * sum_xx) /
                  ArrayLength);
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

unsigned int getmaxpos(const int* array, unsigned int ArrayLength)
{
    int maxtick = -INT_MAX;
    unsigned int postick = 0;
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        if (array[j] > maxtick)
        {
            maxtick = array[j];
            postick = j;
        }
    }
    return postick;
}

void calculateTotalFromFile(unsigned int count,
                            FILE* rawfile,
                            unsigned int ArrayLength,
                            double threshold)
{
    errno = 0;
    if (fseek(rawfile, 0, SEEK_SET) == -1)
    {
        (void)fprintf(stderr, "fseek fauled with %d\n", errno);
        return;
    }
    double* all = calloc(count, sizeof(double));
    unsigned int index = 0;
    if (all)
    {
        size_t bufsize = BUF_SIZE;
        char* buf = malloc(bufsize * sizeof(char));
        while (getline(&buf, &bufsize, rawfile) > 0 && index < count)
        {
            if (buf[0] != '#')
            {
                all[index++] = getDouble(buf);
            }
        }
        free(buf);
        calculateTotal(count, all, ArrayLength, threshold);
        free(all);
    }
}

void calculateTotal(unsigned int count,
                    double* maxpos,
                    unsigned int ArrayLength,
                    double threshold)
{
    double par_b = 0.0;
    double par_a = 0.0;
    double par_s = 0.0;
    double* xarr = calloc(count, sizeof(double));
    if (xarr == NULL)
    {
        exit(EXIT_FAILURE);
    }

    for (unsigned int i = 0; i < count; ++i)
    {
        xarr[i] = (double)i;
    }

    linreg(xarr, maxpos, count, &par_a, &par_b, &par_s);

    /*
       par_a /= ArrayLength*ArrayLength;
       par_b /= ArrayLength;
       par_s /= rate;
     */

    (void)fprintf(stderr,
                  "raw rate: %f s/d, %d samples\n",
                  -par_b * SECS_DAY / ArrayLength,
                  count);
    unsigned int maxIndex = 0;

    double deviation;

    for (unsigned int i = 0; i < count; ++i)
    {
        deviation = fabs((maxpos[i] - (par_a + xarr[i] * par_b)) / par_s);
        if (deviation < threshold)
        {
            maxpos[maxIndex] = maxpos[i];
            xarr[maxIndex] = xarr[i];
            maxIndex++;
        }
    }
    linreg(xarr, maxpos, maxIndex, &par_a, &par_b, &par_s);

    (void)fprintf(stderr,
                  "after %.1fσ removal: %.2f s/d, %d samples\n",
                  threshold,
                  -par_b * SECS_DAY / ArrayLength,
                  maxIndex);
    free(xarr);
}

double getBeatError(const struct myarr* totaltick, double rate, int verbose)
{
    unsigned int ArrayLength = totaltick->ArrayLength;
    int cross[ArrayLength / 2];
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

int shiftHalf(unsigned int value, unsigned int ArrayLength)
{
    return ((int)value + (int)ArrayLength / 2) % (int)(ArrayLength) -
           (int)(ArrayLength / 2);
}

// mods an int with a signed int, but makes sure the result is positive
int modSigned(int value, unsigned int ArrayLength)
{
    return (value % (int)ArrayLength + (int)ArrayLength) % (int)ArrayLength;
}
