#include "capture_helpers.h"

#include "analysis.h"
#include "mydefs.h"
#include "myfft.h"
#include "mymath.h"
#include "mysync.h"
#include "parseargs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
                 double avgPos,
                 unsigned int correlationThreshold)
{
    while (maxpos < (int)mod)
    {
        maxpos += (int)mod;
    }
    while (avgPos < (double)mod)
    {
        avgPos += (double)mod;
    }
    const unsigned int default_columns = 80;
    columns = columns > MAX_COLUMNS ? default_columns : columns;
    size_t width = (size_t)modSigned(maxpos, mod) * columns / mod;
    size_t widtha =
        (size_t)modSigned((int)lround(avgPos), mod) * columns / mod;

    char spaces[MAX_COLUMNS];
    memset(spaces, ' ', width);
    spaces[width] = '\0';
    if (widtha < width)
    {
        spaces[widtha] = '|';
    }

    (void)fprintf(
        stderr,
        "%s%s%X\033[0m",
        spaces,
        (unsigned int)hexvalue < correlationThreshold ? "\033[31m" : "\033[32m",
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

void printFinals(CapConfig* cfg,
                 AppResources* res,
                 unsigned int ArrayLength,
                 unsigned int totalTickTock)
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
                int toothshift = getshift(*res->teethArray[0], cumulativeTick);
                for (unsigned int j = 0; j < ArrayLength; ++j)
                {
                    (void)fprintf(cfg->fptotal,
                                  "%d %d %u %d\n",
                                  (int)j,
                                  cumulativeTick.arr[j],
                                  t,
                                  toothshift);
                }
            }
            else
            {
                break;
            }
            (void)fprintf(cfg->fptotal, "\n\n");
        }
    }
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
                        "not enough values in -D <default peak file>\n 4 columns required, %u samples and %u teeth\n",
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
            reference->arr[reference->ArrayLength / 4 - (unsigned int)peakpos[i]] =
                peakheight[i];
            reference->arr[3 * reference->ArrayLength / 4 - (unsigned int)peakpos[i]] =
                peakheight[i];
        }
    }
}

void shiftBufferData(unsigned int* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals)
{
    memmove(subpos->arrd,
            subpos->arrd + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(double));
    memmove(maxpos->arr,
            maxpos->arr + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(int));
    memmove(maxvals->arrd,
            maxvals->arrd + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(double));
    *ticktock -= ARRAY_BUFFER_SIZE;
}
