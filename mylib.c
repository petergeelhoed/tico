#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "myfft.h"
#include "mylib.h"
#include "mysync.h"
#include "parseargs.h"

#define MAX_COLUMNS 1024
#define MAX_COLUMNS 1024
#define BEAT_WIDTH 5
#define RATE_WIDTH 8
#define DECIMAL 10

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
            line + BEAT_WIDTH - 1, RATE_WIDTH, "ms%+5.1f", fittedRate);
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
                 double corvalue,
                 unsigned int mod,
                 unsigned int columns,
                 double avg_pos,
                 unsigned int correlationThreshold)
{
    while (maxpos < (int)mod)
        maxpos += mod;
    while (avg_pos < (double)mod)
        avg_pos += (double)mod;
    columns = columns > MAX_COLUMNS ? 80 : columns;
    size_t width = (size_t)modSigned(maxpos, mod) * columns / mod;
    size_t widtha = (size_t)modSigned(lround(avg_pos), mod) * columns / mod;

    char spaces[MAX_COLUMNS];
    memset(spaces, ' ', width);
    spaces[width] = '\0';
    if (widtha < width)
        spaces[widtha] = '|';

    (void)fprintf(stderr,
                  "%s%s%X\033[0m",
                  spaces,
                  (unsigned int)corvalue < correlationThreshold ? "\033[31m"
                                                                : "\033[32m",
                  (int)corvalue);

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
            unsigned int NN,
            double* a,
            double* b,
            double* s)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (unsigned int i = 0; i < NN; ++i)
    {
        y += yarr[i];
        xx += xarr[i] * xarr[i];
        x += xarr[i];
        xy += xarr[i] * yarr[i];
        yy += yarr[i] * yarr[i];
    }

    *a = (y * xx - x * xy) / (NN * xx - x * x);
    *b = (NN * xy - x * y) / (NN * xx - x * x);
    *s = sqrt((yy - 2 * (*a) * y - 2 * (*b) * xy + 2 * (*a) * (*b) * x +
               (*a) * (*a) * NN + (*b) * (*b) * xx) /
              NN);
}

void writefile(FILE* fp, int* array, unsigned int NN)
{
    if (fp)
    {
        for (unsigned int j = 0; j < NN; j++)
            (void)fprintf(fp, "%d\n", array[j]);
    }
}

void writefileDouble(FILE* fp, double* array, unsigned int NN)
{
    if (fp)
    {
        for (unsigned int j = 0; j < NN; j++)
            (void)fprintf(fp, "%f\n", array[j]);
    }
}

unsigned int getmaxpos(int* array, unsigned int NN)
{
    int maxtick = -INT_MAX;
    unsigned int postick = 0;
    for (unsigned int j = 0; j < NN; j++)
    {
        if (array[j] > maxtick)
        {
            maxtick = array[j];
            postick = j;
        }
    }
    return postick;
}

void calculateTotalFromFile(unsigned int n,
                            FILE* rawfile,
                            unsigned int NN,
                            double threshold)
{
    fseek(rawfile, 0, SEEK_SET);
    double* all = calloc(n, sizeof(double));
    unsigned int i = 0;
    if (all)
    {
        size_t bufsize = 256;
        char* buf = malloc(bufsize * sizeof(char));
        while (getline(&buf, &bufsize, rawfile) > 0 && i < n)
        {
            if (buf[0] != '#')
            {
                all[i] = atof(buf);
                i++;
            }
        }
        free(buf);
        calculateTotal(n, all, NN, threshold);
        free(all);
    }
}

void calculateTotal(unsigned int n,
                    double* maxpos,
                    unsigned int NN,
                    double threshold)
{
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    double xarr[n];

    for (unsigned int i = 0; i < n; ++i)
    {
        xarr[i] = (double)i;
    }

    linreg(xarr, maxpos, n, &a, &b, &s);

    /*
       a /= NN*NN;
       b /= NN;
       s /= rate;
     */

    (void)fprintf(stderr, "raw rate: %f s/d, %d samples\n", -b * 86400 / NN, n);
    unsigned int m = 0;

    double e;

    for (unsigned int i = 0; i < n; ++i)
    {
        e = fabs(((double)maxpos[i] - (a + xarr[i] * b)) / s);
        if (e < threshold)
        {
            maxpos[m] = maxpos[i];
            xarr[m] = xarr[i];
            m++;
        }
    }
    linreg(xarr, maxpos, m, &a, &b, &s);

    (void)fprintf(stderr,
                  "after %.1fσ removal: %.2f s/d, %d samples\n",
                  threshold,
                  -b * 86400 / NN,
                  m);
}

double getBeatError(const struct myarr* totaltick, double rate, int verbose)
{
    unsigned int NN = totaltick->NN;
    int cross[NN / 2];
    crosscorint(NN / 2, totaltick->arr, totaltick->arr + NN / 2, cross);
    if (verbose)
    {
        syncwrite(cross, NN / 2, "beaterror");
        syncwrite(totaltick->arr, NN / 2, "t1");
        syncwrite(totaltick->arr + NN / 2, NN / 2, "t2");
    }
    unsigned int postick = getmaxpos(cross, NN / 2);
    return shiftHalf(postick, NN / 2) * 1000 / rate;
}

int checkUIntArg(int name, unsigned int* value, char* optarg)
{
    *value = (unsigned int)atoi(optarg);
    if (*value == 0)
    {
        printf("invalid integer argument for -%c: '%s'\n", (char)name, optarg);
        return -1;
    }
    return 0;
}

int checkFileArg(int name, FILE** fp, const char* optarg, const char* mode)
{
    if (*optarg == '-')
    {
        (void)fprintf(
            stderr, "expecting -%c <file>\n got -w %s\n", (char)name, optarg);
        return -1;
    }

    *fp = fopen(optarg, mode);
    if (*fp == NULL)
    {
        (void)fprintf(stderr,
                      "cannot open file -%c '%s' for mode %s\n",
                      (char)name,
                      optarg,
                      mode);
        return -4;
    }
    return 0;
}

void fillReference(FILE* fpDefPeak, struct myarr* reference, unsigned int teeth)
{
    if (fpDefPeak != NULL)
    {
        int dummy0, dummy1;
        int shift = 0;
        for (unsigned int t = 0; t < teeth; t++)
        {
            for (unsigned int j = 0; j < reference->NN; j++)
            {
                int value = 0;
                if (fscanf(fpDefPeak,
                           "%d %d %d %d",
                           &dummy0,
                           &value,
                           &dummy1,
                           &shift) != 4)
                {
                    (void)fprintf(
                        stderr,
                        "not enough values in -D <default peak file>\n 4 "
                        "columns required, %u samples and %u teeth\n",
                        reference->NN,
                        teeth);
                    exit(-5);
                }
                reference
                    ->arr[((int)j + (int)reference->NN) % (int)reference->NN] =
                    value;
            }
        }
        fclose(fpDefPeak);
    }
    else
    {
        reference->arr[reference->NN / 4] = 100000;
        reference->arr[reference->NN / 4 - 400] = 80000;
        reference->arr[reference->NN / 4 - 800] = 60000;
        reference->arr[3 * reference->NN / 4] = 100000;
        reference->arr[3 * reference->NN / 4 - 400] = 80000;
        reference->arr[3 * reference->NN / 4 - 800] = 60000;
    }
}

int shiftHalf(unsigned int value, unsigned int NN)
{
    return ((int)value + (int)NN / 2) % (int)(NN) - (int)(NN / 2);
}

// mods an int with a signed int, but makes sure the result is positive
int modSigned(int value, unsigned int NN)
{
    return (value % (int)NN + (int)NN) % (int)NN;
}

int getInt(char* ptr)
{
    char* endptr;

    errno = 0;
    long val = INT_MIN;
    val = strtol(ptr, &endptr, DECIMAL);

    // If ptr == endptr, no more numbers were found on this line
    if (ptr == endptr || errno == ERANGE || (val < INT_MIN || val > INT_MAX))
    {
        (void)fprintf(stderr, "Invalid long or out of range\n");
        return INT_MIN;
    }

    return (int)val;
}

double getDouble(char* ptr)
{
    char* endptr;

    // Reset errno before the call to accurately catch new range errors
    errno = 0;
    double val = strtod(ptr, &endptr);

    if (ptr == endptr || errno == ERANGE)
    {
        (void)fprintf(stderr, "Invalid double or out of range\n");
        val = NAN;
    }

    return val;
}

#include <ctype.h>

/**
 * Reads one line from stdin and parses up to max_count doubles into array.
 * Returns the number of doubles stored (>=0), or -1 on I/O/memory error.
 */
int getDoublesFromStdin(size_t max_count, double* arr)
{
    if (!arr || max_count == 0)
    {
        return 0; // nothing to do
    }

    char* line = NULL;
    size_t len = 0;
    ssize_t nread = getline(&line, &len, stdin);
    if (nread == -1)
    {
        free(line);
        return -1; // EOF or read error
    }

    const char* ptr = line;
    char* endptr = NULL;
    size_t parsed = 0;

    // Parse up to max_count doubles from THIS line only
    while (*ptr != '\0' && parsed < max_count)
    {
        errno = 0;
        double val = strtod(ptr, &endptr);
        if (endptr == ptr)
        {
            // No number at current position; skip one char
            ptr++;
            continue;
        }
        arr[parsed++] = val;
        ptr = endptr;
    }

    free(line);
    return (int)parsed; // returns 0..max_count
}
