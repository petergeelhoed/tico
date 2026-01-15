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
        maxpos += mod;
    }
    while (avg_pos < (double)mod)
    {
        avg_pos += (double)mod;
    }
    const unsigned int default_columns = 80;
    columns = columns > MAX_COLUMNS ? default_columns : columns;
    size_t width = (size_t)modSigned(maxpos, mod) * columns / mod;
    size_t widtha = (size_t)modSigned(lround(avg_pos), mod) * columns / mod;

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

int checkUIntArg(int name, unsigned int* value, char* optarg)
{
    *value = (unsigned int)getInt(optarg);
    if (*value == 0)
    {
        printf("invalid integer argument for -%c: '%s'\n", (char)name, optarg);
        return -1;
    }
    return 0;
}

int checkFileArg(int name, FILE** filePtr, const char* optarg, const char* mode)
{
    if (*optarg == '-')
    {
        (void)fprintf(
            stderr, "expecting -%c <file>\n got -w %s\n", (char)name, optarg);
        return -1;
    }

    *filePtr = fopen(optarg, mode);
    if (*filePtr == NULL)
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
            reference->arr[reference->ArrayLength / 4 - peakpos[i]] =
                peakheight[i];
            reference->arr[3 * reference->ArrayLength / 4 - peakpos[i]] =
                peakheight[i];
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

int getIntsFromStdin(size_t max_count, int* arr)
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
        int val = strtol(ptr, &endptr, DECIMAL);
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
