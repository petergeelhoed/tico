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

/*Prints header on line or at the top */
void printheader(double b, unsigned int l, double beatError)
{
    if (l)
    {
        char line[14 + 1];
        memset(line, ' ', 14);
        snprintf(line, 5, "%4.2f", beatError);
        snprintf(line + 4, 8, "ms%+5.1f", b);
        sprintf(line + 11, "s/d");
        fprintf(stderr, "%s", line);
    }
    else
    {
        fprintf(stderr,
                "\033[s\033[2;0H\033[0K%8.2fms   %7.1fs/d\033[u",
                beatError,
                b);
    }
}

void printspaces(int maxpos,
                 int hexvalue,
                 int mod,
                 unsigned int columns,
                 double avg_pos,
                 int cvalue)
{
    while (maxpos < mod)
        maxpos += mod;
    while (avg_pos < (double)mod)
        avg_pos += (double)mod;
    columns = columns > 1024 ? 80 : columns;
    int width = (maxpos % mod) * (int)columns / mod;
    int widtha = (((int)avg_pos) % mod) * (int)columns / mod;

    char spaces[1024];
    memset(spaces, ' ', (size_t)width);
    spaces[width] = '\0';
    if (widtha < width)
        spaces[widtha] = '|';

    fprintf(stderr,
            "%s%s%X\033[0m",
            spaces,
            hexvalue < cvalue ? "\033[31m" : "\033[32m",
            hexvalue);

    memset(spaces, ' ', (size_t)width);
    if (widtha > width)
    {
        spaces[widtha - width - 1] = '|';
        spaces[widtha - width] = '\0';
        fprintf(stderr, "%s", spaces);
    }
    fprintf(stderr, "\n");
}

void linregd(const float* xarr,
             const float* yarr,
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

void linreg(const int* xarr,
            const int* yarr,
            unsigned int NN,
            double* a,
            double* b,
            double* s)
{
    // remove the average to pospone overflows
    double mx = (xarr[NN - 1] + xarr[0]) / 2;
    double my = (yarr[NN - 1] + yarr[0]) / 2;
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (unsigned int i = 0; i < NN; ++i)
    {
        y += yarr[i] - my;
        xx += (xarr[i] - mx) * (xarr[i] - mx);
        x += xarr[i] - mx;
        xy += (xarr[i] - mx) * (yarr[i] - my);
        yy += (yarr[i] - my) * (yarr[i] - my);
    }

    *a = (y * xx - x * xy) / (NN * xx - x * x);
    *b = (NN * xy - x * y) / (NN * xx - x * x);
    *s = sqrt((yy - 2 * (*a) * y - 2 * (*b) * xy + 2 * (*a) * (*b) * x +
               (*a) * (*a) * NN + (*b) * (*b) * xx) /
              NN);
    *a -= *b * mx - my;
}

void fit10secs(double* a,
               double* b,
               double* s,
               const unsigned int i,
               const int* maxvals,
               const int* maxes,
               const int cvalue,
               const unsigned int npeaks)
{
    unsigned int m = 0;
    unsigned int fitwindow = i > npeaks ? npeaks : i;

    if (i >= fitwindow)
    {
        int xarr[fitwindow];
        int yarr[fitwindow];
        for (unsigned int k = 0; k < fitwindow; k++)
        {
            if (maxvals[i - k] > cvalue)
            {
                yarr[m] = maxes[i - k];
                xarr[m] = (int)k;
                m++;
            }
        }
        if (m > 1)
        {
            linreg(xarr, yarr, m, a, b, s);
        }
    }
}

void writefile(FILE* fp, int* array, unsigned int NN)
{
    if (fp)
    {
        for (unsigned int j = 0; j < NN; j++)
            fprintf(fp, "%d\n", array[j]);
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
    int* all = calloc(n, sizeof(int));
    unsigned int i = 0;
    if (all)
    {

        unsigned int bufsize = 256;
        char* buf = malloc(bufsize * sizeof(char));
        while (getline(&buf, &bufsize, rawfile) > 0 && i < n)
        {
            if (buf[0] != '#')
            {
                all[i] = atoi(buf);
                i++;
            }
        }
        free(buf);
        calculateTotal(n, all, NN, threshold);
    }
}

void calculateTotal(unsigned int n,
                    int* maxpos,
                    unsigned int NN,
                    double threshold)
{
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int xarr[n];

    for (unsigned int i = 0; i < n; ++i)
    {
        xarr[i] = (int)i;
    }

    linreg(xarr, maxpos, n, &a, &b, &s);

    /*
       a /= NN*NN;
       b /= NN;
       s /= rate;
     */

    fprintf(stderr, "raw rate: %f s/d, %d samples\n", -b * 86400 / NN, n);
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

    fprintf(stderr,
            "after %.1fσ removal: %.2f s/d, %d samples\n",
            threshold,
            -b * 86400 / NN,
            m);
}

double
getBeatError(int* totaltick, unsigned int NN, unsigned int rate, int verbose)
{

    int cross[NN / 2];
    crosscorint(NN / 2, totaltick, totaltick + NN / 2, cross);
    if (verbose)
    {
        syncwrite(cross, NN / 2, "beaterror");
        syncwrite(totaltick, NN / 2, "t1");
        syncwrite(totaltick + NN / 2, NN / 2, "t2");
    }
    unsigned int postick = getmaxpos(cross, NN / 2);
    return shiftHalf(postick, NN / 2) * 1000 / (double)rate;
}

int checkUIntArg(int name, unsigned int* value, char* optarg)
{
    *value = 0;
    *value = (unsigned int)atoi(optarg);
    if (*value == 0)
    {
        printf("invalid integer argument for -%c: '%s'\n", (char)name, optarg);
        return -1;
    }
    return 0;
}

int checkFileArg(int name, FILE** fp, char* optarg, char* mode)
{
    if (*optarg == '-')
    {
        fprintf(
            stderr, "expecting -%c <file>\n got -w %s\n", (char)name, optarg);
        return -1;
    }

    *fp = fopen(optarg, mode);
    if (*fp == 0)
    {
        fprintf(stderr,
                "cannot open file -%c '%s' for mode %s\n",
                (char)name,
                optarg,
                mode);
        return -4;
    }
    return 0;
}

void fillReference(FILE* fpDefPeak, int* reference, unsigned int NN)
{
    // read default peak
    if (fpDefPeak != 0)
    {
        for (unsigned int j = 0; j < NN; j++)
        {

            if (fscanf(fpDefPeak, "%d", reference + j) != 1)
            {
                fprintf(stderr,
                        "not enough values in -D <default peak file>\n");

                exit(-5);
            }
        }
        fclose(fpDefPeak);
    }
    else
    {
        reference[NN / 4] = 100000;
        reference[3 * NN / 4] = 100000;
    }
}

void checkAndFlip(int* totaltick,
                  int* reference,
                  unsigned int NN,
                  unsigned int verbose)
{
    int* cross = malloc(NN * sizeof(int));
    crosscorint(NN, totaltick, reference, cross);
    unsigned int flipmaxp = getmaxpos(cross, NN);
    if (verbose)
    {
        FILE* fp = fopen("flip", "w");
        if (fp)
        {
            for (unsigned int j = 0; j < NN; j++)
            {
                fprintf(fp,
                        "%d %d %d %d\n",
                        j,
                        totaltick[j],
                        reference[j],
                        cross[j]);
            }
            fclose(fp);
        }
    }

    if (flipmaxp > NN / 4 && flipmaxp < NN * 3 / 4)
    {
        fprintf(stderr, "FLIPPING peaks pos %d\n", flipmaxp);

        int tmp = 0;
        for (unsigned int j = 0; j < NN / 2; j++)
        {
            tmp = reference[j + NN / 2];
            reference[j + NN / 2] = reference[j];
            reference[j] = tmp;
        }
    }

    free(cross);
}

int shiftHalf(unsigned int value, unsigned int NN)
{
    return ((int)value + (int)NN / 2) % (int)(NN) - (int)(NN / 2);
}
