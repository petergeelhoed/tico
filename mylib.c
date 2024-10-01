#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mylib.h"

/*Prints header on line or at the top */
void printheader(double b, int NN, int l, float beatError)
{
    if (l)
    {
        char line[14 + 1];
        snprintf(line, 5, "%4.2f", beatError);
        snprintf(line + 4, 8, "ms%+5.1f", b * 86400 / NN);
        sprintf(line + 11, "s/d");
        fprintf(stderr, "%s", line);
    }
    else
    {
        fprintf(stderr,
                "\033[s\033[2;0H\033[K%+8.2fms   %+7.1fs/d\033[u",
                beatError,
                b * 86400 / NN);
    }
}

void printspaces(int maxpos,
                 int hexvalue,
                 char* spaces,
                 int mod,
                 int columns,
                 double avg_pos,
                 int cvalue)
{
    while (maxpos < mod)
        maxpos += mod;
    while (avg_pos < mod)
        avg_pos += (double)mod;
    int width = (maxpos % mod) * columns / mod;
    int widtha = (((int)avg_pos) % mod) * columns / mod;

    memset(spaces, ' ', width);
    spaces[width] = '\0';
    if (widtha < width)
        spaces[widtha] = '|';

    fprintf(stderr,
            "%s%s%X\033[0m",
            spaces,
            hexvalue < cvalue ? "\033[31m" : "\033[32m",
            hexvalue);

    memset(spaces, ' ', width);
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
             int NN,
             double* a,
             double* b,
             double* s)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (int i = 0; i < NN; ++i)
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

void linreg(
    const int* xarr, const int* yarr, int NN, double* a, double* b, double* s)
{
    double x = 0;
    double y = 0;
    double xx = 0;
    double xy = 0;
    double yy = 0;
    for (int i = 0; i < NN; ++i)
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

void fit10secs(double* a,
               double* b,
               double* s,
               int i,
               int* maxvals,
               int* maxes,
               int cvalue,
               int npeaks)
{
    int m = 0;
    int fitwindow = i > npeaks ? npeaks : i;

    if (i >= fitwindow)
    {
        int xarr[fitwindow];
        int yarr[fitwindow];
        for (int k = 0; k < fitwindow; k++)
        {
            if (maxvals[i - k] > cvalue)
            {
                yarr[m] = maxes[i - k];
                xarr[m] = k;
                m++;
            }
        }
        if (m > 1)
        {
            linreg(xarr, yarr, m, a, b, s);
        }
    }
}

void writefile(FILE* fp, int* array, int NN)
{
    if (fp)
    {
        for (int j = 0; j < NN; j++)
            fprintf(fp, "%d\n", array[j]);
    }
}

int getmaxpos(int* array, int NN)
{
    int maxtick = -INT_MAX;
    int postick = 0;
    for (int j = 0; j < NN; j++)
    {
        if (array[j] > maxtick)
        {
            maxtick = array[j];
            postick = j;
        }
    }
    return postick;
}

void calculateTotal(int n, int* maxpos, int NN, double threshold)
{
    double b = 0.0;
    double a = 0.0;
    double s = 0.0;
    int xarr[n];

    for (int i = 0; i < n; ++i)
    {
        xarr[i] = i;
    }

    linreg(xarr, maxpos, n, &a, &b, &s);

    /*
       a /= NN*NN;
       b /= NN;
       s /= rate;
     */

    fprintf(stderr, "raw rate: %f s/d\n", -b * 86400 / NN);
    int m = 0;

    double e;

    for (int i = 0; i < n; ++i)
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

    fprintf(
        stderr, "after %.1fσ removal: %.2f s/d\n", threshold, -b * 86400 / NN);
}

int getBeatError(int* totaltick, int NN, int verbose)
{

    int cross[NN / 2];
    crosscorint(NN / 2, totaltick, totaltick + NN / 2, cross);
    if (verbose)
    {
        syncwrite(cross, NN / 2, "beaterror");
        syncwrite(totaltick, NN / 2, "t1");
        syncwrite(totaltick + NN / 2, NN / 2, "t2");
    }
    int postick = getmaxpos(cross, NN / 2);
    return (postick + NN / 4) % (NN / 2) - NN / 4;
}

void writearray(int* arr, int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    for (int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %d\n", j, arr[j]);
    }
    fclose(fp);
}

void syncappend(int* input, int NN, FILE* file)
{
    struct mystruct
    {
        int* array;
        FILE* file;
        int NN;
    };

    struct mystruct* info = malloc(sizeof *info);

    int* copyarr = malloc(NN * sizeof(int));
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    info->file = file;
    info->NN = NN;

    pthread_t tid;
    pthread_create(&tid, NULL, threadAppend, info);
    pthread_detach(tid);
}

void* threadAppend(void* inStruct)
{
    struct mystruct
    {
        int* array;
        FILE* file;
        int NN;
    } mine = *(struct mystruct*)inStruct;

    int* arrptr = mine.array;
    FILE* file = mine.file;
    int* copyarr = malloc(mine.NN * sizeof(int));
    memcpy(copyarr, arrptr, mine.NN * sizeof(int));
    mine.array = copyarr;
    free(arrptr);
    free(inStruct);

    for (int j = 0; j < mine.NN; j++)
    {
        fprintf(file, "%d\n", mine.array[j]);
    }
    fflush(file);

    free(copyarr);
    pthread_exit(NULL);
}

void syncwrite(int* input, int NN, char* file)
{
    struct mystruct
    {
        int* array;
        char file[20];
        int NN;
    };

    struct mystruct* info = malloc(sizeof *info);

    int* copyarr = malloc(NN * sizeof(int));
    memcpy(copyarr, input, NN * sizeof(int));
    info->array = copyarr;
    strcpy(info->file, file);
    info->NN = NN;

    pthread_t tid;
    pthread_create(&tid, NULL, threadWrite, info);
    pthread_detach(tid);
}

void* threadWrite(void* inStruct)
{
    struct mystruct
    {
        int* array;
        char file[20];
        int NN;
    } mine = *(struct mystruct*)inStruct;

    int* arrptr = mine.array;
    int* copyarr = malloc(mine.NN * sizeof(int));
    memcpy(copyarr, arrptr, mine.NN * sizeof(int));
    mine.array = copyarr;
    free(arrptr);
    free(inStruct);

    writearray(mine.array, mine.NN, mine.file);

    pthread_exit(NULL);
}

void rescale(int* total, int NN)
{
    if (total[NN / 2] > 100000000 || total[0] > 100)
    {

        long int avg = 0;

        for (int j = 0; j < NN; j++)
        {
            avg += total[j];
        }
        avg /= NN;
        int avi = (int)avg;
        if (avi > 100)
        {
            for (int j = 0; j < NN; j++)
            {
                total[j] -= avi;
            }
        }
        else
        {
            for (int j = 0; j < NN; j++)
            {
                total[j] /= 2;
            }
        }
    }
}
