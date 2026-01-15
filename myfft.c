#include <assert.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mysync.h"

fftw_complex* makeFilter(unsigned int evalue, unsigned int ArrayLength)
{
    fftw_complex* filter = fftw_alloc_complex(ArrayLength);
    fftw_complex* filterFFT = fftw_alloc_complex(ArrayLength);
    fftw_plan makefilter = fftw_plan_dft_1d(
        (int)ArrayLength, filter, filterFFT, FFTW_FORWARD, FFTW_ESTIMATE);

    if (filter == NULL || filterFFT == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in_data makeFilter\n");
        return NULL;
    }

    if (evalue != 0)
    {
        for (unsigned int j = 0; j < evalue * GAUSSPOINTS; j++)
        {
            filter[j][0] =
                GAUSSIAN_CONST / evalue *
                exp(-((double)(j * j)) / (double)(evalue * evalue) / 2);
            filter[j][1] = 0.0;
        }
        for (unsigned int j = evalue * GAUSSPOINTS;
             j < ArrayLength - evalue * GAUSSPOINTS;
             j++)
        {
            filter[j][0] = 0.0;
            filter[j][1] = 0.0;
        }
        for (unsigned int j = ArrayLength - (evalue * GAUSSPOINTS);
             j < ArrayLength;
             j++)
        {
            filter[j][0] =
                GAUSSIAN_CONST / evalue *
                exp(-((double)(ArrayLength - j) * (ArrayLength - j)) /
                    (double)(evalue * evalue) / 2);
            filter[j][1] = 0.0;
        }
    }
    else
    {
        filter[0][0] = 1;
        filter[0][1] = 0;
        for (unsigned int j = 1; j < ArrayLength; j++)
        {
            filter[j][0] = 0;
            filter[j][1] = 0;
        }
    }

    fftw_execute(makefilter);
    fftw_destroy_plan(makefilter);
    fftw_free(filter);
    fftw_cleanup();
    return filterFFT;
}

void remove50hz(unsigned int ArrayLength, int* array, unsigned int rate)
{
    const unsigned int freq = 50;
    unsigned int ofj = ArrayLength * freq / rate;
    fftw_complex* in_data = fftw_alloc_complex(ArrayLength);
    fftw_complex* out_data = fftw_alloc_complex(ArrayLength);
    fftw_plan forward = fftw_plan_dft_1d(
        (int)ArrayLength, in_data, out_data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d(
        (int)ArrayLength, out_data, in_data, FFTW_BACKWARD, FFTW_ESTIMATE);

    if (in_data == NULL || out_data == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in_data remove50hz\n");
        return;
    }

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        in_data[j][0] = (double)array[j];
        in_data[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        out_data[j][0] /= ArrayLength;
        out_data[j][1] /= ArrayLength;
    }

    out_data[ofj + 1][1] = 0.0;
    out_data[ofj - 1][0] = 0.0;
    out_data[ofj + 1][1] = 0.0;
    out_data[ofj - 1][0] = 0.0;

    for (unsigned int j = ofj; j < ArrayLength; j += ofj)
    {
        out_data[j][0] = 0.0;
        out_data[j][1] = 0.0;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        array[j] = (int)(in_data[j][0]);
    }
    fftw_free(in_data);
    fftw_free(out_data);
    fftw_cleanup();
}

fftw_complex* convolute(const struct myarr array, fftw_complex* filterFFT)
{
    unsigned int ArrayLength = array.ArrayLength;
    fftw_complex* in_data = fftw_alloc_complex(ArrayLength);
    fftw_complex* out_data = fftw_alloc_complex(ArrayLength);
    fftw_plan forward = fftw_plan_dft_1d(
        (int)ArrayLength, in_data, out_data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d(
        (int)ArrayLength, out_data, in_data, FFTW_BACKWARD, FFTW_ESTIMATE);

    if (in_data == NULL || out_data == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in_data convolute\n");
        return NULL;
    }

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        in_data[j][0] = (double)array.arr[j];
        in_data[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        out_data[j][0] = (out_data[j][0] * filterFFT[j][0] -
                          out_data[j][1] * filterFFT[j][1]) /
                         ArrayLength;
        out_data[j][1] = (out_data[j][0] * filterFFT[j][1] +
                          out_data[j][1] * filterFFT[j][0]) /
                         ArrayLength;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);
    fftw_free(out_data);
    fftw_cleanup();
    return in_data;
}

fftw_complex*
crosscor(unsigned int ArrayLength, fftw_complex* array, fftw_complex* ref)
{
    normalise(ArrayLength, array);
    normalise(ArrayLength, ref);

    fftw_complex* tmparr = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmpref = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmp = fftw_alloc_complex(ArrayLength);
    fftw_complex* corr = fftw_alloc_complex(ArrayLength);

    if (tmparr == NULL || tmpref == NULL || tmp == NULL || corr == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in_data crosscor\n");
        return NULL;
    }

    fftw_plan arrFour = fftw_plan_dft_1d(
        (int)ArrayLength, array, tmparr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan refFour = fftw_plan_dft_1d(
        (int)ArrayLength, ref, tmpref, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan correverse = fftw_plan_dft_1d(
        (int)ArrayLength, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(arrFour);
    fftw_execute(refFour);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        tmp[j][0] = (tmparr[j][0] * tmpref[j][0] + tmparr[j][1] * tmpref[j][1]);
        tmp[j][1] =
            (-tmparr[j][0] * tmpref[j][1] + tmparr[j][1] * tmpref[j][0]);
    }
    fftw_execute(correverse);

    double scale = 1. / ArrayLength / ArrayLength;
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        corr[j][0] *= scale;
    }

    fftw_destroy_plan(correverse);
    fftw_destroy_plan(arrFour);
    fftw_destroy_plan(refFour);

    fftw_free(tmparr);
    fftw_free(tmpref);
    fftw_free(tmp);
    fftw_cleanup();
    return corr;
}

int getshift(const struct myarr xarr, const struct myarr yarr)
{
    unsigned int ArrayLength = xarr.ArrayLength;
    fftw_complex* F_x = fftw_alloc_complex(ArrayLength);
    fftw_complex* F_y = fftw_alloc_complex(ArrayLength);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        F_x[j][0] = (double)xarr.arr[j];
        F_x[j][1] = 0.0;
        F_y[j][0] = (double)yarr.arr[j];
        F_y[j][1] = 0.0;
    }
    fftw_complex* corr = crosscor(ArrayLength, F_x, F_y);

    unsigned int poscor = getmaxfftw(corr, ArrayLength);
    fftw_free(F_x);
    fftw_free(corr);
    fftw_free(F_y);
    fftw_cleanup();
    assert(ArrayLength > 0);
    return (((int)poscor + (int)ArrayLength / 2) % (int)(ArrayLength)) -
           (int)(ArrayLength / 2);
}

unsigned int fftfit(const struct myarr input,
                    int* total,
                    const int* base,
                    double* corvalue,
                    fftw_complex* filterFFT,
                    int verb,
                    double* subpos)
{
    unsigned int ArrayLength = input.ArrayLength;
    fftw_complex* Fbase = fftw_alloc_complex(ArrayLength);
    fftw_complex* filteredinput = convolute(input, filterFFT);

    if (Fbase == NULL || filteredinput == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed in_data fftfit\n");
        return 0;
    }

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        Fbase[j][0] = (double)base[j];
        Fbase[j][1] = 0.0;
    }

    fftw_complex* corr = crosscor(ArrayLength, filteredinput, Fbase);

    if (verb)
    {
        syncwrite(total, ArrayLength, "total");
        syncwrite(input.arr, input.ArrayLength, "input");
        writefftw(filteredinput, ArrayLength, "filteredinput");
        writefftw(Fbase, ArrayLength, "Fbase");
        writefftw(corr, ArrayLength, "crosscor");
    }

    unsigned int poscor = getmaxfftw(corr, ArrayLength);

    double maxcor = corr[poscor][0];
    *corvalue = maxcor;

    if (total)
    {
        rescale(total, ArrayLength);

        for (unsigned int j = 0; j < ArrayLength; j++)
        {
            const int magic = 2000;
            total[j] += (int)(magic * maxcor * maxcor *
                              filteredinput[(j + poscor + ArrayLength) %
                                            ArrayLength][0]);
        }
    }

    const double half = 0.5;
    assert(ArrayLength > 0);
    *subpos = -half *
              (corr[(poscor - 1 + ArrayLength) % ArrayLength][0] -
               corr[(poscor + 1) % ArrayLength][0]) /
              (2 * corr[poscor][0] -
               corr[(poscor - 1 + ArrayLength) % ArrayLength][0] -
               corr[(poscor + 1) % ArrayLength][0]);

    fftw_free(filteredinput);
    fftw_free(Fbase);
    fftw_free(corr);
    fftw_cleanup();
    return poscor;
}

void rescale(int* total, unsigned int ArrayLength)
{
    double avg = 0.0;
    int maxval = -INT_MAX;
    int minval = INT_MAX;

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        avg += (double)total[j];
        maxval = total[j] > maxval ? total[j] : maxval;
        minval = total[j] < minval ? total[j] : minval;
    }
    avg /= ArrayLength;
    int avi = (int)avg;

    const int magic = 100000000;
    const int ten = 10;
    if (maxval > magic || minval < -magic)
    {
        for (unsigned int j = 0; j < ArrayLength; j++)
        {
            total[j] /= 2;
        }
    }
    else if (abs(avi) > ten)
    {
        for (unsigned int j = 0; j < ArrayLength; j++)
        {
            total[j] -= avi;
        }
    }
}

unsigned int getmaxfftw(fftw_complex* array, unsigned int ArrayLength)
{
    double maxtick = -INT_MAX;
    unsigned int postick = 0;
    for (unsigned int j = 0; j < ArrayLength / 3; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    for (unsigned int j = ArrayLength * 2 / 3; j < ArrayLength; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    return postick;
}

void writefftw(fftw_complex* arr, unsigned int ArrayLength, const char* file)
{
    FILE* filePtr = fopen(file, "w");
    if (filePtr == NULL)
    {
        (void)fprintf(stderr, "File opening failed in_data writefftw\n");
        return;
    }
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        (void)fprintf(filePtr, "%d %f %f\n", j, arr[j][0], arr[j][1]);
    }
    (void)fclose(filePtr);
}

void normalise(unsigned int ArrayLength, fftw_complex* in_data)
{
    double Sum_x = 0.0;
    double Sum_xx = 0.0;
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        in_data[j][0] = in_data[j][0];
        Sum_x += in_data[j][0];
        Sum_xx += in_data[j][0] * in_data[j][0];
    }
    double mean = Sum_x / ArrayLength;
    double stdev = sqrt((Sum_xx / ArrayLength) - (mean * mean));
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        in_data[j][0] = (in_data[j][0] - mean) / stdev;
    }
}
