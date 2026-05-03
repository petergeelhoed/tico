#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mysync.h"

fftw_complex* makeFilter(size_t evalue, size_t ArrayLength)
{
    fftw_complex* filter = fftw_alloc_complex(ArrayLength);
    fftw_complex* filterFFT = fftw_alloc_complex(ArrayLength);
    fftw_plan makefilter = fftw_plan_dft_1d((int)ArrayLength,
                                            filter,
                                            filterFFT,
                                            FFTW_FORWARD,
                                            FFTW_ESTIMATE);

    if (filter == NULL || filterFFT == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed inData makeFilter\n");
        return NULL;
    }

    if (evalue != 0)
    {
        for (size_t j = 0; j < evalue * GAUSSPOINTS; j++)
        {
            filter[j][0] =
                GAUSSIAN_CONST / (double)evalue *
                exp(-((double)(j * j)) / ((double)(evalue * evalue) * HALF));
            filter[j][1] = 0.0;
        }
        for (size_t j = evalue * GAUSSPOINTS;
             j < ArrayLength - evalue * GAUSSPOINTS;
             j++)
        {
            filter[j][0] = 0.0;
            filter[j][1] = 0.0;
        }
        for (size_t j = ArrayLength - evalue * GAUSSPOINTS; j < ArrayLength;
             j++)
        {
            filter[j][0] =
                GAUSSIAN_CONST / (double)evalue *
                exp(-((double)((ArrayLength - j) * (ArrayLength - j))) /
                    ((double)(evalue * evalue) * HALF));
            filter[j][1] = 0.0;
        }
    }
    else
    {
        filter[0][0] = 1;
        filter[0][1] = 0;
        for (size_t j = 1; j < ArrayLength; j++)
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

void remove50hz(size_t ArrayLength, int* array, unsigned int rate)
{
    const size_t freq = 50;
    size_t ofj = (size_t)((double)ArrayLength * (double)freq / (double)rate);
    fftw_complex* inData = fftw_alloc_complex(ArrayLength);
    fftw_complex* outData = fftw_alloc_complex(ArrayLength);
    fftw_plan forward = fftw_plan_dft_1d((int)ArrayLength,
                                         inData,
                                         outData,
                                         FFTW_FORWARD,
                                         FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d((int)ArrayLength,
                                         outData,
                                         inData,
                                         FFTW_BACKWARD,
                                         FFTW_ESTIMATE);

    if (inData == NULL || outData == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed inData remove50hz\n");
        return;
    }

    for (size_t j = 0; j < ArrayLength; j++)
    {
        inData[j][0] = (double)array[j];
        inData[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (size_t j = 0; j < ArrayLength; j++)
    {
        outData[j][0] /= (double)ArrayLength;
        outData[j][1] /= (double)ArrayLength;
    }

    if (ofj >= 1 && ofj + 1 < ArrayLength)
    {
        outData[ofj + 1][1] = 0.0;
        outData[ofj - 1][0] = 0.0;
        outData[ofj + 1][1] = 0.0;
        outData[ofj - 1][0] = 0.0;
    }

    if (ofj > 0)
    {
        for (size_t j = ofj; j < ArrayLength; j += ofj)
        {
            outData[j][0] = 0.0;
            outData[j][1] = 0.0;
        }
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        array[j] = (int)(inData[j][0]);
    }
    fftw_free(inData);
    fftw_free(outData);
    fftw_cleanup();
}

fftw_complex* convolute(const struct myarr array, fftw_complex* filterFFT)
{
    size_t ArrayLength = array.ArrayLength;
    fftw_complex* inData = fftw_alloc_complex(ArrayLength);
    fftw_complex* outData = fftw_alloc_complex(ArrayLength);
    fftw_plan forward = fftw_plan_dft_1d((int)ArrayLength,
                                         inData,
                                         outData,
                                         FFTW_FORWARD,
                                         FFTW_ESTIMATE);
    fftw_plan reverse = fftw_plan_dft_1d((int)ArrayLength,
                                         outData,
                                         inData,
                                         FFTW_BACKWARD,
                                         FFTW_ESTIMATE);

    if (inData == NULL || outData == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed inData convolute\n");
        return NULL;
    }

    for (size_t j = 0; j < ArrayLength; j++)
    {
        inData[j][0] = (double)array.arr[j];
        inData[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (size_t j = 0; j < ArrayLength; j++)
    {
        double real = outData[j][0];
        double imag = outData[j][1];
        outData[j][0] = (real * filterFFT[j][0] - imag * filterFFT[j][1]) /
                        (double)ArrayLength;
        outData[j][1] = (real * filterFFT[j][1] + imag * filterFFT[j][0]) /
                        (double)ArrayLength;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);
    fftw_free(outData);
    fftw_cleanup();
    return inData;
}

fftw_complex* crosscor(size_t ArrayLength,
                       fftw_complex* array,
                       fftw_complex* ref)
{
    normalise(ArrayLength, array);
    normalise(ArrayLength, ref);

    fftw_complex* tmparr = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmpref = fftw_alloc_complex(ArrayLength);
    fftw_complex* tmp = fftw_alloc_complex(ArrayLength);
    fftw_complex* corr = fftw_alloc_complex(ArrayLength);

    if (tmparr == NULL || tmpref == NULL || tmp == NULL || corr == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed inData crosscor\n");
        return NULL;
    }

    fftw_plan arrFour = fftw_plan_dft_1d((int)ArrayLength,
                                         array,
                                         tmparr,
                                         FFTW_FORWARD,
                                         FFTW_ESTIMATE);
    fftw_plan refFour = fftw_plan_dft_1d((int)ArrayLength,
                                         ref,
                                         tmpref,
                                         FFTW_FORWARD,
                                         FFTW_ESTIMATE);
    fftw_plan correverse = fftw_plan_dft_1d((int)ArrayLength,
                                            tmp,
                                            corr,
                                            FFTW_BACKWARD,
                                            FFTW_ESTIMATE);

    fftw_execute(arrFour);
    fftw_execute(refFour);

    for (size_t j = 0; j < ArrayLength; j++)
    {
        tmp[j][0] = (tmparr[j][0] * tmpref[j][0] + tmparr[j][1] * tmpref[j][1]);
        tmp[j][1] =
            (-tmparr[j][0] * tmpref[j][1] + tmparr[j][1] * tmpref[j][0]);
    }
    fftw_execute(correverse);

    double scale = 1.0 / ((double)ArrayLength * (double)ArrayLength);
    for (size_t j = 0; j < ArrayLength; j++)
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
    if (xarr.ArrayLength != yarr.ArrayLength || xarr.ArrayLength == 0)
    {
        return 0;
    }

    fftw_complex* F_x = fftw_alloc_complex(xarr.ArrayLength);
    fftw_complex* F_y = fftw_alloc_complex(yarr.ArrayLength);

    for (size_t j = 0; j < xarr.ArrayLength; j++)
    {
        F_x[j][0] = (double)xarr.arr[j];
        F_x[j][1] = 0.0;
        F_y[j][0] = (double)yarr.arr[j];
        F_y[j][1] = 0.0;
    }
    fftw_complex* corr = crosscor(xarr.ArrayLength, F_x, F_y);

    size_t poscor = getmaxfftw(corr, xarr.ArrayLength);
    fftw_free(F_x);
    fftw_free(corr);
    fftw_free(F_y);
    fftw_cleanup();

    return (((int)poscor + (int)xarr.ArrayLength / 2) %
            (int)(xarr.ArrayLength)) -
           (int)(xarr.ArrayLength / 2);
}

size_t fftfit(const struct myarr input,
              int* total,
              const int* base,
              double* corvalue,
              fftw_complex* filterFFT,
              int verb,
              double* subpos)
{
    size_t ArrayLength = input.ArrayLength;
    if (ArrayLength == 0)
    {
        return 0;
    }

    fftw_complex* Fbase = fftw_alloc_complex(ArrayLength);
    fftw_complex* filteredinput = convolute(input, filterFFT);

    if (Fbase == NULL || filteredinput == NULL)
    {
        (void)fprintf(stderr, "Memory allocation failed inData fftfit\n");
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

    size_t poscor = getmaxfftw(corr, ArrayLength);

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

void rescale(int* total, size_t ArrayLength)
{
    double avg = 0.0;
    int maxval = -INT_MAX;
    int minval = INT_MAX;

    for (size_t j = 0; j < ArrayLength; j++)
    {
        avg += (double)total[j];
        maxval = total[j] > maxval ? total[j] : maxval;
        minval = total[j] < minval ? total[j] : minval;
    }
    avg /= (double)ArrayLength;
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

size_t getmaxfftw(fftw_complex* array, size_t ArrayLength)
{
    double maxtick = -INT_MAX;
    size_t postick = 0;
    for (size_t j = 0; j < ArrayLength / 3; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    for (size_t j = ArrayLength * 2 / 3; j < ArrayLength; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    return postick;
}

void writefftw(fftw_complex* arr, size_t ArrayLength, const char* file)
{
    FILE* filePtr = fopen(file, "w");
    if (filePtr == NULL)
    {
        (void)fprintf(stderr, "File opening failed inData writefftw\n");
        return;
    }
    for (size_t j = 0; j < ArrayLength; j++)
    {
        (void)fprintf(filePtr, "%zu %f %f\n", j, arr[j][0], arr[j][1]);
    }
    (void)fclose(filePtr);
}

void normalise(size_t ArrayLength, fftw_complex* inData)
{
    double Sum_x = 0.0;
    double Sum_xx = 0.0;
    for (size_t j = 0; j < ArrayLength; j++)
    {
        inData[j][0] = inData[j][0];
        Sum_x += inData[j][0];
        Sum_xx += inData[j][0] * inData[j][0];
    }
    double mean = Sum_x / (double)ArrayLength;
    double stdev = sqrt((Sum_xx / (double)ArrayLength) - (mean * mean));
    for (size_t j = 0; j < ArrayLength; j++)
    {
        inData[j][0] = (inData[j][0] - mean) / stdev;
    }
}
