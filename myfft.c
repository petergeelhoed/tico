#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"
#include "mysync.h"

#define MAX_COLUMNS 1024
#define GAUSSIAN_CONST .398942280401

fftw_complex* makeFilter(unsigned int evalue, unsigned int NN)
{
    fftw_complex* in2 = fftw_alloc_complex(NN);
    fftw_complex* filterFFT = fftw_alloc_complex(NN);
    fftw_plan makefilter =
        fftw_plan_dft_1d((int)NN, in2, filterFFT, FFTW_FORWARD, FFTW_ESTIMATE);

    if (in2 == NULL || filterFFT == NULL)
    {
        fprintf(stderr, "Memory allocation failed in makeFilter\n");
        return NULL;
    }

    if (evalue != 0)
    {
        for (unsigned int j = 0; j < evalue * 5; j++)
        {
            in2[j][0] = GAUSSIAN_CONST / evalue *
                        exp(-((double)(j * j)) / (double)(evalue * evalue) / 2);
            in2[j][1] = 0.0;
        }
        for (unsigned int j = evalue * 5; j < NN - evalue * 5; j++)
        {
            in2[j][0] = 0.0;
            in2[j][1] = 0.0;
        }
        for (unsigned int j = NN - evalue * 5; j < NN; j++)
        {
            in2[j][0] = GAUSSIAN_CONST / evalue *
                        exp(-((double)(NN - j) * (NN - j)) /
                            (double)(evalue * evalue) / 2);
            in2[j][1] = 0.0;
        }
    }
    else
    {
        in2[0][0] = 100;
        in2[0][1] = 0;
        for (unsigned int j = 1; j < NN; j++)
        {
            in2[j][0] = 0;
            in2[j][1] = 0;
        }
    }

    fftw_execute(makefilter);
    fftw_destroy_plan(makefilter);
    fftw_free(in2);
    fftw_cleanup();
    return filterFFT;
}

void remove50hz(unsigned int NN, int* array, unsigned int rate)
{
    unsigned int freq = 50;
    unsigned int ofj = NN * freq / rate;
    fftw_complex* in = fftw_alloc_complex(NN);
    fftw_complex* out = fftw_alloc_complex(NN);
    fftw_plan forward =
        fftw_plan_dft_1d((int)NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse =
        fftw_plan_dft_1d((int)NN, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    if (in == NULL || out == NULL)
    {
        fprintf(stderr, "Memory allocation failed in remove50hz\n");
        return;
    }

    for (unsigned int j = 0; j < NN; j++)
    {
        in[j][0] = (double)array[j];
        in[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (unsigned int j = 0; j < NN; j++)
    {
        out[j][0] /= NN;
        out[j][1] /= NN;
    }

    out[ofj + 1][1] = 0.0;
    out[ofj - 1][0] = 0.0;
    out[ofj + 1][1] = 0.0;
    out[ofj - 1][0] = 0.0;

    for (unsigned int j = ofj; j < NN; j += ofj)
    {
        out[j][0] = 0.0;
        out[j][1] = 0.0;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);

    for (unsigned int j = 0; j < NN; j++)
    {
        array[j] = (int)(in[j][0]);
    }
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
}

fftw_complex* convolute(const struct myarr array, fftw_complex* filterFFT)
{
    unsigned int NN = array.NN;
    fftw_complex* in = fftw_alloc_complex(NN);
    fftw_complex* out = fftw_alloc_complex(NN);
    fftw_plan forward =
        fftw_plan_dft_1d((int)NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse =
        fftw_plan_dft_1d((int)NN, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    if (in == NULL || out == NULL)
    {
        fprintf(stderr, "Memory allocation failed in convolute\n");
        return NULL;
    }

    for (unsigned int j = 0; j < NN; j++)
    {
        in[j][0] = (double)array.arr[j];
        in[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (unsigned int j = 0; j < NN; j++)
    {
        out[j][0] =
            (out[j][0] * filterFFT[j][0] - out[j][1] * filterFFT[j][1]) / NN;
        out[j][1] =
            (out[j][0] * filterFFT[j][1] + out[j][1] * filterFFT[j][0]) / NN;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);
    fftw_free(out);
    fftw_cleanup();
    return in;
}

fftw_complex* crosscor(unsigned int NN, fftw_complex* array, fftw_complex* ref)
{
    normalise(NN, array);
    normalise(NN, ref);

    fftw_complex* tmparr = fftw_alloc_complex(NN);
    fftw_complex* tmpref = fftw_alloc_complex(NN);
    fftw_complex* tmp = fftw_alloc_complex(NN);
    fftw_complex* corr = fftw_alloc_complex(NN);

    if (tmparr == NULL || tmpref == NULL || tmp == NULL || corr == NULL)
    {
        fprintf(stderr, "Memory allocation failed in crosscor\n");
        return NULL;
    }

    fftw_plan arrFour =
        fftw_plan_dft_1d((int)NN, array, tmparr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan refFour =
        fftw_plan_dft_1d((int)NN, ref, tmpref, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan correverse =
        fftw_plan_dft_1d((int)NN, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(arrFour);
    fftw_execute(refFour);

    for (unsigned int j = 0; j < NN; j++)
    {
        tmp[j][0] = (tmparr[j][0] * tmpref[j][0] + tmparr[j][1] * tmpref[j][1]);
        tmp[j][1] =
            (-tmparr[j][0] * tmpref[j][1] + tmparr[j][1] * tmpref[j][0]);
    }
    fftw_execute(correverse);

    double scale = 1. / NN / NN;
    for (unsigned int j = 0; j < NN; j++)
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

void applyFilter(const struct myarr input, fftw_complex* filterFFT, double* out)
{
    fftw_complex* filteredinput = convolute(input, filterFFT);
    if (filteredinput == NULL)
    {
        fprintf(stderr, "Convolution failed in applyFilter\n");
        return;
    }
    for (unsigned int j = 0; j < input.NN; j++)
    {
        out[j] = filteredinput[j][0];
    }
    fftw_free(filteredinput);
    fftw_cleanup();
}

int getshift(const struct myarr x, const struct myarr y)
{
    unsigned int NN = x.NN;
    fftw_complex* Fx = fftw_alloc_complex(NN);
    fftw_complex* Fy = fftw_alloc_complex(NN);

    for (unsigned int j = 0; j < NN; j++)
    {
        Fx[j][0] = (double)x.arr[j];
        Fx[j][1] = 0.0;
        Fy[j][0] = (double)y.arr[j];
        Fy[j][1] = 0.0;
    }
    fftw_complex* corr = crosscor(NN, Fx, Fy);

    unsigned int poscor = getmaxfftw(corr, NN);
    fftw_free(Fx);
    fftw_free(corr);
    fftw_free(Fy);
    fftw_cleanup();
    return ((int)poscor + (int)NN / 2) % (int)(NN) - (int)(NN / 2);
}

unsigned int fftfit(const struct myarr input,
                    int* total,
                    int* base,
                    double* corvalue,
                    fftw_complex* filterFFT,
                    int verb)
{
    unsigned int NN = input.NN;
    fftw_complex* Fbase = fftw_alloc_complex(NN);
    fftw_complex* filteredinput = convolute(input, filterFFT);

    if (Fbase == NULL || filteredinput == NULL)
    {
        fprintf(stderr, "Memory allocation failed in fftfit\n");
        return 0;
    }

    for (unsigned int j = 0; j < NN; j++)
    {
        Fbase[j][0] = (double)base[j];
        Fbase[j][1] = 0.0;
    }

    fftw_complex* corr = crosscor(NN, filteredinput, Fbase);

    if (verb)
    {
        syncwrite(total, NN, "total");
        syncwrite(input.arr, input.NN, "input");
        writefftw(filteredinput, NN, "filteredinput");
        writefftw(Fbase, NN, "Fbase");
        writefftw(corr, NN, "crosscor");
    }

    unsigned int poscor = getmaxfftw(corr, NN);

    double maxcor = corr[poscor][0];
    *corvalue = maxcor;

    if (total)
    {
        rescale(total, NN);

        for (unsigned int j = 0; j < NN; j++)
        {
            total[j] += (int)(2000 * maxcor * maxcor *
                              filteredinput[(j + poscor + NN) % NN][0]);
        }
    }
    fftw_free(filteredinput);
    fftw_free(Fbase);
    fftw_free(corr);
    fftw_cleanup();
    return poscor;
}

void rescale(int* total, unsigned int NN)
{
    double avg = 0.0;
    int maxval = -INT_MAX;
    int minval = INT_MAX;

    for (unsigned int j = 0; j < NN; j++)
    {
        avg += (double)total[j];
        maxval = total[j] > maxval ? total[j] : maxval;
        minval = total[j] < minval ? total[j] : minval;
    }
    avg /= NN;
    int avi = (int)avg;

    if (maxval > 100000000 || minval < -100000000)
    {
        for (unsigned int j = 0; j < NN; j++)
        {
            total[j] /= 2;
        }
    }
    else if (abs(avi) > 10)
    {
        for (unsigned int j = 0; j < NN; j++)
        {
            total[j] -= avi;
        }
    }
}

unsigned int getmaxfftw(fftw_complex* array, unsigned int NN)
{
    double maxtick = -INT_MAX;
    unsigned int postick = 0;
    for (unsigned int j = 0; j < NN / 3; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    for (unsigned int j = NN * 2 / 3; j < NN; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    return postick;
}

void writefftw(fftw_complex* arr, unsigned int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "File opening failed in writefftw\n");
        return;
    }
    for (unsigned int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %f %f\n", j, arr[j][0], arr[j][1]);
    }
    fclose(fp);
}

void normalise(unsigned int NN, fftw_complex* in)
{
    double ix = 0.0;
    double ixx = 0.0;
    for (unsigned int j = 0; j < NN; j++)
    {
        in[j][0] = in[j][0];
        ix += in[j][0];
        ixx += in[j][0] * in[j][0];
    }
    double m = ix / NN;
    double s = sqrt(ixx / NN - m * m);
    for (unsigned int j = 0; j < NN; j++)
    {
        in[j][0] = (in[j][0] - m) / s;
    }
}
