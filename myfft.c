#include <limits.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"


fftw_complex* makeFilter(int evalue, int NN)
{
    fftw_complex* in2 = fftw_alloc_complex(NN);
    fftw_complex* filterFFT = fftw_alloc_complex(NN);
    fftw_plan makefilter =
        fftw_plan_dft_1d(NN, in2, filterFFT, FFTW_FORWARD, FFTW_ESTIMATE);

    if (evalue != 0)
    {
        // make filter array
        for (int j = 0; j < evalue * 5; j++)
        {
            in2[j][0] =
                .398942280401 / evalue *
                (exp(-((double)(j * j)) / (double)(evalue * evalue) / 2));
            in2[j][1] = 0.0;
        }
        for (int j = evalue * 5; j < NN - evalue * 5; j++)
        {
            in2[j][0] = 0.0;
            in2[j][1] = 0.0;
        }

        for (int j = NN - evalue * 5; j < NN; j++)
        {
            in2[j][0] =
                .398942280401 / evalue * (exp(-((double)(NN - j) * (NN - j)) /
                                              (double)(evalue * evalue) / 2));
            in2[j][1] = 0.0;
        }
    }
    else
    {
        in2[0][0] = 100;
        in2[0][1] = 0;
        for (int j = 1; j < NN; j++)
        {
            in2[j][0] = 0;
            in2[j][1] = 0;
        }
    }

    fftw_execute(makefilter);
    fftw_destroy_plan(makefilter);
    fftw_free(in2);
    return filterFFT;
}

void remove50hz(int NN, int* array, int rate)
{
    int freq = 50;
    int ofj = NN * freq / rate;
    if (ofj > 0)
    {
        fftw_complex* in = fftw_alloc_complex(NN);
        fftw_complex* out = fftw_alloc_complex(NN);
        fftw_plan forward =
            fftw_plan_dft_1d(NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan reverse =
            fftw_plan_dft_1d(NN, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
        for (int j = 0; j < NN; j++)
        {
            in[j][0] = (double)array[j];
            in[j][1] = 0.0;
        }
        fftw_execute(forward);

        for (int j = 0; j < NN; j++)
        {
            out[j][0] = out[j][0] / NN;
            out[j][1] = out[j][1] / NN;
        }

        out[ofj + 1][1] = 0.0;
        out[ofj - 1][0] = 0.0;
        out[ofj + 1][1] = 0.0;
        out[ofj - 1][0] = 0.0;

        for (int j = ofj; j < NN; j += ofj)
        {
            out[j][0] = 0.0;
            out[j][1] = 0.0;
        }

        fftw_execute(reverse);
        fftw_destroy_plan(forward);
        fftw_destroy_plan(reverse);
        for (int j = 0; j < NN; j++)
        {
            array[j] = (int)(in[j][0]);
        }
        fftw_free(*in);
        fftw_free(*out);
    }
}

fftw_complex* convolute(int NN, int* array, fftw_complex* filterFFT)
{
    fftw_complex* in = fftw_alloc_complex(NN);
    fftw_complex* out = fftw_alloc_complex(NN);
    fftw_plan forward =
        fftw_plan_dft_1d(NN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan reverse =
        fftw_plan_dft_1d(NN, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int j = 0; j < NN; j++)
    {
        in[j][0] = (double)array[j];
        in[j][1] = 0.0;
    }
    fftw_execute(forward);

    for (int j = 0; j < NN; j++)
    {
        out[j][0] =
            (out[j][0] * filterFFT[j][0] - out[j][1] * filterFFT[j][1]) / NN;
        out[j][1] =
            (out[j][0] * filterFFT[j][1] + out[j][1] * filterFFT[j][0]) / NN;
    }

    fftw_execute(reverse);
    fftw_destroy_plan(forward);
    fftw_destroy_plan(reverse);
    fftw_free(*out);
    return in;
}

fftw_complex* crosscor(int NN, fftw_complex* array, fftw_complex* ref)
{
    normalise(NN, array);
    normalise(NN, ref);
    fftw_complex* tmparr = fftw_alloc_complex(NN);
    fftw_complex* tmpref = fftw_alloc_complex(NN);
    fftw_complex* tmp = fftw_alloc_complex(NN);
    fftw_complex* corr = fftw_alloc_complex(NN);

    fftw_plan arrFour =
        fftw_plan_dft_1d(NN, array, tmparr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan refFour =
        fftw_plan_dft_1d(NN, ref, tmpref, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan correverse =
        fftw_plan_dft_1d(NN, tmp, corr, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(arrFour);
    fftw_execute(refFour);

    for (int j = 0; j < NN; j++)
    {
        tmp[j][0] = (tmparr[j][0] * tmpref[j][0] + tmparr[j][1] * tmpref[j][1]);
        tmp[j][1] =
            (-tmparr[j][0] * tmpref[j][1] + tmparr[j][1] * tmpref[j][0]);
    }
    fftw_execute(correverse);

    double scale = 1. / NN / NN;
    for (int j = 0; j < NN; j++)
    {
        corr[j][0] *= scale;
    }

    fftw_destroy_plan(correverse);
    fftw_destroy_plan(arrFour);
    fftw_destroy_plan(refFour);

    fftw_free(*tmparr);
    fftw_free(*tmpref);
    fftw_free(*tmp);
    return corr;
}

void applyFilter(int* input, int NN, fftw_complex* filterFFT, double* out)
{
    fftw_complex* filteredinput = convolute(NN, input, filterFFT);
    for (int j = 0; j < NN; j++)
    {
        out[j] = filteredinput[j][0];
    }
    fftw_free(filteredinput);
}

int fftfit(int* input,
           int* total,
           int* base,
           int* hexvalue,
           fftw_complex* filterFFT,
           int NN,
           int verb)
{
    fftw_complex* Fbase = fftw_alloc_complex(NN);
    fftw_complex* filteredinput = convolute(NN, input, filterFFT);

    for (int j = 0; j < NN; j++)
    {
        Fbase[j][0] = (double)base[j];
        Fbase[j][1] = 0.0;
    }

    fftw_complex* corr = crosscor(NN, filteredinput, Fbase);
    if (verb)
    {
        syncwrite(total, NN, "total");
        syncwrite(input, NN, "input");
    }

        if (verb) writefftw(filteredinput,NN,"filteredinput");
        if (verb) writefftw(Fbase,NN,"Fbase");
        if (verb) writefftw(corr,NN,"crosscor");

    int poscor = getmaxfftw(corr, NN);

    // for hexadecimal print
    double maxcor = corr[poscor][0];
    *hexvalue = (int)(maxcor * 16);

    // needed to add to total
    poscor -= NN / 2;

    // rescale if large
    if (total)
    {
        rescale(total, NN);

        // weigh with square of correlation
        for (int j = 0; j < NN; j++)
        {
            total[j] = (total[j] +
                        (int)(2000 * maxcor * maxcor) *
                            filteredinput[(j + poscor + NN / 2 + NN) % NN][0]);
        }
    }
    fftw_free(*filteredinput);
    fftw_free(*Fbase);
    fftw_free(*corr);

    return poscor;
}

int getmaxfftw(fftw_complex* array, int NN)
{
    double maxtick = -INT_MAX;
    int postick = 0;
    for (int j = NN / 3; j < 2 * NN / 3; j++)
    {
        if (array[j][0] > maxtick)
        {
            maxtick = array[j][0];
            postick = j;
        }
    }
    return postick;
}

void crosscorint(int NN, int* array, int* ref, int* cross)
{

    fftw_complex* tmparr = fftw_alloc_complex(NN);
    fftw_complex* tmpref = fftw_alloc_complex(NN);
    for (int j = 0; j < NN; j++)
    {
        tmparr[j][0] = array[j];
        tmparr[j][1] = 0.0;
        tmpref[j][0] = ref[j];
        tmpref[j][1] = 0.0;
    }
    fftw_complex* coor = crosscor(NN, tmparr, tmpref);

    for (int j = 0; j < NN; j++)
    {
        cross[j] = (int)(coor[j][0] * NN);
    }
    fftw_free(tmpref);
    fftw_free(tmparr);
    fftw_free(coor);
}

void writefftw(fftw_complex* arr, int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    for (int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %f %f\n", j, arr[j][0], arr[j][1]);
    }
    fclose(fp);
}


void normalise(int NN, fftw_complex* in)
{
    double ix = 0.0;
    double ixx = 0.0;
    for (int j = 0; j < NN; j++)
    {
        in[j][0] = in[j][0];
        ix += in[j][0];
        ixx += in[j][0] * in[j][0];
    }
    double m = ix / NN;
    double s = sqrt(ixx / NN - m * m);
    for (int j = 0; j < NN; j++)
    {
        in[j][0] = (in[j][0] - m) / s;
    }
}

