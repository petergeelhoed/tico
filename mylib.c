#include <limits.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mylib.h"

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

snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate)
{
    int err;
    snd_pcm_t* capture_handle;
    snd_pcm_hw_params_t* hw_params;

    if ((err = snd_pcm_open(
             &capture_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0)
    {
        fprintf(stderr,
                "cannot open audio device %s (%s)\n",
                device,
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_malloc(&hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot allocate hardware parameter structure (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_any(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot initialize hardware parameter structure (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_access(
             capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0)
    {
        fprintf(stderr, "cannot set access type (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_format(
             capture_handle, hw_params, format)) < 0)
    {
        fprintf(stderr, "cannot set sample format (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_rate_near(
             capture_handle, hw_params, &rate, 0)) < 0)
    {
        fprintf(stderr, "cannot set sample rate (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1)) <
        0)
    {
        fprintf(stderr, "cannot set channel count (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr, "cannot set parameters (%s)\n", snd_strerror(err));
        exit(1);
    }

    snd_pcm_hw_params_free(hw_params);

    if ((err = snd_pcm_prepare(capture_handle)) < 0)
    {
        fprintf(stderr,
                "cannot prepare audio interface for use (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    return capture_handle;
}

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

fftw_complex* convolute(int NN, int* array, const fftw_complex* filterFFT)
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

void readBufferRaw(snd_pcm_t* capture_handle, int NN, char* buffer, int* in)
{
    unsigned char lsb;
    signed char msb;
    int err;
    if ((err = snd_pcm_readi(capture_handle, buffer, NN)) != NN)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        exit(-1);
    }
    for (int j = 0; j < NN * 2; j += 2)
    {
        msb = *(buffer + j + 1);
        lsb = *(buffer + j);
        in[j / 2] = (msb << 8) | lsb;
    }
}

int readBuffer(snd_pcm_t* capture_handle, int NN, char* buffer, int* derivative)
{
    unsigned char lsb;
    signed char msb;
    int err;
    if ((err = snd_pcm_readi(capture_handle, buffer, NN)) != NN)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        return err;
    }
    for (int j = 0; j < NN * 2; j += 2)
    {
        msb = *(buffer + j + 1);
        lsb = *(buffer + j);
        derivative[j / 2] = (msb << 8) | lsb;
    }
    //       remove50hz(NN,in,48000);

    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
    derivative[NN] = 0;
    return err;
}

int readShiftedBuffer(int* derivative,
                      snd_pcm_t* capture_handle,
                      int NN,
                      char* buffer,
                      int maxpos,
                      int* totalshift,
                      FILE* fpInput)
{
    int ret;
    int shift = (int)sqrt(abs(maxpos));

    
    if (maxpos < 0)
    {
        *totalshift -= shift;
        memcpy(derivative + NN - shift, derivative, shift * sizeof(int));
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<NN-shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else
        {
        ret = readBuffer(capture_handle, NN - shift, buffer, derivative);
        }
    }
    else if (maxpos > 0)
    {
        *totalshift += shift;
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
            for (int j=0 ; j<NN; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else
        {
            ret = readBuffer(capture_handle, shift, buffer, derivative);

            if (ret != -32)
            {
                ret = readBuffer(capture_handle, NN, buffer, derivative);
            }
        }
    }
    else
    {
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else 
        {
            ret = readBuffer(capture_handle, NN, buffer, derivative);
        }
    }
    return ret;
}

void printspaces(int maxpos,
                 int hexvalue,
                 char* spaces,
                 int mod,
                 int columns,
                 double a,
                 double b,
                 int NN,
                 int cvalue,
                 float beatError)
{
    while (maxpos < mod)
        maxpos += mod;
    while (a < mod)
        a += (double)mod;

    int width = (maxpos % mod) * columns / mod;
    int widtha = (((int)a) % mod) * columns / mod;
    fprintf(stderr, "%5.2fms%6.1fs/d", beatError, b * 86400 / NN);
    memset(spaces, ' ', columns);
    spaces[widtha] = '|';
    spaces[width] = '\0';
    fprintf(stderr,
            "%s%s%X\e[0m",
            spaces,
            hexvalue < cvalue ? "\e[31m" : "\e[32m",
            hexvalue);
    memset(spaces, ' ', columns);
    if (widtha > width)
    {
        spaces[widtha - width - 1] = '|';
        spaces[widtha - width - 1 + 1] = '\0';
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

void writefiles(
    FILE* fptotal, FILE* rawfile, int* totaltick, int* maxpos, int n, int NN)
{
    if (fptotal)
    {
        for (int j = 0; j < NN; j++)
            fprintf(fptotal, "%d\n", totaltick[j]);
        fclose(fptotal);
    }
    if (rawfile)
    {
        for (int i = 0; i < n; ++i)
            fprintf(rawfile, "%d %d\n", i, maxpos[i]);
        fclose(rawfile);
    }
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

int getmaxposscaled(int* array, int NN)
{
    int maxtick = -INT_MAX;
    int postick = 0;
    int half = NN/2;

    for (int j = 0; j < NN; j++)
    {
        int scaled = (j<half)?j:NN-j;

        if (array[j] * (half - scaled) > maxtick * half )

        {
            maxtick = array[j];
            postick = j;
        }
    }
    return postick;
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
    //int postick = getmaxposscaled(cross, NN / 2);
    return (postick + NN / 4) % (NN / 2) - NN / 4;
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

void writearraydouble(double* arr, int NN, const char* file)
{
    FILE* fp = fopen(file, "w");
    for (int j = 0; j < NN; j++)
    {
        fprintf(fp, "%d %f\n", j, arr[j]);
    }
    fclose(fp);
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
