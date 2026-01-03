#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mymath.h"
#include "mysignal.h"
#include "mysound.h"
#include "mysync.h"

#define ARR_BUFF 512
#define DEFAULT_RATE 48000
#define DEFAULT_BPH 21600
#define DEFAULT_ZOOM 10
#define DEFAULT_EVALUE 4
#define DEFAULT_FITN 30
#define DEFAULT_TIME 30
#define DEFAULT_TICKTOCK_WRITE 30
#define DEFAULT_TEETH 1
#define DEFAULT_SDTHRESHOLD 3.0
#define DEFAULT_CVALUE 7
#define PRESHIFT_THRESHOLD 100
#define PRESHIFT_THRESHOLD_ROOT 10
#define AUTOCOR_LIMIT 1
#define HALF 0.5
#define DEFAULT_COLUMNS 80
#define PCM_WIDTH 8

#define ERROR_ALLOCATE_MEM -5
#define ERROR_NO_SOURCE -6

volatile int keepRunning = 1;
volatile unsigned int columns = DEFAULT_COLUMNS;

int main(int argc, char* argv[])
{
    unsigned int everyline = 0;
    unsigned int cvalue = DEFAULT_CVALUE;
    unsigned int verbose = 0;
    unsigned int writeinterval = DEFAULT_TICKTOCK_WRITE;
    unsigned int lastWrite = 0;
    unsigned int fitN = DEFAULT_FITN;
    unsigned int teeth = DEFAULT_TEETH;
    unsigned int ticktockBuffer = ARR_BUFF * 2;
    double SDthreshold = DEFAULT_SDTHRESHOLD;
    char* device = NULL;
    FILE* fpposition = NULL;
    FILE* fpmaxcor = NULL;
    FILE* fptotal = NULL;
    FILE* fpDefPeak = NULL;
    FILE* fpInput = NULL;
    double a = 0.0;
    double b = 0.0;
    int maxposition = 0;
    int totalshift = 0;
    struct winsize w;

    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    columns = w.ws_col;
    set_signal_action();

    CapConfig cfg;
    cfg.rate = DEFAULT_RATE;
    cfg.bph = DEFAULT_BPH;
    cfg.evalue = DEFAULT_EVALUE;
    cfg.zoom = DEFAULT_ZOOM;
    cfg.time = 0;
    cfg.everyline = everyline;
    cfg.cvalue = cvalue;
    cfg.verbose = verbose;
    cfg.fitN = fitN;
    cfg.teeth = teeth;
    cfg.SDthreshold = SDthreshold;
    cfg.device = device;
    cfg.fpposition = fpposition;
    cfg.fpmaxcor = fpmaxcor;
    cfg.fptotal = fptotal;
    cfg.fpDefPeak = fpDefPeak;
    cfg.fpInput = fpInput;

    parse_arguments(argc, argv, &cfg);

    everyline = cfg.everyline;
    cvalue = cfg.cvalue;
    verbose = cfg.verbose;
    fitN = cfg.fitN;
    teeth = cfg.teeth;
    SDthreshold = cfg.SDthreshold;
    device = cfg.device;
    fpposition = cfg.fpposition;
    fpmaxcor = cfg.fpmaxcor;
    fptotal = cfg.fptotal;
    fpDefPeak = cfg.fpDefPeak;
    fpInput = cfg.fpInput;

    unsigned int actualRate = (unsigned int)(cfg.rate + HALF);

    if (fitN > ticktockBuffer / 2)
    {
        printf("Local fit N(%d) cannot be larger than %d\n",
               fitN,
               ticktockBuffer / 2);
    }

    device = (device == NULL) ? "default:2" : device;

    // initialize sound source
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = NULL;

    if (fpInput == NULL)
    {
        // rate could change, if not available
        printf("Casting inputrate %f(double) to soundcard rate %d(int)\n",
               cfg.rate,
               actualRate);
        capture_handle = initAudio(format, device, &actualRate);
        printf("Actual rate %d, calculating with %f\n", actualRate, cfg.rate);
    }

    if (fpInput == NULL && capture_handle == NULL)
    {
        (void)fprintf(stderr, "No inputfile or soundcard");
        return ERROR_NO_SOURCE;
    }

    unsigned int NN = actualRate * 2 * SECS_HOUR / cfg.bph;
    NN = (NN + NN % 2);
    unsigned int mod = NN / cfg.zoom;
    unsigned int maxtime =
        (unsigned int)cfg.rate * (cfg.time ? cfg.time : DEFAULT_TIME) / NN;

    fftw_complex* filterFFT = makeFilter(cfg.evalue, NN);

    struct myarr teethArray[teeth];
    for (unsigned int t = 0; t < teeth; t++)
    {
        teethArray[t].arr = calloc(NN, sizeof(int));
        if (teethArray[t].arr == NULL)
        {
            for (unsigned int idx = 0; idx < t; idx++)
            {
                free(teethArray[idx].arr);
            }

            (void)fprintf(stderr, "Could not allocate memory");
            return ERROR_ALLOCATE_MEM;
        }

        teethArray[t].arrd = NULL;
        teethArray[t].NN = NN;
    }

    struct myarr subpos = {
        0, calloc(ticktockBuffer, sizeof(double)), ticktockBuffer};
    struct myarr maxpos = {
        calloc(ticktockBuffer, sizeof(int)), 0, ticktockBuffer};
    struct myarr maxvals = {
        0, calloc(ticktockBuffer, sizeof(double)), ticktockBuffer};
    struct myarr derivative = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr tmpder = {calloc(NN, sizeof(int)), 0, NN};
    struct myarr reference = {calloc(NN, sizeof(int)), 0, NN};

    char* buffer =
        calloc(NN, (unsigned int)snd_pcm_format_width(format) / PCM_WIDTH);
    if (buffer == NULL || reference.arr == NULL || maxvals.arrd == NULL ||
        maxpos.arr == NULL || subpos.arrd == NULL || filterFFT == NULL ||
        derivative.arr == NULL || tmpder.arr == NULL)
    {
        free(buffer);
        freemyarr(&reference);
        freemyarr(&maxvals);
        freemyarr(&maxpos);
        freemyarr(&subpos);
        fftw_free(&filterFFT);
        freemyarr(&derivative);
        freemyarr(&tmpder);

        (void)fprintf(stderr, "Could not allocate memory");
        return ERROR_ALLOCATE_MEM;
    }

    (void)fprintf(stderr,
                  "Found COLUMNS=%d, width = %.3fms / "
                  "%.1fμs/character\n",
                  columns,
                  mod * 1000. / cfg.rate,
                  mod * 1000000. / cfg.rate / (columns - everyline));

    fillReference(fpDefPeak, &reference, teeth);

    sigset_t block;
    sigset_t non_block;
    setup_block_signals(&block);
    int toothshift = 0;
    unsigned int ticktock = 0;
    unsigned int totalTickTock = 0;
    while (keepRunning && !(totalTickTock > maxtime && cfg.time))
    {
        if (ticktock == ticktockBuffer)
        {
            // shift data back, has been written already
            memcpy(
                subpos.arrd, subpos.arrd + ARR_BUFF, ARR_BUFF * sizeof(double));
            memcpy(maxpos.arr, maxpos.arr + ARR_BUFF, ARR_BUFF * sizeof(int));
            memcpy(maxvals.arrd,
                   maxvals.arrd + ARR_BUFF,
                   ARR_BUFF * sizeof(double));
            ticktock -= ARR_BUFF;
        }

        block_signal(&block, &non_block);
        int err = getData(fpposition,
                          fpInput,
                          capture_handle,
                          format,
                          device,
                          cfg.rate,
                          buffer,
                          derivative);
        unblock_signal(&non_block);

        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }

        struct myarr* cumulativeTick = &teethArray[totalTickTock % teeth];
        if (totalTickTock >= AUTOCOR_LIMIT * teeth)
        {
            // make sure this is only done after the j teethArray are filled at
            // least once
            if (teeth > 1)
            {
                toothshift =
                    getshift(teethArray[0], teethArray[totalTickTock % teeth]);
            }

            if (totalTickTock == AUTOCOR_LIMIT * teeth)
            {
                free(reference.arr);
            }
            // use the appropriate tick as a reference
            reference.arr = cumulativeTick->arr;
        }

        for (unsigned int j = 0; j < NN; ++j)
        {
            // preshift the derivative
            tmpder.arr[j] =
                derivative.arr[modSigned(totalshift + j + toothshift, NN)];
        }

        maxposition =
            shiftHalf(fftfit(tmpder,
                             cumulativeTick->arr,
                             reference.arr,
                             maxvals.arrd + ticktock,
                             filterFFT,
                             totalTickTock > 0 && totalTickTock == verbose,
                             subpos.arrd + ticktock),
                      NN);

        maxpos.arr[ticktock] = totalshift + maxposition;

        if (totalTickTock > AUTOCOR_LIMIT &&
            *(maxvals.arrd + ticktock) > (double)cvalue / 16 &&
            totalTickTock % teeth == 0)
        {
            if (abs(maxposition) > PRESHIFT_THRESHOLD)
            {
                // shift with √ outside threshold
                maxposition = (int)(PRESHIFT_THRESHOLD_ROOT * maxposition /
                                    sqrt(abs(maxposition)));
            }

            totalshift += maxposition;
        }

        if (ticktock > 0 && totalTickTock % writeinterval == 0)
        {
            lastWrite = totalTickTock;

            if (fpposition)
            {
                struct myarr syncarr = {
                    0, calloc(writeinterval, sizeof(double)), writeinterval};
                if (syncarr.arrd != NULL)
                {
                    for (unsigned int k = 0; k < writeinterval; ++k)
                    {
                        syncarr.arrd[k] =
                            subpos.arrd[ticktock - writeinterval + k] +
                            (double)maxpos.arr[ticktock - writeinterval + k];
                    }
                    syncAppendMyarr(&syncarr, fpposition);
                    free(syncarr.arrd);
                }
            }
            if (fpmaxcor != NULL)
            {
                struct myarr syncarr = {
                    NULL, calloc(writeinterval, sizeof(double)), writeinterval};
                if (syncarr.arrd != NULL)
                {
                    memcpy(syncarr.arrd,
                           maxvals.arrd + ticktock - writeinterval,
                           writeinterval * sizeof(double));
                    syncAppendMyarr(&syncarr, fpmaxcor);
                    free(syncarr.arrd);
                }
            }
            syncwrite(teethArray->arr, NN, "/home/peter/tmp/livepeak");
        }

        fitNpeaks(&a, &b, ticktock, &maxvals, &maxpos, &subpos, fitN);

        printheader(b * 86400 / NN,
                    everyline,
                    getBeatError(cumulativeTick, cfg.rate, 0),
                    (double)totalTickTock * NN / cfg.rate);
        printspaces(maxpos.arr[ticktock],
                    maxvals.arrd[ticktock] * 16,
                    mod,
                    columns - everyline,
                    a,
                    cvalue);

        ticktock++;
        totalTickTock++;
    }

    if (teeth > 1)
    {
        printf("peak   shift beaterr\n");
        for (unsigned int k = 0; k < teeth; ++k)
        {
            double beaterr = getBeatError(&teethArray[k], cfg.rate, 0);
            printf("%6d%6d%6.2f\n",
                   k,
                   getshift(teethArray[0], teethArray[k]),
                   beaterr);
        }
    }

    free(buffer);
    freemyarr(&derivative);
    freemyarr(&tmpder);
    fftw_free(filterFFT);

    freemyarr(&reference);
    wait();
    if (fpmaxcor)
    {
        printTOD(fpmaxcor);
        writefileDouble(fpmaxcor,
                        maxvals.arrd + ticktock - (totalTickTock - lastWrite),
                        totalTickTock - lastWrite);
        (void)fclose(fpmaxcor);
    }
    if (fpposition)
    {
        thread_lock();
        unsigned int writelength = totalTickTock - lastWrite;
        struct myarr syncarr = {
            0, calloc(writelength, sizeof(double)), writelength};
        if (syncarr.arrd != NULL)
        {
            for (unsigned int k = 0; k < writelength; ++k)
            {
                syncarr.arrd[k] =
                    subpos.arrd[ticktock - writelength + k] +
                    (double)maxpos.arr[ticktock - writelength + k];
            }
            syncAppendMyarr(&syncarr, fpposition);
            freemyarr(&syncarr);
        }

        calculateTotalFromFile(totalTickTock, fpposition, NN, SDthreshold);
        thread_unlock();
        wait();
        (void)fclose(fpposition);
    }
    freemyarr(&subpos);

    freemyarr(&maxvals);
    freemyarr(&maxpos);

    for (unsigned int t = 0; t < teeth; ++t)
    {
        struct myarr* cumulativeTick = &teethArray[t];
        if (fptotal)
        {
            int toothshift = getshift(teethArray[0], *cumulativeTick);
            for (unsigned int j = 0; j < NN; ++j)
            {
                (void)fprintf(fptotal,
                              "%d %d %u %d\n",
                              shiftHalf(j + toothshift, NN),
                              cumulativeTick->arr[j],
                              t,
                              toothshift);
            }
            (void)fprintf(fptotal, "\n\n");
        }
    }
    for (unsigned int t = 0; t < teeth; ++t)
    {
        free(teethArray[t].arr);
    }
    if (fpInput)
    {
        (void)fclose(fpInput);
    }
    if (fptotal)
    {
        (void)fclose(fptotal);
    }

    if (capture_handle != NULL)
    {
        snd_pcm_close(capture_handle);
        snd_pcm_hw_free(capture_handle);
    }

    (void)fprintf(stderr,
                  "width = %.3fms / %.1fμs/character\n",
                  mod * 1000. / cfg.rate,
                  mod * 1000000. / cfg.rate / (columns - everyline));
    return 0;
}
