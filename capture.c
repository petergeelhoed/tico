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

#define KILO 1000.
#define MEGA 1000000.
#define ARR_BUFF 512
#define SECS_DAY 86400
#define HEXDEC 16
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

void shift_buffer_data(unsigned int* ticktock,
                       struct myarr* subpos,
                       struct myarr* maxpos,
                       struct myarr* maxvals)
{
    memmove(subpos->arrd, subpos->arrd + ARR_BUFF, ARR_BUFF * sizeof(double));
    memmove(maxpos->arr, maxpos->arr + ARR_BUFF, ARR_BUFF * sizeof(int));
    memmove(maxvals->arrd, maxvals->arrd + ARR_BUFF, ARR_BUFF * sizeof(double));

    *ticktock -= ARR_BUFF;
}

int main(int argc, char* argv[])
{
    unsigned int writeinterval = DEFAULT_TICKTOCK_WRITE;
    unsigned int lastWrite = 0;
    unsigned int ticktockBuffer = ARR_BUFF * 2;
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
    cfg.everyline = 0;
    cfg.cvalue = DEFAULT_CVALUE;
    cfg.verbose = 0;
    cfg.fitN = DEFAULT_FITN;
    cfg.teeth = DEFAULT_TEETH;
    cfg.SDthreshold = DEFAULT_SDTHRESHOLD;
    cfg.device = NULL;
    cfg.fpposition = NULL;
    cfg.fpmaxcor = NULL;
    cfg.fptotal = NULL;
    cfg.fpDefPeak = NULL;
    cfg.fpInput = NULL;

    parse_arguments(argc, argv, &cfg);

    unsigned int actualRate = (unsigned int)(cfg.rate + HALF);

    if (cfg.fitN > ticktockBuffer / 2)
    {
        printf("Local fit N(%d) cannot be larger than %d\n",
               cfg.fitN,
               ticktockBuffer / 2);
    }

    cfg.device = (cfg.device == NULL) ? "default:2" : cfg.device;

    // initialize sound source
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = NULL;

    if (cfg.fpInput == NULL)
    {
        // rate could change, if not available
        printf("Casting inputrate %f(double) to soundcard rate %d(int)\n",
               cfg.rate,
               actualRate);
        capture_handle = initAudio(format, cfg.device, &actualRate);
        printf("Actual rate %d, calculating with %f\n", actualRate, cfg.rate);
    }

    if (cfg.fpInput == NULL && capture_handle == NULL)
    {
        (void)fprintf(stderr, "No inputfile or soundcard");
        return ERROR_NO_SOURCE;
    }

    unsigned int NN = actualRate * 2 * SECS_HOUR / cfg.bph;
    NN = (NN + NN % 2);
    unsigned int mod = NN / cfg.zoom;
    const unsigned int maxtime =
        (unsigned int)cfg.rate * (cfg.time ? cfg.time : DEFAULT_TIME) / NN;

    fftw_complex* filterFFT = makeFilter(cfg.evalue, NN);

    struct myarr* teethArray[cfg.teeth];
    for (unsigned int t = 0; t < cfg.teeth; t++)
    {
        teethArray[t] = makemyarr(NN);
    }

    struct myarr* subpos = makemyarrd(ticktockBuffer);
    struct myarr* maxpos = makemyarr(ticktockBuffer);
    struct myarr* maxvals = makemyarrd(ticktockBuffer);
    struct myarr* derivative = makemyarr(NN);
    struct myarr* tmpder = makemyarr(NN);
    struct myarr* reference = makemyarr(NN);
    // struct myarr reference = {calloc(NN, sizeof(int)), 0, NN};

    char* buffer =
        calloc(NN, (unsigned int)snd_pcm_format_width(format) / PCM_WIDTH);
    if (buffer == NULL || reference->arr == NULL || maxvals->arrd == NULL ||
        maxpos->arr == NULL || subpos->arrd == NULL || filterFFT == NULL ||
        derivative->arr == NULL || tmpder->arr == NULL)
    {
        free(buffer);
        freemyarr(subpos);
        freemyarr(reference);
        freemyarr(maxvals);
        freemyarr(maxpos);
        fftw_free(&filterFFT);
        freemyarr(derivative);
        freemyarr(tmpder);

        for (unsigned int idx = 0; idx < cfg.teeth; idx++)
        {
            freemyarr(teethArray[idx]);
        }
        (void)fprintf(stderr, "Could not allocate memory");
        return ERROR_ALLOCATE_MEM;
    }

    (void)fprintf(stderr,
                  "Found COLUMNS=%d, width = %.3fms / "
                  "%.1fμs/character\n",
                  columns,
                  mod * KILO / cfg.rate,
                  mod * MEGA / cfg.rate / (columns - cfg.everyline));

    fillReference(cfg.fpDefPeak, reference, cfg.teeth);

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
            shift_buffer_data(&ticktock, subpos, maxpos, maxvals);
        }

        block_signal(&block, &non_block);
        int err = getData(cfg.fpposition,
                          cfg.fpInput,
                          capture_handle,
                          format,
                          cfg.device,
                          cfg.rate,
                          buffer,
                          *derivative);
        unblock_signal(&non_block);

        if (err < 0)
        {
            printf("capture error %d\n", err);
            break;
        }

        struct myarr* cumulativeTick = teethArray[totalTickTock % cfg.teeth];
        if (cfg.teeth > 1 && totalTickTock >= AUTOCOR_LIMIT * cfg.teeth)
        {
            // make sure this is only done after the j teethArray are filled at
            // least once
            toothshift = getshift(*teethArray[0],
                                  *teethArray[totalTickTock % cfg.teeth]);
        }

        for (unsigned int j = 0; j < NN; ++j)
        {
            // preshift the derivative
            tmpder->arr[j] =
                derivative->arr[modSigned(totalshift + j + toothshift, NN)];
        }

        maxposition =
            shiftHalf(fftfit(*tmpder,
                             cumulativeTick->arr,
                             (totalTickTock < AUTOCOR_LIMIT * cfg.teeth)
                                 ? reference->arr
                                 : cumulativeTick->arr,
                             maxvals->arrd + ticktock,
                             filterFFT,
                             totalTickTock > 0 && totalTickTock == cfg.verbose,
                             subpos->arrd + ticktock),
                      NN);

        maxpos->arr[ticktock] = totalshift + maxposition;

        if (totalTickTock > AUTOCOR_LIMIT &&
            *(maxvals->arrd + ticktock) > (double)cfg.cvalue / HEXDEC &&
            totalTickTock % cfg.teeth == 0)
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

            if (cfg.fpposition)
            {
                struct myarr* syncarr = makemyarrd(writeinterval);
                if (syncarr != NULL)
                {
                    for (unsigned int k = 0; k < writeinterval; ++k)
                    {
                        syncarr->arrd[k] =
                            subpos->arrd[ticktock - writeinterval + k] +
                            (double)maxpos->arr[ticktock - writeinterval + k];
                    }
                    syncAppendMyarr(syncarr, cfg.fpposition);
                    freemyarr(syncarr);
                }
            }
            if (cfg.fpmaxcor != NULL)
            {
                struct myarr* syncarr = makemyarrd(writeinterval);
                if (syncarr != NULL)
                {
                    memcpy(syncarr->arrd,
                           maxvals->arrd + ticktock - writeinterval,
                           writeinterval * sizeof(double));
                    syncAppendMyarr(syncarr, cfg.fpmaxcor);
                    freemyarr(syncarr);
                }
            }
            //  syncwrite(teethArray->arr, NN, "/home/peter/tmp/livepeak");
        }

        fitNpeaks(&a, &b, ticktock, maxvals, maxpos, subpos, cfg.fitN);

        printheader(b * SECS_DAY / NN,
                    cfg.everyline,
                    getBeatError(cumulativeTick, cfg.rate, 0),
                    (double)totalTickTock * NN / cfg.rate);
        printspaces(maxpos->arr[ticktock],
                    maxvals->arrd[ticktock] * HEXDEC,
                    mod,
                    columns - cfg.everyline,
                    a,
                    cfg.cvalue);

        ticktock++;
        totalTickTock++;
    }

    if (cfg.teeth > 1)
    {
        printf("peak   shift beaterr\n");
        for (unsigned int k = 0; k < cfg.teeth; ++k)
        {
            double beaterr = getBeatError(teethArray[k], cfg.rate, 0);
            printf("%6d%6d%6.2f\n",
                   k,
                   getshift(*teethArray[0], *teethArray[k]),
                   beaterr);
        }
    }

    free(buffer);
    freemyarr(derivative);
    freemyarr(tmpder);
    fftw_free(filterFFT);

    wait();
    if (cfg.fpmaxcor)
    {
        printTOD(cfg.fpmaxcor);
        writefileDouble(cfg.fpmaxcor,
                        maxvals->arrd + ticktock - (totalTickTock - lastWrite),
                        totalTickTock - lastWrite);
        (void)fclose(cfg.fpmaxcor);
    }
    if (cfg.fpposition)
    {
        thread_lock();
        unsigned int writelength = totalTickTock - lastWrite;
        struct myarr* syncarr = makemyarrd(writelength);
        if (syncarr != NULL)
        {
            for (unsigned int k = 0; k < writelength; ++k)
            {
                syncarr->arrd[k] =
                    subpos->arrd[ticktock - writelength + k] +
                    (double)maxpos->arr[ticktock - writelength + k];
            }
            syncAppendMyarr(syncarr, cfg.fpposition);
            freemyarr(syncarr);
        }

        calculateTotalFromFile(
            totalTickTock, cfg.fpposition, NN, cfg.SDthreshold);
        thread_unlock();
        wait();
        (void)fclose(cfg.fpposition);
    }
    freemyarr(subpos);
    freemyarr(reference);
    freemyarr(maxvals);
    freemyarr(maxpos);
    for (unsigned int t = 0; t < cfg.teeth; ++t)
    {
        struct myarr* cumulativeTick = teethArray[t];
        if (cfg.fptotal)
        {
            int toothshift = getshift(*teethArray[0], *cumulativeTick);
            for (unsigned int j = 0; j < NN; ++j)
            {
                (void)fprintf(cfg.fptotal,
                              "%d %d %u %d\n",
                              shiftHalf(j + toothshift, NN),
                              cumulativeTick->arr[j],
                              t,
                              toothshift);
            }
            (void)fprintf(cfg.fptotal, "\n\n");
        }
    }

    for (unsigned int t = 0; t < cfg.teeth; ++t)
    {
        freemyarr(teethArray[t]);
    }
    if (cfg.fpInput)
    {
        (void)fclose(cfg.fpInput);
    }
    if (cfg.fptotal)
    {
        (void)fclose(cfg.fptotal);
    }

    if (capture_handle != NULL)
    {
        snd_pcm_close(capture_handle);
        snd_pcm_hw_free(capture_handle);
    }

    (void)fprintf(stderr,
                  "width = %.3fms / %.1fμs/character\n",
                  mod * KILO / cfg.rate,
                  mod * MEGA / cfg.rate / (columns - cfg.everyline));
    return 0;
}
