#include "config.h"
#include "myarr.h"
#include "myfft.h"
#include "mylib.h"
#include "mymath.h"
#include "mysignal.h"
#include "mysound.h"
#include "mysync.h"
#include "parseargs.h"
#include "resources.h"

#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>

volatile int keepRunning = 1;
volatile unsigned int columns = DEFAULT_COLUMNS;
static int init_audio_source(CapConfig* cfg,
                             snd_pcm_t** handle,
                             unsigned int* actualRate)
{
    if (cfg->fpInput == NULL && *cfg->device == '\0')
    {
        printf("device and file are NULL\n");
        exit(EXIT_FAILURE);
    }

    *actualRate = (unsigned int)(cfg->rate + HALF);
    if (cfg->fpInput == NULL)
    {
        printf("Casting inputrate %f to soundcard(%s) rate %d\n",
               cfg->rate,
               cfg->device,
               *actualRate);
        *handle = initAudio(SND_PCM_FORMAT_S16_LE, cfg->device, actualRate);
        printf("Actual rate %d, calculating with %f\n", *actualRate, cfg->rate);
    }
    return (*handle == NULL && cfg->fpInput == NULL) ? ERROR_NO_SOURCE : 0;
}

static AppResources allocate_resources(unsigned int ArrayLength,
                                       unsigned int ticktockBuffer,
                                       unsigned int teeth,
                                       unsigned int evalue)
{
    AppResources res = {0};
    res.subpos = makemyarrd(ticktockBuffer);
    res.maxpos = makemyarr(ticktockBuffer);
    res.maxvals = makemyarrd(ticktockBuffer);
    res.derivative = makemyarr(ArrayLength);
    res.tmpder = makemyarr(ArrayLength);
    res.reference = makemyarr(ArrayLength);
    res.filterFFT = makeFilter(evalue, ArrayLength);
    res.audioBuffer =
        calloc(ArrayLength, 2 * sizeof(*res.audioBuffer)); // 16bit
    res.teethArray = calloc(teeth, sizeof(*res.teethArray));
    if (res.audioBuffer == NULL || res.teethArray == NULL)
    {
        free(res.audioBuffer);
        free(res.teethArray);
        (void)fprintf(stderr, "Failed memory allocation\n");
        exit(EXIT_FAILURE);
    }

    for (unsigned int t = 0; t < teeth; t++)
    {
        res.teethArray[t] = makemyarr(ArrayLength);
    }
    return res;
}

static void cleanup_resources(AppResources* res, CapConfig* cfg)
{
    if (cfg->fpposition)
    {
        (void)fclose(cfg->fpposition);
    }
    if (cfg->fpInput)
    {
        (void)fclose(cfg->fpInput);
    }
    if (cfg->fpmaxcor)
    {
        (void)fclose(cfg->fpmaxcor);
    }
    if (cfg->fptotal)
    {
        (void)fclose(cfg->fptotal);
    }
    free(res->audioBuffer);
    freemyarr(res->subpos);
    freemyarr(res->maxpos);
    freemyarr(res->maxvals);
    freemyarr(res->derivative);
    freemyarr(res->tmpder);
    freemyarr(res->reference);
    fftw_free(res->filterFFT);
    for (unsigned int t = 0; t < cfg->teeth; t++)
    {
        freemyarr(res->teethArray[t]);
    }
    free(res->teethArray);
}

static void process_logging(CapConfig* cfg,
                            AppResources* res,
                            unsigned int totalTime,
                            unsigned int writeinterval)
{
    if (totalTime > 0 && totalTime % writeinterval == 0)
    {
        if (cfg->fpposition)
        {
            struct myarr* sync = makemyarrd(writeinterval);
            for (unsigned int k = 0; k < writeinterval; ++k)
            {
                sync->arrd[k] =
                    res->subpos->arrd[totalTime - writeinterval + k] +
                    (double)res->maxpos->arr[totalTime - writeinterval + k];
            }
            syncAppendMyarr(sync, cfg->fpposition);
            freemyarr(sync);
        }
        if (cfg->fpmaxcor)
        {
            struct myarr* sync = makemyarrd(writeinterval);
            memcpy(sync->arrd,
                   res->maxvals->arrd + totalTime - writeinterval,
                   writeinterval * sizeof(double));
            syncAppendMyarr(sync, cfg->fpmaxcor);
            freemyarr(sync);
        }
    }
}

int main(int argc, char* argv[])
{
    CapConfig cfg = {.rate = DEFAULT_RATE,
                     .bph = DEFAULT_BPH,
                     .evalue = DEFAULT_EVALUE,
                     .zoom = DEFAULT_ZOOM,
                     .fitN = DEFAULT_FITN,
                     .teeth = DEFAULT_TEETH,
                     .SDthreshold = DEFAULT_SDTHRESHOLD,
                     .device = "default:2",
                     .cvalue = DEFAULT_CVALUE,
                     .fpposition = NULL,
                     .fpmaxcor = NULL,
                     .fptotal = NULL,
                     .fpDefPeak = NULL,
                     .fpInput = NULL};

    parse_arguments(argc, argv, &cfg);

    struct winsize windowSize;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &windowSize);
    columns = windowSize.ws_col;
    set_signal_action();

    snd_pcm_t* capture_handle = NULL;
    unsigned int actualRate;
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    if (init_audio_source(&cfg, &capture_handle, &actualRate))
    {
        return EXIT_FAILURE;
    }
    CaptureCtx ctx;
    unsigned int countdown = 0; // show countdown concurrently with capture
    if (capture_setup(&ctx,
                      capture_handle,
                      actualRate,
                      cfg.bph,
                      countdown,
                      format) < 0)
    {
        (void)fprintf(stderr, "capture_setup failed\n");
        snd_pcm_close(capture_handle);
        return EXIT_FAILURE;
    }

    unsigned int ArrayLength = (actualRate * 2 * SECS_HOUR / cfg.bph);
    ArrayLength += (ArrayLength % 2);
    unsigned int mod = ArrayLength / cfg.zoom;
    const unsigned int maxtime = (unsigned int)cfg.rate *
                                 (cfg.time ? cfg.time : DEFAULT_TIME) /
                                 ArrayLength;

    AppResources res =
        allocate_resources(ArrayLength, ARR_BUFF * 2, cfg.teeth, cfg.evalue);
    fillReference(cfg.fpDefPeak, res.reference, cfg.teeth);

    sigset_t block;
    // sigset_t non_block;
    setup_block_signals(&block);

    int totalshift = 0;
    int toothshift = 0;
    unsigned int ticktock = 0;
    unsigned int totalTickTock = 0;

    while (keepRunning && !(totalTickTock > maxtime && cfg.time))
    {
        if (ticktock == ARR_BUFF * 2)
        {
            shift_buffer_data(&ticktock, res.subpos, res.maxpos, res.maxvals);
        }
        //  block_signal(&block, &non_block);
        int err = getData(cfg.fpInput, *res.derivative, &ctx);
        //  unblock_signal(&non_block);
        if (err < 0)
        {
            break;
        }
        struct myarr* cumulativeTick =
            res.teethArray[totalTickTock % cfg.teeth];
        if (cfg.teeth > 1 && totalTickTock >= AUTOCOR_LIMIT * cfg.teeth)
        {
            toothshift = getshift(*res.teethArray[0],
                                  *res.teethArray[totalTickTock % cfg.teeth]);
        }
        for (int j = 0; j < (int)ArrayLength; ++j)
        {
            res.tmpder->arr[j] =
                res.derivative
                    ->arr[modSigned(totalshift + j + toothshift, ArrayLength)];
        }
        int maxposition =
            shiftHalf(fftfit(*res.tmpder,
                             cumulativeTick->arr,
                             (totalTickTock < AUTOCOR_LIMIT * cfg.teeth)
                                 ? res.reference->arr
                                 : cumulativeTick->arr,
                             res.maxvals->arrd + ticktock,
                             res.filterFFT,
                             totalTickTock == cfg.verbose,
                             res.subpos->arrd + ticktock),
                      ArrayLength);

        res.maxpos->arr[ticktock] = totalshift + maxposition;
        if (totalTickTock > AUTOCOR_LIMIT &&
            res.maxvals->arrd[ticktock] > (double)cfg.cvalue / HEXDEC &&
            totalTickTock % cfg.teeth == 0)
        {
            if (abs(maxposition) > PRESHIFT_THRESHOLD)
            {
                maxposition = (int)(PRESHIFT_THRESHOLD_ROOT * maxposition /
                                    sqrt(abs(maxposition)));
            }
            totalshift += maxposition;
        }

        process_logging(&cfg, &res, ticktock, DEFAULT_TICKTOCK_WRITE);

        double intercept = 0.0;
        double slope = 0.0;
        fitNpeaks(&intercept,
                  &slope,
                  ticktock,
                  res.maxvals,
                  res.maxpos,
                  res.subpos,
                  cfg.fitN,
                  cfg.SDthreshold);

        printheader(slope * SECS_DAY / ArrayLength,
                    cfg.everyline,
                    getBeatError(cumulativeTick, cfg.rate, 0),
                    (double)totalTickTock * ArrayLength / cfg.rate);

        printspaces(res.maxpos->arr[ticktock],
                    res.maxvals->arrd[ticktock] * HEXDEC,
                    mod,
                    columns - cfg.everyline,
                    intercept,
                    cfg.cvalue);
        ticktock++;
        totalTickTock++;
    }

    print_finals(&cfg, &res, ArrayLength, totalTickTock, toothshift);
    cleanup_resources(&res, &cfg);
    if (capture_handle)
    {
        snd_pcm_close(capture_handle);
    }
    return 0;
}
