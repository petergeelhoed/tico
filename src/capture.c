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

typedef struct
{
    unsigned int arrayLength;
    unsigned int mod;
    unsigned int maxTime;
} RuntimeParams;

typedef struct
{
    int totalshift;
    unsigned int ticktock;
    unsigned int totalTickTock;
} LoopState;

static int initAudioSource(CapConfig* cfg, unsigned int* actualRate)
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
        cfg->captureHandle =
            initAudio(SND_PCM_FORMAT_S16_LE, cfg->device, actualRate);
        printf("Actual rate %d, calculating with %f\n", *actualRate, cfg->rate);
    }
    return (cfg->captureHandle == NULL && cfg->fpInput == NULL)
               ? ERROR_NO_SOURCE
               : 0;
}

static AppResources allocateResources(unsigned int ArrayLength,
                                      unsigned int ticktockBuffer,
                                      CapConfig* cfg)
{
    AppResources res = {0};
    res.subpos = makemyarrd(ticktockBuffer);
    res.maxpos = makemyarr(ticktockBuffer);
    res.maxvals = makemyarrd(ticktockBuffer);
    res.derivative = makemyarr(ArrayLength);
    res.tmpder = makemyarr(ArrayLength);
    res.reference = makemyarr(ArrayLength);
    res.filterFFT = makeFilter(cfg->evalue, ArrayLength);
    res.audioBuffer16 =
        calloc(ArrayLength, sizeof(*res.audioBuffer16)); // 16bit
    res.teethArray = calloc(cfg->teeth, sizeof(*res.teethArray));
    if (res.teethArray == NULL || res.audioBuffer16 == NULL)
    {
        free(res.audioBuffer16);
        free(res.teethArray);
        (void)fprintf(stderr, "Failed memory allocation\n");
        exit(EXIT_FAILURE);
    }

    for (unsigned int t = 0; t < cfg->teeth; t++)
    {
        res.teethArray[t] = makemyarr(ArrayLength);
    }
    return res;
}

static void cleanupResources(AppResources* res, CapConfig* cfg, CaptureCtx* ctx)
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
    if (cfg->captureHandle)
    {
        snd_pcm_close(cfg->captureHandle);
    }
    free(res->audioBuffer16);
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
    captureTeardown(ctx);
}

static void processLogging(CapConfig* cfg,
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

static void fitAndPrint(unsigned int ticktock,
                        unsigned int totalTickTock,
                        struct myarr* cumulativeTick,
                        AppResources* res,
                        CapConfig* cfg,
                        unsigned int ArrayLength,
                        unsigned int mod)
{

    double intercept = 0.0;
    double slope = 0.0;
    fitNpeaks(&intercept,
              &slope,
              ticktock,
              res->maxvals,
              res->maxpos,
              res->subpos,
              cfg->fitN,
              cfg->SDthreshold);

    printheader(slope * SECS_DAY / ArrayLength,
                cfg->everyline,
                getBeatError(cumulativeTick, cfg->rate, 0),
                (double)totalTickTock * ArrayLength / cfg->rate);

    printspaces(res->maxpos->arr[ticktock],
                res->maxvals->arrd[ticktock] * HEXDEC,
                mod,
                columns - cfg->everyline,
                intercept,
                cfg->cvalue);
}

static void rotateDerivativeWindow(AppResources* res,
                                   unsigned int arrayLength,
                                   int totalshift)
{
    for (int j = 0; j < (int)arrayLength; ++j)
    {
        res->tmpder->arr[j] =
            res->derivative->arr[modSigned(totalshift + j, arrayLength)];
    }
}

static int findMaxPosition(AppResources* res,
                           struct myarr* cumulativeTick,
                           unsigned int totalTickTock,
                           unsigned int ticktock,
                           unsigned int arrayLength,
                           CapConfig* cfg)
{
    const int useReference = (totalTickTock < AUTOCOR_LIMIT * cfg->teeth);
    return shiftHalf(
        fftfit(*res->tmpder,
               cumulativeTick->arr,
               useReference ? res->reference->arr : cumulativeTick->arr,
               res->maxvals->arrd + ticktock,
               res->filterFFT,
               totalTickTock == cfg->verbose,
               res->subpos->arrd + ticktock),
        arrayLength);
}

static int updateTotalShiftIfNeeded(int totalshift,
                                    int maxposition,
                                    unsigned int totalTickTock,
                                    unsigned int ticktock,
                                    AppResources* res,
                                    CapConfig* cfg)
{
    if (totalTickTock > AUTOCOR_LIMIT &&
        res->maxvals->arrd[ticktock] > (double)cfg->cvalue / HEXDEC &&
        totalTickTock % cfg->teeth == 0)
    {
        int delta = maxposition;
        if (abs(delta) > PRESHIFT_THRESHOLD)
        {
            delta = (int)(PRESHIFT_THRESHOLD_ROOT * delta / sqrt(abs(delta)));
        }
        totalshift += delta;
    }
    return totalshift;
}

static RuntimeParams computeRuntimeParams(const CapConfig* cfg,
                                          unsigned int actualRate)
{
    RuntimeParams params = {0};
    params.arrayLength = (actualRate * 2 * SECS_HOUR / cfg->bph);
    params.arrayLength += (params.arrayLength % 2);
    params.mod = params.arrayLength / cfg->zoom;
    params.maxTime = (unsigned int)cfg->rate *
                     (cfg->time ? cfg->time : DEFAULT_TIME) /
                     params.arrayLength;
    return params;
}

static int processCaptureTick(CapConfig* cfg,
                              AppResources* res,
                              CaptureCtx* ctx,
                              const RuntimeParams* params,
                              LoopState* state)
{
    if (state->ticktock == ARR_BUFF * 2)
    {
        shiftBufferData(&state->ticktock,
                          res->subpos,
                          res->maxpos,
                          res->maxvals);
    }

    int err = getData(cfg->fpInput, res->derivative, ctx, res->audioBuffer16);
    if (err < 0)
    {
        return -1;
    }

    struct myarr* cumulativeTick =
        res->teethArray[state->totalTickTock % cfg->teeth];
    rotateDerivativeWindow(res, params->arrayLength, state->totalshift);
    int maxposition = findMaxPosition(res,
                                      cumulativeTick,
                                      state->totalTickTock,
                                      state->ticktock,
                                      params->arrayLength,
                                      cfg);

    res->maxpos->arr[state->ticktock] = state->totalshift + maxposition;
    state->totalshift = updateTotalShiftIfNeeded(state->totalshift,
                                                 maxposition,
                                                 state->totalTickTock,
                                                 state->ticktock,
                                                 res,
                                                 cfg);

    processLogging(cfg, res, state->ticktock, ARR_BUFF / DEFAULT_WRITE_FACTOR);

    fitAndPrint(state->ticktock,
                state->totalTickTock,
                cumulativeTick,
                res,
                cfg,
                params->arrayLength,
                params->mod);

    state->ticktock++;
    state->totalTickTock++;
    return 0;
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
                     .fpInput = NULL,
                     .captureHandle = NULL};

    parseArguments(argc, argv, &cfg);

    struct winsize windowSize;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &windowSize);
    columns = windowSize.ws_col;
    setSignalAction();

    unsigned int actualRate;
    if (initAudioSource(&cfg, &actualRate))
    {
        return EXIT_FAILURE;
    }
    CaptureCtx ctx;
    if (captureSetup(&ctx, &cfg, actualRate) < 0)
    {
        (void)fprintf(stderr, "captureSetup failed\n");
        snd_pcm_close(cfg.captureHandle);
        return EXIT_FAILURE;
    }

    RuntimeParams params = computeRuntimeParams(&cfg, actualRate);
    AppResources res =
        allocateResources(params.arrayLength, ARR_BUFF * 2, &cfg);
    fillReference(cfg.fpDefPeak, res.reference, cfg.teeth);

    sigset_t block;
    // sigset_t non_block;
    setupBlockSignals(&block);

    LoopState state = {0};

    while (keepRunning && !(state.totalTickTock > params.maxTime && cfg.time))
    {
        if (processCaptureTick(&cfg, &res, &ctx, &params, &state) < 0)
        {
            break;
        }
    }

    printFinals(&cfg, &res, params.arrayLength, state.totalTickTock);
    cleanupResources(&res, &cfg, &ctx);

    return 0;
}
