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
    int cumulativeShift;
    unsigned int tickIndex;
    unsigned int globalTickIndex;
} LoopState;

static int initAudioSource(CapConfig* cfg, unsigned int* actualRate)
{
    if (cfg->fpInput == NULL && *cfg->device == '\0')
    {
        printf("device and file are NULL\n");
        exit(EXIT_FAILURE);
    }

    *actualRate = (unsigned int)(cfg->rate + ROUNDING_HALF);
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

static AppResources allocateResources(unsigned int arrayLength,
                                      unsigned int ticktockBuffer,
                                      CapConfig* cfg)
{
    AppResources res = {0};
    res.subpos = makemyarrd(ticktockBuffer);
    res.maxpos = makemyarr(ticktockBuffer);
    res.maxvals = makemyarrd(ticktockBuffer);
    res.derivative = makemyarr(arrayLength);
    res.tmpder = makemyarr(arrayLength);
    res.reference = makemyarr(arrayLength);
    res.filterFFT = makeFilter(cfg->evalue, arrayLength);
    res.audioBuffer16 =
        calloc(arrayLength, sizeof(*res.audioBuffer16)); // 16bit
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
        res.teethArray[t] = makemyarr(arrayLength);
    }
    return res;
}

static void cleanupResources(AppResources* res, CapConfig* cfg, CaptureCtx* ctx)
{
    if (cfg->fpposition)
    {
        waitClose(cfg->fpposition);
    }
    if (cfg->fpInput)
    {
        (void)fclose(cfg->fpInput);
    }
    if (cfg->fpmaxcor)
    {
        waitClose(cfg->fpmaxcor);
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
                           unsigned int writeInterval)
{
    if (totalTime > 0 && totalTime % writeInterval == 0)
    {
        if (cfg->fpposition)
        {
            struct myarr* positionBatch = makemyarrd(writeInterval);
            for (unsigned int k = 0; k < writeInterval; ++k)
            {
                positionBatch->arrd[k] =
                    res->subpos->arrd[totalTime - writeInterval + k] +
                    (double)res->maxpos->arr[totalTime - writeInterval + k];
            }
            syncAppendMyarr(positionBatch, cfg->fpposition);
            freemyarr(positionBatch);
        }
        if (cfg->fpmaxcor)
        {
            struct myarr* correlationBatch = makemyarrd(writeInterval);
            memcpy(correlationBatch->arrd,
                   res->maxvals->arrd + totalTime - writeInterval,
                   writeInterval * sizeof(double));
            syncAppendMyarr(correlationBatch, cfg->fpmaxcor);
            freemyarr(correlationBatch);
        }
    }
}

static void fitAndPrint(unsigned int tickIndex,
                        unsigned int globalTickIndex,
                        struct myarr* cumulativeTick,
                        AppResources* res,
                        CapConfig* cfg,
                        unsigned int arrayLength,
                        unsigned int mod)
{

    double intercept = 0.0;
    double slope = 0.0;
    fitNpeaks(&intercept,
              &slope,
              tickIndex,
              res->maxvals,
              res->maxpos,
              res->subpos,
              cfg->fitN,
              cfg->SDthreshold);

    printheader(slope * SECS_DAY / arrayLength,
                cfg->everyline,
                getBeatError(cumulativeTick, cfg->rate, 0),
                (double)globalTickIndex * arrayLength / cfg->rate);

    printspaces(res->maxpos->arr[tickIndex],
                res->maxvals->arrd[tickIndex] * HEX_BASE,
                mod,
                columns - cfg->everyline,
                intercept,
                cfg->cvalue);
}

static void rotateDerivativeWindow(AppResources* res,
                                   unsigned int arrayLength,
                                   int cumulativeShift)
{
    for (int j = 0; j < (int)arrayLength; ++j)
    {
        res->tmpder->arr[j] =
            res->derivative->arr[modSigned(cumulativeShift + j, arrayLength)];
    }
}

static int findMaxPosition(AppResources* res,
                           struct myarr* cumulativeTick,
                           unsigned int globalTickIndex,
                           unsigned int tickIndex,
                           unsigned int arrayLength,
                           CapConfig* cfg)
{
    const int useReference = (globalTickIndex < AUTOCOR_LIMIT * cfg->teeth);
    return shiftHalf(
        fftfit(*res->tmpder,
               cumulativeTick->arr,
               useReference ? res->reference->arr : cumulativeTick->arr,
                   res->maxvals->arrd + tickIndex,
               res->filterFFT,
               globalTickIndex == cfg->verbose,
                   res->subpos->arrd + tickIndex),
        arrayLength);
}

static int updateTotalShiftIfNeeded(int cumulativeShift,
                                    int peakOffset,
                                    unsigned int globalTickIndex,
                                    unsigned int tickIndex,
                                    AppResources* res,
                                    CapConfig* cfg)
{
    if (globalTickIndex > AUTOCOR_LIMIT &&
        res->maxvals->arrd[tickIndex] > (double)cfg->cvalue / HEX_BASE &&
        globalTickIndex % cfg->teeth == 0)
    {
        int delta = peakOffset;
        if (abs(delta) > PRESHIFT_THRESHOLD)
        {
            delta = (int)(PRESHIFT_THRESHOLD_ROOT * delta / sqrt(abs(delta)));
        }
        cumulativeShift += delta;
    }
    return cumulativeShift;
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
    if (state->tickIndex == ARRAY_BUFFER_SIZE * 2)
    {
        shiftBufferData(&state->tickIndex,
                          res->subpos,
                          res->maxpos,
                          res->maxvals);
    }

    int readStatus =
        getData(cfg->fpInput, res->derivative, ctx, res->audioBuffer16);
    if (readStatus < 0)
    {
        return -1;
    }

    struct myarr* cumulativeTick =
        res->teethArray[state->globalTickIndex % cfg->teeth];
    rotateDerivativeWindow(res, params->arrayLength, state->cumulativeShift);
    int peakOffset = findMaxPosition(res,
                                      cumulativeTick,
                                      state->globalTickIndex,
                                      state->tickIndex,
                                      params->arrayLength,
                                      cfg);

    res->maxpos->arr[state->tickIndex] =
        state->cumulativeShift + peakOffset;
    state->cumulativeShift = updateTotalShiftIfNeeded(state->cumulativeShift,
                                                 peakOffset,
                                                 state->globalTickIndex,
                                                 state->tickIndex,
                                                 res,
                                                 cfg);

    processLogging(
        cfg,
        res,
        state->tickIndex,
        ARRAY_BUFFER_SIZE / DEFAULT_WRITE_FACTOR);

    fitAndPrint(state->tickIndex,
                state->globalTickIndex,
                cumulativeTick,
                res,
                cfg,
                params->arrayLength,
                params->mod);

    state->tickIndex++;
    state->globalTickIndex++;
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
        allocateResources(params.arrayLength, ARRAY_BUFFER_SIZE * 2, &cfg);
    fillReference(cfg.fpDefPeak, res.reference, cfg.teeth);

    sigset_t block;
    // sigset_t non_block;
    setupBlockSignals(&block);

    LoopState state = {0};

    while (keepRunning && !(state.globalTickIndex > params.maxTime && cfg.time))
    {
        if (processCaptureTick(&cfg, &res, &ctx, &params, &state) < 0)
        {
            break;
        }
    }

    printFinals(&cfg, &res, params.arrayLength, state.globalTickIndex);
    cleanupResources(&res, &cfg, &ctx);

    return 0;
}
