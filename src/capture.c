#include "analysis.h"
#include "capture_helpers.h"
#include "config.h"
#include "myarr.h"
#include "myfft.h"
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
    size_t tickIndex;
    unsigned int globalTickIndex;
} LoopState;

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

    if (getData(cfg->fpInput, res->derivative, ctx, res->audioBuffer16) < 0)
    {
        return -1;
    }

    struct myarr* cumulativeTick =
        res->teethArray[state->globalTickIndex % cfg->teeth];
    rotateDerivativeWindow(res, params->arrayLength, state->cumulativeShift);
    int peakOffset = findMaxPosition(res,
                                     cumulativeTick,
                                     state->globalTickIndex,
                                     (unsigned int)state->tickIndex,
                                     params->arrayLength,
                                     cfg);

    res->maxpos->arr[state->tickIndex] = state->cumulativeShift + peakOffset;
    state->cumulativeShift = updateTotalShiftIfNeeded(state->cumulativeShift,
                                                      peakOffset,
                                                      state->globalTickIndex,
                                                      state->tickIndex,
                                                      res,
                                                      cfg);

    processLogging(cfg,
                   res,
                   state->tickIndex,
                   ARRAY_BUFFER_SIZE / DEFAULT_WRITE_FACTOR);

    fitAndPrint(state->tickIndex,
                state->globalTickIndex,
                cumulativeTick,
                res,
                cfg,
                params->arrayLength,
                params->mod,
                columns);

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
