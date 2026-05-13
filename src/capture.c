
#include "printing.h"

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

/**
 * @brief Processes a single ticktock of audio capture and analysis.
 *
 * @param cfg Pointer to the configuration structure containing user-defined
 * settings.
 * @param res Pointer to the application resources structure for managing
 * buffers and state.
 * @param ctx Pointer to the capture context structure for managing audio
 * @param state Pointer to the loop state structure for tracking cumulative
 * shifts and tick indices.
 * @return 0 on success, -1 on failure (e.g., if audio data cannot be read).
 */
static int processTickTock(CapConfig* cfg,
                           AppResources* res,
                           CaptureCtx* ctx,
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
    rotateDerivativeWindow(res, state->cumulativeShift);
    int peakOffset = findMaxPosition(res, cumulativeTick, state, cfg);

    res->maxpos->arr[state->tickIndex] = state->cumulativeShift + peakOffset;
    updateTotalShiftIfNeeded(state, peakOffset, res, cfg);

    processLogging(cfg,
                   res,
                   state->tickIndex,
                   ARRAY_BUFFER_SIZE / DEFAULT_WRITE_FACTOR);

    fitAndPrint(state, cumulativeTick, res, cfg, columns);

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
                     .device = "",
                     .cvalue = DEFAULT_CVALUE,
                     .fpposition = NULL,
                     .fpmaxcor = NULL,
                     .fptotal = NULL,
                     .fpDefPeak = NULL,
                     .fpInput = NULL,
                     .captureHandle = NULL};

    parseArguments(argc, argv, &cfg);

    if (strlen(cfg.device) == 0)
    {
        const char* defaultDevice = get_default_device();
        if (defaultDevice)
        {
            strncpy(cfg.device, defaultDevice, sizeof(cfg.device) - 1);
            cfg.device[sizeof(cfg.device) - 1] = '\0';
        }
        else
        {
            print("No default audio device found.\n");
            return EXIT_FAILURE;
        }
    }

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
        print("captureSetup failed\n");
        snd_pcm_close(cfg.captureHandle);
        return EXIT_FAILURE;
    }

    unsigned int arrayLength = (actualRate * 2 * SECS_HOUR / cfg.bph);
    arrayLength += (arrayLength % 2);
    AppResources res =
        allocateResources(arrayLength, ARRAY_BUFFER_SIZE * 2, &cfg);
    fillReference(cfg.fpDefPeak, res.reference, cfg.teeth);

    sigset_t block;
    // sigset_t non_block;
    setupBlockSignals(&block);

    LoopState state = {0};

    while (keepRunning && !(state.globalTickIndex > res.maxTime && cfg.time))
    {
        if (processTickTock(&cfg, &res, &ctx, &state) < 0)
        {
            break;
        }
    }

    printFinals(&cfg, &res, state.globalTickIndex);
    cleanupResources(&res, &cfg, &ctx);

    return 0;
}
