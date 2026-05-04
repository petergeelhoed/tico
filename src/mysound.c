
#include <stdio.h>
#include <stdlib.h>
#include "mysound.h"


// Print all ALSA logical capture devices (suggested input devices)
void get_suggested_device(void) {
    printf("[Suggested ALSA Input Devices] (from arecord -L):\n");
    FILE *fa = popen("arecord -L", "r");
    if (!fa) {
        perror("arecord -L");
        return;
    }
    char line[256];
    while (fgets(line, sizeof(line), fa)) {
        // Only print lines that look like device names (not indented)
        if (line[0] != '\t' && line[0] != '\n') {
            printf("  %s", line);
        }
    }
    pclose(fa);
}
#include "mysound.h"
#include "config.h"
#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mysync.h"
#include "parseargs.h"

#include <alsa/asoundlib.h>
#include <errno.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <poll.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // strlen, strncpy
#include <sys/timerfd.h>
#include <time.h>
#include <unistd.h> // getopt, read

static int derived(int* derivative, size_t ArrayLength, int16_t* samples)
{

    int clipCount = 0;

    for (size_t k = 0; k < ArrayLength - 1; k++)
    {
        if (samples[k] == INT16_MAX || samples[k] == INT16_MIN)
        {
            ++clipCount;
        }
        derivative[k] = abs(samples[k] - samples[k + 1]);
    }

    if (clipCount > 1)
    {
        (void)fprintf(stderr, "%d audio 16-bit clipping event(s)\n", clipCount);
    }
    derivative[ArrayLength - 1] = 0;
    return (int)ArrayLength;
}

// Helper to handle repetitive ALSA parameter setting and error reporting
static void checkAlsaErr(int err,
                         const char* device,
                         const char* msg,
                         snd_pcm_t* handle,
                         snd_pcm_hw_params_t* params)
{
    if (err < 0)
    {
        (void)fprintf(stderr,
                      "Device %s: %s (%s)\n",
                      device,
                      msg,
                      snd_strerror(err));
        if (params)
        {
            snd_pcm_hw_params_free(params);
        }
        if (handle)
        {
            snd_pcm_close(handle);
        }
        exit(INIT_ERROR);
    }
}

// Initialize audio capture
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate)
{
    snd_pcm_t* captureHandle = NULL;
    snd_pcm_hw_params_t* hw_params = NULL;
    unsigned int requestedRate = *rate;

    // 1. Open Device
    checkAlsaErr(
        snd_pcm_open(&captureHandle, device, SND_PCM_STREAM_CAPTURE, 0),
        device,
        "cannot open audio device",
        NULL,
        NULL);

    // 2. Allocate and Init Params
    checkAlsaErr(snd_pcm_hw_params_malloc(&hw_params),
                 device,
                 "cannot allocate hardware parameter structure",
                 captureHandle,
                 NULL);

    checkAlsaErr(snd_pcm_hw_params_any(captureHandle, hw_params),
                 device,
                 "cannot initialize hardware parameter structure",
                 captureHandle,
                 hw_params);

    // 3. Set Hardware Configurations
    checkAlsaErr(snd_pcm_hw_params_set_access(captureHandle,
                                              hw_params,
                                              SND_PCM_ACCESS_RW_INTERLEAVED),
                 device,
                 "cannot set access type",
                 captureHandle,
                 hw_params);

    checkAlsaErr(snd_pcm_hw_params_set_format(captureHandle, hw_params, format),
                 device,
                 "cannot set sample format",
                 captureHandle,
                 hw_params);

    checkAlsaErr(
        snd_pcm_hw_params_set_rate_near(captureHandle, hw_params, rate, 0),
        device,
        "cannot set sample rate",
        captureHandle,
        hw_params);

    checkAlsaErr(snd_pcm_hw_params_set_channels(captureHandle, hw_params, 1),
                 device,
                 "cannot set channel count",
                 captureHandle,
                 hw_params);

    // 4. Apply Params and Prepare
    checkAlsaErr(snd_pcm_hw_params(captureHandle, hw_params),
                 device,
                 "cannot set parameters",
                 captureHandle,
                 hw_params);

    snd_pcm_hw_params_free(hw_params);

    checkAlsaErr(snd_pcm_prepare(captureHandle),
                 device,
                 "cannot prepare audio interface",
                 captureHandle,
                 NULL);

    // Minor logic check
    if (*rate != requestedRate)
    {
        (void)fprintf(stderr,
                      "Requested audiorate %u unavailable, using %u\n",
                      requestedRate,
                      *rate);
    }

    return captureHandle;
}

int initAudioSource(CapConfig* cfg, unsigned int* actualRate)
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

int readBufferOrFile(int* derivative,
                     size_t ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx,
                     int16_t* buffer16)
{
    int ret = READ_FAILED;

    if (fpInput)
    {
        char* line = NULL;
        size_t len = 0;
        unsigned int index = 0;

        // Read entire lines until we have ArrayLength numbers
        while (index < ArrayLength && getline(&line, &len, fpInput) != -1)
        {
            char* ptr = line;
            char* endptr;

            while (index < ArrayLength)
            {
                errno = 0;
                // Use strtol for %d equivalent; use strtod for floating point
                long val = strtol(ptr, &endptr, DECIMAL);

                // If ptr == endptr, no more numbers were found on this line
                if (ptr == endptr)
                {
                    break;
                }

                // Error checking (optional but recommended)
                if (errno == ERANGE)
                { /* Handle overflow */
                    free(line);
                    return INPUT_OVERFLOW;
                }

                buffer16[index++] = (int16_t)val;
                ptr = endptr; // Advance to the rest of the string
            }
        }

        free(line); // getline allocates memory that must be freed

        if (index < ArrayLength)
        {
            return INPUT_FILE_ERROR;
        }
    }
    else
    {
        ret = readSamples(ctx->cap, ArrayLength, buffer16);
        if (ret < 0)
        {
            return ret;
        }
        if ((unsigned)ret != ArrayLength)
        {
            return INPUT_SOUND_ERROR;
        }
    }

    return derived(derivative, ArrayLength, buffer16);
}

// Get data from audio capture
int getData(FILE* fpInput,
            struct myarr* derivative,
            CaptureCtx* ctx,
            int16_t* out)
{
    int err = readBufferOrFile(derivative->arr,
                               derivative->ArrayLength,
                               fpInput,
                               ctx,
                               out);
    if (err == INPUT_FILE_ERROR)
    {
        (void)fprintf(stderr,
                      "Could not read integer from inputfile or audio\n");
    }
    return err;
}

/* -------------------- API: setup / next_block / teardown --------------------
 */

/**
 * captureSetup
 * - Accepts an already-opened/initialized ALSA handle (from initAudio).
 * - Computes ArrayLength = rate * SECS_HOUR * 2 / bph
 * - Queries periodSize and builds poll descriptors.
 */
int captureSetup(CaptureCtx* ctx, CapConfig* cfg, unsigned int rate)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->cap = cfg->captureHandle;
    ctx->rate = rate;

    // Geometry
    ctx->ArrayLength =
        rate * SECS_HOUR * 2 / cfg->bph; // e.g., ~16000 @ 48k/21600bph

    // Query ALSA buffer/period sizes (read per period)
    if (cfg->fpInput == 0 &&
        snd_pcm_get_params(ctx->cap, &ctx->bufferSize, &ctx->periodSize) < 0)
    {
        (void)fprintf(
            stderr,
            "ALSA: snd_pcm_get_params failed; defaulting periodSize=%d\n",
            DEFAULT_PERIOD);
        ctx->periodSize = DEFAULT_PERIOD;
    }
    if (ctx->periodSize == 0)
    {
        ctx->periodSize = DEFAULT_PERIOD;
    }
    if (ctx->periodSize > ctx->ArrayLength)
    {
        ctx->periodSize = ctx->ArrayLength;
    }

    return 0;
}

void captureTeardown(CaptureCtx* ctx)
{
    if (!ctx)
    {
        return;
    }
    memset(ctx, 0, sizeof(*ctx));
}

int readSamples(snd_pcm_t* cap, size_t ArrayLength, int16_t* out)
{
    const size_t TARGET = ArrayLength;
    size_t collected = 0;

    while (collected < TARGET)
    {
        snd_pcm_sframes_t got =
            snd_pcm_readi(cap,
                          out + collected,
                          (snd_pcm_uframes_t)(TARGET - collected));

        if (got == -EAGAIN)
        {
            /* No data available right now */
            /* Yield or sleep very briefly */
            const unsigned int SLEEP_US = 1000;
            usleep(SLEEP_US); /* 1 ms */
            continue;
        }

        if (got == -EPIPE || got == -ESTRPIPE)
        {
            /* XRUN or suspend */
            if (snd_pcm_recover(cap, (int)got, 1) < 0)
            {
                return -1;
            }
            continue;
        }

        if (got < 0)
        {
            (void)fprintf(stderr,
                          "ALSA read failed: %s\n",
                          snd_strerror((int)got));
            return -1;
        }

        collected += (unsigned)got;
    }
    return (int)collected;
}
