#include "mysound.h"
#include "config.h"
#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
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

int derived(int* derivative, unsigned int ArrayLength, int16_t* samples);
// Helper to handle repetitive ALSA parameter setting and error reporting
static void check_alsa_err(int err,
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
    snd_pcm_t* capture_handle = NULL;
    snd_pcm_hw_params_t* hw_params = NULL;
    unsigned int requestedRate = *rate;

    // 1. Open Device
    check_alsa_err(
        snd_pcm_open(&capture_handle, device, SND_PCM_STREAM_CAPTURE, 0),
        device,
        "cannot open audio device",
        NULL,
        NULL);

    // 2. Allocate and Init Params
    check_alsa_err(snd_pcm_hw_params_malloc(&hw_params),
                   device,
                   "cannot allocate hardware parameter structure",
                   capture_handle,
                   NULL);

    check_alsa_err(snd_pcm_hw_params_any(capture_handle, hw_params),
                   device,
                   "cannot initialize hardware parameter structure",
                   capture_handle,
                   hw_params);

    // 3. Set Hardware Configurations
    check_alsa_err(snd_pcm_hw_params_set_access(capture_handle,
                                                hw_params,
                                                SND_PCM_ACCESS_RW_INTERLEAVED),
                   device,
                   "cannot set access type",
                   capture_handle,
                   hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_format(capture_handle, hw_params, format),
        device,
        "cannot set sample format",
        capture_handle,
        hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_rate_near(capture_handle, hw_params, rate, 0),
        device,
        "cannot set sample rate",
        capture_handle,
        hw_params);

    check_alsa_err(snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1),
                   device,
                   "cannot set channel count",
                   capture_handle,
                   hw_params);

    // 4. Apply Params and Prepare
    check_alsa_err(snd_pcm_hw_params(capture_handle, hw_params),
                   device,
                   "cannot set parameters",
                   capture_handle,
                   hw_params);

    snd_pcm_hw_params_free(hw_params);

    check_alsa_err(snd_pcm_prepare(capture_handle),
                   device,
                   "cannot prepare audio interface",
                   capture_handle,
                   NULL);

    // Minor logic check
    if (*rate != requestedRate)
    {
        (void)fprintf(stderr,
                      "Requested audiorate %u unavailable, using %u\n",
                      requestedRate,
                      *rate);
    }

    return capture_handle;
}

int derived(int* derivative, unsigned int ArrayLength, int16_t* samples)
{

    int clip_count = 0;

    for (unsigned int k = 0; k < ArrayLength - 1; k++)
    {
        if (samples[k] == INT16_MAX || samples[k] == INT16_MIN)
        {
            ++clip_count;
        }
        derivative[k] = abs(samples[k] - samples[k + 1]);
    }

    if (clip_count > 1)
    {
        (void)fprintf(stderr,
                      "%d audio 16-bit clipping event(s)\n",
                      clip_count);
    }
    derivative[ArrayLength - 1] = 0;
    return (int)ArrayLength;
}

int readBufferOrFile(int* derivative,
                     unsigned int ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx,
                     int16_t* samples)
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

                samples[index++] = (int16_t)val;
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
        ret = read_samples(ctx->cap, ArrayLength, samples);
        if (ret < 0)
        {
            return ret;
        }
        if ((unsigned)ret != ArrayLength)
        {
            return INPUT_SOUND_ERROR;
        }
    }

    return derived(derivative, ArrayLength, samples);
}

// Get data from audio capture
int getData(FILE* fpInput,
            struct myarr derivative,
            CaptureCtx* ctx,
            int16_t* out)
{
    int err = readBufferOrFile(derivative.arr,
                               derivative.ArrayLength,
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
 * capture_setup
 * - Accepts an already-opened/initialized ALSA handle (from initAudio).
 * - Computes ArrayLength = rate * SECS_HOUR * 2 / bph
 * - Queries period_size and builds poll descriptors.
 */
int capture_setup(CaptureCtx* ctx,
                  CapConfig* cfg,
                  unsigned int rate,
                  unsigned int bph)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->cap = cfg->capture_handle;
    ctx->rate = rate;

    // Geometry
    ctx->ArrayLength =
        rate * SECS_HOUR * 2 / bph; // e.g., ~16000 @ 48k/21600bph

    // Query ALSA buffer/period sizes (read per period)
    if (cfg->fpInput == 0 &&
        snd_pcm_get_params(ctx->cap, &ctx->buffer_size, &ctx->period_size) < 0)
    {
        (void)fprintf(
            stderr,
            "ALSA: snd_pcm_get_params failed; defaulting period_size=%d\n",
            DEFAULT_PERIOD);
        ctx->period_size = DEFAULT_PERIOD;
    }
    if (ctx->period_size == 0)
    {
        ctx->period_size = DEFAULT_PERIOD;
    }
    if (ctx->period_size > ctx->ArrayLength)
    {
        ctx->period_size = ctx->ArrayLength;
    }

    return 0;
}

void capture_teardown(CaptureCtx* ctx)
{
    if (!ctx)
    {
        return;
    }
    memset(ctx, 0, sizeof(*ctx));
}

int read_samples(snd_pcm_t* cap, unsigned int ArrayLength, int16_t* out)
{
    const unsigned TARGET = ArrayLength;
    unsigned collected = 0;

    while (collected < TARGET)
    {
        snd_pcm_sframes_t got =
            snd_pcm_readi(cap, out + collected, TARGET - collected);

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
