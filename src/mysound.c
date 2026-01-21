#include <alsa/asoundlib.h>
#include <errno.h>
#include <fftw3.h>
#include <limits.h> // INT16_MAX/INT16_MIN
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "mylib.h"
#include "mysound.h"

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

void readBufferRaw(snd_pcm_t* capture_handle,
                   char* buffer,
                   struct myarr* data_in)
{
    unsigned char lsb;
    signed char msb;
    long err = snd_pcm_readi(capture_handle, buffer, data_in->ArrayLength);
    if (err != (long)data_in->ArrayLength)
    {
        (void)fprintf(stderr,
                      "read from audio interface failed %ld (%s)\n",
                      err,
                      snd_strerror((int)err));
        exit(READ_FAILED);
    }
    int overflow = 0;
    for (unsigned int index = 0; index < 2 * data_in->ArrayLength; index += 2)
    {
        msb = (signed char)buffer[index + 1];
        lsb = *(buffer + index);
        data_in->arr[index / 2] = (msb << BITS_IN_BYTE) | lsb;
        overflow += (data_in->arr[index / 2] == SHRT_MAX);
        overflow += (data_in->arr[index / 2] == SHRT_MIN);
    }
    if (overflow > 1)
    {
        (void)fprintf(stderr, "%d audio 16bit overflows\n", overflow);
    }
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative)
{
    // --- Make stderr unbuffered once so prints are visible even if redirected
    // ---
    static int stderr_unbuffered = 0;
    if (!stderr_unbuffered)
    {
        setvbuf(stderr, NULL, _IONBF, 0);
        stderr_unbuffered = 1;
    }

    if (ArrayLength == 0)
    {
        return 0;
    }

    // --- Format assumptions: S16_LE mono ---
    const unsigned int BYTES_PER_SAMPLE = 2; // 16-bit
    const unsigned int CHANNELS = 1;         // mono
    const unsigned int BYTES_PER_FRAME = BYTES_PER_SAMPLE * CHANNELS;

    // --- Read exactly ArrayLength frames, handling xruns and short reads ---
    snd_pcm_uframes_t remaining = (snd_pcm_uframes_t)ArrayLength;
    char* write_ptr = buffer;
    while (remaining > 0)
    {
        snd_pcm_sframes_t r =
            snd_pcm_readi(capture_handle, write_ptr, remaining);

        if (r == -EAGAIN)
        {
            // Non-blocking mode: try again
            continue;
        }
        else if (r == -EPIPE || r == -ESTRPIPE)
        {
            // Overrun (capture) or suspend: recover and continue
            fprintf(stderr,
                    "ALSA: overrun/suspend detected (%s)\n",
                    snd_strerror((int)r));
            int rr = snd_pcm_recover(capture_handle, (int)r, /*silent=*/1);
            if (rr < 0)
            {
                fprintf(stderr, "ALSA: recover failed: %s\n", snd_strerror(rr));
                return rr;
            }
            // After recover, try reading again
            continue;
        }
        else if (r < 0)
        {
            // Other fatal error
            fprintf(stderr,
                    "ALSA: snd_pcm_readi failed: %s\n",
                    snd_strerror((int)r));
            return (int)r;
        }

        // Successful short/long read
        write_ptr += (size_t)r * BYTES_PER_FRAME;
        remaining -= (snd_pcm_uframes_t)r;
    }

    // --- Post-read ALSA status snapshot: XRUN/overrange reporting ---
    {
        snd_pcm_status_t* status;
        snd_pcm_status_alloca(&status);
        if (snd_pcm_status(capture_handle, status) == 0)
        {
            if (snd_pcm_status_get_state(status) == SND_PCM_STATE_XRUN)
            {
                fprintf(stderr,
                        "ALSA: XRUN detected in status (post-read). "
                        "Recovering...\n");
                (void)snd_pcm_prepare(capture_handle);
            }
            unsigned long over = snd_pcm_status_get_overrange(status);
            if (over > 0)
            {
                fprintf(stderr,
                        "ALSA: ADC overrange reported %lu time(s)\n",
                        over);
            }
        }
    }

    // --- Timestamp-based discontinuity detection (no signature change) ---
    // We query the current hardware rate/period once and then use
    // snd_pcm_htimestamp() to check for unexpected gaps.
    {
        static int rate_inited = 0;
        static unsigned int sample_rate = 0;
        static int period_inited = 0;
        static snd_pcm_uframes_t period_size = 0;
        static int have_prev_ts = 0;
        static snd_htimestamp_t prev_ts = {0};

        if (!rate_inited || !period_inited)
        {
            snd_pcm_hw_params_t* hw;
            snd_pcm_hw_params_alloca(&hw);
            if (snd_pcm_hw_params_current(capture_handle, hw) == 0)
            {
                int dir = 0;
                (void)snd_pcm_hw_params_get_rate(hw, &sample_rate, &dir);
                (void)snd_pcm_hw_params_get_period_size(hw, &period_size, &dir);
            }
            rate_inited = period_inited = 1;
        }

        if (sample_rate > 0)
        {
            snd_htimestamp_t ts;
            snd_pcm_uframes_t avail_hw = 0;
            if (snd_pcm_htimestamp(capture_handle, &avail_hw, &ts) == 0)
            {
                if (have_prev_ts)
                {
                    double dt = (double)(ts.tv_sec - prev_ts.tv_sec) +
                                (double)(ts.tv_nsec - prev_ts.tv_nsec) * 1e-9;
                    double expected_frames = dt * (double)sample_rate;
                    double actual_frames = (double)ArrayLength;

                    // Tolerance: one ALSA period if known, else ~1 ms of audio
                    double tol_frames = (period_size > 0)
                                            ? (double)period_size
                                            : (double)sample_rate * 0.001;

                    double gap = expected_frames - actual_frames;
                    if (gap > tol_frames)
                    {
                        fprintf(stderr,
                                "ALSA: discontinuity suspected: ts advanced by "
                                "~%.0f frames, "
                                "but read %.0f (gap ~%.0f, tol ~%.0f)\n",
                                expected_frames,
                                actual_frames,
                                gap,
                                tol_frames);
                    }
                }
                prev_ts = ts;
                have_prev_ts = 1;
            }
        }
    }

    // --- Convert bytes -> int16 samples; detect clipping; compute derivative
    // ---
    int clip_overflow = 0;

    // First store samples into derivative[] (int) then reuse for |x[n]-x[n+1]|
    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        const uint8_t lsb = (uint8_t)buffer[2 * i + 0];
        const int8_t msb = (int8_t)buffer[2 * i + 1];
        const int16_t sample = (int16_t)(((int)msb << 8) | lsb); // S16_LE

        derivative[i] = (int)sample;

        // Clipping check (ADC full-scale, not XRUN)
        if (sample == INT16_MAX || sample == INT16_MIN)
        {
            clip_overflow++;
        }
    }

    if (clip_overflow > 1)
    {
        fprintf(stderr, "%d audio 16-bit clipping event(s)\n", clip_overflow);
    }

    // Compute derivative: |x[n] - x[n+1]|; last sample = 0
    for (unsigned int i = 0; i + 1 < ArrayLength; ++i)
    {
        int d = derivative[i] - derivative[i + 1];
        derivative[i] = (d < 0) ? -d : d;
    }
    derivative[ArrayLength - 1] = 0;

    // Success: exactly ArrayLength frames were read
    return (int)ArrayLength;
}

int readBufferOrFile(int* derivative,
                     snd_pcm_t* capture_handle,
                     unsigned int ArrayLength,
                     char* buffer,
                     FILE* fpInput)
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
                    return INPUT_OVERFLOW;
                }

                derivative[index++] = (int)val;
                ptr = endptr; // Advance to the rest of the string
            }
        }

        free(line); // getline allocates memory that must be freed

        if (index < ArrayLength)
        {
            return INPUT_FILE_ERROR;
        }
        ret = (int)ArrayLength;

        // Perform derivative calculation
        for (unsigned int k = 0; k < ArrayLength - 1; k++)
        {
            derivative[k] = abs(derivative[k] - derivative[k + 1]);
        }
        derivative[ArrayLength - 1] = 0;
    }
    else
    {
        ret = readBuffer(capture_handle, ArrayLength, buffer, derivative);
    }
    return ret;
}

// Get data from audio capture
int getData(FILE* rawfile,
            FILE* fpInput,
            snd_pcm_t* capture_handle,
            snd_pcm_format_t format,
            char* device,
            unsigned int rate,
            char* buffer,
            struct myarr derivative)
{
    int err = REINIT_ERROR;
    while (err == REINIT_ERROR)
    {
        err = readBufferOrFile(derivative.arr,
                               capture_handle,
                               derivative.ArrayLength,
                               buffer,
                               fpInput);
        if (err == REINIT_ERROR)
        {
            (void)fprintf(stderr, "Reinitializing capture_handle");
            if (rawfile)
            {
                (void)fprintf(rawfile, "# Reinitializing capture_handle");
            }
            snd_pcm_close(capture_handle);
            capture_handle = initAudio(format, device, &rate);
            err = readBufferOrFile(derivative.arr,
                                   capture_handle,
                                   derivative.ArrayLength,
                                   buffer,
                                   fpInput);
        }
        if (err == INPUT_FILE_ERROR)
        {
            (void)fprintf(stderr,
                          "Could not read integer from inputfile or audio\n");
        }
    }
    return err;
}
