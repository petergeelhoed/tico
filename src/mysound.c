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

// Complete function with a small helper to read exactly N frames.

// Helper: read exactly `frames` frames into `buf` or recover from xruns.
// Returns 0 on success, negative ALSA error on failure.
static int read_exact_frames(snd_pcm_t* h,
                             char* buf,
                             snd_pcm_uframes_t frames,
                             unsigned bytes_per_frame)
{
    snd_pcm_uframes_t total = 0;

    while (total < frames)
    {
        snd_pcm_sframes_t r =
            snd_pcm_readi(h, buf + (total * bytes_per_frame), frames - total);

        if (r == -EAGAIN)
        {
            // Non-blocking mode: try again.
            continue;
        }
        else if (r == -EPIPE || r == -ESTRPIPE)
        {
            // Overrun or stream suspended: attempt recovery.
            int rr = snd_pcm_recover(h, (int)r, /*silent=*/1);
            if (rr < 0)
            {
                fprintf(stderr,
                        "snd_pcm_recover failed: %s\n",
                        snd_strerror(rr));
                return rr;
            }
            // After successful recover, try again.
            continue;
        }
        else if (r < 0)
        {
            // Other error
            fprintf(stderr, "snd_pcm_readi failed: %s\n", snd_strerror((int)r));
            return (int)r;
        }

        total += (snd_pcm_uframes_t)r;
    }

    return 0;
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength, // frames
               char* buffer, // raw PCM bytes: must hold ArrayLength * 2 bytes
                             // (S16_LE mono)
               int* derivative) // length >= ArrayLength
{
    if (ArrayLength == 0)
    {
        return 0;
    }

    // --- Format assumptions: S16_LE, mono ---
    const unsigned int BYTES_PER_SAMPLE = 2; // 16-bit
    const unsigned int CHANNELS = 1;         // mono
    const unsigned int BYTES_PER_FRAME = BYTES_PER_SAMPLE * CHANNELS;

    // 1) Read exactly ArrayLength frames (blocking until filled or error)
    int rc = read_exact_frames(capture_handle,
                               buffer,
                               (snd_pcm_uframes_t)ArrayLength,
                               BYTES_PER_FRAME);
    if (rc < 0)
    {
        // Already logged; return ALSA error code.
        return rc;
    }

    // 2) Optional: Check ALSA status (XRUN state, overrange)
    {
        snd_pcm_status_t* status;
        snd_pcm_status_alloca(&status);
        if (snd_pcm_status(capture_handle, status) == 0)
        {
            snd_pcm_state_t st = snd_pcm_status_get_state(status);
            if (st == SND_PCM_STATE_XRUN)
            {
                fprintf(stderr,
                        "ALSA: XRUN detected (post-read); recovering...\n");
                (void)snd_pcm_prepare(capture_handle);
            }

            // Some drivers/plugins report ADC overrange count here
            unsigned long over = snd_pcm_status_get_overrange(status);
            if (over > 0)
            {
                fprintf(stderr,
                        "ALSA: ADC overrange reported %lu time(s)\n",
                        over);
            }
        }
    }

    // 3) Convert bytes -> int16_t (S16_LE) and detect clipping.
    //    We first store the sample values into derivative[]; then reuse the
    //    same array to hold |x[n] - x[n+1]|.
    int clip_overflow_count = 0;

    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        const uint8_t lsb = (uint8_t)buffer[2 * i + 0];
        const int8_t msb = (int8_t)buffer[2 * i + 1];
        // Combine as little-endian 16-bit
        const int16_t sample = (int16_t)(((int)msb << 8) | lsb);

        derivative[i] = (int)sample;

        // Clipping detection
        if (sample == INT16_MAX || sample == INT16_MIN)
        {
            clip_overflow_count++;
        }
    }

    if (clip_overflow_count > 1)
    {
        fprintf(stderr,
                "%d audio 16-bit clipping event(s)\n",
                clip_overflow_count);
    }

    // 4) Compute derivative: |x[n] - x[n+1]|, last sample = 0
    for (unsigned int i = 0; i + 1 < ArrayLength; ++i)
    {
        int d = derivative[i] - derivative[i + 1];
        derivative[i] = (d < 0) ? -d : d;
    }
    derivative[ArrayLength - 1] = 0;

    // Success: we read exactly ArrayLength frames
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
