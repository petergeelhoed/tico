#include <alsa/asoundlib.h>
#include <errno.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "mylib.h"
#include "mysound.h"

// Helper to handle repetitive ALSA parameter setting and error reporting
static void check_alsa_err(int err,
                           const char* msg,
                           snd_pcm_t* handle,
                           snd_pcm_hw_params_t* params)
{
    if (err < 0)
    {
        (void)fprintf(stderr, "%s (%s)\n", msg, snd_strerror(err));
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
        "cannot open audio device",
        NULL,
        NULL);

    // 2. Allocate and Init Params
    check_alsa_err(snd_pcm_hw_params_malloc(&hw_params),
                   "cannot allocate hardware parameter structure",
                   capture_handle,
                   NULL);

    check_alsa_err(snd_pcm_hw_params_any(capture_handle, hw_params),
                   "cannot initialize hardware parameter structure",
                   capture_handle,
                   hw_params);

    // 3. Set Hardware Configurations
    check_alsa_err(snd_pcm_hw_params_set_access(capture_handle,
                                                hw_params,
                                                SND_PCM_ACCESS_RW_INTERLEAVED),
                   "cannot set access type",
                   capture_handle,
                   hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_format(capture_handle, hw_params, format),
        "cannot set sample format",
        capture_handle,
        hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_rate_near(capture_handle, hw_params, rate, 0),
        "cannot set sample rate",
        capture_handle,
        hw_params);

    check_alsa_err(snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1),
                   "cannot set channel count",
                   capture_handle,
                   hw_params);

    // 4. Apply Params and Prepare
    check_alsa_err(snd_pcm_hw_params(capture_handle, hw_params),
                   "cannot set parameters",
                   capture_handle,
                   hw_params);

    snd_pcm_hw_params_free(hw_params);

    check_alsa_err(snd_pcm_prepare(capture_handle),
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
    int err = snd_pcm_readi(capture_handle, buffer, data_in->ArrayLength);
    if (err != (int)data_in->ArrayLength)
    {
        (void)fprintf(stderr,
                      "read from audio interface failed %d (%s)\n",
                      err,
                      snd_strerror(err));
        exit(READ_FAILED);
    }
    for (unsigned int index = 0; index < 2 * data_in->ArrayLength; index += 2)
    {
        msb = (signed char)buffer[index + 1];
        lsb = *(buffer + index);
        data_in->arr[index / 2] = (msb << BITS_IN_BYTE) | lsb;
    }
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative)
{
    unsigned char lsb;
    signed char msb;
    int err =
        snd_pcm_readi(capture_handle, buffer, (long unsigned int)ArrayLength);
    if (err < 0)
    {
        (void)fprintf(stderr,
                      "read from audio interface failed %d (%s)\n",
                      err,
                      snd_strerror(err));
        return err;
    }
    if (err != (int)ArrayLength)
    {
        (void)fprintf(stderr, "reread from audio interface  %d \n", err);
        err = snd_pcm_readi(
            capture_handle, buffer + err, (long unsigned int)ArrayLength - err);
        if (err < 0)
        {
            (void)fprintf(stderr,
                          "reread from audio interface failed %d (%s)\n",
                          err,
                          snd_strerror(err));
        }
        else
        {
            err = (int)ArrayLength;
        }
    }
    for (unsigned int index = 0; index < ArrayLength * 2; index += 2)
    {
        msb = (signed char)buffer[index + 1];
        lsb = *(buffer + index);
        derivative[index / 2] = (msb << BITS_IN_BYTE) | lsb;
    }
    //       remove50hz(ArrayLength,data_in,48000);

    for (unsigned int index = 0; index < ArrayLength - 1; index++)
    {
        derivative[index] = abs(derivative[index] - derivative[index + 1]);
    }
    derivative[ArrayLength - 1] = 0;
    return err;
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
