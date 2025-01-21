#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "mylib.h"
#include "mysound.h"

#define INIT_ERROR -2
#define READ_FAILED -1
#define REINIT_ERROR -32
#define INPUT_FILE_ERROR -33

// Initialize audio capture
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate)
{
    int err;
    snd_pcm_t* capture_handle = NULL;
    snd_pcm_hw_params_t* hw_params = NULL;

    if ((err = snd_pcm_open(
             &capture_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0)
    {
        fprintf(stderr,
                "cannot open audio device %s (%s)\n",
                device,
                snd_strerror(err));
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_malloc(&hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot allocate hardware parameter structure (%s)\n",
                snd_strerror(err));
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_any(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot initialize hardware parameter structure (%s)\n",
                snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_set_access(
             capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0)
    {
        fprintf(stderr, "cannot set access type (%s)\n", snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_set_format(
             capture_handle, hw_params, format)) < 0)
    {
        fprintf(stderr, "cannot set sample format (%s)\n", snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_set_rate_near(
             capture_handle, hw_params, &rate, 0)) < 0)
    {
        fprintf(stderr, "cannot set sample rate (%s)\n", snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1)) <
        0)
    {
        fprintf(stderr, "cannot set channel count (%s)\n", snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    if ((err = snd_pcm_hw_params(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr, "cannot set parameters (%s)\n", snd_strerror(err));
        if (hw_params)
            snd_pcm_hw_params_free(hw_params);
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    snd_pcm_hw_params_free(hw_params);

    if ((err = snd_pcm_prepare(capture_handle)) < 0)
    {
        fprintf(stderr,
                "cannot prepare audio interface for use (%s)\n",
                snd_strerror(err));
        if (capture_handle)
            snd_pcm_close(capture_handle);
        exit(INIT_ERROR);
    }

    return capture_handle;
}

void readBufferRaw(snd_pcm_t* capture_handle,
                   unsigned int NN,
                   char* buffer,
                   int* in)
{
    unsigned char lsb;
    signed char msb;
    int err = snd_pcm_readi(capture_handle, buffer, (long unsigned int)NN);
    if (err != (int)NN)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        exit(READ_FAILED);
    }
    for (unsigned int j = 0; j < 2 * NN; j += 2)
    {
        msb = (signed char)buffer[j + 1];
        lsb = *(buffer + j);
        in[j / 2] = (msb << 8) | lsb;
    }
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int NN,
               char* buffer,
               int* derivative)
{
    unsigned char lsb;
    signed char msb;
    int err = snd_pcm_readi(capture_handle, buffer, (long unsigned int)NN);
    if (err < 0)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        return err;
    }
    if (err != (int)NN)
    {
        fprintf(stderr, "reread from audio interface  %d \n", err);
        err = snd_pcm_readi(
            capture_handle, buffer + err, (long unsigned int)NN - err);
        if (err < 0)
        {
            fprintf(stderr,
                    "reread from audio interface failed %d (%s)\n",
                    err,
                    snd_strerror(err));
        }
        else
        {
            err = (int)NN;
        }
    }
    for (unsigned int j = 0; j < NN * 2; j += 2)
    {
        msb = (signed char)buffer[j + 1];
        lsb = *(buffer + j);
        derivative[j / 2] = (msb << 8) | lsb;
    }
    //       remove50hz(NN,in,48000);

    for (unsigned int j = 0; j < NN - 1; j++)
    {
        derivative[j] = abs(derivative[j] - derivative[j + 1]);
    }
    derivative[NN - 1] = 0;
    return err;
}

int readBufferOrFile(int* derivative,
                     snd_pcm_t* capture_handle,
                     unsigned int NN,
                     char* buffer,
                     FILE* fpInput)
{
    int ret = READ_FAILED;
    if (fpInput)
    {
        ret = (int)NN;
        for (unsigned int j = 0; j < NN; ++j)
        {
            if (fscanf(fpInput, "%d", derivative + j) != 1)
            {
                ret = INPUT_FILE_ERROR;
                break;
            }
        }
        for (unsigned int j = 0; j < NN - 1; j++)
        {
            derivative[j] = abs(derivative[j] - derivative[j + 1]);
        }
        derivative[NN - 1] = 0;
    }
    else
    {
        ret = readBuffer(capture_handle, NN, buffer, derivative);
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
        err = readBufferOrFile(
            derivative.arr, capture_handle, derivative.NN, buffer, fpInput);
        if (err == REINIT_ERROR)
        {
            fprintf(stderr, "Reinitializing capture_handle");
            if (rawfile)
            {
                fprintf(rawfile, "# Reinitializing capture_handle");
            }
            snd_pcm_close(capture_handle);
            capture_handle = initAudio(format, device, rate);
            err = readBufferOrFile(
                derivative.arr, capture_handle, derivative.NN, buffer, fpInput);
        }
        if (err == INPUT_FILE_ERROR)
        {
            fprintf(stderr, "Could not read integer from inputfile or audio\n");
        }
    }
    return err;
}
