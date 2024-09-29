#include <limits.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mysound.h"

snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int rate)
{
    int err;
    snd_pcm_t* capture_handle;
    snd_pcm_hw_params_t* hw_params;

    if ((err = snd_pcm_open(
             &capture_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0)
    {
        fprintf(stderr,
                "cannot open audio device %s (%s)\n",
                device,
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_malloc(&hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot allocate hardware parameter structure (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_any(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr,
                "cannot initialize hardware parameter structure (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_access(
             capture_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0)
    {
        fprintf(stderr, "cannot set access type (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_format(
             capture_handle, hw_params, format)) < 0)
    {
        fprintf(stderr, "cannot set sample format (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_rate_near(
             capture_handle, hw_params, &rate, 0)) < 0)
    {
        fprintf(stderr, "cannot set sample rate (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1)) <
        0)
    {
        fprintf(stderr, "cannot set channel count (%s)\n", snd_strerror(err));
        exit(1);
    }

    if ((err = snd_pcm_hw_params(capture_handle, hw_params)) < 0)
    {
        fprintf(stderr, "cannot set parameters (%s)\n", snd_strerror(err));
        exit(1);
    }

    snd_pcm_hw_params_free(hw_params);

    if ((err = snd_pcm_prepare(capture_handle)) < 0)
    {
        fprintf(stderr,
                "cannot prepare audio interface for use (%s)\n",
                snd_strerror(err));
        exit(1);
    }

    return capture_handle;
}

void readBufferRaw(snd_pcm_t* capture_handle, int NN, char* buffer, int* in)
{
    unsigned char lsb;
    signed char msb;
    int err;
    if ((err = snd_pcm_readi(capture_handle, buffer, NN)) != NN)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        exit(-1);
    }
    for (int j = 0; j < NN * 2; j += 2)
    {
        msb = *(buffer + j + 1);
        lsb = *(buffer + j);
        in[j / 2] = (msb << 8) | lsb;
    }
}

int readBuffer(snd_pcm_t* capture_handle, int NN, char* buffer, int* derivative)
{
    unsigned char lsb;
    signed char msb;
    int err;
    if ((err = snd_pcm_readi(capture_handle, buffer, NN)) != NN)
    {
        fprintf(stderr,
                "read from audio interface failed %d (%s)\n",
                err,
                snd_strerror(err));
        return err;
    }
    for (int j = 0; j < NN * 2; j += 2)
    {
        msb = *(buffer + j + 1);
        lsb = *(buffer + j);
        derivative[j / 2] = (msb << 8) | lsb;
    }
    //       remove50hz(NN,in,48000);

    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
    derivative[NN] = 0;
    return err;
}

int readShiftedBuffer(int* derivative,
                      snd_pcm_t* capture_handle,
                      int NN,
                      char* buffer,
                      int maxpos,
                      int* totalshift,
                      FILE* fpInput)
{
    int ret;
    int shift = (int)sqrt(abs(maxpos));

    
    if (maxpos < 0)
    {
        *totalshift -= shift;
        memcpy(derivative + NN - shift, derivative, shift * sizeof(int));
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<NN-shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else
        {
        ret = readBuffer(capture_handle, NN - shift, buffer, derivative);
        }
    }
    else if (maxpos > 0)
    {
        *totalshift += shift;
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
            for (int j=0 ; j<NN; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else
        {
            ret = readBuffer(capture_handle, shift, buffer, derivative);

            if (ret != -32)
            {
                ret = readBuffer(capture_handle, NN, buffer, derivative);
            }
        }
    }
    else
    {
        if (fpInput)
        {
            ret = 0;
            for (int j=0 ; j<shift; ++j)
            {
                if (fscanf(fpInput,"%d",derivative+j) != 1)
                {
                    ret = -32;
                    break;
                }
            }
    for (int j = 0; j < NN - 1; j++)
    {
        derivative[j] = fabs(derivative[j] - derivative[j + 1]);
    }
        }
        else 
        {
            ret = readBuffer(capture_handle, NN, buffer, derivative);
        }
    }
    return ret;
}
