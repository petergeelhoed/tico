#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>

#include "mylib.h"


int main (int argc, char *argv[])
{
    unsigned int rate = 48000;
    char *device = 0;
    // declarations
    int NN = 48000*1;

    device = device==0?"default:1":device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);
    char *buffer = malloc(NN * snd_pcm_format_width(format) / 8);

    int rawin[NN];
    readBufferRaw(capture_handle, 8000, buffer, rawin);
    readBufferRaw(capture_handle, NN, buffer, rawin);
    double out[NN];
    for (int j=0; j <NN ; j++)
    {
        out[j] = rawin[j];
    }
    remove50hz(NN,rawin,rate);

    for (int j=0; j <NN ; j++)
    {
        printf("%d %f %d\n",j,out[j],rawin[j]);
    }
    exit (0);
}

