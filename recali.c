#include <alsa/asoundlib.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"

int main()
{
    unsigned int rate = 48000;
    char* device = 0;
    // declarations
    unsigned int NN = 48000 * 1;

    device = device == 0 ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, rate);
    char* buffer = malloc(NN * (unsigned int)snd_pcm_format_width(format) / 8);

    int rawin[NN];
    readBufferRaw(capture_handle, 8000, buffer, rawin);
    readBufferRaw(capture_handle, NN, buffer, rawin);
    double out[NN];
    for (unsigned int j = 0; j < NN; j++)
    {
        out[j] = rawin[j];
    }
    remove50hz(NN, rawin, rate);

    for (unsigned int j = 0; j < NN; j++)
    {
        printf("%d %f %d\n", j, out[j], rawin[j]);
    }
    exit(0);
}
