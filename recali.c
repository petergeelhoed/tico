#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

#include "mylib.h"


int main (int argc, char *argv[])
{
    unsigned int rate = 48000;
    int evalue = 4;
    char *device = 0;
    // declarations
    int NN = 48000*2;

    fftw_complex *filterFFT = makeFilter(evalue, NN);

    device = device==0?"default:1":device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);
    char *buffer = malloc(NN * snd_pcm_format_width(format) / 8);

    int rawin[NN];
    readBufferRaw(capture_handle, 8000, buffer, rawin);
    readBufferRaw(capture_handle, NN, buffer, rawin);
    double out[NN];
    applyFilter50(rawin,NN,filterFFT,out);

    for (int j=0; j <NN ; j++)
    {
        printf("%d %f %d\n",j,out[j],rawin[j]);
    }
    fftw_free(filterFFT);
    exit (0);
}

