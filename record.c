#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

#include "mylib.h"
#include "mysound.h"

int main(int argc, char* argv[])
{
    unsigned int rate = 48000;
    int bph = 21600;
    int time = 30;
    int c;
    int evalue = 4;
    char* device = 0;
    // declarations
    int NN = rate * 3600 / bph * 2 ;

    while ((c = getopt(argc, argv, "b:r:ht:d:e:")) != -1)
    {
        switch (c)
        {
        case 'e':
            evalue = atoi(optarg);
            break;
        case 'd':
            device = optarg;
            break;
        case 't':
            time = atoi(optarg);
            break;
        case 'b':
            bph = atoi(optarg);
            break;
        case 'r':
            rate = atoi(optarg);
            break;
        case 'h':
        default:
            fprintf(
                stderr,
                "usage: capture \n"
                "capture reads from the microphone and timegraphs your watch\n"
                "options:\n"
                " -d <capture device> (default: 'default:1')\n"
                " -b bph of the watch (default: 21600/h) \n"
                " -r sampling rate (default: 48000Hz)\n"
                " -t time to record (default: 30s )\n");
            exit(0);
            break;
        }
    }
    int length = time * bph / 7200;

    fftw_complex* filterFFT = makeFilter(evalue, NN);

    device = device == 0 ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, rate);
    char* buffer = malloc(NN * snd_pcm_format_width(format) / 8);
    int rawread[NN];
    FILE* fp = fopen("recorded", "w");
    readBufferRaw(capture_handle, 8000, buffer, rawread);
    readBufferRaw(capture_handle, 8000, buffer, rawread);
    readBufferRaw(capture_handle, 8000, buffer, rawread);
    while (length)
    {
        length--;
        readBufferRaw(capture_handle, NN, buffer, rawread);

        //double out[NN];
        //applyFilter(rawread, NN, filterFFT, out);
 /*       for (int j = 0; j < NN; j++)
        {
            fprintf(fp, "%d\n", rawread[j]);
        }
   */
        syncappend(rawread, NN , fp);
        printf("%d\n", length);
    }

    fftw_free(filterFFT);
    exit(0);
}
