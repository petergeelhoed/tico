#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>

#include "mylib.h"


int main (int argc, char *argv[])
{
    unsigned int rate = 48000;
    int bph = 21600;
    int time = 30;
    int c;
    char *device = 0;
    // declarations
    int NN = rate*3600/bph;

    while ((c = getopt (argc, argv, "b:r:ht:d:w:p:")) != -1)
    {
        switch (c)
        {
            case 'd':
                device = optarg;
                break;
            case 't':
                time = atoi(optarg);
                NN=time*NN;
                break;
            case 'b':
                bph = atoi(optarg);
                break;
            case 'r':
                rate = atoi(optarg);
                break;
            case 'h':
            default:
                fprintf (stderr,
                        "usage: capture \n"\
                        "capture reads from the microphone and timegraphs your watch\n" 
                        "options:\n"\
                        " -d <capture device> (default: 'default:1')\n"\
                        " -b bph of the watch (default: 21600/h) \n"\
                        " -r sampling rate (default: 48000Hz)\n"\
                        " -t <measurment time> (default: 30s)\n");
                exit(0);
                break;
        }
    }


    device = device==0?"default:1":device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t *capture_handle = initAudio(format, device, rate);
    char *buffer = malloc(NN * snd_pcm_format_width(format) / 8);
    int derivative[NN];
    readBuffer(capture_handle, NN, buffer, derivative);

    for (int j=0; j <NN ; j++)
    {
        printf("%d %d\n",j,derivative[j]);
    }
    exit (0);
}

