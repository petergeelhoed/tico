#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"

int main(int argc, char* argv[])
{
    unsigned int rate = DEFAULT_RATE;
    unsigned int bph = DEFAULT_BPH;
    unsigned int time = 3;
    int c;
    unsigned int evalue = 4;
    char* device = 0;
    // declarations

    while ((c = getopt(argc, argv, "b:r:ht:d:e:")) != -1)
    {
        int retVal = 0;
        switch (c)
        {
        case 'e':
            retVal = checkUIntArg(c, &evalue, optarg);
            if (evalue == 0)
            {
                printf("invalid integer argument for -e '%s'\n", optarg);
                return -1;
            }
            break;
        case 'd':
            device = optarg;
            break;
        case 't':
            retVal = checkUIntArg(c, &time, optarg);
            break;
        case 'b':
            retVal = checkUIntArg(c, &bph, optarg);
            break;
        case 'r':
            retVal = checkUIntArg(c, &rate, optarg);
            break;
        case 'h':
        default:
            (void)fprintf(
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
        if (retVal != 0)
        {
            return retVal;
        }
    }

    device = device == 0 ? "default:1" : device;

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device, &rate);
    unsigned int ArrayLength = rate * SECS_HOUR * 2 / bph;
    unsigned int length = time * bph / 2 / SECS_HOUR;

    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);
    char* buffer =
        malloc(ArrayLength * (unsigned int)snd_pcm_format_width(format) / BITS_IN_BYTE);
    struct myarr rawread = {calloc(ArrayLength, sizeof(int)), 0, ArrayLength};

    FILE* fp = fopen("recorded", "w");
    readBufferRaw(capture_handle, buffer, &rawread);
    readBufferRaw(capture_handle, buffer, &rawread);
    while (length)
    {
        length--;
        readBufferRaw(capture_handle, buffer, &rawread);

        syncAppendMyarr(&rawread, fp);
        printf("%d\n", length);
    }

    fftw_free(filterFFT);
    wait();
    thread_lock();
    (void)fclose(fp);
    thread_unlock();
    exit(0);
}
