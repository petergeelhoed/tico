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
#include "parseargs.h"

int main(int argc, char* argv[])
{
    unsigned int rate = DEFAULT_RATE;
    unsigned int bph = DEFAULT_BPH;
    unsigned int time = 3;
    int flag;
    unsigned int evalue = 4;
    const char* device = NULL;
    // declarations

    while ((flag = getopt(argc, argv, "b:r:ht:d:e:")) != -1)
    {
        int retVal = 0;
        switch (flag)
        {
        case 'e':
            retVal = checkUIntArg(flag, &evalue, optarg);
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
            retVal = checkUIntArg(flag, &time, optarg);
            break;
        case 'b':
            retVal = checkUIntArg(flag, &bph, optarg);
            break;
        case 'r':
            retVal = checkUIntArg(flag, &rate, optarg);
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

    if (device == 0)
    {
        device = "default:1";
    }

    size_t device_len = strlen(device);
    char* device_mutable = (char*)malloc(device_len);
    if (device_mutable == NULL)
    {
        exit(EXIT_FAILURE);
    }
    memcpy(device_mutable, device, device_len);

    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device_mutable, &rate);
    unsigned int ArrayLength = rate * SECS_HOUR * 2 / bph;
    unsigned int length = time * bph / 2 / SECS_HOUR;

    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);
    char* buffer =
        malloc(ArrayLength * (unsigned int)snd_pcm_format_width(format) /
               BITS_IN_BYTE);
    struct myarr rawread = {calloc(ArrayLength, sizeof(int)), 0, ArrayLength};

    FILE* filePtr = fopen("recorded", "w");
    readBufferRaw(capture_handle, buffer, &rawread);
    readBufferRaw(capture_handle, buffer, &rawread);
    while (length)
    {
        length--;
        readBufferRaw(capture_handle, buffer, &rawread);

        syncAppendMyarr(&rawread, filePtr);
        printf("%d\n", length);
    }

    fftw_free(filterFFT);
    wait();
    thread_lock();
    (void)fclose(filePtr);
    thread_unlock();
    exit(0);
}
