
#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // OPTION A: for sleep()

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

    unsigned int countdown =
        0; // OPTION A: seconds to count down before capture

    // Parse args
    while ((flag = getopt(argc, argv, "b:r:ht:d:e:c:")) !=
           -1) // OPTION A: added 'c:'
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
        case 'c': // OPTION A: countdown seconds
            retVal = checkUIntArg(flag, &countdown, optarg);
            break;
        case 'h':
        default:
            (void)fprintf(stderr,
                          "usage: capture \n"
                          "capture reads from the microphone and timegraphs "
                          "your watch\n"
                          "options:\n"
                          " -d <capture device> (default: 'default:1')\n"
                          " -b bph of the watch (default: 21600/h)\n"
                          " -r sampling rate (default: 48000Hz)\n"
                          " -t time to record (default: 30s)\n"
                          " -e envelope level (default: 4)\n"
                          " -c countdown seconds before starting capture "
                          "(default: 0)\n" // OPTION A
            );
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
        device = "default:2";
    }

    size_t device_len = strlen(device);
    char* device_mutable = (char*)malloc(device_len + 1);
    if (device_mutable == NULL)
    {
        (void)fprintf(stderr, "devive memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    strncpy(device_mutable, device, device_len + 1);
    device_mutable[device_len] = '\0';

    // Init ALSA
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device_mutable, &rate);

    // OPTION A: Do a blocking countdown BEFORE any reads so capture isn't
    // starved
    if (countdown > 0)
    {
        for (int i = (int)countdown; i > 0; --i)
        {
            printf("%d\n", i);
            fflush(stdout);
            sleep(1);
        }
        // Ensure fresh start after countdown (discard any pending data if the
        // device auto-ran)
        if (snd_pcm_prepare(capture_handle) < 0)
        {
            fprintf(stderr, "ALSA: re-prepare after countdown failed\n");
            // Not fatal—continue, recover will happen later if needed
        }
    }

    unsigned int ArrayLength = rate * SECS_HOUR * 2 / bph;
    unsigned int length = time * bph / 2 / SECS_HOUR;

    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);

    char* buffer =
        (char*)malloc(ArrayLength * (unsigned int)snd_pcm_format_width(format) /
                      BITS_IN_BYTE);
    if (buffer == NULL)
    {
        free(device_mutable);
        (void)fprintf(stderr, "buffer mamory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    struct myarr* rawread = makemyarr(ArrayLength);

    FILE* filePtr = fopen("recorded", "w");

    // Prime the stream (now that countdown is done)
    readBufferRaw(capture_handle, buffer, rawread);
    readBufferRaw(capture_handle, buffer, rawread);

    while (length)
    {
        length--;
        readBuffer(capture_handle, ArrayLength, buffer, rawread->arr);

        syncAppendMyarr(rawread, filePtr);
        printf("%d\n", length);
        fflush(stdout);
    }

    free(buffer);
    free(device_mutable);
    fftw_free(filterFFT);
    fftw_cleanup();
    wait_close(filePtr);
    if (capture_handle)
    {
        snd_pcm_close(capture_handle);
    }
    freemyarr(rawread);
    exit(0);
}
