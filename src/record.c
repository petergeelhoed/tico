#include <alsa/asoundlib.h>
#include <errno.h>
#include <poll.h>
#include <stdint.h> // uint64_t
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // strlen, strncpy
#include <sys/timerfd.h>
#include <time.h>
#include <unistd.h> // getopt, read

#include "config.h"
#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mysound.h"
#include "mysync.h"
#include "parseargs.h"

/* -------------------- Capture Context -------------------- */

/* -------------------- Main -------------------- */

int main(int argc, char* argv[])
{
    unsigned int rate = DEFAULT_RATE; // e.g., 48000
    unsigned int bph = DEFAULT_BPH;   // e.g., 21600
    unsigned int time = 3;            // seconds to record (as blocks)
    unsigned int evalue = 4;          // filter param
    const char* device = NULL;

    int flag;
    while ((flag = getopt(argc, argv, "b:r:ht:d:e:c:")) != -1)
    {
        int retVal = 0;
        switch (flag)
        {
        case 'e':
            retVal = checkUIntArg(flag, &evalue, optarg);
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
                "usage: capture\n"
                "capture reads from the microphone and timegraphs your watch\n"
                "options:\n"
                " -d <capture device> (default: 'default:2')\n"
                " -b bph of the watch (default: %u/h)\n"
                " -r sampling rate (default: %u Hz)\n"
                " -t time to record in seconds (default: 3)\n"
                " -e envelope level (default: 4)\n"
                "(default: 0)\n",
                DEFAULT_BPH,
                DEFAULT_RATE);
            return 0;
        }
        if (retVal != 0)
        {
            return retVal;
        }
    }

    if (device == NULL)
    {
        device = "default:2";
    }

    // Mutable device string
    size_t deviceLen = strlen(device);
    char* deviceMutable = (char*)malloc(deviceLen + 1);
    if (!deviceMutable)
    {
        (void)fprintf(stderr, "device memory allocation failed\n");
        return EXIT_FAILURE;
    }
    strncpy(deviceMutable, device, deviceLen + 1);
    deviceMutable[deviceLen] = '\0';

    // Init ALSA (your function)
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* cap = initAudio(format, deviceMutable, &rate);
    if (!cap)
    {
        free(deviceMutable);
        return EXIT_FAILURE;
    }

    // Setup capture context (poll + timerfd, buffers, etc.)
    CaptureCtx ctx;
    CapConfig cfg = {.rate = rate,
                     .bph = bph,
                     .evalue = DEFAULT_EVALUE,
                     .zoom = DEFAULT_ZOOM,
                     .fitN = DEFAULT_FITN,
                     .teeth = DEFAULT_TEETH,
                     .SDthreshold = DEFAULT_SDTHRESHOLD,
                     .device = "default:2",
                     .cvalue = DEFAULT_CVALUE,
                     .fpposition = NULL,
                     .fpmaxcor = NULL,
                     .fptotal = NULL,
                     .fpDefPeak = NULL,
                     .fpInput = NULL,
                     .captureHandle = cap};

    if (captureSetup(&ctx, &cfg, rate) < 0)
    {
        (void)fprintf(stderr, "captureSetup failed\n");
        free(deviceMutable);
        snd_pcm_close(cap);
        return EXIT_FAILURE;
    }

    // Open output file (buffered) — optional
    FILE* filePtr = fopen("recorded", "w");

    // Number of blocks to capture in total
    unsigned int blocksLeft = time * bph / 2 / SECS_HOUR;
    const unsigned int ArrayLength = 16000;

    int16_t* out = calloc(ArrayLength, sizeof(int16_t)); // 16bit
    while (out != NULL && blocksLeft > 0)
    {

        struct myarr* filled = makemyarr(ArrayLength);
        int read = readSamples(ctx.cap, ArrayLength, out);
        if (!read)
        {
            (void)fprintf(stderr, "capture_next_block failed; stopping\n");
            break;
        }

        for (unsigned int i = 0; i < ArrayLength; ++i)
        {
            filled->arr[i] = out[i];
        }

        // Your existing persist/processing step
        syncAppendMyarr(filled, filePtr);

        // Progress print (once per block)
        printf("%u\n", blocksLeft - 1);
        --blocksLeft;
    }

    // Cleanup
    if (filePtr)
    {
        waitClose(filePtr);
        captureTeardown(&ctx);
    }
    if (cap)
    {
        snd_pcm_close(cap);
    }
    if (deviceMutable)
    {
        free(deviceMutable);
    }
    free(out);
    captureTeardown(&ctx);
    return 0;
}
