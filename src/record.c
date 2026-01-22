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

#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
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
    size_t device_len = strlen(device);
    char* device_mutable = (char*)malloc(device_len + 1);
    if (!device_mutable)
    {
        (void)fprintf(stderr, "device memory allocation failed\n");
        return EXIT_FAILURE;
    }
    strncpy(device_mutable, device, device_len + 1);
    device_mutable[device_len] = '\0';

    // Init ALSA (your function)
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* cap = initAudio(format, device_mutable, &rate);
    if (!cap)
    {
        free(device_mutable);
        return EXIT_FAILURE;
    }

    // Setup capture context (poll + timerfd, buffers, etc.)
    CaptureCtx ctx;
    if (capture_setup(&ctx, cap, rate, bph, format) < 0)
    {
        (void)fprintf(stderr, "capture_setup failed\n");
        free(device_mutable);
        snd_pcm_close(cap);
        return EXIT_FAILURE;
    }

    // Open output file (buffered) — optional
    FILE* filePtr = fopen("recorded", "w");
    if (!filePtr)
    {
        (void)fprintf(stderr,
                      "warning: could not open output file 'recorded'\n");
    }

    // Number of blocks to capture in total
    unsigned int blocks_left = time * bph / 2 / SECS_HOUR;
    const int POLL_TIMEOUT_MS = 2000;

    while (blocks_left > 0)
    {
        //       printTOD(filePtr);
        struct myarr* filled = capture_next_block(&ctx, POLL_TIMEOUT_MS);
        if (!filled)
        {
            (void)fprintf(stderr, "capture_next_block failed; stopping\n");
            break;
        }

        // Your existing persist/processing step
        syncAppendMyarr(filled, filePtr);

        // Progress print (once per block)
        printf("%u\n", blocks_left - 1);
        --blocks_left;
    }

    // Cleanup
    if (filePtr)
    {
        wait_close(filePtr);
        capture_teardown(&ctx);
    }
    if (cap)
    {
        snd_pcm_close(cap);
    }
    if (device_mutable)
    {
        free(device_mutable);
    }
    return 0;
}
