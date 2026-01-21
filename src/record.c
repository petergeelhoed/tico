
// main.c
#include <alsa/asoundlib.h>
#include <errno.h>
#include <fftw3.h>
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

// ---- Helper: build ALSA poll descriptors ----
static int build_alsa_pollfds(snd_pcm_t* h,
                              struct pollfd** fds_out,
                              nfds_t* nfds_out)
{
    int count = snd_pcm_poll_descriptors_count(h);
    if (count <= 0)
    {
        fprintf(stderr, "ALSA: invalid poll descriptors count: %d\n", count);
        return -1;
    }
    struct pollfd* fds = (struct pollfd*)calloc((size_t)count, sizeof(*fds));
    if (!fds)
    {
        fprintf(stderr, "alloc failed for alsa fds\n");
        return -1;
    }
    int err = snd_pcm_poll_descriptors(h, fds, (unsigned)count);
    if (err < 0)
    {
        fprintf(stderr,
                "ALSA: snd_pcm_poll_descriptors failed: %s\n",
                snd_strerror(err));
        free(fds);
        return -1;
    }
    *fds_out = fds;
    *nfds_out = (nfds_t)count;
    return 0;
}

// ---- Helper: create/rearm a timerfd ----
static int make_timerfd_ms(unsigned initial_ms, unsigned interval_ms)
{
    int tfd = timerfd_create(CLOCK_MONOTONIC, TFD_NONBLOCK | TFD_CLOEXEC);
    if (tfd < 0)
    {
        perror("timerfd_create");
        return -1;
    }
    struct itimerspec its;
    memset(&its, 0, sizeof(its));
    its.it_value.tv_sec = initial_ms / 1000;
    its.it_value.tv_nsec = (long)(initial_ms % 1000) * 1000000L;
    its.it_interval.tv_sec = interval_ms / 1000;
    its.it_interval.tv_nsec = (long)(interval_ms % 1000) * 1000000L;

    if (timerfd_settime(tfd, 0, &its, NULL) < 0)
    {
        perror("timerfd_settime");
        close(tfd);
        return -1;
    }
    return tfd;
}

int main(int argc, char* argv[])
{
    unsigned int rate = DEFAULT_RATE; // e.g., 48000
    unsigned int bph = DEFAULT_BPH;   // e.g., 21600
    unsigned int time = 3;            // seconds to record (as blocks)
    unsigned int evalue = 4;          // filter param
    const char* device = NULL;

    unsigned int countdown = 0; // seconds to show concurrently with capture

    // ---- Parse args (added -c for countdown) ----
    int flag;
    while ((flag = getopt(argc, argv, "b:r:ht:d:e:c:")) != -1)
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
        case 'c': // countdown seconds to show concurrently with capture
            retVal = checkUIntArg(flag, &countdown, optarg);
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
                " -c countdown seconds shown concurrently with capture "
                "(default: 0)\n",
                DEFAULT_BPH,
                DEFAULT_RATE);
            exit(0);
        }
        if (retVal != 0)
            return retVal;
    }

    if (device == NULL)
        device = "default:2";

    // ---- Mutable copy of device ----
    size_t device_len = strlen(device);
    char* device_mutable = (char*)malloc(device_len + 1);
    if (!device_mutable)
    {
        (void)fprintf(stderr, "device memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    strncpy(device_mutable, device, device_len + 1);
    device_mutable[device_len] = '\0';

    // ---- Init ALSA (assumes initAudio sets period_event=1 and non-blocking)
    // ----
    snd_pcm_format_t format = SND_PCM_FORMAT_S16_LE;
    snd_pcm_t* capture_handle = initAudio(format, device_mutable, &rate);

    // ---- Compute your processing block size ----
    // ArrayLength (frames per block) and number of blocks to capture
    unsigned int ArrayLength =
        rate * SECS_HOUR * 2 / bph; // ~16000 @48k/21600bph
    unsigned int blocks_left = time * bph / 2 / SECS_HOUR; // number of blocks

    // ---- Build filter (kept for parity with your original) ----
    fftw_complex* filterFFT = makeFilter(evalue, ArrayLength);

    // ---- Allocate buffers ----
    const unsigned bytes_per_sample =
        (unsigned)snd_pcm_format_width(format) / 8; // 2
    const unsigned channels = 1;
    const unsigned bytes_per_frame = bytes_per_sample * channels;

    char* block_buf = (char*)malloc((size_t)ArrayLength * bytes_per_frame);
    if (!block_buf)
    {
        free(device_mutable);
        (void)fprintf(stderr, "buffer memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    struct myarr* rawread = makemyarr(ArrayLength);

    FILE* filePtr = fopen("recorded", "w");
    if (filePtr)
    {
        static char file_buf[1 << 20]; // 1 MiB buffer to reduce I/O stalls
        setvbuf(filePtr, file_buf, _IOFBF, sizeof(file_buf));
    }
    else
    {
        fprintf(stderr, "warning: could not open output file 'recorded'\n");
    }

    // ---- Query ALSA buffer/period sizes (read per period) ----
    snd_pcm_uframes_t period_size = 0, buffer_size = 0;
    if (snd_pcm_get_params(capture_handle, &buffer_size, &period_size) < 0)
    {
        fprintf(
            stderr,
            "ALSA: snd_pcm_get_params failed; defaulting period_size=1024\n");
        period_size = 1024;
    }
    if (period_size == 0)
        period_size = 1024;
    if (period_size > ArrayLength)
        period_size = ArrayLength;

    // ---- Optional warm-up reads ----
    readBufferRaw(capture_handle, block_buf, rawread);
    readBufferRaw(capture_handle, block_buf, rawread);

    // ---- Build ALSA poll descriptors ----
    struct pollfd* alsa_fds = NULL;
    nfds_t alsa_nfds = 0;
    if (build_alsa_pollfds(capture_handle, &alsa_fds, &alsa_nfds) < 0)
    {
        fprintf(stderr, "ALSA: failed to build pollfds\n");
        goto cleanup_all;
    }

    // ---- Create a 1 Hz timerfd for countdown (only if requested) ----
    int tfd = -1;
    if (countdown > 0)
    {
        tfd = make_timerfd_ms(/*initial*/ 1000, /*interval*/ 1000);
        if (tfd < 0)
        {
            fprintf(stderr,
                    "warning: timerfd not available; countdown disabled\n");
            countdown = 0;
        }
    }

    // ---- Combine ALSA fds + timerfd into a single poll array ----
    nfds_t nfds = alsa_nfds + ((countdown > 0) ? 1 : 0);
    struct pollfd* fds = (struct pollfd*)calloc(nfds, sizeof(*fds));
    if (!fds)
    {
        fprintf(stderr, "alloc failed for combined fds\n");
        goto cleanup_all;
    }
    // Copy ALSA fds
    for (nfds_t i = 0; i < alsa_nfds; ++i)
        fds[i] = alsa_fds[i];
    // Append timerfd (last entry)
    if (countdown > 0)
    {
        fds[alsa_nfds].fd = tfd;
        fds[alsa_nfds].events = POLLIN;
    }
    free(alsa_fds);
    alsa_fds = NULL;

    // ---- Per-block state ----
    unsigned frames_collected_in_block = 0; // progress inside the current block
    size_t byte_off = 0;                    // write offset in block_buf
    size_t samp_off = 0;                    // write offset in rawread->arr

    // ---- Main single-thread event loop ----
    const int POLL_TIMEOUT_MS = 2000; // guard timeout in case device stalls

    while (blocks_left > 0)
    {
        int pr = poll(fds, nfds, POLL_TIMEOUT_MS);
        if (pr == 0)
        {
            fprintf(stderr, "poll: timeout\n");
            continue;
        }
        else if (pr < 0)
        {
            if (errno == EINTR)
                continue;
            perror("poll");
            break;
        }

        // 1) Handle ALSA events (only pass ALSA fds to revents helper)
        unsigned short revents = 0;
        if (snd_pcm_poll_descriptors_revents(capture_handle,
                                             fds,
                                             (unsigned)alsa_nfds,
                                             &revents) < 0)
        {
            fprintf(stderr, "ALSA: revents failed\n");
        }
        else
        {
            if (revents & POLLERR)
            {
                // Recover and rebuild ALSA descriptors in the combined array
                fprintf(stderr, "ALSA: POLLERR → recover\n");
                if (snd_pcm_recover(capture_handle, -EPIPE, 1) < 0)
                {
                    fprintf(stderr, "ALSA: recover failed\n");
                    break;
                }
                struct pollfd* new_alsa = NULL;
                nfds_t new_n = 0;
                if (build_alsa_pollfds(capture_handle, &new_alsa, &new_n) < 0)
                {
                    fprintf(stderr, "ALSA: rebuild pollfds failed\n");
                    break;
                }
                // Rebuild combined fds
                nfds_t new_total = new_n + ((countdown > 0) ? 1 : 0);
                struct pollfd* tmp =
                    (struct pollfd*)calloc(new_total, sizeof(*tmp));
                if (!tmp)
                {
                    free(new_alsa);
                    break;
                }
                for (nfds_t i = 0; i < new_n; ++i)
                    tmp[i] = new_alsa[i];
                if (countdown > 0)
                {
                    tmp[new_n] =
                        fds[alsa_nfds]; // carry timerfd over (old index)
                }
                free(new_alsa);
                free(fds);
                fds = tmp;
                alsa_nfds = new_n;
                nfds = new_total;
            }
            if (revents & POLLIN)
            {
                // Read exactly one period (or final remainder of the block)
                unsigned frames_left_in_block =
                    ArrayLength - frames_collected_in_block;
                unsigned chunk = (frames_left_in_block > period_size)
                                     ? (unsigned)period_size
                                     : frames_left_in_block;

                int got = readBuffer(capture_handle,
                                     chunk,
                                     block_buf + byte_off,
                                     rawread->arr + samp_off);
                if (got < 0)
                {
                    // Try additional recovery and continue
                    int rr = snd_pcm_recover(capture_handle, got, 1);
                    if (rr < 0)
                    {
                        fprintf(
                            stderr,
                            "ALSA: readBuffer failed & recover failed: %s\n",
                            snd_strerror(rr));
                        break;
                    }
                    // After recover, ALSA fds may change; rebuild as above
                    struct pollfd* new_alsa = NULL;
                    nfds_t new_n = 0;
                    if (build_alsa_pollfds(capture_handle, &new_alsa, &new_n) ==
                        0)
                    {
                        nfds_t new_total = new_n + ((countdown > 0) ? 1 : 0);
                        struct pollfd* tmp =
                            (struct pollfd*)calloc(new_total, sizeof(*tmp));
                        if (tmp)
                        {
                            for (nfds_t i = 0; i < new_n; ++i)
                                tmp[i] = new_alsa[i];
                            if (countdown > 0)
                                tmp[new_n] = fds[alsa_nfds]; // carry timerfd
                            free(fds);
                            fds = tmp;
                            alsa_nfds = new_n;
                            nfds = new_total;
                        }
                        free(new_alsa);
                    }
                    continue; // retry next poll
                }

                // Advance write offsets within the current block
                byte_off += (size_t)got * bytes_per_frame;
                samp_off += (size_t)got;
                frames_collected_in_block += (unsigned)got;

                // If we’ve completed one full block → process/persist
                if (frames_collected_in_block >= ArrayLength)
                {
                    syncAppendMyarr(rawread, filePtr);

                    // Optional: progress print (1 per block)
                    printf("%u\n", blocks_left - 1);
                    fflush(stdout);

                    // Reset for next block
                    frames_collected_in_block = 0;
                    byte_off = 0;
                    samp_off = 0;
                    if (blocks_left > 0)
                    {
                        --blocks_left;
                    }
                }
            }
        }

        // 2) Handle timerfd (countdown) — only if requested
        if (countdown > 0)
        {
            struct pollfd* t =
                &fds[alsa_nfds]; // timerfd is placed after ALSA fds
            if (t->revents & POLLIN)
            {
                uint64_t expirations = 0;
                ssize_t rd = read(t->fd, &expirations, sizeof(expirations));
                if (rd == (ssize_t)sizeof(expirations))
                {
                    while (expirations-- && countdown > 0)
                    {
                        --countdown;
                        printf("%u\n", countdown);
                        fflush(stdout);
                    }
                    if (countdown == 0)
                    {
                        // Disarm timer; keep fd in array but with events=0
                        struct itimerspec its = {0};
                        timerfd_settime(t->fd, 0, &its, NULL);
                        t->events = 0;
                    }
                }
            }
        }

        // Clear revents
        for (nfds_t i = 0; i < nfds; ++i)
            fds[i].revents = 0;
    }

    // ---- Cleanup & exit ----
cleanup_all:
    if (filePtr)
        wait_close(filePtr);
    if (capture_handle)
        snd_pcm_close(capture_handle);
    if (block_buf)
        free(block_buf);
    if (device_mutable)
        free(device_mutable);
    if (fds)
        free(fds);
    if (rawread)
        freemyarr(rawread);
    if (filterFFT)
        fftw_free(filterFFT);
    fftw_cleanup();
    return 0;
}
