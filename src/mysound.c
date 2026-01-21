#include "mysound.h"
#include "config.h"
#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
#include "mysound.h"
#include "mysync.h"
#include "parseargs.h"

#include <alsa/asoundlib.h>
#include <errno.h>
#include <fftw3.h>
#include <limits.h>
#include <math.h>
#include <poll.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // strlen, strncpy
#include <sys/timerfd.h>
#include <time.h>
#include <unistd.h> // getopt, read

// Helper to handle repetitive ALSA parameter setting and error reporting
static void check_alsa_err(int err,
                           const char* device,
                           const char* msg,
                           snd_pcm_t* handle,
                           snd_pcm_hw_params_t* params)
{
    if (err < 0)
    {
        (void)fprintf(stderr,
                      "Device %s: %s (%s)\n",
                      device,
                      msg,
                      snd_strerror(err));
        if (params)
        {
            snd_pcm_hw_params_free(params);
        }
        if (handle)
        {
            snd_pcm_close(handle);
        }
        exit(INIT_ERROR);
    }
}

// Initialize audio capture
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate)
{
    snd_pcm_t* capture_handle = NULL;
    snd_pcm_hw_params_t* hw_params = NULL;
    unsigned int requestedRate = *rate;

    // 1. Open Device
    check_alsa_err(
        snd_pcm_open(&capture_handle, device, SND_PCM_STREAM_CAPTURE, 0),
        device,
        "cannot open audio device",
        NULL,
        NULL);

    // 2. Allocate and Init Params
    check_alsa_err(snd_pcm_hw_params_malloc(&hw_params),
                   device,
                   "cannot allocate hardware parameter structure",
                   capture_handle,
                   NULL);

    check_alsa_err(snd_pcm_hw_params_any(capture_handle, hw_params),
                   device,
                   "cannot initialize hardware parameter structure",
                   capture_handle,
                   hw_params);

    // 3. Set Hardware Configurations
    check_alsa_err(snd_pcm_hw_params_set_access(capture_handle,
                                                hw_params,
                                                SND_PCM_ACCESS_RW_INTERLEAVED),
                   device,
                   "cannot set access type",
                   capture_handle,
                   hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_format(capture_handle, hw_params, format),
        device,
        "cannot set sample format",
        capture_handle,
        hw_params);

    check_alsa_err(
        snd_pcm_hw_params_set_rate_near(capture_handle, hw_params, rate, 0),
        device,
        "cannot set sample rate",
        capture_handle,
        hw_params);

    check_alsa_err(snd_pcm_hw_params_set_channels(capture_handle, hw_params, 1),
                   device,
                   "cannot set channel count",
                   capture_handle,
                   hw_params);

    // 4. Apply Params and Prepare
    check_alsa_err(snd_pcm_hw_params(capture_handle, hw_params),
                   device,
                   "cannot set parameters",
                   capture_handle,
                   hw_params);

    snd_pcm_hw_params_free(hw_params);

    check_alsa_err(snd_pcm_prepare(capture_handle),
                   device,
                   "cannot prepare audio interface",
                   capture_handle,
                   NULL);

    // Minor logic check
    if (*rate != requestedRate)
    {
        (void)fprintf(stderr,
                      "Requested audiorate %u unavailable, using %u\n",
                      requestedRate,
                      *rate);
    }

    return capture_handle;
}

void readBufferRaw(snd_pcm_t* capture_handle,
                   char* buffer,
                   struct myarr* data_in)
{
    unsigned char lsb;
    signed char msb;
    long err = snd_pcm_readi(capture_handle, buffer, data_in->ArrayLength);
    if (err != (long)data_in->ArrayLength)
    {
        (void)fprintf(stderr,
                      "read from audio interface failed %ld (%s)\n",
                      err,
                      snd_strerror((int)err));
        exit(READ_FAILED);
    }
    int overflow = 0;
    for (unsigned int index = 0; index < 2 * data_in->ArrayLength; index += 2)
    {
        msb = (signed char)buffer[index + 1];
        lsb = *(buffer + index);
        data_in->arr[index / 2] = (msb << BITS_IN_BYTE) | lsb;
        overflow += (data_in->arr[index / 2] == SHRT_MAX);
        overflow += (data_in->arr[index / 2] == SHRT_MIN);
    }
    if (overflow > 1)
    {
        (void)fprintf(stderr, "%d audio 16bit overflows\n", overflow);
    }
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative)
{
    if (ArrayLength == 0)
    {
        return 0;
    }

    // --- Format assumptions: S16_LE mono ---
    const unsigned int BYTES_PER_SAMPLE = 2; // 16-bit
    const unsigned int CHANNELS = 1;         // mono
    const unsigned int BYTES_PER_FRAME = BYTES_PER_SAMPLE * CHANNELS;

    // --- Read exactly ArrayLength frames, handling xruns and short reads ---
    snd_pcm_uframes_t remaining = (snd_pcm_uframes_t)ArrayLength;
    char* write_ptr = buffer;
    while (remaining > 0)
    {
        snd_pcm_sframes_t r =
            snd_pcm_readi(capture_handle, write_ptr, remaining);

        if (r == -EAGAIN)
        {
            (void)fprintf(stderr,
                          "EAGAIN%d %s %s\n",
                          __LINE__,
                          __func__,
                          __FILE__);
            // Non-blocking mode: try again
            continue;
        }
        else if (r == -EPIPE || r == -ESTRPIPE)
        {
            // Overrun (capture) or suspend: recover and continue
            fprintf(stderr,
                    "ALSA: overrun/suspend detected (%s)\n",
                    snd_strerror((int)r));
            int rr = snd_pcm_recover(capture_handle, (int)r, /*silent=*/1);
            if (rr < 0)
            {
                fprintf(stderr, "ALSA: recover failed: %s\n", snd_strerror(rr));
                return rr;
            }
            // After recover, try reading again
            continue;
        }
        else if (r < 0)
        {
            // Other fatal error
            fprintf(stderr,
                    "ALSA: snd_pcm_readi failed: %s\n",
                    snd_strerror((int)r));
            return (int)r;
        }

        // Successful short/long read
        write_ptr += (size_t)r * BYTES_PER_FRAME;
        remaining -= (snd_pcm_uframes_t)r;
    }

    // --- Post-read ALSA status snapshot: XRUN/overrange reporting ---
    {
        snd_pcm_status_t* status;
        snd_pcm_status_alloca(&status);
        if (snd_pcm_status(capture_handle, status) == 0)
        {
            if (snd_pcm_status_get_state(status) == SND_PCM_STATE_XRUN)
            {
                fprintf(stderr,
                        "ALSA: XRUN detected in status (post-read). "
                        "Recovering...\n");
                (void)snd_pcm_prepare(capture_handle);
            }
            unsigned long over = snd_pcm_status_get_overrange(status);
            if (over > 0)
            {
                fprintf(stderr,
                        "ALSA: ADC overrange reported %lu time(s)\n",
                        over);
            }
        }
    }

    // --- Timestamp-based discontinuity detection (no signature change) ---
    // We query the current hardware rate/period once and then use
    // snd_pcm_htimestamp() to check for unexpected gaps.
    {
        static int rate_inited = 0;
        static unsigned int sample_rate = 0;
        static int period_inited = 0;
        static snd_pcm_uframes_t period_size = 0;
        static int have_prev_ts = 0;
        static snd_htimestamp_t prev_ts = {0};

        if (!rate_inited || !period_inited)
        {
            snd_pcm_hw_params_t* hw;
            snd_pcm_hw_params_alloca(&hw);
            if (snd_pcm_hw_params_current(capture_handle, hw) == 0)
            {
                int dir = 0;
                (void)snd_pcm_hw_params_get_rate(hw, &sample_rate, &dir);
                (void)snd_pcm_hw_params_get_period_size(hw, &period_size, &dir);
            }
            rate_inited = period_inited = 1;
        }

        if (sample_rate > 0)
        {
            snd_htimestamp_t ts;
            snd_pcm_uframes_t avail_hw = 0;
            if (snd_pcm_htimestamp(capture_handle, &avail_hw, &ts) == 0)
            {
                if (have_prev_ts)
                {
                    double dt = (double)(ts.tv_sec - prev_ts.tv_sec) +
                                (double)(ts.tv_nsec - prev_ts.tv_nsec) * 1e-9;
                    double expected_frames = dt * (double)sample_rate;
                    double actual_frames = (double)ArrayLength;

                    // Tolerance: one ALSA period if known, else ~1 ms of audio
                    double tol_frames = (period_size > 0)
                                            ? (double)period_size
                                            : (double)sample_rate * 0.001;

                    double gap = expected_frames - actual_frames;
                    if (gap > tol_frames)
                    {
                        fprintf(stderr,
                                "ALSA: discontinuity suspected: ts advanced by "
                                "~%.0f frames, "
                                "but read %.0f (gap ~%.0f, tol ~%.0f)\n",
                                expected_frames,
                                actual_frames,
                                gap,
                                tol_frames);
                    }
                }
                prev_ts = ts;
                have_prev_ts = 1;
            }
        }
    }

    // --- Convert bytes -> int16 samples; detect clipping; compute derivative
    // ---
    int clip_overflow = 0;

    // First store samples into derivative[] (int) then reuse for |x[n]-x[n+1]|
    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        const uint8_t lsb = (uint8_t)buffer[2 * i + 0];
        const int8_t msb = (int8_t)buffer[2 * i + 1];
        const int16_t sample = (int16_t)(((int)msb << 8) | lsb); // S16_LE

        derivative[i] = (int)sample;

        // Clipping check (ADC full-scale, not XRUN)
        if (sample == INT16_MAX || sample == INT16_MIN)
        {
            clip_overflow++;
        }
    }

    if (clip_overflow > 1)
    {
        fprintf(stderr, "%d audio 16-bit clipping event(s)\n", clip_overflow);
    }

    // Compute derivative: |x[n] - x[n+1]|; last sample = 0
    for (unsigned int i = 0; i + 1 < ArrayLength; ++i)
    {
        int d = derivative[i] - derivative[i + 1];
        derivative[i] = (d < 0) ? -d : d;
    }
    derivative[ArrayLength - 1] = 0;

    // Success: exactly ArrayLength frames were read
    return (int)ArrayLength;
}

int readBufferOrFile(int* derivative,
                     unsigned int ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx)
{
    int ret = READ_FAILED;

    if (fpInput)
    {
        char* line = NULL;
        size_t len = 0;
        unsigned int index = 0;

        // Read entire lines until we have ArrayLength numbers
        while (index < ArrayLength && getline(&line, &len, fpInput) != -1)
        {
            char* ptr = line;
            char* endptr;

            while (index < ArrayLength)
            {
                errno = 0;
                // Use strtol for %d equivalent; use strtod for floating point
                long val = strtol(ptr, &endptr, DECIMAL);

                // If ptr == endptr, no more numbers were found on this line
                if (ptr == endptr)
                {
                    break;
                }

                // Error checking (optional but recommended)
                if (errno == ERANGE)
                { /* Handle overflow */
                    return INPUT_OVERFLOW;
                }

                derivative[index++] = (int)val;
                ptr = endptr; // Advance to the rest of the string
            }
        }

        free(line); // getline allocates memory that must be freed

        if (index < ArrayLength)
        {
            return INPUT_FILE_ERROR;
        }
        ret = (int)ArrayLength;

        // Perform derivative calculation
        for (unsigned int k = 0; k < ArrayLength - 1; k++)
        {
            derivative[k] = abs(derivative[k] - derivative[k + 1]);
        }
        derivative[ArrayLength - 1] = 0;
    }
    else
    {
        // ret = readBuffer(capture_handle, ArrayLength, buffer, derivative);
        //  change to new read.
        const int POLL_TIMEOUT_MS = 2000;
        struct myarr* filled = capture_next_block(ctx, POLL_TIMEOUT_MS);
        if (!filled)
        {
            (void)fprintf(stderr, "capture_next_block failed; stopping\n");
            return ret;
        }
        for (unsigned int k = 0; k < filled->ArrayLength - 1; k++)
        {
            derivative[k] = abs(filled->arr[k] - filled->arr[k + 1]);
        }
        derivative[ArrayLength - 1] = 0;
        ret = (int)ArrayLength;
    }
    return ret;
}

// Get data from audio capture
int getData(FILE* fpInput, struct myarr derivative, CaptureCtx* ctx)
{
    int err =
        readBufferOrFile(derivative.arr, derivative.ArrayLength, fpInput, ctx);
    if (err == INPUT_FILE_ERROR)
    {
        (void)fprintf(stderr,
                      "Could not read integer from inputfile or audio\n");
    }
    return err;
}

/* -------------------- Helpers -------------------- */

int build_alsa_pollfds(snd_pcm_t* h, struct pollfd** fds_out, nfds_t* nfds_out)
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

int make_timerfd_ms(unsigned initial_ms, unsigned interval_ms)
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

/* -------------------- API: setup / next_block / teardown --------------------
 */

/**
 * capture_setup
 * - Accepts an already-opened/initialized ALSA handle (from initAudio).
 * - Computes ArrayLength = rate * SECS_HOUR * 2 / bph
 * - Queries period_size and builds poll descriptors.
 */
int capture_setup(CaptureCtx* ctx,
                  snd_pcm_t* cap,
                  unsigned int rate,
                  unsigned int bph,
                  snd_pcm_format_t fmt /* expect SND_PCM_FORMAT_S16_LE */)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->cap = cap;
    ctx->rate = rate;
    ctx->tfd = -1;

    // Geometry
    ctx->ArrayLength =
        rate * SECS_HOUR * 2 / bph; // e.g., ~16000 @ 48k/21600bph
    const unsigned bytes_per_sample =
        (unsigned)snd_pcm_format_width(fmt) / 8; // 2
    const unsigned channels = 1; // you configure mono in initAudio
    ctx->bytes_per_frame = bytes_per_sample * channels;

    // Query ALSA buffer/period sizes (read per period)
    if (snd_pcm_get_params(cap, &ctx->buffer_size, &ctx->period_size) < 0)
    {
        fprintf(
            stderr,
            "ALSA: snd_pcm_get_params failed; defaulting period_size=1024\n");
        ctx->period_size = 1024;
    }
    if (ctx->period_size == 0)
        ctx->period_size = 1024;
    if (ctx->period_size > ctx->ArrayLength)
        ctx->period_size = ctx->ArrayLength;

    // Buffers
    ctx->block_buf =
        (char*)malloc((size_t)ctx->ArrayLength * ctx->bytes_per_frame);
    if (!ctx->block_buf)
    {
        fprintf(stderr, "alloc failed for block_buf\n");
        return -1;
    }
    ctx->rawread = makemyarr(ctx->ArrayLength);
    if (!ctx->rawread)
    {
        fprintf(stderr, "alloc failed for rawread\n");
        free(ctx->block_buf);
        ctx->block_buf = NULL;
        return -1;
    }

    // Warm-up reads to settle the pipeline
    readBufferRaw(ctx->cap, ctx->block_buf, ctx->rawread);
    readBufferRaw(ctx->cap, ctx->block_buf, ctx->rawread);

    // Build ALSA poll fds
    if (build_alsa_pollfds(ctx->cap, &ctx->fds, &ctx->alsa_nfds) < 0)
    {
        fprintf(stderr, "failed to build ALSA pollfds\n");
        freemyarr(ctx->rawread);
        ctx->rawread = NULL;
        free(ctx->block_buf);
        ctx->block_buf = NULL;
        return -1;
    }

    // Combine fds (ALSA + optional timerfd)
    ctx->nfds = ctx->alsa_nfds;
    if (ctx->nfds == ctx->alsa_nfds)
    {
        // ALSA-only; ctx->fds is already set
        return 0;
    }

    struct pollfd* combined =
        (struct pollfd*)calloc(ctx->nfds, sizeof(*combined));
    if (!combined)
    {
        fprintf(stderr, "alloc failed for combined fds\n");
        free(ctx->fds);
        ctx->fds = NULL;
        ctx->alsa_nfds = 0;
        freemyarr(ctx->rawread);
        ctx->rawread = NULL;
        free(ctx->block_buf);
        ctx->block_buf = NULL;
        return -1;
    }
    for (nfds_t i = 0; i < ctx->alsa_nfds; ++i)
        combined[i] = ctx->fds[i];
    combined[ctx->alsa_nfds].fd = ctx->tfd;
    combined[ctx->alsa_nfds].events = POLLIN;

    free(ctx->fds);
    ctx->fds = combined;

    return 0;
}

/**
 * capture_next_block
 * - Single-threaded event loop step: polls ALSA (and timerfd if enabled),
 *   reads per period until a full block is collected, then returns
 * ctx->rawread.
 * - Returns NULL on unrecoverable error.
 * - Does NOT write to file — caller owns persistence/processing.
 */
struct myarr* capture_next_block(CaptureCtx* ctx, int poll_timeout_ms)
{
    ctx->frames_collected = 0;
    ctx->byte_off = 0;
    ctx->samp_off = 0;

    for (;;)
    {
        int pr = poll(ctx->fds, ctx->nfds, poll_timeout_ms);
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
            return NULL;
        }

        // 1) Handle ALSA events
        unsigned short revents = 0;
        if (snd_pcm_poll_descriptors_revents(ctx->cap,
                                             ctx->fds,
                                             (unsigned)ctx->alsa_nfds,
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
                if (snd_pcm_recover(ctx->cap, -EPIPE, 1) < 0)
                {
                    fprintf(stderr, "ALSA: recover failed\n");
                    return NULL;
                }
                struct pollfd* new_alsa = NULL;
                nfds_t new_n = 0;
                if (build_alsa_pollfds(ctx->cap, &new_alsa, &new_n) < 0)
                {
                    fprintf(stderr, "ALSA: rebuild pollfds failed\n");
                    return NULL;
                }
                nfds_t new_total = new_n;
                struct pollfd* tmp =
                    (struct pollfd*)calloc(new_total, sizeof(*tmp));
                if (!tmp)
                {
                    free(new_alsa);
                    return NULL;
                }
                for (nfds_t i = 0; i < new_n; ++i)
                    tmp[i] = new_alsa[i];
                free(new_alsa);
                free(ctx->fds);
                ctx->fds = tmp;
                ctx->alsa_nfds = new_n;
                ctx->nfds = new_total;
            }
            if (revents & POLLIN)
            {
                for (;;)
                {
                    unsigned frames_left_in_block =
                        ctx->ArrayLength - ctx->frames_collected;
                    if (frames_left_in_block == 0)
                        break;

                    snd_pcm_sframes_t avail = snd_pcm_avail_update(ctx->cap);
                    if (avail < 0)
                    {
                        int rr = snd_pcm_recover(ctx->cap, (int)avail, 1);
                        if (rr < 0)
                        {
                            fprintf(stderr,
                                    "ALSA: avail_update recover failed: %s\n",
                                    snd_strerror(rr));
                            return NULL;
                        }
                        // After recover, let poll() fetch the next event
                        break;
                    }

                    if ((snd_pcm_uframes_t)avail < ctx->period_size)
                    {
                        // Not enough frames yet to read another period — go
                        // back to poll()
                        break;
                    }

                    // Read up to what's available, but not more than the
                    // remainder of the block.
                    unsigned chunk =
                        (unsigned)((avail <
                                    (snd_pcm_sframes_t)frames_left_in_block)
                                       ? avail
                                       : (snd_pcm_sframes_t)
                                             frames_left_in_block);
                    // Optionally cap chunk to multiples of period_size to keep
                    // cadence:
                    if (chunk > ctx->period_size)
                    {
                        chunk -= (chunk % (unsigned)ctx->period_size);
                        if (chunk == 0)
                            chunk = (unsigned)ctx->period_size;
                    }

                    int got = readBuffer(ctx->cap,
                                         chunk,
                                         ctx->block_buf + ctx->byte_off,
                                         ctx->rawread->arr + ctx->samp_off);
                    if (got < 0)
                    {
                        int rr = snd_pcm_recover(ctx->cap, got, 1);
                        if (rr < 0)
                        {
                            fprintf(stderr,
                                    "ALSA: readBuffer failed & recover failed: "
                                    "%s\n",
                                    snd_strerror(rr));
                            return NULL;
                        }
                        // After recover, descriptors may change; rebuild
                        struct pollfd* new_alsa = NULL;
                        nfds_t new_n = 0;
                        if (build_alsa_pollfds(ctx->cap, &new_alsa, &new_n) ==
                            0)
                        {
                            nfds_t new_total = new_n;
                            struct pollfd* tmp =
                                (struct pollfd*)calloc(new_total, sizeof(*tmp));
                            if (tmp)
                            {
                                for (nfds_t i = 0; i < new_n; ++i)
                                    tmp[i] = new_alsa[i];
                                free(ctx->fds);
                                ctx->fds = tmp;
                                ctx->alsa_nfds = new_n;
                                ctx->nfds = new_total;
                            }
                            free(new_alsa);
                        }
                        // Go back to poll() and continue
                        break;
                    }

                    // Advance offsets
                    ctx->byte_off += (size_t)got * ctx->bytes_per_frame;
                    ctx->samp_off += (size_t)got;
                    ctx->frames_collected += (unsigned)got;

                    if (ctx->frames_collected >= ctx->ArrayLength)
                    {
                        // Completed a full block
                        for (nfds_t i = 0; i < ctx->nfds; ++i)
                            ctx->fds[i].revents = 0;
                        return ctx->rawread;
                    }
                }
            }
        }

        // Clear revents
        for (nfds_t i = 0; i < ctx->nfds; ++i)
        {
            ctx->fds[i].revents = 0;
        }
    }
}

void capture_teardown(CaptureCtx* ctx)
{
    if (!ctx)
        return;
    if (ctx->fds)
        free(ctx->fds);
    if (ctx->tfd >= 0)
        close(ctx->tfd);
    if (ctx->block_buf)
        free(ctx->block_buf);
    if (ctx->rawread)
        freemyarr(ctx->rawread);
    memset(ctx, 0, sizeof(*ctx));
}
