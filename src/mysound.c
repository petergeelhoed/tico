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
    // Choose period target (e.g., 1024) and larger buffer (e.g., 8192 = 8×)
    snd_pcm_uframes_t period_size = 1024;
    snd_pcm_uframes_t buffer_size = period_size * 12;

    snd_pcm_hw_params_set_period_size_near(capture_handle,
                                           hw_params,
                                           &period_size,
                                           0);
    snd_pcm_hw_params_set_buffer_size_near(capture_handle,
                                           hw_params,
                                           &buffer_size);

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
                                            ? 1.5 * (double)period_size
                                            : (double)sample_rate * 0.003;

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

        // Report any continuity warning for this block
        if (ctx->last_block_had_gap)
        {
            double ms = 1000.0 * ctx->last_deficit_frames / (double)ctx->rate;
            fprintf(stderr,
                    "CONTINUITY (block): gap ~%.0f frames (~%.2f ms)\n",
                    ctx->last_deficit_frames,
                    ms);
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

static void rotate_int_block(int* a, unsigned n, int shift)
{
    if (n == 0)
        return;
    // Normalize to [-n+1, n-1]
    int s = shift % (int)n;
    if (s == 0)
        return;
    if (s < 0)
        s += (int)n; // convert negative to positive left-rotate amount

    // Rotate-left by s: [s..n-1] + [0..s-1]
    // Use temporary buffer for simplicity (n is ~16000 -> OK)
    int* tmp = (int*)malloc(n * sizeof(int));
    if (!tmp)
        return; // fail silently (no snap) if OOM
    memcpy(tmp, a + s, (n - (unsigned)s) * sizeof(int));
    memcpy(tmp + (n - (unsigned)s), a, (unsigned)s * sizeof(int));
    memcpy(a, tmp, n * sizeof(int));
    free(tmp);
}

/**
 * capture_setup
 * - Accepts an already-opened/initialized ALSA handle (from initAudio).
 * - Computes ArrayLength = rate * SECS_HOUR * 2 / bph
 * - Queries period_size and builds poll descriptors.
 * - Allocates buffers and creates a 1Hz timerfd if countdown > 0.
 */
int capture_setup(CaptureCtx* ctx,
                  snd_pcm_t* cap,
                  unsigned int rate,
                  unsigned int bph,
                  unsigned int countdown_seconds,
                  snd_pcm_format_t fmt /* expect SND_PCM_FORMAT_S16_LE */)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->cap = cap;
    ctx->rate = rate;
    ctx->tfd = -1;
    ctx->countdown = countdown_seconds;

    ctx->frames_total = 0;
    ctx->have_t0 = 0;
    ctx->rate_est = (double)ctx->rate; // start from nominal
    ctx->rate_alpha = 0.02;            // ~2% update per block

    ctx->last_boundary_lag = 0;
    ctx->apply_boundary_snap = 1; // turn on "snap" by default

    snd_htimestamp_t ts;
    snd_pcm_uframes_t dummy = 0;
    if (snd_pcm_htimestamp(ctx->cap, &dummy, &ts) == 0)
    {
        ctx->t0_pcm = ts;
        ctx->have_t0 = 1;
    }

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

    // After warm-up reads
    ctx->frames_total = 0;
    ctx->have_t0 = 0;
    ctx->last_deficit_frames = 0.0;
    ctx->last_block_had_gap = 0;

    // Set t0 from the device's hardware timestamp
    {
        snd_htimestamp_t ts1;
        snd_pcm_uframes_t avail_dummy = 0;
        if (snd_pcm_htimestamp(ctx->cap, &avail_dummy, &ts1) == 0)
        {
            ctx->t0_pcm = ts1;
            ctx->have_t0 = 1;
        }
    }

    // Optional: allocate a tail window for correlation (e.g., 512 samples)
    ctx->tail_len =
        (ctx->period_size >= 512) ? 512u : (unsigned)ctx->period_size;
    ctx->prev_tail = (int*)calloc(ctx->tail_len, sizeof(int));

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

    // Add timerfd if a countdown is requested
    if (ctx->countdown > 0)
    {
        ctx->tfd = make_timerfd_ms(/*initial*/ 1000, /*interval*/ 1000);
        if (ctx->tfd < 0)
        {
            fprintf(stderr,
                    "warning: timerfd not available; countdown disabled\n");
            ctx->countdown = 0;
        }
    }

    // Combine fds (ALSA + optional timerfd)
    ctx->nfds = ctx->alsa_nfds + ((ctx->countdown > 0) ? 1 : 0);
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

struct myarr* capture_next_block(CaptureCtx* ctx, int poll_timeout_ms)
{
    // Per-block reset
    ctx->frames_collected = 0;
    ctx->byte_off = 0;
    ctx->samp_off = 0;
    ctx->last_block_had_gap = 0;
    ctx->last_deficit_frames = 0.0;
    ctx->last_boundary_lag = 0;

    // Tunables
    const double MIN_CORR_MS = 2.0; // ignore corr shifts < 2 ms
    const double GAP_MIN_MS = 3.0; // block gap must exceed 3 ms and 1.5×period
    const double GAP_MULT_PER = 1.5; // block gap must exceed 1.5×period
    const double NCC_MIN =
        0.80; // zero-lag NCC must be below this to consider shift

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

        // ALSA events
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
                nfds_t new_total = new_n + ((ctx->countdown > 0) ? 1 : 0);
                struct pollfd* tmp =
                    (struct pollfd*)calloc(new_total, sizeof(*tmp));
                if (!tmp)
                {
                    free(new_alsa);
                    return NULL;
                }
                for (nfds_t i = 0; i < new_n; ++i)
                    tmp[i] = new_alsa[i];
                if (ctx->countdown > 0)
                    tmp[new_n] = ctx->fds[ctx->alsa_nfds];
                free(new_alsa);
                free(ctx->fds);
                ctx->fds = tmp;
                ctx->alsa_nfds = new_n;
                ctx->nfds = new_total;
            }

            if (revents & POLLIN)
            {
                // Drain as many whole periods as fit the block remainder
                for (;;)
                {
                    const unsigned frames_left_in_block =
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
                        break; // back to poll
                    }
                    if ((snd_pcm_uframes_t)avail < ctx->period_size)
                    {
                        break; // wait for next period
                    }

                    unsigned chunk =
                        (unsigned)((avail <
                                    (snd_pcm_sframes_t)frames_left_in_block)
                                       ? avail
                                       : (snd_pcm_sframes_t)
                                             frames_left_in_block);
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
                        struct pollfd* new_alsa = NULL;
                        nfds_t new_n = 0;
                        if (build_alsa_pollfds(ctx->cap, &new_alsa, &new_n) ==
                            0)
                        {
                            nfds_t new_total =
                                new_n + ((ctx->countdown > 0) ? 1 : 0);
                            struct pollfd* tmp =
                                (struct pollfd*)calloc(new_total, sizeof(*tmp));
                            if (tmp)
                            {
                                for (nfds_t i = 0; i < new_n; ++i)
                                    tmp[i] = new_alsa[i];
                                if (ctx->countdown > 0)
                                    tmp[new_n] = ctx->fds[ctx->alsa_nfds];
                                free(ctx->fds);
                                ctx->fds = tmp;
                                ctx->alsa_nfds = new_n;
                                ctx->nfds = new_total;
                            }
                            free(new_alsa);
                        }
                        break; // back to poll
                    }

                    // Advance counters
                    ctx->byte_off += (size_t)got * ctx->bytes_per_frame;
                    ctx->samp_off += (size_t)got;
                    ctx->frames_collected += (unsigned)got;
                    ctx->frames_total += (uint64_t)got;

                    if (ctx->frames_collected >= ctx->ArrayLength)
                    {
                        // ---- Block boundary ----
                        // (A) Block-level continuity (PLL-smoothed)
                        if (ctx->have_t0)
                        {
                            snd_htimestamp_t ts_now;
                            snd_pcm_uframes_t dummy2 = 0;
                            if (snd_pcm_htimestamp(ctx->cap,
                                                   &dummy2,
                                                   &ts_now) == 0)
                            {
                                const double dt =
                                    (double)(ts_now.tv_sec -
                                             ctx->t0_pcm.tv_sec) +
                                    1e-9 * (double)(ts_now.tv_nsec -
                                                    ctx->t0_pcm.tv_nsec);
                                if (dt > 0.0)
                                {
                                    const double rate_meas =
                                        (double)ctx->frames_total / dt;
                                    ctx->rate_est = (1.0 - ctx->rate_alpha) *
                                                        ctx->rate_est +
                                                    ctx->rate_alpha * rate_meas;

                                    const double expected_frames =
                                        dt * ctx->rate_est;
                                    const double deficit =
                                        expected_frames -
                                        (double)ctx->frames_total;

                                    const double ms_min_frames =
                                        (GAP_MIN_MS / 1000.0) * ctx->rate_est;
                                    const double gap_thr = fmax(
                                        GAP_MULT_PER * (double)ctx->period_size,
                                        ms_min_frames);

                                    ctx->last_deficit_frames = deficit;
                                    ctx->last_block_had_gap =
                                        (deficit > gap_thr) ? 1 : 0;

                                    if (ctx->last_block_had_gap)
                                    {
                                        const double ms =
                                            1000.0 * deficit / ctx->rate_est;
                                        fprintf(stderr,
                                                "CONTINUITY (block): deficit "
                                                "~%.0f frames (~%.2f ms)\n",
                                                deficit,
                                                ms);
                                    }
                                }
                            }
                        }
                        else
                        {
                            snd_pcm_htimestamp(ctx->cap,
                                               (snd_pcm_uframes_t[]){0},
                                               &ctx->t0_pcm);
                            ctx->have_t0 = 1;
                        }

                        // (B) Boundary correlation (normalized) + optional
                        // "snap"
                        if (ctx->prev_tail && ctx->tail_len > 0)
                        {
                            const unsigned W = ctx->tail_len;
                            const int* tail = ctx->prev_tail;
                            const int* head = ctx->rawread->arr;

                            // Zero-lag NCC
                            long long sumt2 = 0, sumh2_0 = 0, dot0 = 0;
                            for (unsigned i = 0; i < W; ++i)
                            {
                                sumt2 += (long long)tail[i] * tail[i];
                                sumh2_0 += (long long)head[i] * head[i];
                                dot0 += (long long)tail[i] * head[i];
                            }
                            const double normt = sqrt((double)sumt2);
                            const double normh0 = sqrt((double)sumh2_0);
                            const double ncc0 =
                                (normt > 0.0 && normh0 > 0.0)
                                    ? ((double)dot0 / (normt * normh0))
                                    : 0.0;

                            // Best lag NCC over ±period/2
                            const int maxlag = (int)(ctx->period_size / 2);
                            double best_ncc = -1.0;
                            int best_lag = 0;
                            for (int lag = -maxlag; lag <= maxlag; ++lag)
                            {
                                long long dot = 0, sumh2 = 0;
                                for (unsigned i = 0; i < W; ++i)
                                {
                                    const int hi = (int)i + lag;
                                    if (hi < 0 || hi >= (int)W)
                                        continue;
                                    dot += (long long)tail[i] * head[hi];
                                    sumh2 += (long long)head[hi] * head[hi];
                                }
                                const double normh = sqrt((double)sumh2);
                                const double ncc =
                                    (normt > 0.0 && normh > 0.0)
                                        ? ((double)dot / (normt * normh))
                                        : -1.0;
                                if (ncc > best_ncc)
                                {
                                    best_ncc = ncc;
                                    best_lag = lag;
                                }
                            }

                            // Convert min lag in ms to frames, and combine with
                            // period fraction
                            double lag_ms_frames_d =
                                (MIN_CORR_MS / 1000.0) * ctx->rate_est;
                            if (lag_ms_frames_d < 0.0)
                                lag_ms_frames_d = 0.0;
                            if (lag_ms_frames_d > (double)INT_MAX)
                                lag_ms_frames_d = (double)INT_MAX;
                            const int lag_min_frames_ms = (int)lag_ms_frames_d;
                            const int lag_min_frames_per =
                                (int)(ctx->period_size / 8);
                            const int LAG_MIN =
                                (lag_min_frames_ms > lag_min_frames_per)
                                    ? lag_min_frames_ms
                                    : lag_min_frames_per;

                            // Decide: only if zero-lag similarity is low AND
                            // lag is significant
                            if (ncc0 < NCC_MIN && abs(best_lag) >= LAG_MIN)
                            {
                                ctx->last_boundary_lag = best_lag;

                                if (ctx->apply_boundary_snap)
                                {
                                    // Rotate current block by -best_lag so head
                                    // aligns to prev tail
                                    rotate_int_block(ctx->rawread->arr,
                                                     ctx->ArrayLength,
                                                     -best_lag);
                                    // Also rotate PCM bytes if you need them
                                    // phase-aligned elsewhere:
                                    // rotate_bytes(ctx->block_buf,
                                    // ctx->ArrayLength, ctx->bytes_per_frame,
                                    // -best_lag);

                                    // Update prev_tail from the *aligned*
                                    // block's tail
                                    memcpy(ctx->prev_tail,
                                           ctx->rawread->arr +
                                               (ctx->ArrayLength - W),
                                           W * sizeof(int));
                                }
                                else
                                {
                                    // No snap: still update prev_tail from
                                    // current block
                                    memcpy(ctx->prev_tail,
                                           ctx->rawread->arr +
                                               (ctx->ArrayLength - W),
                                           W * sizeof(int));
                                    // And log once
                                    const double ms = 1000.0 *
                                                      (double)best_lag /
                                                      ctx->rate_est;
                                    fprintf(stderr,
                                            "CONTINUITY (corr): boundary "
                                            "offset ~%d frames (~%.2f ms), "
                                            "ncc0=%.2f best=%.2f\n",
                                            best_lag,
                                            ms,
                                            ncc0,
                                            best_ncc);
                                }
                            }
                            else
                            {
                                // Already aligned; update prev_tail normally
                                memcpy(
                                    ctx->prev_tail,
                                    ctx->rawread->arr + (ctx->ArrayLength - W),
                                    W * sizeof(int));
                            }
                        }

                        // Return the (possibly snapped) block
                        for (nfds_t i = 0; i < ctx->nfds; ++i)
                            ctx->fds[i].revents = 0;
                        return ctx->rawread;
                    }
                } // drain loop
            }
        }

        // timerfd countdown (optional)
        if (ctx->countdown > 0)
        {
            struct pollfd* t = &ctx->fds[ctx->alsa_nfds];
            if (t->revents & POLLIN)
            {
                uint64_t expirations = 0;
                ssize_t rd = read(t->fd, &expirations, sizeof(expirations));
                if (rd == (ssize_t)sizeof(expirations))
                {
                    while (expirations-- && ctx->countdown > 0)
                    {
                        --ctx->countdown;
                        printf("%u\n", ctx->countdown);
                        fflush(stdout);
                    }
                    if (ctx->countdown == 0)
                    {
                        struct itimerspec its = {0};
                        timerfd_settime(t->fd, 0, &its, NULL);
                        t->events = 0;
                    }
                }
            }
        }

        for (nfds_t i = 0; i < ctx->nfds; ++i)
            ctx->fds[i].revents = 0;
    }
}
