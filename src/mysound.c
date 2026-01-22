#include "mysound.h"
#include "config.h"
#include "myarr.h"
#include "mydefs.h"
#include "myfft.h"
#include "mylib.h"
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

/* Initialize from current hw_param params (safe even if not all fields are
 * available)
 */
static inline void missdet_init(MissDet* missdet, snd_pcm_t* pcm)
{
    memset(missdet, 0, sizeof(*missdet));

    snd_pcm_hw_params_t* hw_param = NULL;
    int err = snd_pcm_hw_params_malloc(&hw_param);

    if (err < 0 || !hw_param)
    {
        // handle allocation failure
        (void)fprintf(stderr,
                      "snd_pcm_hw_params_malloc failed: %s\n",
                      snd_strerror(err));
        exit(err ? err : -ENOMEM);
    }

    if (snd_pcm_hw_params_current(pcm, hw_param) == 0)
    {
        int dir = 0;
        (void)snd_pcm_hw_params_get_rate(hw_param, &missdet->sample_rate, &dir);
        (void)snd_pcm_hw_params_get_period_size(hw_param,
                                                &missdet->period_size,
                                                &dir);
    }
    snd_pcm_hw_params_free(hw_param);

    /* Build a practical tolerance:
       - If we know period_size: allow one period worth of jitter
       - Else: allow ~1 ms worth of frames
    */
    if (missdet->period_size > 0)
    {
        missdet->tol_frames = (double)missdet->period_size;
    }
    else if (missdet->sample_rate > 0)
    {
        missdet->tol_frames = MILLI * (double)missdet->sample_rate; /* ~1 ms */
    }
    else
    {
        missdet->tol_frames = FRAME_FALLBACK; /* conservative fallback */
    }
}

/* Call AFTER each successful read() that returned 'got' frames.
   Returns:
     >0  suspected missing data (gap larger than tolerance)
      0  OK (no gap)
     <0  not enough info (e.g., sample_rate unknown or timestamp unavailable)
*/
static inline int missdet_check(MissDet* missdet, snd_pcm_t* pcm, unsigned got)
{
    if (!missdet || missdet->sample_rate == 0)
    {
        return -1;
    }

    snd_htimestamp_t timestmp = {0, 0};
    snd_pcm_uframes_t avail_hw = 0; /* not used, but required by API */
    if (snd_pcm_htimestamp(pcm, &avail_hw, &timestmp) != 0)
    {
        return -1;
    }

    if (missdet->inited)
    {
        /* Compute elapsed wallclock time between reads */
        double delta_t =
            (double)(timestmp.tv_sec - missdet->prev_ts.tv_sec) +
            (double)(timestmp.tv_nsec - missdet->prev_ts.tv_nsec) * NANO;
        if (delta_t > 0.0)
        {
            double expected = delta_t * (double)missdet->sample_rate;
            double actual = (double)got;
            double gap = expected - actual;

            if (gap > missdet->tol_frames)
            {
                /* Missing data suspected */
                missdet->prev_ts = timestmp;
                return 1;
            }
        }
    }

    missdet->prev_ts = timestmp;
    missdet->inited = 1;
    return 0;
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

static int read_exact_frames(snd_pcm_t* pcm,
                             char* buf,
                             snd_pcm_uframes_t frames_to_read,
                             unsigned int bytes_per_frame)
{
    snd_pcm_uframes_t remaining = frames_to_read;
    char* write_ptr = buf;

    while (remaining > 0)
    {
        snd_pcm_sframes_t got = snd_pcm_readi(pcm, write_ptr, remaining);

        if (got == -EAGAIN)
        {
            // Non-blocking: try again
            continue;
        }
        if (got == -EPIPE || got == -ESTRPIPE)
        {
            // Overrun / suspend: recover
            int retval = snd_pcm_recover(pcm, (int)got, /*silent=*/1);
            if (retval < 0)
            {
                (void)fprintf(stderr,
                              "ALSA recover failed: %s\n",
                              snd_strerror(retval));
                return retval;
            }
            continue;
        }
        if (got < 0)
        {
            // Other error
            (void)fprintf(stderr,
                          "ALSA read failed: %s\n",
                          snd_strerror((int)got));
            return (int)got;
        }

        write_ptr += (size_t)got * bytes_per_frame;
        remaining -= (snd_pcm_uframes_t)got;
    }
    return 0; // Success
}

int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative)
{
    if (!capture_handle || !buffer || !derivative)
    {
        return -EINVAL;
    }
    if (ArrayLength == 0)
    {
        return 0;
    }

    // S16_LE mono
    const unsigned int BYTES_PER_SAMPLE = 2;
    const unsigned int CHANNELS = 1;
    const unsigned int BYTES_PER_FRAME = BYTES_PER_SAMPLE * CHANNELS;

    // Read exactly ArrayLength frames
    int retval = read_exact_frames(capture_handle,
                                   buffer,
                                   (snd_pcm_uframes_t)ArrayLength,
                                   BYTES_PER_FRAME);
    if (retval < 0)
    {
        return retval;
    }

    // Convert to samples, detect clipping, store in derivative[]
    int clip_count = 0;
    for (unsigned int i = 0; i < ArrayLength; ++i)
    {
        const uint8_t lsb = (uint8_t)buffer[2 * i + 0];
        const int8_t msb = (int8_t)buffer[2 * i + 1];
        const int16_t samp = (int16_t)(((int)msb << 8) | lsb); // S16_LE

        derivative[i] = (int)samp;

        if (samp == INT16_MAX || samp == INT16_MIN)
        {
            ++clip_count;
        }
    }

    if (clip_count > 1)
    {
        (void)fprintf(stderr,
                      "%d audio 16-bit clipping event(s)\n",
                      clip_count);
    }

    // Compute |x[n] - x[n+1]|, last = 0
    for (unsigned int i = 0; i + 1 < ArrayLength; ++i)
    {
        derivative[i] = abs(derivative[i] - derivative[i + 1]);
    }
    derivative[ArrayLength - 1] = 0;

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

int build_alsa_pollfds(snd_pcm_t* handle,
                       struct pollfd** fds_out,
                       nfds_t* nfds_out)
{
    int count = snd_pcm_poll_descriptors_count(handle);
    if (count <= 0)
    {
        (void)fprintf(stderr,
                      "ALSA: invalid poll descriptors count: %d\n",
                      count);
        return -1;
    }
    struct pollfd* fds = (struct pollfd*)calloc((size_t)count, sizeof(*fds));
    if (!fds)
    {
        (void)fprintf(stderr, "alloc failed for alsa fds\n");
        return -1;
    }
    int err = snd_pcm_poll_descriptors(handle, fds, (unsigned)count);
    if (err < 0)
    {
        (void)fprintf(stderr,
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
    its.it_value.tv_sec = initial_ms / THOUSAND;
    its.it_value.tv_nsec = (long)(initial_ms % THOUSAND) * MILLION;
    its.it_interval.tv_sec = interval_ms / THOUSAND;
    its.it_interval.tv_nsec = (long)(interval_ms % THOUSAND) * MILLION;

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
        (unsigned)snd_pcm_format_width(fmt) / BITS_IN_BYTE; // 2
    const unsigned channels = 1; // you configure mono in initAudio
    ctx->bytes_per_frame = bytes_per_sample * channels;

    // Query ALSA buffer/period sizes (read per period)
    if (snd_pcm_get_params(cap, &ctx->buffer_size, &ctx->period_size) < 0)
    {
        (void)fprintf(
            stderr,
            "ALSA: snd_pcm_get_params failed; defaulting period_size=%d\n",
            DEFAULT_PERIOD);
        ctx->period_size = DEFAULT_PERIOD;
    }
    if (ctx->period_size == 0)
    {
        ctx->period_size = DEFAULT_PERIOD;
    }
    if (ctx->period_size > ctx->ArrayLength)
    {
        ctx->period_size = ctx->ArrayLength;
    }

    // Buffers
    ctx->block_buf =
        (char*)malloc((size_t)ctx->ArrayLength * ctx->bytes_per_frame);
    if (!ctx->block_buf)
    {
        (void)fprintf(stderr, "alloc failed for block_buf\n");
        return -1;
    }
    ctx->rawread = makemyarr(ctx->ArrayLength);
    if (!ctx->rawread)
    {
        (void)fprintf(stderr, "alloc failed for rawread\n");
        free(ctx->block_buf);
        ctx->block_buf = NULL;
        return -1;
    }
    missdet_init(&ctx->missdet, ctx->cap);

    // Warm-up reads to settle the pipeline
    readBufferRaw(ctx->cap, ctx->block_buf, ctx->rawread);
    readBufferRaw(ctx->cap, ctx->block_buf, ctx->rawread);

    // Build ALSA poll fds
    if (build_alsa_pollfds(ctx->cap, &ctx->fds, &ctx->alsa_nfds) < 0)
    {
        (void)fprintf(stderr, "failed to build ALSA pollfds\n");
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
        (void)fprintf(stderr, "alloc failed for combined fds\n");
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
    {
        combined[i] = ctx->fds[i];
    }
    combined[ctx->alsa_nfds].fd = ctx->tfd;
    combined[ctx->alsa_nfds].events = POLLIN;

    free(ctx->fds);
    ctx->fds = combined;

    return 0;
}

void capture_teardown(CaptureCtx* ctx)
{
    if (!ctx)
    {
        return;
    }
    if (ctx->fds)
    {
        free(ctx->fds);
    }
    if (ctx->tfd >= 0)
    {
        close(ctx->tfd);
    }
    if (ctx->block_buf)
    {
        free(ctx->block_buf);
    }
    if (ctx->rawread)
    {
        freemyarr(ctx->rawread);
    }
    memset(ctx, 0, sizeof(*ctx));
}

static inline void clear_revents(struct pollfd* fds, nfds_t n)
{
    for (nfds_t i = 0; i < n; ++i)
    {
        fds[i].revents = 0;
    }
}

static int rebuild_alsa_pollfds(CaptureCtx* ctx)
{
    struct pollfd* new_alsa = NULL;
    nfds_t new_n = 0;
    if (build_alsa_pollfds(ctx->cap, &new_alsa, &new_n) < 0)
    {
        (void)fprintf(stderr, "ALSA: rebuild pollfds failed\n");
        return -1;
    }

    struct pollfd* tmp = (struct pollfd*)calloc(new_n, sizeof(*tmp));
    if (!tmp)
    {
        free(new_alsa);
        return -1;
    }

    memcpy(tmp, new_alsa, new_n * sizeof(*tmp));
    free(new_alsa);

    free(ctx->fds);
    ctx->fds = tmp;
    ctx->alsa_nfds = new_n;
    ctx->nfds = new_n; // Keeping ALSA-only set for simplicity
    return 0;
}

static int recover_and_rebuild(CaptureCtx* ctx, int err)
{
    int recover = snd_pcm_recover(ctx->cap, err, /*silent=*/1);
    if (recover < 0)
    {
        (void)fprintf(stderr,
                      "ALSA: recover failed: %s\n",
                      snd_strerror(recover));
        return recover;
    }
    return rebuild_alsa_pollfds(ctx);
}

static inline void log_xrun_if_any(snd_pcm_t* cap)
{
    snd_pcm_status_t* status = NULL;
    int err = snd_pcm_status_malloc(&status);
    if (err < 0 || !status)
    {
        (void)fprintf(stderr,
                      "snd_pcm_status_malloc failed: %s\n",
                      snd_strerror(err));
        return;
    }

    if (snd_pcm_status(cap, status) == 0)
    {
        if (snd_pcm_status_get_state(status) == SND_PCM_STATE_XRUN)
        {
            (void)fprintf(stderr,
                          "MISSING DATA: XRUN state reported (post-read)\n");
            (void)snd_pcm_prepare(cap); // best effort; ignore errors here
        }
    }
    snd_pcm_status_free(status);
}

/* Read a single burst worth of frames (<= period_size, also <= remaining).
   Returns:
     1  = block complete after this read
     0  = more data needed
    <0  = fatal error
*/
static int read_one_burst(CaptureCtx* ctx)
{
    unsigned int remaining = ctx->ArrayLength - ctx->frames_collected;
    if (remaining == 0)
    {
        return 1;
    }

    unsigned int want =
        (unsigned int)((ctx->period_size > 0 &&
                        (snd_pcm_uframes_t)remaining > ctx->period_size)
                           ? ctx->period_size
                           : remaining);

    int got = readBuffer(ctx->cap,
                         want,
                         ctx->block_buf + ctx->byte_off,
                         /* reuse/target int buffer inside rawread->arr */
                         ctx->rawread->arr + ctx->samp_off);
    if (got < 0)
    {
        int recover = recover_and_rebuild(ctx, got);
        if (recover < 0)
        {
            return recover;
        }
        return 0; // recovered; let poll drive the next attempt
    }

    log_xrun_if_any(ctx->cap);

    if (got > 0)
    {
        int md_flag = missdet_check(&ctx->missdet, ctx->cap, (unsigned)got);
        if (md_flag > 0)
        {
            ctx->lost_events++;
            (void)fprintf(
                stderr,
                "MISSING DATA suspected: time advanced >> frames read "
                "(event #%lu)\n",
                ctx->lost_events);
        }
    }

    // Advance offsets
    ctx->byte_off += (size_t)got * ctx->bytes_per_frame;
    ctx->samp_off += (size_t)got;
    ctx->frames_collected += (unsigned)got;

    return (ctx->frames_collected >= ctx->ArrayLength) ? 1 : 0;
}

struct myarr* capture_next_block(CaptureCtx* ctx, int poll_timeout_ms)
{
    if (!ctx || !ctx->cap || !ctx->fds || ctx->nfds == 0 || !ctx->block_buf ||
        !ctx->rawread || ctx->bytes_per_frame == 0)
    {
        errno = EINVAL;
        return NULL;
    }

    ctx->frames_collected = 0;
    ctx->byte_off = 0;
    ctx->samp_off = 0;

    for (;;)
    {
        int poll_ret = poll(ctx->fds, ctx->nfds, poll_timeout_ms);
        if (poll_ret == 0)
        {
            // Timeout: keep polling; upper layer can decide if it wants to bail
            // out. (You could return NULL with a custom error code instead.)
            // fprintf(stderr, "poll: timeout\n");
            continue;
        }
        if (poll_ret < 0)
        {
            if (errno == EINTR)
            {
                continue; // interrupted -> retry
            }
            perror("poll");
            return NULL;
        }

        unsigned short revents = 0;
        if (snd_pcm_poll_descriptors_revents(ctx->cap,
                                             ctx->fds,
                                             (unsigned)ctx->alsa_nfds,
                                             &revents) < 0)
        {
            (void)fprintf(stderr, "ALSA: revents failed\n");
            clear_revents(ctx->fds, ctx->nfds);
            continue;
        }

        if (revents & POLLERR)
        {
            if (recover_and_rebuild(ctx, -EPIPE) < 0)
            {
                return NULL;
            }
            clear_revents(ctx->fds, ctx->nfds);
            continue;
        }

        if (revents & POLLIN)
        {
            int burst = read_one_burst(ctx);
            if (burst < 0)
            {
                // Fatal
                return NULL;
            }
            if (burst > 0)
            {
                // Completed full block
                clear_revents(ctx->fds, ctx->nfds);
                return ctx->rawread;
            }
            // else need more data; fall through to clear and poll again
        }

        clear_revents(ctx->fds, ctx->nfds);
    }
}
