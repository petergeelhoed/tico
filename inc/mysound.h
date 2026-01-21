#pragma once
#include "myarr.h"
#include <alsa/asoundlib.h> // IWYU pragma: export
#include <poll.h>
#include <stdio.h>
#include <sys/timerfd.h>

/* -------------------- Capture Context -------------------- */

typedef struct CaptureCtx
{
    /* ALSA */
    snd_pcm_t* cap;
    unsigned int rate;
    snd_pcm_uframes_t period_size;
    snd_pcm_uframes_t buffer_size;

    /* Block geometry */
    unsigned int ArrayLength;     // frames per processing block
    unsigned int bytes_per_frame; // for S16_LE mono -> 2

    /* Accumulation state for current block */
    unsigned frames_collected;
    size_t byte_off;
    size_t samp_off;

    /* Buffers */
    char* block_buf;       // raw PCM bytes (ArrayLength * bytes_per_frame)
    struct myarr* rawread; // holds int samples (ArrayLength)

    /* Poll */
    struct pollfd* fds; // ALSA fds + optional timerfd
    nfds_t nfds;        // total fds
    nfds_t alsa_nfds;   // number of ALSA fds at head of fds[]

    /* Countdown (optional, via timerfd) */
    int tfd;                // timerfd or -1
    unsigned int countdown; // seconds remaining; 0 disables
} CaptureCtx;

// NOLINTNEXTLINE(misc-include-cleaner)
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate);
int readBuffer(snd_pcm_t* capture_handle,
               unsigned int ArrayLength,
               char* buffer,
               int* derivative);
void readBufferRaw(snd_pcm_t* capture_handle,
                   char* buffer,
                   struct myarr* data_in);

int readBufferOrFile(int* derivative,
                     snd_pcm_t* capture_handle,
                     unsigned int ArrayLength,
                     char* buffer,
                     FILE* fpInput);
int getData(FILE* rawfile,
            FILE* fpInput,
            snd_pcm_t* capture_handle,
            snd_pcm_format_t format,
            char* device,
            unsigned int rate,
            char* buffer,
            struct myarr derivative);
//=====
void capture_teardown(CaptureCtx* ctx);
struct myarr* capture_next_block(CaptureCtx* ctx, int poll_timeout_ms);
int capture_setup(CaptureCtx* ctx,
                  snd_pcm_t* cap,
                  unsigned int rate,
                  unsigned int bph,
                  unsigned int countdown_seconds,
                  snd_pcm_format_t fmt /* expect SND_PCM_FORMAT_S16_LE */);
int make_timerfd_ms(unsigned initial_ms, unsigned interval_ms);
int build_alsa_pollfds(snd_pcm_t* h, struct pollfd** fds_out, nfds_t* nfds_out);
