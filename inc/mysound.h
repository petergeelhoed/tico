#pragma once
#include "config.h"
#include "myarr.h"
#include <alsa/asoundlib.h> // IWYU pragma: export
#include <stdio.h>

/* -------------------- Capture Context -------------------- */

typedef struct CaptureCtx
{
    /* ALSA */
    snd_pcm_t* cap;
    unsigned int rate;
    snd_pcm_uframes_t period_size;
    snd_pcm_uframes_t buffer_size;
    unsigned int ArrayLength; // frames per processing block
} CaptureCtx;

// NOLINTNEXTLINE(misc-include-cleaner)
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate);

int readBufferOrFile(int* derivative,
                     unsigned int ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx,
                     int16_t* samples);

int getData(FILE* fpInput,
            struct myarr* derivative,
            CaptureCtx* ctx,
            int16_t* out);

void capture_teardown(CaptureCtx* ctx);

int capture_setup(CaptureCtx* ctx,
                  CapConfig* cfg,
                  unsigned int rate,
                  unsigned int bph);

int read_samples(snd_pcm_t* cap, unsigned int ArrayLength, int16_t* out);
