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
    snd_pcm_uframes_t periodSize;
    snd_pcm_uframes_t bufferSize;
    unsigned int ArrayLength; // frames per processing block
} CaptureCtx;

// NOLINTNEXTLINE(misc-include-cleaner)
snd_pcm_t* initAudio(snd_pcm_format_t format, char* device, unsigned int* rate);

int readBufferOrFile(int* derivative,
                     unsigned int ArrayLength,
                     FILE* fpInput,
                     CaptureCtx* ctx,
                     int16_t* buffer16);

int getData(FILE* fpInput,
            struct myarr* derivative,
            CaptureCtx* ctx,
            int16_t* out);

void captureTeardown(CaptureCtx* ctx);

int captureSetup(CaptureCtx* ctx, CapConfig* cfg, unsigned int rate);

int readSamples(snd_pcm_t* cap, unsigned int ArrayLength, int16_t* out);
