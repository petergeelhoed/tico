#pragma once
#include "config.h"
#include "mysound.h"
#include "myarr.h"

#include "fftw3.h"

#include <stdint.h>

typedef struct
{
    struct myarr* subpos;
    struct myarr* maxpos;
    struct myarr* maxvals;
    struct myarr* derivative;
    struct myarr* tmpder;
    struct myarr* reference;
    struct myarr** teethArray;
    fftw_complex* filterFFT;
    int16_t* audioBuffer16;
} AppResources;

AppResources allocateResources(unsigned int arrayLength,
                               unsigned int ticktockBuffer,
                               CapConfig* cfg);

void cleanupResources(AppResources* res, CapConfig* cfg, CaptureCtx* ctx);
