#include "resources.h"
#include "printing.h"

#include "myfft.h"
#include "mysync.h"

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

AppResources allocateResources(size_t arrayLength,
                               size_t ticktockBuffer,
                               CapConfig* cfg)
{
    AppResources res = {
        .arrayLength = arrayLength,
        .subpos = makemyarrd(ticktockBuffer),
        .maxpos = makemyarr(ticktockBuffer),
        .maxvals = makemyarrd(ticktockBuffer),
        .derivative = makemyarr(arrayLength),
        .tmpder = makemyarr(arrayLength),
        .reference = makemyarr(arrayLength),
        .teethArray = calloc(cfg->teeth, sizeof(*res.teethArray)),
        .filterFFT = makeFilter(cfg->evalue, arrayLength),
        .audioBuffer16 = calloc(arrayLength, sizeof(*res.audioBuffer16))};

    if (res.teethArray == NULL || res.audioBuffer16 == NULL)
    {
        free(res.audioBuffer16);
        free(res.teethArray);
        print("Failed memory allocation\n");
        exit(EXIT_FAILURE);
    }

    for (size_t t = 0; t < cfg->teeth; t++)
    {
        res.teethArray[t] = makemyarr(arrayLength);
    }

    res.mod = arrayLength / cfg->zoom;
    res.maxTime = (unsigned int)cfg->rate *
                  (cfg->time ? cfg->time : DEFAULT_TIME) / arrayLength;
    return res;
}

void cleanupResources(AppResources* res, CapConfig* cfg, CaptureCtx* ctx)
{
    if (cfg->fpposition)
    {
        waitClose(cfg->fpposition);
    }
    if (cfg->fpInput)
    {
        (void)fclose(cfg->fpInput);
    }
    if (cfg->fpmaxcor)
    {
        waitClose(cfg->fpmaxcor);
    }
    if (cfg->fptotal)
    {
        (void)fclose(cfg->fptotal);
    }
    if (cfg->captureHandle)
    {
        snd_pcm_close(cfg->captureHandle);
    }
    free(res->audioBuffer16);
    freemyarr(res->subpos);
    freemyarr(res->maxpos);
    freemyarr(res->maxvals);
    freemyarr(res->derivative);
    freemyarr(res->tmpder);
    freemyarr(res->reference);
    fftw_free(res->filterFFT);
    for (unsigned int t = 0; t < cfg->teeth; t++)
    {
        freemyarr(res->teethArray[t]);
    }
    free(res->teethArray);
    captureTeardown(ctx);
}
