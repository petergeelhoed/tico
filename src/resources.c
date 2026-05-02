#include "resources.h"

#include "myfft.h"
#include "mysync.h"

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

AppResources allocateResources(size_t arrayLength,
                               size_t ticktockBuffer,
                               CapConfig* cfg)
{
    AppResources res = {0};
    res.subpos = makemyarrd(ticktockBuffer);
    res.maxpos = makemyarr(ticktockBuffer);
    res.maxvals = makemyarrd(ticktockBuffer);
    res.derivative = makemyarr(arrayLength);
    res.tmpder = makemyarr(arrayLength);
    res.reference = makemyarr(arrayLength);
    res.filterFFT = makeFilter(cfg->evalue, arrayLength);
    res.audioBuffer16 = calloc(arrayLength, sizeof(*res.audioBuffer16));
    res.teethArray = calloc(cfg->teeth, sizeof(*res.teethArray));
    if (res.teethArray == NULL || res.audioBuffer16 == NULL)
    {
        free(res.audioBuffer16);
        free(res.teethArray);
        (void)fprintf(stderr, "Failed memory allocation\n");
        exit(EXIT_FAILURE);
    }

    for (size_t t = 0; t < cfg->teeth; t++)
    {
        res.teethArray[t] = makemyarr(arrayLength);
    }
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
