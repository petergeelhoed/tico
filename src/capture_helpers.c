#include "capture_helpers.h"

#include "analysis.h"
#include "modsigned.h"
#include "mydefs.h"
#include "myfft.h"
#include "mymath.h"
#include "mysync.h"
#include "parseargs.h"
#include "printing.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdarg.h>

void printFinals(CapConfig* cfg, AppResources* res, size_t totalTickTock)
{
    if (cfg->fpposition)
    {
        calculateTotalFromFile(totalTickTock,
                               cfg->fpposition,
                               res->arrayLength,
                               cfg->SDthreshold,
                               cfg->rate);
    }
    if (cfg->fpmaxcor)
    {
        printTOD(cfg->fpmaxcor);
    }

    if (cfg->fptotal)
    {
        for (size_t t = 0; t < cfg->teeth; ++t)
        {
            struct myarr* tmp = res->teethArray[t];
            if (tmp != NULL)
            {
                struct myarr cumulativeTick = *tmp;
                int toothshift = getshift(*res->teethArray[0], cumulativeTick);
                for (size_t j = 0; j < res->arrayLength; ++j)
                {
                    (void)fprintf(cfg->fptotal,
                                  "%zu %d %zu %d\n",
                                  j,
                                  cumulativeTick.arr[j],
                                  t,
                                  toothshift);
                }
            }
            else
            {
                break;
            }
            (void)fprintf(cfg->fptotal, "\n\n");
        }
    }
}

void fillReference(FILE* fpDefPeak, struct myarr* reference, size_t teeth)
{
    if (fpDefPeak != NULL)
    {
        int arr[4];
        for (size_t t = 0; t < teeth; t++)
        {
            for (size_t j = 0; j < reference->ArrayLength; j++)
            {
                int int_count = getIntsFromStdin(4, arr);
                if (int_count < 0)
                {
                    break;
                }
                if (int_count < 3)
                {
                    (void)fprintf(
                        stderr,
                        "not enough values in -D <default peak file>\n 4 "
                        "columns required, %zu samples and %zu teeth\n",
                        reference->ArrayLength,
                        teeth);
                    exit(EXIT_FAILURE);
                }
                int value = arr[1];

                reference->arr[((int)j + (int)reference->ArrayLength) %
                               (int)reference->ArrayLength] = value;
            }
        }
        (void)fclose(fpDefPeak);
        fpDefPeak = NULL;
    }
    else
    {
        const int peakheight[3] = {100000, 80000, 60000};
        const size_t peakpos[3] = {0, 400, 800};

        for (int i = 0; i < 3; i++)
        {
            reference->arr[reference->ArrayLength / 4 - peakpos[i]] =
                peakheight[i];
            reference->arr[3 * reference->ArrayLength / 4 - peakpos[i]] =
                peakheight[i];
        }
    }
}

void shiftBufferData(size_t* ticktock,
                     struct myarr* subpos,
                     struct myarr* maxpos,
                     struct myarr* maxvals)
{
    memmove(subpos->arrd,
            subpos->arrd + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(double));
    memmove(maxpos->arr,
            maxpos->arr + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(int));
    memmove(maxvals->arrd,
            maxvals->arrd + ARRAY_BUFFER_SIZE,
            ARRAY_BUFFER_SIZE * sizeof(double));
    *ticktock -= ARRAY_BUFFER_SIZE;
}

void processLogging(CapConfig* cfg,
                    AppResources* res,
                    size_t totalTime,
                    size_t writeInterval)
{
    if (totalTime > 0 && totalTime % writeInterval == 0)
    {
        if (cfg->fpposition)
        {
            struct myarr* positionBatch = makemyarrd(writeInterval);
            for (size_t k = 0; k < writeInterval; ++k)
            {
                positionBatch->arrd[k] =
                    res->subpos->arrd[totalTime - writeInterval + k] +
                    (double)res->maxpos->arr[totalTime - writeInterval + k];
            }
            syncAppendMyarr(positionBatch, cfg->fpposition);
            freemyarr(positionBatch);
        }
        if (cfg->fpmaxcor)
        {
            struct myarr* correlationBatch = makemyarrd(writeInterval);
            memcpy(correlationBatch->arrd,
                   res->maxvals->arrd + totalTime - writeInterval,
                   writeInterval * sizeof(double));
            syncAppendMyarr(correlationBatch, cfg->fpmaxcor);
            freemyarr(correlationBatch);
        }
    }
}

void fitAndPrint(const LoopState* state,
                 struct myarr* cumulativeTick,
                 AppResources* res,
                 CapConfig* cfg,
                 size_t currentColumns)
{
    double intercept = 0.0;
    double slope = 0.0;
    size_t arrayLength = res->arrayLength;
    fitNpeaks(&intercept,
              &slope,
              (unsigned int)state->tickIndex,
              res->maxvals,
              res->maxpos,
              res->subpos,
              cfg->fitN,
              cfg->SDthreshold);

    printheader(
        slope * SECS_DAY / (double)arrayLength,
        cfg->everyline,
        getBeatError(cumulativeTick, cfg->rate, 0),
        (double)state->globalTickIndex * (double)arrayLength / cfg->rate);

    printspaces(
        res->maxpos->arr[(unsigned int)state->tickIndex],
        (size_t)(res->maxvals->arrd[(unsigned int)state->tickIndex] * HEX_BASE),
        res->mod,
        currentColumns - cfg->everyline,
        intercept,
        cfg->cvalue);
}

void rotateDerivativeWindow(AppResources* res, int cumulativeShift)
{
    const size_t arrayLength = res->arrayLength;
    memmove(res->tmpder->arr,
            res->derivative->arr + modSigned(cumulativeShift, arrayLength),
            arrayLength * sizeof(int));
    memmove(res->tmpder->arr + arrayLength -
                modSigned(cumulativeShift, arrayLength),
            res->derivative->arr,
            modSigned(cumulativeShift, arrayLength) * sizeof(int));
}

int findMaxPosition(AppResources* res,
                    struct myarr* cumulativeTick,
                    const LoopState* state,
                    CapConfig* cfg)
{
    const int useReference =
        (state->globalTickIndex < AUTOCOR_LIMIT * cfg->teeth);
    return shiftHalf(
        fftfit(*res->tmpder,
               cumulativeTick->arr,
               useReference ? res->reference->arr : cumulativeTick->arr,
               res->maxvals->arrd + state->tickIndex,
               res->filterFFT,
               state->globalTickIndex == cfg->verbose,
               res->subpos->arrd + state->tickIndex),
        cumulativeTick->ArrayLength);
}

void updateTotalShiftIfNeeded(LoopState* state,
                              int peakOffset,
                              AppResources* res,
                              CapConfig* cfg)
{
    if (state->globalTickIndex > AUTOCOR_LIMIT &&
        res->maxvals->arrd[state->tickIndex] > (double)cfg->cvalue / HEX_BASE &&
        state->globalTickIndex % cfg->teeth == 0)
    {
        int delta = peakOffset;
        if (abs(delta) > PRESHIFT_THRESHOLD)
        {
            delta = (int)(PRESHIFT_THRESHOLD_ROOT * delta / sqrt(abs(delta)));
        }
        state->cumulativeShift += delta;
    }
}
