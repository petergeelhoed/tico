#include <math.h>
#include <stdio.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

// Magic number constants
#define SHIFTBUFFER_EXTRA_FACTOR 2
#define SHIFTBUFFER_OFFSET 0.25
#define SHIFTBUFFER_LARGE_OFFSET 1000.5
#define SHIFTBUFFER_TOLERANCE 1e-12
#define RETURN_TICKTOCK_MISMATCH 2

int main(void)
{
    const unsigned int totalLength =
        ARRAY_BUFFER_SIZE * SHIFTBUFFER_EXTRA_FACTOR;
    size_t ticktockCounter = totalLength;

    struct myarr* subpos = makemyarrd(totalLength);
    struct myarr* maxpos = makemyarr(totalLength);
    struct myarr* maxvals = makemyarrd(totalLength);

    if (subpos == NULL || maxpos == NULL || maxvals == NULL)
    {
        return 1;
    }

    for (unsigned int index = 0; index < totalLength; ++index)
    {
        subpos->arrd[index] = (double)index + SHIFTBUFFER_OFFSET;
        maxpos->arr[index] = (int)index;
        maxvals->arrd[index] = (double)index + SHIFTBUFFER_LARGE_OFFSET;
    }

    shiftBufferData(&ticktockCounter, subpos, maxpos, maxvals);

    if (ticktockCounter != ARRAY_BUFFER_SIZE)
    {
        (void)fprintf(stderr, "ticktock mismatch: %zu\n", ticktockCounter);
        return RETURN_TICKTOCK_MISMATCH;
    }

    for (unsigned int i = 0; i < ARRAY_BUFFER_SIZE; ++i)
    {
        if (fabs(subpos->arrd[i] -
                 ((double)(i + ARRAY_BUFFER_SIZE) + SHIFTBUFFER_OFFSET)) >
            SHIFTBUFFER_TOLERANCE)
        {
            (void)fprintf(stderr, "subpos mismatch at %u\n", i);
            return 3;
        }
        if (maxpos->arr[i] != (int)(i + ARRAY_BUFFER_SIZE))
        {
            (void)fprintf(stderr, "maxpos mismatch at %u\n", i);
            return 4;
        }
        if (fabs(maxvals->arrd[i] -
                 ((double)(i + ARRAY_BUFFER_SIZE) + SHIFTBUFFER_LARGE_OFFSET)) >
            SHIFTBUFFER_TOLERANCE)
        {
            (void)fprintf(stderr, "maxvals mismatch at %u\n", i);
            return 5; // NOLINT(readability-magic-numbers)
        }
    }

    freemyarr(subpos);
    freemyarr(maxpos);
    freemyarr(maxvals);
    return 0;
}
