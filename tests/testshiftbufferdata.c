#include <math.h>
#include <stdio.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

int main(void)
{
    const unsigned int totalLength = ARRAY_BUFFER_SIZE * 2;
    unsigned int ticktockCounter = totalLength;

    struct myarr* subPositionArray = makemyarrd(totalLength);
    struct myarr* maxPositionArray = makemyarr(totalLength);
    struct myarr* maxValueArray = makemyarrd(totalLength);

    if (subPositionArray == NULL || maxPositionArray == NULL || maxValueArray == NULL)
    {
        return 1;
    }

    for (unsigned int index = 0; index < totalLength; ++index)
    {
        subPositionArray->arrd[index] = (double)index + 0.25; // NOLINT(readability-magic-numbers)
        maxPositionArray->arr[index] = (int)index;
        maxValueArray->arrd[index] =
            (double)index + 1000.5; // NOLINT(readability-magic-numbers)
    }

    shiftBufferData(&ticktockCounter, subPositionArray, maxPositionArray, maxValueArray);

    if (ticktockCounter != ARRAY_BUFFER_SIZE)
    {
        (void)fprintf(stderr, "ticktock mismatch: %u\n", ticktockCounter);
        return 2;
    }

    for (unsigned int index = 0; index < ARRAY_BUFFER_SIZE; ++index)
    {
        if (fabs(subPositionArray->arrd[index] -
                 ((double)(index + ARRAY_BUFFER_SIZE) +
                  0.25)) > // NOLINT(readability-magic-numbers)
            1e-12)         // NOLINT(readability-magic-numbers)
        {
            (void)fprintf(stderr, "subpos mismatch at %u\n", index);
            return 3;
        }
        if (maxpos->arr[i] != (int)(i + ARRAY_BUFFER_SIZE))
        {
            (void)fprintf(stderr, "maxpos mismatch at %u\n", i);
            return 4;
        }
        if (fabs(maxvals->arrd[i] -
                 ((double)(i + ARRAY_BUFFER_SIZE) + 1000.5)) >
            1e-12) // NOLINT(readability-magic-numbers)
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
