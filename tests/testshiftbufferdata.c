#include <stdio.h>
#include <math.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

int main(void)
{
    const unsigned int len = ARRAY_BUFFER_SIZE * 2;
    unsigned int ticktock = len;

    struct myarr* subpos = makemyarrd(len);
    struct myarr* maxpos = makemyarr(len);
    struct myarr* maxvals = makemyarrd(len);

    if (subpos == NULL || maxpos == NULL || maxvals == NULL)
    {
        return 1;
    }

    for (unsigned int i = 0; i < len; ++i)
    {
        subpos->arrd[i] = (double)i + 0.25;
        maxpos->arr[i] = (int)i;
        maxvals->arrd[i] = (double)i + 1000.5;
    }

    shiftBufferData(&ticktock, subpos, maxpos, maxvals);

    if (ticktock != ARRAY_BUFFER_SIZE)
    {
        (void)fprintf(stderr, "ticktock mismatch: %u\n", ticktock);
        return 2;
    }

    for (unsigned int i = 0; i < ARRAY_BUFFER_SIZE; ++i)
    {
        if (fabs(subpos->arrd[i] -
             ((double)(i + ARRAY_BUFFER_SIZE) + 0.25)) > 1e-12)
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
             ((double)(i + ARRAY_BUFFER_SIZE) + 1000.5)) > 1e-12)
        {
            (void)fprintf(stderr, "maxvals mismatch at %u\n", i);
            return 5;
        }
    }

    freemyarr(subpos);
    freemyarr(maxpos);
    freemyarr(maxvals);
    return 0;
}
