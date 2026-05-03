#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "crosscorint.h"
#include "myfft.h"
// Magic number constants
#define CROSSCOR_ARRAY_LENGTH 20
#define CROSSCOR_INDEX_1 4
#define CROSSCOR_INDEX_2 7

int main(void)
{
    int peak[CROSSCOR_ARRAY_LENGTH];
    int peak2[CROSSCOR_ARRAY_LENGTH];
    int cross[CROSSCOR_ARRAY_LENGTH];
    for (unsigned int j = 0; j < CROSSCOR_ARRAY_LENGTH; j++)
    {
        peak2[(j + 1) % CROSSCOR_ARRAY_LENGTH] =
            ((j == CROSSCOR_INDEX_1) + (j == CROSSCOR_INDEX_2));
        peak[j] = ((j == CROSSCOR_INDEX_1) + (j == CROSSCOR_INDEX_2));
    }

    crosscorint(CROSSCOR_ARRAY_LENGTH, peak, peak2, cross);
    for (unsigned int i = 0; i < CROSSCOR_ARRAY_LENGTH; i++)
    {
        printf("%d %d %d %d\n", i, peak[i], peak2[i], cross[i]);
    }

    exit(0);
}
