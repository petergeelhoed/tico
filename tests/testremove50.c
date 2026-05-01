#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

// Magic number constants
#define SIGNAL_ARRAY_LENGTH 48000
#define SINE_FREQ1 50
#define SINE_FREQ2 52
#define SINE_AMP1 1000
#define SINE_AMP2 800
// Magic number constants
#define PI 3.1415926
#define TWO 2.0
int main(void)
{
    int signalArray[SIGNAL_ARRAY_LENGTH];
    int originalArray[SIGNAL_ARRAY_LENGTH];
    for (unsigned int sampleIndex = 0; sampleIndex < SIGNAL_ARRAY_LENGTH;
         sampleIndex++)
    {
        signalArray[sampleIndex] =
            (int)(SINE_AMP1 * sin(PI * TWO * sampleIndex * SINE_FREQ1 /
                                  SIGNAL_ARRAY_LENGTH)) +
            (int)(SINE_AMP2 * sin(PI * TWO * sampleIndex * SINE_FREQ2 /
                                  SIGNAL_ARRAY_LENGTH));
        originalArray[sampleIndex] = signalArray[sampleIndex];
    }

    remove50hz(SIGNAL_ARRAY_LENGTH, signalArray, SIGNAL_ARRAY_LENGTH);

    for (unsigned int sampleIndex = 0; sampleIndex < SIGNAL_ARRAY_LENGTH;
         sampleIndex++)
    {
        printf("%3u %12d %12d\n",
               sampleIndex,
               signalArray[sampleIndex],
               originalArray[sampleIndex]);
    }

    exit(0);
}
