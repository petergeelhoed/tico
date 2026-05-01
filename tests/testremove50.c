#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

    #define kArrayLength 48000
    int main(void)
    {
        int signalArray[kArrayLength];
        int originalArray[kArrayLength];
        for (unsigned int sampleIndex = 0; sampleIndex < kArrayLength; sampleIndex++)
        {
            signalArray[sampleIndex] =
                (int)(1000 *
                      sin((double)(3.1415926 * 2 * sampleIndex * 50. /
                                   48000.))) + // NOLINT(readability-magic-numbers)
                (int)(800 *
                      sin((double)(3.1415926 * 2 * sampleIndex * 52. /
                                   48000.))); // NOLINT(readability-magic-numbers)
            originalArray[sampleIndex] = signalArray[sampleIndex];
        }

        remove50hz(kArrayLength, signalArray, 48000); // NOLINT(readability-magic-numbers)

        for (unsigned int sampleIndex = 0; sampleIndex < kArrayLength; sampleIndex++)
    {
        printf("%3d %12d %12d\n", j, blah[j], orig[j]);
    }

    exit(0);
}
