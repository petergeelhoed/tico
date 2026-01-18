#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

#define ArrayLength 48000
int main(void)
{

    int blah[ArrayLength];
    int orig[ArrayLength];
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        blah[j] = (int)(1000 * sin((float)(3.1415926 * 2 * j * 50. / 48000.))) +
                  (int)(800 * sin((float)(3.1415926 * 2 * j * 52. / 48000.)));
        orig[j] = blah[j];
    }

    remove50hz(ArrayLength, blah, 48000);

    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        printf("%3d %12d %12d\n", j, blah[j], orig[j]);
    }

    exit(0);
}
