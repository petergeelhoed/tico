#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"

int main()
{
    unsigned int NN = 48000;

    int blah[NN];
    int orig[NN];
    for (unsigned int j = 0; j < NN; j++)
    {
        blah[j] = (int)(1000 * sin((float)(3.1415926 * 2 * j * 50. / 48000.))) +
                  (int)(800 * sin((float)(3.1415926 * 2 * j * 52. / 48000.)));
        orig[j] = blah[j];
    }

    remove50hz(NN, blah, 48000);

    for (unsigned int j = 0; j < NN; j++)
    {
        printf("%3d %12d %12d\n", j, blah[j], orig[j]);
    }

    exit(0);
}
