#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"
#define ArrayLength 20
int main(void)
{
    int peak[ArrayLength];
    int peak2[ArrayLength];
    int cross[ArrayLength];
    for (unsigned int j = 0; j < ArrayLength; j++)
    {
        peak2[(j + 1) % ArrayLength] = ((j == 4) + (j == 7));
        peak[j] = ((j == 4) + (j == 7));
    }

    crosscorint(ArrayLength, peak, peak2, cross);
    for (unsigned int i = 0; i < ArrayLength; i++)
    {
        printf("%d %d %d %d\n", i, peak[i], peak2[i], cross[i]);
    }

    exit(0);
}
