#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mylib.h"
#include "myfft.h"

int main()
{
    int NN = 20;
    int peak[NN];
    int peak2[NN];
    int cross[NN];
    for (int j = 0; j < NN; j++)
    {
        peak2[(j + 1) % NN] = ((j == 4) + (j == 7));
        peak[j] = ((j == 4) + (j == 7));
    }

    crosscorint(NN, peak, peak2, cross);
    for (int i = 0; i < NN; i++)
    {
        printf("%d %d %d %d\n",i,peak[i],peak2[i],cross[i]);
    }

    exit(0);
}
