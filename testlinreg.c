#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <alsa/asoundlib.h>
#include <fftw3.h>

#include "mylib.h"

int main (int argc, char *argv[])
{
    int fitwindow = 14;
    int xarr[fitwindow];
    int yarr[fitwindow];
    double a = 0;
    double b = 0;
    double s = 0;
    int n =0;
    for (int k = 0; k < fitwindow;k++)
    {
        n++;
        xarr[k] = k;
        yarr[k] = 4+k*2 + k%2;

    }

    linreg(xarr,yarr, n, &a, &b, &s);

    fprintf(stderr,"%d %f %f %f\n",n,a,b,s);

    exit (0);
}
