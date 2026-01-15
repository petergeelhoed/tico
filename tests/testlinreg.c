#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mylib.h"
#include "mysync.h"

int main(void)
{
    unsigned int fitwindow = 14;
    double xarr[fitwindow];
    double yarr[fitwindow];
    double a = 0;
    double b = 0;
    double s = 0;
    unsigned int n = 0;
    for (int k = 0; k < (int)fitwindow; k++)
    {
        n++;
        xarr[k] = k;
        yarr[k] = 4 + k * 2 + k % 2;
    }

    linreg(xarr, yarr, n, &a, &b, &s);

    fprintf(stderr, "%d %f %f %f\n", n, a, b, s);
    yarr[0] = 600;
    xarr[0] = 1;
    yarr[1] = 592;
    xarr[1] = 2;
    yarr[1] = 552;
    xarr[2] = 6;
    yarr[2] = 458;
    xarr[3] = 7;
    yarr[3] = 120;
    xarr[4] = 10;
    yarr[4] = 398;
    xarr[5] = 11;
    yarr[5] = 87;
    xarr[6] = 12;
    yarr[6] = 394;
    xarr[7] = 14;
    yarr[7] = 358;
    xarr[8] = 16;
    yarr[8] = 332;
    xarr[9] = 19;
    yarr[9] = 16;
    n = 10;
    linreg(xarr, yarr, n, &a, &b, &s);

    fprintf(stderr, "%d %f %f %f\n", n, a, b, s);

    exit(0);
}
