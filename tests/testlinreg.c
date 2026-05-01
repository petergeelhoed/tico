#include <alsa/asoundlib.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myfft.h"
#include "mymath.h"
#include "mysync.h"

// Magic number constants
#define FITWINDOW 14
#define YARR_0 600
#define XARR_0 1
#define YARR_1A 592
#define LINREG_IDX_7 7
#define LINREG_IDX_8 8
#define LINREG_IDX_9 9
int main(void)
{
    // NOLINTBEGIN[readability-magic-numbers]
    double xarr[FITWINDOW] =
        {1, 2, 6, 7, 10, 11, 12, 14, 16, 19, 10, 11, 12, 13};
    double yarr[FITWINDOW] =
        {600, 592, 458, 120, 398, 87, 394, 358, 332, 16, 100, 200, 300, 400};
    // NOLINTEND[readability-magic-numbers]
    double slope = 0;
    double intercept = 0;
    double stddev = 0;
    unsigned int sampleCount = FITWINDOW;

    linreg(xarr, yarr, sampleCount, &slope, &intercept, &stddev);

    (void)
        fprintf(stderr, "%d %f %f %f\n", sampleCount, slope, intercept, stddev);
    linreg(xarr, yarr, sampleCount, &slope, &intercept, &stddev);

    (void)
        fprintf(stderr, "%d %f %f %f\n", sampleCount, slope, intercept, stddev);

    exit(0);
}
