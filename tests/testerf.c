#include "erf.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define EXPECTED_MEAN 2936.000
#define EXPECTED_STDEV 146.106
#define EPSILON 0.001

int main(void)
{
    // NOLINTBEGIN(readability-magic-numbers)
    double samples[] = {2699,
                        2746,
                        2799,
                        2901,
                        2914,
                        2921,
                        2922,
                        2978,
                        3018,
                        3077,
                        3107,
                        3150};
    // NOLINTEND(readability-magic-numbers)

    struct stats result =
        fit_erf(samples, sizeof(samples) / sizeof(samples[0]));

    return ((fabs(result.mean - EXPECTED_MEAN) >= EPSILON) ||
            (fabs(result.stdev - EXPECTED_STDEV) >= EPSILON));
}
