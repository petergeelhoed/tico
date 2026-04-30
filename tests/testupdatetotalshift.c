#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

int main(void)
{
    CapConfig cfg = {0};
    cfg.cvalue = 8;
    cfg.teeth = 2;

    AppResources res = {0};
    res.maxvals = makemyarrd(4);
    if (res.maxvals == NULL)
    {
        return 1;
    }

    const double threshold = (double)cfg.cvalue / HEX_BASE;
    res.maxvals->arrd[0] = threshold + 0.01;
    res.maxvals->arrd[1] = threshold - 0.01;

    int shift = 100;

    int out = updateTotalShiftIfNeeded(shift, 9, AUTOCOR_LIMIT, 0, &res, &cfg);
    if (out != shift)
    {
        (void)fprintf(stderr, "unexpected shift update at AUTOCOR_LIMIT\n");
        return 2;
    }

    out = updateTotalShiftIfNeeded(shift, 9, AUTOCOR_LIMIT + 1, 1, &res, &cfg);
    if (out != shift)
    {
        (void)fprintf(stderr, "unexpected shift update below threshold\n");
        return 3;
    }

    res.maxvals->arrd[0] = threshold + 0.02;
    out = updateTotalShiftIfNeeded(shift, 9, AUTOCOR_LIMIT + 1, 0, &res, &cfg);
    if (out != shift + 9)
    {
        (void)fprintf(stderr, "small-delta shift mismatch: got %d\n", out);
        return 4;
    }

    int largeDelta = 250;
    int absDelta = abs(largeDelta);
    int expectedDelta =
        (int)(PRESHIFT_THRESHOLD_ROOT * largeDelta / sqrt((double)absDelta));
    out = updateTotalShiftIfNeeded(shift,
                                   largeDelta,
                                   AUTOCOR_LIMIT + 3,
                                   0,
                                   &res,
                                   &cfg);
    if (out != shift + expectedDelta)
    {
        (void)fprintf(stderr,
                      "large-delta shift mismatch: got %d expected %d\n",
                      out,
                      shift + expectedDelta);
        return 5;
    }

    freemyarr(res.maxvals);
    return 0;
}
