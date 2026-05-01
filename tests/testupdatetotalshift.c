#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

int main(void)
{
    CapConfig captureConfig = {0};
    captureConfig.cvalue = 8; // NOLINT(readability-magic-numbers)
    captureConfig.teeth = 2;

    AppResources resources = {0};
    resources.maxvals = makemyarrd(4);
    if (resources.maxvals == NULL)
    {
        return 1;
    }

    const double threshold = (double)captureConfig.cvalue / HEX_BASE;
    resources.maxvals->arrd[0] =
        threshold + 0.01; // NOLINT(readability-magic-numbers)
    resources.maxvals->arrd[1] =
        threshold - 0.01; // NOLINT(readability-magic-numbers)

    int shiftValue = 100;

    int shiftResult =
        updateTotalShiftIfNeeded(shiftValue,
                                 9,
                                 AUTOCOR_LIMIT,
                                 0,
                                 &resources,
                                 &captureConfig); // NOLINT(readability-magic-numbers)
    if (shiftResult != shiftValue)
    {
        (void)fprintf(stderr, "unexpected shift update at AUTOCOR_LIMIT\n");
        return 2;
    }

    shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                   9,
                                   AUTOCOR_LIMIT + 1,
                                   1,
                                   &resources,
                                   &captureConfig); // NOLINT(readability-magic-numbers)
    if (shiftResult != shiftValue)
    {
        (void)fprintf(stderr, "unexpected shift update below threshold\n");
        return 3;
    }

    resources.maxvals->arrd[0] =
        threshold + 0.02; // NOLINT(readability-magic-numbers)
    shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                   9,
                                   AUTOCOR_LIMIT + 1,
                                   0,
                                   &resources,
                                   &captureConfig); // NOLINT(readability-magic-numbers)
    if (shiftResult != shiftValue + 9)                 // NOLINT(readability-magic-numbers)
    {
        (void)fprintf(stderr, "small-delta shift mismatch: got %d\n", shiftResult);
        return 4;
    }

    int largeDeltaValue = 250; // NOLINT(readability-magic-numbers)
    int absDeltaValue = abs(largeDeltaValue);
    int expectedDeltaValue =
        (int)(PRESHIFT_THRESHOLD_ROOT * largeDeltaValue / sqrt((double)absDeltaValue));
    shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                   largeDeltaValue,
                                   AUTOCOR_LIMIT + 3,
                                   0,
                                   &resources,
                                   &captureConfig);
    if (shiftResult != shiftValue + expectedDeltaValue)
    {
        (void)fprintf(stderr,
                      "large-delta shift mismatch: got %d expected %d\n",
                      shiftResult,
                      shiftValue + expectedDeltaValue);
        return 5; // NOLINT(readability-magic-numbers)
    }

    freemyarr(resources.maxvals);
    return 0;
}
