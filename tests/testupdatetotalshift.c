#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "capture_helpers.h"
#include "myarr.h"
#include "mydefs.h"

// Meaningful constants for magic numbers
#define CAPTURE_CVALUE 8
#define CAPTURE_TEETH 2
#define MAXVALS_SIZE 4
#define SHIFT_INIT 100
#define UPDATE_INDEX 9
#define DELTA_SMALL 0.01
#define DELTA_LARGER 0.02
#define DELTA_LARGE 250
#define RETURN_UNEXPECTED_UPDATE 2
#define RETURN_BELOW_THRESHOLD 3
#define RETURN_SMALL_DELTA_MISMATCH 4
#define RETURN_LARGE_DELTA_MISMATCH 5

int main(void)
{
    CapConfig captureConfig = {0};
    captureConfig.cvalue = CAPTURE_CVALUE;
    captureConfig.teeth = CAPTURE_TEETH;

    AppResources resources = {0};
    resources.maxvals = makemyarrd(MAXVALS_SIZE);
    if (resources.maxvals == NULL)
    {
        return 1;
    }

    const double threshold = (double)captureConfig.cvalue / HEX_BASE;
    resources.maxvals->arrd[0] = threshold + DELTA_SMALL;
    resources.maxvals->arrd[1] = threshold - DELTA_SMALL;

    int shiftValue = SHIFT_INIT;

    int shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                               UPDATE_INDEX,
                                               AUTOCOR_LIMIT,
                                               0,
                                               &resources,
                                               &captureConfig);
    if (shiftResult != shiftValue)
    {
        (void)fprintf(stderr, "unexpected shift update at AUTOCOR_LIMIT\n");
        return RETURN_UNEXPECTED_UPDATE;
    }

    shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                           UPDATE_INDEX,
                                           AUTOCOR_LIMIT + 1,
                                           1,
                                           &resources,
                                           &captureConfig);
    if (shiftResult != shiftValue)
    {
        (void)fprintf(stderr, "unexpected shift update below threshold\n");
        return RETURN_BELOW_THRESHOLD;
    }

    resources.maxvals->arrd[0] = threshold + DELTA_LARGER;
    shiftResult = updateTotalShiftIfNeeded(shiftValue,
                                           UPDATE_INDEX,
                                           AUTOCOR_LIMIT + 1,
                                           0,
                                           &resources,
                                           &captureConfig);
    if (shiftResult != shiftValue + UPDATE_INDEX)
    {
        (void)fprintf(stderr,
                      "small-delta shift mismatch: got %d\n",
                      shiftResult);
        return RETURN_SMALL_DELTA_MISMATCH;
    }

    int largeDeltaValue = DELTA_LARGE;
    int absDeltaValue = abs(largeDeltaValue);
    int expectedDeltaValue = (int)(PRESHIFT_THRESHOLD_ROOT * largeDeltaValue /
                                   sqrt((double)absDeltaValue));
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
        return RETURN_LARGE_DELTA_MISMATCH;
    }

    freemyarr(resources.maxvals);
    return 0;
}
