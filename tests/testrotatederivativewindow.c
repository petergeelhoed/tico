#include <stdio.h>

#include "capture_helpers.h"
#include "myarr.h"

// Magic number constants
#define DERIVATIVE_ARRAY_LENGTH 5
#define DERIVATIVE_INIT_OFFSET 10
#define ROTATE_MINUS_ONE -1
#define ROTATE_PLUS_TWO 2
#define RETURN_ROTATE_MISMATCH 2
static const int kExpectedRotated[DERIVATIVE_ARRAY_LENGTH] = {14,
                                                              10,
                                                              11,
                                                              12,
                                                              13};

int main(void)
{
    AppResources resources = {0};
    resources.derivative = makemyarr(DERIVATIVE_ARRAY_LENGTH);
    resources.tmpder = makemyarr(DERIVATIVE_ARRAY_LENGTH);

    if (resources.derivative == NULL || resources.tmpder == NULL)
    {
        return 1;
    }

    for (unsigned int index = 0; index < DERIVATIVE_ARRAY_LENGTH; ++index)
    {
        resources.derivative->arr[index] =
            (int)(DERIVATIVE_INIT_OFFSET + index);
    }

    rotateDerivativeWindow(&resources,
                           DERIVATIVE_ARRAY_LENGTH,
                           ROTATE_MINUS_ONE);

    for (unsigned int index = 0; index < DERIVATIVE_ARRAY_LENGTH; ++index)
    {
        if (resources.tmpder->arr[index] != kExpectedRotated[index])
        {
            (void)fprintf(stderr,
                          "rotate(-1) mismatch at %u: got %d expected %d\n",
                          index,
                          resources.tmpder->arr[index],
                          kExpectedRotated[index]);
            return RETURN_ROTATE_MISMATCH;
        }
    }

    rotateDerivativeWindow(&resources,
                           DERIVATIVE_ARRAY_LENGTH,
                           ROTATE_PLUS_TWO);

    const int expected1[DERIVATIVE_ARRAY_LENGTH] = {12, 13, 14, 10, 11};
    for (unsigned int i = 0; i < DERIVATIVE_ARRAY_LENGTH; ++i)
    {
        if (resources.tmpder->arr[i] != expected1[i])
        {
            (void)fprintf(stderr,
                          "rotate(2) mismatch at %u: got %d expected %d\n",
                          i,
                          resources.tmpder->arr[i],
                          expected1[i]);
            return 3;
        }
    }

    freemyarr(resources.derivative);
    freemyarr(resources.tmpder);
    return 0;
}
