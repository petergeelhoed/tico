#include <stdio.h>

#include "capture_helpers.h"
#include "myarr.h"

int main(void)
{
    const unsigned int arrayLength = 5; // NOLINT(readability-magic-numbers)
    AppResources resources = {0};
    resources.derivative = makemyarr(arrayLength);
    resources.tmpder = makemyarr(arrayLength);

    if (resources.derivative == NULL || resources.tmpder == NULL)
    {
        return 1;
    }

    for (unsigned int index = 0; index < arrayLength; ++index)
    {
        resources.derivative->arr[index] =
            (int)(10 + index); // NOLINT(readability-magic-numbers)
    }

    rotateDerivativeWindow(&resources, arrayLength, -1);

    const int expectedRotated[5] = {14, 10, 11, 12, 13}; // NOLINT(readability-magic-numbers)
    for (unsigned int index = 0; index < arrayLength; ++index)
    {
        if (resources.tmpder->arr[index] != expectedRotated[index])
        {
            (void)fprintf(stderr,
                          "rotate(-1) mismatch at %u: got %d expected %d\n",
                          index,
                          resources.tmpder->arr[index],
                          expectedRotated[index]);
            return 2;
        }
    }

    rotateDerivativeWindow(&resources, arrayLength, 2);

    const int expected1[5] = {12, 13, 14, 10, 11};
    for (unsigned int i = 0; i < len; ++i)
    {
        if (res.tmpder->arr[i] != expected1[i])
        {
            (void)fprintf(stderr,
                          "rotate(2) mismatch at %u: got %d expected %d\n",
                          i,
                          res.tmpder->arr[i],
                          expected1[i]);
            return 3;
        }
    }

    freemyarr(res.derivative);
    freemyarr(res.tmpder);
    return 0;
}
