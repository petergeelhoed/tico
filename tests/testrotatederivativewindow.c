#include <stdio.h>

#include "capture_helpers.h"
#include "myarr.h"

int main(void)
{
    const unsigned int len = 5;
    AppResources res = {0};
    res.derivative = makemyarr(len);
    res.tmpder = makemyarr(len);

    if (res.derivative == NULL || res.tmpder == NULL)
    {
        return 1;
    }

    for (unsigned int i = 0; i < len; ++i)
    {
        res.derivative->arr[i] = (int)(10 + i);
    }

    rotateDerivativeWindow(&res, len, -1);

    const int expected0[5] = {14, 10, 11, 12, 13};
    for (unsigned int i = 0; i < len; ++i)
    {
        if (res.tmpder->arr[i] != expected0[i])
        {
            (void)fprintf(stderr,
                          "rotate(-1) mismatch at %u: got %d expected %d\n",
                          i,
                          res.tmpder->arr[i],
                          expected0[i]);
            return 2;
        }
    }

    rotateDerivativeWindow(&res, len, 2);

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
