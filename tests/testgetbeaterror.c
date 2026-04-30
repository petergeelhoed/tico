#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "analysis.h"
#include "myarr.h"

int main(void)
{
    struct myarr* total = makemyarr(16);
    if (total == NULL)
    {
        return 1;
    }

    const int pattern[8] = {3, 1, 4, 1, 5, 9, 2, 6};
    for (unsigned int i = 0; i < 8; ++i)
    {
        total->arr[i] = pattern[i];
        total->arr[8 + i] = pattern[i];
    }

    double beatError = getBeatError(total, 48000.0, 0);
    if (fabs(beatError) > 1e-12)
    {
        (void)fprintf(stderr, "expected beat error 0, got %g\n", beatError);
        freemyarr(total);
        return 2;
    }

    freemyarr(total);
    return 0;
}
