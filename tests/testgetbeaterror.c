#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "analysis.h"
#include "myarr.h"

int main(void)
{
// Magic number constants
#define TESTGETBEATERROR_TOTAL_SIZE 16
#define TESTGETBEATERROR_PATTERN_SIZE 8
#define TESTGETBEATERROR_SAMPLE_RATE 48000.0
#define TESTGETBEATERROR_TOLERANCE 1e-12

    struct myarr* total = makemyarr(TESTGETBEATERROR_TOTAL_SIZE);
    if (total == NULL)
    {
        return 1;
    }

    const int pattern[TESTGETBEATERROR_PATTERN_SIZE] = {3, 1, 4, 1, 5, 9, 2, 6};
    for (unsigned int i = 0; i < TESTGETBEATERROR_PATTERN_SIZE; ++i)
    {
        total->arr[i] = pattern[i];
        total->arr[TESTGETBEATERROR_PATTERN_SIZE + i] = pattern[i];
    }

    double beatError = getBeatError(total, TESTGETBEATERROR_SAMPLE_RATE, 0);
    if (fabs(beatError) > TESTGETBEATERROR_TOLERANCE)
    {
        (void)fprintf(stderr, "expected beat error 0, got %g\n", beatError);
        freemyarr(total);
        return 2;
    }

    freemyarr(total);
    return 0;
}
