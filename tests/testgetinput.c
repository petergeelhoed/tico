#include "mylib.h"
#include "mymath.h"
#include "parseargs.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct exp
{
    char input[BUFFER_SIZE];
    double dbl;
    int integer;
};

int main(void)
{
    struct exp testCases[NR_TESTS];
    // NOLINTBEGIN(readability-magic-numbers)
    testCases[0] = (struct exp){.input = " 34", .integer = 34, .dbl = 34};
    testCases[1] = (struct exp){.input = " 34.1", .integer = 34, .dbl = 34.1};
    testCases[2] =
        (struct exp){.input = "", .integer = INT_MIN, .dbl = (double)NAN};
    testCases[3] = (struct exp){.input = " -1.7 ", .integer = -1, .dbl = -1.7};
    testCases[4] = (struct exp){.input = " 1111111111111111",
                                .integer = INT_MIN,
                                .dbl = 1111111111111111.};
    // NOLINTEND(readability-magic-numbers)

    for (size_t testIndex = 0; testIndex < NR_TESTS; testIndex++)
    {
        struct exp testCase = testCases[testIndex];
        int intResult = getInt(testCase.input);
        double b = getDouble(test1.input);
        if (!nearlyEqual(b, test1.dbl))
        {
            if (isnan(b) && isnan(test1.dbl))
            {
                continue;
            }
            printf("test %lu:\n"
                   "input '%s' => (dbl)%lf : expected %lf\n",
                   i,
                   test1.input,
                   b,
                   test1.dbl);
            return -1;
        }
        if (a != test1.integer)
        {
            printf("test %lu:\n"
                   "input '%s' => (int)%d :expected %d \n",
                   i,
                   test1.input,
                   a,
                   test1.integer);
            return -1;
        }
    }

    puts("all fine");
    return 0;
}
