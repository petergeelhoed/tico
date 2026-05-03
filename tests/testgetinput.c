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

// Magic number constants
#define TEST_INT_1 34
#define TEST_DBL_1 34.0
#define TEST_INPUT_2 " 34.1"
#define TEST_INT_2 34
#define TEST_DBL_2 34.1
#define TEST_INPUT_3 ""
#define TEST_INT_3 INT_MIN
#define TEST_DBL_3 ((double)NAN)
#define TEST_INPUT_4 " -1.7 "
#define TEST_INT_4 -1
#define TEST_DBL_4 -1.7
#define TEST_INPUT_5 " 1111111111111111"
#define TEST_INT_5 INT_MIN
#define TEST_DBL_5 1111111111111111.0

int main(void)
{
    struct exp testCases[NR_TESTS];
    // NOLINTBEGIN(readability-magic-numbers)
    testCases[0] =
        (struct exp){.input = " 34", .integer = TEST_INT_1, .dbl = TEST_DBL_1};
    testCases[1] = (struct exp){.input = TEST_INPUT_2,
                                .integer = TEST_INT_2,
                                .dbl = TEST_DBL_2};
    testCases[2] = (struct exp){.input = TEST_INPUT_3,
                                .integer = TEST_INT_3,
                                .dbl = TEST_DBL_3};
    testCases[3] = (struct exp){.input = TEST_INPUT_4,
                                .integer = TEST_INT_4,
                                .dbl = TEST_DBL_4};
    testCases[4] = (struct exp){.input = TEST_INPUT_5,
                                .integer = TEST_INT_5,
                                .dbl = TEST_DBL_5};
    // NOLINTEND(readability-magic-numbers)

    for (size_t testIndex = 0; testIndex < NR_TESTS; testIndex++)
    {
        struct exp testCase = testCases[testIndex];
        int intResult = getInt(testCase.input);
        double dblResult = getDouble(testCase.input);
        if (!nearlyEqual(dblResult, testCase.dbl))
        {
            if (isnan(dblResult) && isnan(testCase.dbl))
            {
                continue;
            }
            printf("test %zu:\n"
                   "input '%s' => (dbl)%lf : expected %lf\n",
                   testIndex,
                   testCase.input,
                   dblResult,
                   testCase.dbl);
            return -1;
        }
        if (intResult != testCase.integer)
        {
            printf("test %zu:\n"
                   "input '%s' => (int)%d :expected %d \n",
                   testIndex,
                   testCase.input,
                   intResult,
                   testCase.integer);
            return -1;
        }
    }

    puts("all fine");
    return 0;
}
