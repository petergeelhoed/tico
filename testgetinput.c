#include "mylib.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF_SIZE 32
#define NR_TESTS 4

struct exp
{
    char input[BUF_SIZE];
    double dbl;
    int integer;
};

int main(void)
{
    struct exp array[NR_TESTS];
    // NOLINTBEGIN(readability-magic-numbers)
    array[0] = (struct exp){.input = " 34", .integer = 34, .dbl = 34};
    array[1] = (struct exp){.input = " 34.1", .integer = 34, .dbl = 34.1};
    array[2] = (struct exp){.input = "", .integer = INT_MIN, .dbl = NAN};
    array[3] = (struct exp){.input = " -1.7 ", .integer = -1, .dbl = -1.7};
    // NOLINTEND(readability-magic-numbers)

    for (size_t i = 0; i < NR_TESTS; i++)
    {
        struct exp test1 = array[i];
        int a = getInt(test1.input);
        double b = getDouble(test1.input);
        if (b != test1.dbl)
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
