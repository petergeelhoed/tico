
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compare.h"

int main(void)
{
    double doublearray[3] = {1., -1., 0.};
    qsort(doublearray, 3, sizeof(double), compare_double);
    if (doublearray[0] > -1. || doublearray[1] > 0. || doublearray[2] > 1. ||
        doublearray[0] < -1. || doublearray[1] < 0. || doublearray[2] < 1.)
    {
        printf("%f %f %f\n", doublearray[0], doublearray[1], doublearray[2]);
        return (EXIT_FAILURE);
    }
    return (0);
}
