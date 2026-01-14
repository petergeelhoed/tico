
#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>
// this gives a liniar weight , not a squared one like gnuplot or matlinreg from
// mymath
int main(void)
{
    double Sum_w = 0.0;
    double Sum_wx = 0.0;
    double Sum_wy = 0.0;
    double Sum_wxx = 0.0;
    double Sum_wxy = 0.0;
    double triplet[3];
    size_t count = 0;

    for (;;)
    {
        int k = getDoublesFromStdin(3, triplet);
        if (k < 0)
        {
            break;
        }
        if (k < 3)
        {
            continue;
        }

        double x = triplet[0];
        double y = triplet[1];
        double w = triplet[2];
        Sum_w += w;
        Sum_wx += w * x;
        Sum_wy += w * y;
        Sum_wxx += w * x * x;
        Sum_wxy += w * x * y;
        count++;
    }

    // Sum_olve weighted simple linear regression y = a + b x
    double denom = Sum_w * Sum_wxx - Sum_wx * Sum_wx;
    double a = 0.0;
    double b = 0.0;
    if (denom != 0.0)
    {
        b = (Sum_w * Sum_wxy - Sum_wx * Sum_wy) / denom;
        a = (Sum_wy - b * Sum_wx) / Sum_w;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }

    printf("%10.6g %10.6g\n", a, b);
    return 0;
}
