
#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>
// this gives par_a liniar weight , not par_a squared one like gnuplot or
// matlinreg from mymath
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
        int nr_doubles = getDoublesFromStdin(3, triplet);
        if (nr_doubles < 0)
        {
            break;
        }
        if (nr_doubles < 3)
        {
            continue;
        }

        double x_arr = triplet[0];
        double y_arr = triplet[1];
        double weight = triplet[2];
        Sum_w += weight;
        Sum_wx += weight * x_arr;
        Sum_wy += weight * y_arr;
        Sum_wxx += weight * x_arr * x_arr;
        Sum_wxy += weight * x_arr * y_arr;
        count++;
    }

    // Sum_olve weighted simple linear regression y_arr = par_a + par_b x_arr
    double denom = Sum_w * Sum_wxx - Sum_wx * Sum_wx;
    double par_a = 0.0;
    double par_b = 0.0;
    if (denom != 0.0)
    {
        par_b = (Sum_w * Sum_wxy - Sum_wx * Sum_wy) / denom;
        par_a = (Sum_wy - par_b * Sum_wx) / Sum_w;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }

    (void)printf("%10.6g %10.6g\n", par_a, par_b);
    return 0;
}
