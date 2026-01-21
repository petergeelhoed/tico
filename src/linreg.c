#include "mydefs.h"
#include "mylib.h"
#include "parseargs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// this gives intercept liniar weight , not intercept squared one like gnuplot
// or matlinreg from mymath
int main(void)
{
    double Sum_w = 0.0;
    double Sum_wx = 0.0;
    double Sum_wy = 0.0;
    double Sum_wxx = 0.0;
    double Sum_wxy = 0.0;
    double triplet[3];

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
    }

    // Solve weighted simple linear regression y_arr = intercept + slope x_arr
    double denom = Sum_w * Sum_wxx - Sum_wx * Sum_wx;
    double intercept = 0.0;
    double slope = 0.0;
    if (fabs(denom) > DOUBLE_LIMIT)
    {
        slope = (Sum_w * Sum_wxy - Sum_wx * Sum_wy) / denom;
        intercept = (Sum_wy - slope * Sum_wx) / Sum_w;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }

    (void)printf("%10.6g %10.6g\n", intercept, slope);
    return 0;
}
