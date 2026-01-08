
#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    double Sw = 0.0;
    double Swx = 0.0;
    double Swy = 0.0;
    double Swxx = 0.0;
    double Swxy = 0.0;
    double triplet[3];
    size_t N = 0;

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
        Sw += w;
        Swx += w * x;
        Swy += w * y;
        Swxx += w * x * x;
        Swxy += w * x * y;
        N++;
    }

    // Solve weighted simple linear regression y = a + b x
    double denom = Sw * Swxx - Swx * Swx;
    double a = 0.0;
    double b = 0.0;
    if (denom != 0.0)
    {
        b = (Sw * Swxy - Swx * Swy) / denom;
        a = (Swy - b * Swx) / Sw;
    }
    else
    {
        (void)fprintf(stderr,
                      "Degenerate data: cannot regress (denominator zero)\n");
    }

    printf("%10.6g %10.6g\n", a, b);
    return 0;
}
