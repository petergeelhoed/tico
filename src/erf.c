#include "erf.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

// NOLINTBEGIN(readability-identifier-length)
static int cmp_double(const void* a, const void* b)
{
    double ua = *(const double*)a;
    double ub = *(const double*)b;

    if (ua < ub)
    {
        return -1;
    }
    if (ua > ub)
    {
        return 1;
    }
    return 0;
}

/* Acklam inverse normal CDF approximation */
double invnorm(double p)
{
    static const double a[] = {-3.969683028665376e+01,
                               2.209460984245205e+02,
                               -2.759285104469687e+02,
                               1.383577518672690e+02,
                               -3.066479806614716e+01,
                               2.506628277459239e+00};

    static const double b[] = {-5.447609879822406e+01,
                               1.615858368580409e+02,
                               -1.556989798598866e+02,
                               6.680131188771972e+01,
                               -1.328068155288572e+01};

    static const double c[] = {-7.784894002430293e-03,
                               -3.223964580411365e-01,
                               -2.400758277161838e+00,
                               -2.549732539343734e+00,
                               4.374664141464968e+00,
                               2.938163982698783e+00};

    static const double d[] = {7.784695709041462e-03,
                               3.224671290700398e-01,
                               2.445134137142996e+00,
                               3.754408661907416e+00};

    double q = 0.0;
    double r = 0.0;

    // NOLINTBEGIN(readability-magic-numbers)
    if (p < 0.02425)
    {
        q = sqrt(-2.0 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
                c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    if (p > 1.0 - 0.02425)
    {
        q = sqrt(-2.0 * log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
                 c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    q = p - 0.5;
    r = q * q;

    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r +
            a[5]) *
           q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}
// NOLINTEND(readability-magic-numbers)
// NOLINTEND(readability-identifier-length)

struct stats fit_erf(double* samples, size_t size)
{

    struct stats stats = {.stdev = 0.0, .mean = 0.0};

    qsort(samples, size, sizeof(samples[0]), cmp_double);

    double sumx = 0.0;
    double sumz = 0.0;
    double sumzz = 0.0;
    double sumzx = 0.0;

    for (size_t i = 0; i < size; i++)
    {
        /* Blom plotting position */
        // NOLINTNEXTLINE(readability-magic-numbers)
        double prob = ((double)i + 0.625) / ((double)size + 0.25);

        double zval = invnorm(prob);
        double xval = samples[i];

        sumx += xval;
        sumz += zval;
        sumzz += zval * zval;
        sumzx += zval * xval;
    }

    stats.stdev = ((double)size * sumzx - sumz * sumx) /
                  ((double)size * sumzz - sumz * sumz);
    stats.mean = (sumx - stats.stdev * sumz) / (double)size;

    if (stats.stdev < 0)
    {
        stats.stdev = -stats.stdev;
    }
    return stats;
}
