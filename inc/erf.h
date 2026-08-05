#include <math.h>
#include <stdint.h>
#include <stdlib.h>

// NOLINTBEGIN(readability-identifier-length)
int cmp_u64(const void* a, const void* b);
/* Acklam inverse normal CDF approximation */
double invnorm(double p);
// NOLINTEND(readability-identifier-length)

struct stats
{
    double stdev;
    double mean;
};

struct stats fit_erf(double* samples, size_t size);
