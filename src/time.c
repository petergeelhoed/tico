#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 12
#define ISO_LENGTH 40

// NOLINTBEGIN(readability-identifier-length)
static int cmp_u64(const void* a, const void* b)
{
    uint64_t ua = *(const uint64_t*)a;
    uint64_t ub = *(const uint64_t*)b;

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
static double invnorm(double p)
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

struct stats
{
    double stdev;
    double mean;
};

static struct stats fit_erf(uint64_t* samples)
{

    struct stats stats = {.stdev = 0.0, .mean = 0.0};

    qsort(samples, N, sizeof(samples[0]), cmp_u64);

    double sumx = 0.0;
    double sumz = 0.0;
    double sumzz = 0.0;
    double sumzx = 0.0;

    for (int i = 0; i < N; i++)
    {

        /* Blom plotting position */
        // NOLINTNEXTLINE(readability-magic-numbers)
        double prob = (i + 0.625) / ((double)N + 0.25);

        double zval = invnorm(prob);
        double xval = (double)samples[i];

        sumx += xval;
        sumz += zval;
        sumzz += zval * zval;
        sumzx += zval * xval;
    }

    stats.stdev = (N * sumzx - sumz * sumx) / (N * sumzz - sumz * sumz);
    stats.mean = (sumx - stats.stdev * sumz) / N;

    if (stats.stdev < 0)
    {
        stats.stdev = -stats.stdev;
    }
    return stats;
}

int main(void)
{
    uint64_t samples[N];
    struct timespec tspec;
    const uint64_t mod = 5000;
    const uint64_t kilo = 1000ULL;
    const uint64_t mega = 1000000ULL;
    const long MEGA = 1000000UL;
    printf("Press ENTER %d times\n\n", N);

    for (int i = 0; i < N; i++)
    {
        printf("[%d/%d] ", i + 1, N);

        while (getchar() != '\n')
        {
        }

        int err = clock_gettime(CLOCK_REALTIME, &tspec);
        if (err)
        {
            exit(-1);
        }

        uint64_t msec =
            (uint64_t)tspec.tv_sec * kilo + (uint64_t)tspec.tv_nsec / mega;

        samples[i] = msec % mod;
        struct tm time;
        char iso8601[ISO_LENGTH];
        if (NULL == gmtime_r(&tspec.tv_sec, &time))
        {
            exit(-1);
        }
        if (0 == strftime(iso8601, sizeof(iso8601), "%Y-%m-%dT%H:%M:%S", &time))
        {
            exit(-1);
        }
        printf("%llu ms  %s.%03ldZ\n",
               (unsigned long long)samples[i],
               iso8601,
               tspec.tv_nsec / MEGA);
    }

    struct stats stats = fit_erf(samples);
    printf("\nSorted samples:\n");
    for (int i = 0; i < N; i++)
    {
        printf("%4llu ", (unsigned long long)samples[i]);
    }
    printf("\n\nQQ Gaussian fit:\n");
    printf("Mean   = %.3f ms\n", stats.mean);
    printf("Stddev = %.3f ms\n", stats.stdev);

    return 0;
}
